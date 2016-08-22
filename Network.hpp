// Author : Joose Rajamäki
// Date: 22.8.2016

#ifndef NETWORK_HPP_
#define NETWORK_HPP_

#include <map>
#include <list>
#include <memory>
#include <cmath>
#include <complex>
#include <Dense>

#include "Impedance.hpp"
#include "Node.hpp"

template<typename Scalar>
class Network{

public:

	typedef Node<Scalar>* Node_ptr;
	typedef std::pair<Node_ptr,Node_ptr > Node_ptr_pair;
	typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> Matrix;
	typedef Eigen::Matrix<std::complex<Scalar>, Eigen::Dynamic, Eigen::Dynamic> Complex_matrix;
	Eigen::Matrix<Scalar,Eigen::Dynamic,1> Vector;

	//These are the impedances, connecting the nodes.
	std::map<Node_ptr_pair, Impedance<Scalar>* > connections_;

	//The storage classes use unique pointers so that we get the performance benefits of using pointers but they are automatically deleted.
	//The nodes between which the impedances are.
	std::list<std::unique_ptr<Node<Scalar> > > nodes_;
	//These are the impedances.
	std::list<std::unique_ptr<Impedance<Scalar> > > impedances_;

	int reference_node_index_; //TODO Implement everything such that the reference node can be changed. Now the system can "leak" in resistance computation.

	//The admittance matrix of the system. The admittances are relative to the zero node.
	Complex_matrix admittance_matrix_;

	Network(){
		connections_ = std::map<Node_ptr_pair, Impedance<Scalar>* >();
		nodes_ = std::list<std::unique_ptr<Node<Scalar> > >();
		impedances_ = std::list<std::unique_ptr<Impedance<Scalar> > >();
		admittance_matrix_ = Complex_matrix::Zero(0,0);
		reference_node_index_ = 0;
	}

	Network(Network& other){

		reference_node_index_ = other.reference_node_index_;
		admittance_matrix_ = other.admittance_matrix_;
		connections_ = std::map<Node_ptr_pair, Impedance<Scalar>* >();
		nodes_ = std::list<std::unique_ptr<Node<Scalar> > >();
		impedances_ = std::list<std::unique_ptr<Impedance<Scalar> > >();

		for (unsigned i = 0; i < other.nodes_.size(); i++){
			add_node();
		}

		auto node_iterator = nodes_.begin();
		for (const std::unique_ptr<Node<Scalar> >& node : other.nodes_){
			**(node_iterator++) = *node;
		}

		int i = 0;
		for (const std::unique_ptr<Node<Scalar> >& node_i : other.nodes_){
			int j = 0;
			for (const std::unique_ptr<Node<Scalar> >& node_j : other.nodes_){

				if (i == j){
					continue;
				}

				const Node_ptr node_one = node_i.get();
				const Node_ptr node_two = node_j.get();

				connect_nodes(i,j,*(other.get_impedance(node_one,node_two)));
				j++;
			}
			i++;
		}

	}

	Network<Scalar> operator=(Network<Scalar>& other){
		reference_node_index_ = other.reference_node_index_;
			admittance_matrix_ = other.admittance_matrix_;
			connections_ = std::map<Node_ptr_pair, Impedance<Scalar>* >();
			nodes_ = std::list<std::unique_ptr<Node<Scalar> > >();
			impedances_ = std::list<std::unique_ptr<Impedance<Scalar> > >();

			for (unsigned i = 0; i < other.nodes_.size(); i++){
				add_node();
			}

			auto node_iterator = nodes_.begin();
			for (const std::unique_ptr<Node<Scalar> >& node : other.nodes_){
				**(node_iterator++) = *node;
			}

			int i = 0;
			for (const std::unique_ptr<Node<Scalar> >& node_i : other.nodes_){
				int j = 0;
				for (const std::unique_ptr<Node<Scalar> >& node_j : other.nodes_){

					if (i == j){
						continue;
					}

					const Node_ptr node_one = node_i.get();
					const Node_ptr node_two = node_j.get();

					connect_nodes(i,j,*(other.get_impedance(node_one,node_two)));
					j++;
				}
				i++;
			}

			return *this;
	}


	void index_nodes(void){
		int i = 0;
		for (std::unique_ptr<Node<Scalar> >& node : nodes_){
			node->index_ = i++;
		}
	}

	void add_node(void){
		nodes_.push_back(std::move(std::unique_ptr<Node<Scalar> >(new Node<Scalar>())));
		index_nodes();

		int last_node_idx = nodes_.size() - 1;
		for (int ii = 0; ii < last_node_idx; ii++){
			connect_nodes(ii,last_node_idx,std::numeric_limits<Scalar>::infinity());
		}
	}

	void set_reference_node(int new_reference_node){
		reference_node_index_ = new_reference_node;
	}

	Node_ptr find_node(int node_index){
		Node_ptr node_ptr = nullptr;
		for (std::unique_ptr<Node<Scalar> >& node : nodes_){
			if (node->index_ == node_index){
				node_ptr = node.get();
				break;
			}
		}
		return node_ptr;
	}



	void remove_unconnected_impedances(void){
		auto list_iter = impedances_.begin();

		while(list_iter != impedances_.end()){

			auto map_iterator = connections_.begin();
			bool found = false;
			while (map_iterator != connections_.end()){

				if (list_iter->get() == map_iterator->second){
					found = true;
					break;
				}
				map_iterator++;

			}

			if (!found){
				list_iter = impedances_.erase(list_iter);
			}
			else{
				list_iter++;
			}


		}

	}

	//Connects nodes. Default resistor with unit resistance.
	void connect_nodes(int node_one_index, int node_two_index, Scalar resistance = 1, Scalar inductance = 0, Scalar capacitance = 0){

		Node_ptr node_one = find_node(node_one_index);
		Node_ptr node_two = find_node(node_two_index);

		std::unique_ptr<Impedance<Scalar> > impedance = std::unique_ptr<Impedance<Scalar> >(new Impedance<Scalar>(resistance,inductance,capacitance));


		if (node_one < node_two){
			connections_[std::make_pair(node_one,node_two)] = impedance.get();
		}
		else{
			connections_[std::make_pair(node_two,node_one)] = impedance.get();
		}


		impedances_.push_back(std::move(impedance));

		remove_unconnected_impedances();

	}

	//Connects nodes. Default resistor with unit resistance.
	void connect_nodes(int node_one_index, int node_two_index, const Impedance<Scalar>& new_impedance){

		Node_ptr node_one = find_node(node_one_index);
		Node_ptr node_two = find_node(node_two_index);

		std::unique_ptr<Impedance<Scalar> > impedance = std::unique_ptr<Impedance<Scalar> >(new Impedance<Scalar>(new_impedance));

		if (node_one < node_two){
			connections_[std::make_pair(node_one,node_two)] = impedance.get();
		}
		else{
			connections_[std::make_pair(node_two,node_one)] = impedance.get();
		}


		impedances_.push_back(std::move(impedance));

		remove_unconnected_impedances();

	}


	Impedance<Scalar>* get_impedance(const Node_ptr node_one, const Node_ptr node_two){

		if (node_one < node_two){
			return connections_[std::make_pair(node_one,node_two)];
		}
		else{
			return connections_[std::make_pair(node_two,node_one)];
		}

	}


	void compute_impedances(const Scalar frequency){

		for (std::unique_ptr<Impedance<Scalar> >& impedance : impedances_){
			impedance->compute_impedance(frequency);
		}

	}

	//Returns the index that corresponds to the row and column of the impedance/admittance matrix.
	int get_node_matrix_index(int node_idx){

		if (node_idx < reference_node_index_){
			return node_idx;
		}
		else{
			if (node_idx == reference_node_index_){
				return -1;
			}
			else{
				return node_idx - 1;
			}
		}

	}

	void compute_admittance_matrix(const Scalar frequency){

		compute_impedances(frequency);

		admittance_matrix_.resize(nodes_.size()-1,nodes_.size()-1);
		admittance_matrix_.setZero();

		//Diagonal elements
		for (std::unique_ptr<Node<Scalar> >& ref_node : nodes_){
			int diagonal_index = get_node_matrix_index(ref_node->index_);
			if (diagonal_index < 0){
				continue;
			}

			for (std::unique_ptr<Node<Scalar> >& connected_node : nodes_){

				if (ref_node.get() == connected_node.get()){
					continue;
				}

				Impedance<Scalar>* impedance_ptr = get_impedance(ref_node.get(),connected_node.get());
				std::complex<Scalar> impedance  = impedance_ptr->get_impedance();

				if (impedance.real() < std::numeric_limits<Scalar>::infinity()){
					admittance_matrix_(diagonal_index,diagonal_index) += 1.0/impedance;
				}

			}
		}


		for (std::unique_ptr<Node<Scalar> >& ref_node : nodes_){
			int row = get_node_matrix_index(ref_node->index_);
			if (row < 0){
				continue;
			}

			for (std::unique_ptr<Node<Scalar> >& connected_node : nodes_){

				if (ref_node.get() == connected_node.get()){
					continue;
				}

				int col = get_node_matrix_index(connected_node->index_);

				if (col < 0){
					continue;
				}

				Impedance<Scalar>* impedance_ptr = get_impedance(ref_node.get(),connected_node.get());
				std::complex<Scalar> impedance  = impedance_ptr->get_impedance();

				if (impedance.real() < std::numeric_limits<Scalar>::infinity()){
					admittance_matrix_(row,col) -= 1.0/impedance;
				}

			}
		}


	}


	//Compute the voltages of the nodes when unit current is fed to the excited node.
	Eigen::Matrix<std::complex<Scalar>,Eigen::Dynamic,1> compute_nodal_voltages(Node_ptr excited_node, Scalar frequency){

		compute_admittance_matrix(frequency);

		int excitation_idx = get_node_matrix_index(excited_node->index_);

		if (excitation_idx < 0){
			return Eigen::Matrix<std::complex<Scalar>,Eigen::Dynamic,1>::Zero(nodes_.size());
		}

		Eigen::Matrix<std::complex<Scalar>,Eigen::Dynamic,1> excitation = Eigen::Matrix<std::complex<Scalar>,Eigen::Dynamic,1>::Zero(nodes_.size()-1);

		excitation(excitation_idx) = std::complex<Scalar>(1.0,0.0);

		Eigen::Matrix<std::complex<Scalar>,Eigen::Dynamic,1> voltages = admittance_matrix_.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(excitation);

		voltages.conservativeResize(voltages.size()+1);
		for(int ii = voltages.size()-1; ii > reference_node_index_; ii--){
			voltages(ii) = voltages(ii-1);
		}
		voltages(reference_node_index_) = 0;

		return voltages;


	}

	//Compute the voltages of the nodes when unit current is fed to the excited node.
	Eigen::Matrix<std::complex<Scalar>,Eigen::Dynamic,1> compute_nodal_voltages(int excited_node, Scalar frequency){
		return compute_nodal_voltages(find_node(excited_node), frequency);
	}



};


#endif /* NETWORK_HPP_ */
