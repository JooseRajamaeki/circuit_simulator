// Author : Joose Rajamäki
// Date: 22.8.2016

#include <iostream>
#include <Dense> //Includes Eigen-library's dense matrices.
#include <memory>
#include "Network.hpp"

int main() {

	//First test, testing three 2 mF capacitors in parallel

	Network<double> parallel_network;

	parallel_network.add_node();
	parallel_network.add_node();
	parallel_network.add_node();
	parallel_network.add_node();

	int reference_node_idx = 1;

	parallel_network.set_reference_node(reference_node_idx);

	double frequency = 1e6;

	std::cout << "Frequency is : " << frequency/(1e6) << " MHz" << std::endl;

	double resistance = 1e-15;
	double capacitance = 0;
	double inductance = 0;

	//Short circuit nodes 0, 1 and 2.
	parallel_network.connect_nodes(0,1,resistance,inductance,capacitance);
	parallel_network.connect_nodes(0,2,resistance,inductance,capacitance);
	parallel_network.connect_nodes(1,2,resistance,inductance,capacitance);


	resistance = 0;
	capacitance = 0.002;
	inductance = 0;

	//Connect nodes 0, 1 and 2 to node 3 with a capacitor.
	parallel_network.connect_nodes(0,3,resistance,inductance,capacitance);
	parallel_network.connect_nodes(1,3,resistance,inductance,capacitance);
	parallel_network.connect_nodes(2,3,resistance,inductance,capacitance);

	int excited_node_index = 3;
	std::cout << "Nodal voltages when node "  << excited_node_index << " is fed with unit current and node " << parallel_network.reference_node_index_ <<" grounded:" << std::endl;
	Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> voltages = parallel_network.compute_nodal_voltages(excited_node_index,frequency);
	std::cout << voltages << std::endl << std::endl;

	std::complex<double> capacitance_tmp = voltages(excited_node_index); //Impedance
	capacitance_tmp = 1.0/capacitance_tmp; //Admittance

	std::cout << "Capacitance between nodes " << parallel_network.reference_node_index_ << " and " << excited_node_index << " in parallel network is: " << capacitance_tmp.imag()/(2*M_PI*frequency) << std::endl;

	for (int i = 0; i < 3 ; i++){
		std::cout << std::endl;
	}



	//Second test testing three 2 mF capacitors in series.

	std::unique_ptr<Network<double> > series_network = std::unique_ptr<Network<double> >(new Network<double>());

	series_network->add_node();
	series_network->add_node();
	series_network->add_node();
	series_network->add_node();

	reference_node_idx = 0;

	series_network->set_reference_node(reference_node_idx);

	resistance = 0;
	capacitance = 0.002;
	inductance = 0;

	series_network->connect_nodes(0,1,resistance,inductance,capacitance);
	series_network->connect_nodes(1,2,resistance,inductance,capacitance);
	series_network->connect_nodes(2,3,resistance,inductance,capacitance);

	excited_node_index = 3;
	std::cout << "Nodal voltages when node "  << excited_node_index << " is fed with unit current and node " << series_network->reference_node_index_ <<" grounded:" << std::endl;
	voltages = series_network->compute_nodal_voltages(excited_node_index,frequency);
	std::cout << voltages << std::endl << std::endl;

	capacitance_tmp = voltages(excited_node_index); //Impedance
	capacitance_tmp = 1.0/capacitance_tmp; //Admittance

	std::cout << "Capacitance between nodes " << series_network->reference_node_index_ << " and " << excited_node_index << " in series network is: " << capacitance_tmp.imag()/(2*M_PI*frequency) << std::endl;

	for (int i = 0; i < 3 ; i++){
		std::cout << std::endl;
	}


	frequency = 1e4;

	//Third test cloning the series network.

	//Network<double> clone_network(*series_network);

	Network<double> clone_network;
	clone_network = *series_network;

	series_network.reset(nullptr);

	excited_node_index = 3;
	std::cout << "Nodal voltages when node "  << excited_node_index << " is fed with unit current and node " << clone_network.reference_node_index_ <<" grounded:" << std::endl;
	voltages = clone_network.compute_nodal_voltages(excited_node_index,frequency);
	std::cout << voltages << std::endl << std::endl;

	capacitance_tmp = voltages(excited_node_index); //Impedance
	capacitance_tmp = 1.0/capacitance_tmp; //Admittance

	std::cout << "Capacitance between nodes " << clone_network.reference_node_index_ << " and " << excited_node_index << " in cloned network is: " << capacitance_tmp.imag()/(2*M_PI*frequency) << std::endl;






	return 0;
}
