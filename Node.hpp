// Author : Joose Rajamäki
// Date: 22.8.2016

#ifndef NODE_HPP_
#define NODE_HPP_

#include <vector>

template<typename Scalar>
class Node{

public:

	std::vector<Node<Scalar>* > connected_nodes_;
	int index_;

	Node(){
		connected_nodes_ = std::vector<Node<Scalar>* >();
		index_ = 0;
	}

};


#endif
