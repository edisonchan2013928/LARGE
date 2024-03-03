#pragma once
#ifndef TREE_H
#define TREE_H

#include "SCAN.h"
#include "Plane.h"

const int leafCapacity = 10;
class Node;
class Tree;

class Node
{
public:
	vector<int> idList; //list of line segment ids
	vector<Node*> childVector;
	double** boundary; //boundary of this node

	virtual Node* createNode() = 0;
	//virtual double lb(double* q) = 0;
	//virtual void update_Aug(Node* node, Tree* t) = 0;
};

class Tree
{
public:
	//int dim;
	int n;
	line_segment* ls_dataset;
	//Point* point_dataset;
	int leafCapacity; //Set to be 10 by default
	Node* rootNode;
};

#endif