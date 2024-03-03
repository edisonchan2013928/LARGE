#pragma once
#ifndef R_TREE_H
#define R_TREE_H

#include "tree.h"

struct R_node : Node
{
	double* center;

	//Used for the line segment
	//int id;

	R_node* createNode();
};

struct R_tree : Tree
{
	vector<R_node*> ls_rect_dataset;
	double bandwidth;

	R_tree(int n, line_segment* ls_dataset, int leafCapacity, double bandwidth);
	void construct_ls_rect_data();
	void build_RTree();
	double refinement(Pixel& p, R_node* r_node);
	double compute_LDF(Pixel& p, R_node* r_node);
	void generate_LDV(statistics& stat);

	void R_Tree_Recur(R_node* r_node, int split_Dim);
	void update_R_node_info(R_node* r_node);
	double obtain_SplitValue(R_node* r_node, int split_Dim);
};

#endif