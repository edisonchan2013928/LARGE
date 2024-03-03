#pragma once
#ifndef PMR_QUAD_TREE_H
#define PMR_QUAD_TREE_H

#include "tree.h"

struct pmr_quad_node : Node
{
	pmr_quad_node();
	int search(double* q, double bandwidth);
	bool check_line_intersection(line_segment* ls_dataset, int lid);
	pmr_quad_node* createNode();
};

struct pmr_quad_tree : Tree
{
	int ls_num_split; //number of line segments for splitting the grid
	vector<int> result_idList;
	bool* sel_id_vec;

	pmr_quad_tree(line_segment* ls_dataset, int n, int max_ls_num);
	void insert_one_line_segment(pmr_quad_node* node, int lid);
	void construct_pmr_quad_tree();

	void obtain_idList(double* q, double bandwidth, pmr_quad_node* node);
	double compute_LDF(double* q, Pixel& p, double bandwidth);
	void generate_LDV(statistics& stat);
};

#endif