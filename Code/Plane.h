#pragma once
#ifndef PLANE_H
#define PLANE_H

#include "Library.h"

const double inf = 9999999999;
const double eps = 0.0001;

//#define STATISTICS

struct Pixel
{
	double x_min;
	double x_max;
	double y_min;
	double y_max;

	double x_center;
	double y_center;

	double density_value;

	//Used in LARGE
	double total_length;
	int ER_region_x_index;
	int ER_region_y_index;
};

struct line_segment
{
	double x_start;
	double y_start;
	double x_end;
	double y_end;
	double m;
	double k;
};

struct Point
{
	double* p;
	int line_index;
};

struct statistics
{
	int n; //number of line segments
	int X; //number of pixels in a row
	int Y; //number of pixels in a column
	char* input_data_fileName; //The name of the input file
	char* output_visual_fileName; //The name of the output file
	int method; //method = 0: SCAN, method = 1: R* tree, method = 2: PMR Quad-tree, method = 3: LARGE
	//int query_type; //query_type = 0: epsilon-LDV, query_type = 1: tau-LDV
	Pixel** plane; //original plane
	line_segment* ls_dataset; //line segment dataset
	double bandwidth; //bandwidth parameter
	//Uses this variable if query_type = 0
	double epsilon;
	//Uses these variables if query_type = 1
	//int num_tau; //the number of thresholds
	//vector<double> tau_vector; //the vector of thresholds

	//boundary
	double boundary_x_min;
	double boundary_x_max;
	double boundary_y_min;
	double boundary_y_max;
	
	//pixel
	double incr_x;
	double incr_y;

	//indexing-based methods (kd-tree, R+ tree, and PMR QUAD-tree)
	//double* q;
	//Point* point_dataset;
	//vector<int> line_id_set;
	//bool* ls_set_select_vec;

	//Used in LARGE and SCAN_Line
	int cur_line_index;

	//LARGE
	/*int cur_x_index;
	int cur_y_index;
	double cur_x_value;
	double cur_y_value;*/
	//The extended region
	int ER_X;
	int ER_Y;
	Pixel** ER_plane;
	Pixel** ER_prefix_plane;
	Pixel** ER_one_D_prefix_plane; //Used in arbitrary-shaped LB and UB
	double ER_boundary_x_min;
	double ER_boundary_y_min;
	double ER_boundary_x_max;
	double ER_boundary_y_max;
	int os_appended_x_grid_num;
	int os_appended_y_grid_num;
	//intersection point
	double intersect_x;
	double intersect_y;
	int ER_cur_index_x;
	int ER_cur_index_y;
	//rectangular-shaped LB and UB
	int LB_expanded_index;
	int UB_expanded_index;
	//arbitrary-shaped LB and UB
	vector<int> LB_arbit_expanded_index_vec;
	vector<int> UB_arbit_expanded_index_vec;

	//SCAN_Line
	bool** is_selected_pixel_matrix;
	int influence_x_grid;
	int influence_y_grid;
	vector<pair<int, int>> influence_pixel_index_vec;
};

#endif