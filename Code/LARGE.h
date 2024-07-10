#pragma once
#ifndef LARGE_H
#define LARGE_H

#include "Plane.h"
#include "SCAN.h"
#include "R_tree.h"

void create_extended_region(statistics& stat); //Step 1
void obtain_accumulated_length(statistics& stat); //Step 2
void construct_prefix_structure(statistics& stat); //Step 3

double intersect_ED_determine(double line_x_start, double line_y_start, double x_coord, double y_coord, double threshold, statistics& stat);
void find_initial_intersection_point(statistics& stat);
void construct_LARGE(statistics& stat);

//Used for LB_{rectangle} and UB_{rectangle}
void obtain_rectangular_mask(statistics& stat);
double bound_rectangle(Pixel& p, statistics& stat, bool is_LB);

//Used for LB_{arbitrary} and UB_{arbitrary}
void obtain_arbitrary_mask(statistics& stat);
double bound_arbit(Pixel& p, statistics& stat, vector<int>& expanded_index_vec);

//void filter_and_refinement(statistics& stat);
void filter_and_refinement(statistics& stat, R_tree& R_tree);
void bound_visual(statistics& stat, R_tree& R_tree);

//For statistical purpose
#ifdef STATISTICS
void filter_and_refinement_STAT(statistics& stat, R_tree& R_tree);
#endif

//Used for debuging
void output_LARGE(statistics& stat);
#endif