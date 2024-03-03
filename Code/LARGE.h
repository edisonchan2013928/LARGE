#pragma once
#ifndef LARGE_H
#define LARGE_H

#include "Plane.h"
#include "SCAN.h"
#include "R_tree.h"

void create_extended_region(statistics& stat); //Step 1
void obtain_accumulated_length(statistics& stat); //Step 2
void construct_prefix_structure(statistics& stat); //Step 3
void construct_LARGE(statistics& stat);

//Used for LB_{rectangle} and UB_{rectangle}
void obtain_rectangular_mask(statistics& stat);
double bound_rectangle(Pixel& p, statistics& stat, bool is_LB);

//Used for LB_{arbitrary} and UB_{arbitrary}
void obtain_arbitrary_mask(statistics& stat);
double bound_arbit(Pixel& p, statistics& stat, vector<int>& expanded_index_vec);

//void filter_and_refinement(statistics& stat);
void filter_and_refinement(statistics& stat, R_tree& R_tree);

//Used for debuging
void output_LARGE(statistics& stat);
#endif