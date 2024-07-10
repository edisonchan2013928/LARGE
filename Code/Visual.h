#pragma once
#ifndef VISUAL_H
#define VISUAL_H

#include "Library.h"
#include "Plane.h"
#include "SCAN.h"
//#include "kd_tree.h"
#include "R_tree.h"
#include "PMR_QUAD_Tree.h"
#include "LARGE.h"
#include "SCAN_line.h"

void init(statistics& stat, int argc, char** argv);
void visual_algorithm(statistics& stat);
void output_visual(statistics& stat);

#endif