#pragma once
#ifndef SCAN_H
#define SCAN_H

#include "Plane.h"

double Euclidean_distance(double x_1, double y_1, double x_2, double y_2);
double compute_length(Pixel& p, line_segment& l, double bandwidth);
void SCAN_one_pixel(statistics& stat, int x_index, int y_index);
void SCAN(statistics& stat);

#endif