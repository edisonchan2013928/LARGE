#pragma once
#ifndef SCAN_LINE_H
#define SCAN_LINE_H

#include "Plane.h"
#include "SCAN.h"
#include "LARGE.h"

void init_influence_pixels(statistics& stat);
void find_influence_pixels(statistics& stat);
void SCAN_line(statistics& stat);

#endif