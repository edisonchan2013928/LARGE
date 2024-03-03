#include "SCAN.h"

double Euclidean_distance(double x_1, double y_1, double x_2, double y_2)
{
	double dist;
	dist = sqrt((x_1 - x_2) * (x_1 - x_2) + (y_1 - y_2) * (y_1 - y_2));

	return dist;
}

double compute_length(Pixel& p, line_segment& l, double bandwidth)
{
	double dist_left;
	double dist_right;
	double length = 0;
	double A, B, C; //Used for the quadratic equation
	double ell, u;

	dist_left = Euclidean_distance(p.x_center, p.y_center, l.x_start, l.y_start);
	dist_right = Euclidean_distance(p.x_center, p.y_center, l.x_end, l.y_end);

	if (dist_left <= bandwidth && dist_right <= bandwidth) //Case 1
		length = Euclidean_distance(l.x_start, l.y_start, l.x_end, l.y_end);
	
	A = 1 + l.m * l.m;
	B = 2 * l.m * l.k - 2 * p.x_center - 2 * l.m * p.y_center;
	C = p.x_center * p.x_center + l.k * l.k - 2 * l.k * p.y_center
		+ p.y_center * p.y_center - bandwidth * bandwidth;

	if (dist_left <= bandwidth && dist_right > bandwidth) //Case 2
	{
		u = (-B + sqrt(B * B - 4 * A * C)) / (2.0 * A);
		length = sqrt(1 + l.m * l.m) * fabs(u - l.x_start);
	}

	if (dist_left > bandwidth && dist_right <= bandwidth) //Case 3
	{
		ell = (-B - sqrt(B * B - 4 * A * C)) / (2.0 * A);
		length = sqrt(1 + l.m * l.m) * fabs(l.x_end - ell);
	}

	if (dist_left > bandwidth && dist_right > bandwidth) //Case 4
	{
		if (B * B - 4 * A * C < 0) //Case 4a
			length = 0;

		if (B * B - 4 * A * C >= 0) 
		{
			ell = (-B - sqrt(B * B - 4 * A * C)) / (2.0 * A);
			u = (-B + sqrt(B * B - 4 * A * C)) / (2.0 * A);

			if (u < l.x_start || ell > l.x_end) //Case 4b
				length = 0;
			else
				length = sqrt(1 + l.m * l.m) * fabs(u - ell); //Case 4c
		}
	}

	return length;
}

void SCAN_one_pixel(statistics& stat, int x_index, int y_index)
{
	Pixel& p = stat.plane[x_index][y_index];
	for (int i = 0;i < stat.n;i++)
	{
		line_segment& l = stat.ls_dataset[i];
		p.density_value += compute_length(p, l, stat.bandwidth);
	}
}

void SCAN(statistics& stat)
{
	for (int x = 0;x < stat.X;x++)
		for (int y = 0;y < stat.Y;y++)
			SCAN_one_pixel(stat, x, y);
}