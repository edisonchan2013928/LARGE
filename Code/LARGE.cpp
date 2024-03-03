#include "LARGE.h"

//Step 1
void create_extended_region(statistics& stat)
{
	double cur_x;
	double cur_y;
	stat.os_appended_x_grid_num = (int)ceil((stat.bandwidth - (stat.incr_x / 2.0)) / stat.incr_x);
	stat.os_appended_y_grid_num = (int)ceil((stat.bandwidth - (stat.incr_y / 2.0)) / stat.incr_y);
	stat.ER_boundary_x_min = stat.boundary_x_min - stat.os_appended_x_grid_num * stat.incr_x;
	stat.ER_boundary_x_max = stat.boundary_x_max + stat.os_appended_x_grid_num * stat.incr_x;
	stat.ER_boundary_y_min = stat.boundary_y_min - stat.os_appended_y_grid_num * stat.incr_y;
	stat.ER_boundary_y_max = stat.boundary_y_max + stat.os_appended_y_grid_num * stat.incr_y;

	stat.ER_X = stat.X + 2 * stat.os_appended_x_grid_num;
	stat.ER_Y = stat.Y + 2 * stat.os_appended_y_grid_num;
	stat.ER_plane = new Pixel*[stat.ER_X];
	stat.ER_prefix_plane = new Pixel*[stat.ER_X];
	stat.ER_one_D_prefix_plane = new Pixel*[stat.ER_X];
	for (int x = 0;x < stat.ER_X;x++)
	{
		stat.ER_plane[x] = new Pixel[stat.ER_Y];
		stat.ER_prefix_plane[x] = new Pixel[stat.ER_Y];
		stat.ER_one_D_prefix_plane[x] = new Pixel[stat.ER_Y];
	}

	for (int x = 0;x < stat.ER_X;x++)
	{
		cur_x = stat.ER_boundary_x_min + x * stat.incr_x;
		for (int y = 0;y < stat.ER_Y;y++)
		{
			cur_y = stat.ER_boundary_y_min + y * stat.incr_y;
			stat.ER_plane[x][y].x_min = cur_x;
			stat.ER_plane[x][y].x_max = cur_x + stat.incr_x;
			stat.ER_plane[x][y].y_min = cur_y;
			stat.ER_plane[x][y].y_max = cur_y + stat.incr_y;
			stat.ER_plane[x][y].x_center = cur_x + stat.incr_x / 2.0;
			stat.ER_plane[x][y].y_center = cur_y + stat.incr_y / 2.0;
			stat.ER_plane[x][y].density_value = 0;
			stat.ER_plane[x][y].total_length = 0;

			stat.ER_prefix_plane[x][y].x_min = cur_x;
			stat.ER_prefix_plane[x][y].x_max = cur_x + stat.incr_x;
			stat.ER_prefix_plane[x][y].y_min = cur_y;
			stat.ER_prefix_plane[x][y].y_max = cur_y + stat.incr_y;
			stat.ER_prefix_plane[x][y].x_center = cur_x + stat.incr_x / 2.0;
			stat.ER_prefix_plane[x][y].y_center = cur_y + stat.incr_y / 2.0;
			stat.ER_prefix_plane[x][y].density_value = 0;
			stat.ER_prefix_plane[x][y].total_length = 0;

			stat.ER_one_D_prefix_plane[x][y].x_min = cur_x;
			stat.ER_one_D_prefix_plane[x][y].x_max = cur_x + stat.incr_x;
			stat.ER_one_D_prefix_plane[x][y].y_min = cur_y;
			stat.ER_one_D_prefix_plane[x][y].y_max = cur_y + stat.incr_y;
			stat.ER_one_D_prefix_plane[x][y].x_center = cur_x + stat.incr_x / 2.0;
			stat.ER_one_D_prefix_plane[x][y].y_center = cur_y + stat.incr_y / 2.0;
			stat.ER_one_D_prefix_plane[x][y].density_value = 0;
			stat.ER_one_D_prefix_plane[x][y].total_length = 0;
		}
	}

	for (int x = 0;x < stat.X;x++)
	{
		for (int y = 0;y < stat.Y;y++)
		{
			stat.plane[x][y].ER_region_x_index = x + stat.os_appended_x_grid_num;
			stat.plane[x][y].ER_region_y_index = y + stat.os_appended_y_grid_num;
		}
	}
}

double intersect_ED_determine(double line_x_start, double line_y_start, double x_coord, double y_coord, double threshold, statistics& stat)
{
	double ED = Euclidean_distance(line_x_start, line_y_start, x_coord, y_coord);
	if (ED > threshold) //Cannot intersect the boundary
	{
		stat.intersect_x = -inf - 1;
		stat.intersect_y = -inf - 1;
	}
	else
	{
		stat.intersect_x = x_coord;
		stat.intersect_y = y_coord;
	}

	return ED;
}

void find_initial_intersection_point(statistics& stat)
{
	int cur_line_index = stat.cur_line_index;
	line_segment* ls_dataset = stat.ls_dataset;
	double x_coord, y_coord;
	double m = ls_dataset[cur_line_index].m;
	double k = ls_dataset[cur_line_index].k;

	//Region V
	if ((ls_dataset[cur_line_index].x_start >= stat.ER_boundary_x_min && ls_dataset[cur_line_index].x_start <= stat.ER_boundary_x_max) &&
		(ls_dataset[cur_line_index].y_start >= stat.ER_boundary_y_min && ls_dataset[cur_line_index].y_start <= stat.ER_boundary_y_max))
	{
		stat.intersect_x = ls_dataset[cur_line_index].x_start;
		stat.intersect_y = ls_dataset[cur_line_index].y_start;
		return;
	}
	
	//Region I
	if (ls_dataset[cur_line_index].x_start < stat.ER_boundary_x_min && ls_dataset[cur_line_index].y_start > stat.ER_boundary_y_max)
	{
		x_coord = (stat.ER_boundary_y_max - k) / m;
		y_coord = stat.ER_boundary_y_max;
		if (x_coord >= stat.ER_boundary_x_min && x_coord <= stat.ER_boundary_x_max)
		{
			intersect_ED_determine(ls_dataset[cur_line_index].x_start, ls_dataset[cur_line_index].y_start, x_coord, y_coord, stat.bandwidth, stat);
			return;
		}
		
		x_coord = stat.ER_boundary_x_min;
		y_coord = m * x_coord + k;
		if (y_coord >= stat.ER_boundary_y_min && y_coord <= stat.ER_boundary_y_max)
		{
			intersect_ED_determine(ls_dataset[cur_line_index].x_start, ls_dataset[cur_line_index].y_start, x_coord, y_coord, stat.bandwidth, stat);
			return;
		}

		stat.intersect_x = -inf - 1;
		stat.intersect_y = -inf - 1;
		return;
	}

	//Region II
	if (ls_dataset[cur_line_index].x_start < stat.ER_boundary_x_min &&
		(ls_dataset[cur_line_index].y_start >= stat.ER_boundary_y_min && ls_dataset[cur_line_index].y_start <= stat.ER_boundary_y_max))
	{
		x_coord = stat.ER_boundary_x_min;
		y_coord = m * x_coord + k;
		if (y_coord >= stat.ER_boundary_y_min && y_coord <= stat.ER_boundary_y_max)
		{
			intersect_ED_determine(ls_dataset[cur_line_index].x_start, ls_dataset[cur_line_index].y_start, x_coord, y_coord, stat.bandwidth, stat);
			return;
		}
		
		stat.intersect_x = -inf - 1;
		stat.intersect_y = -inf - 1;
		return;
	}

	//Region III
	if (ls_dataset[cur_line_index].x_start < stat.ER_boundary_x_min && ls_dataset[cur_line_index].y_start < stat.ER_boundary_y_min)
	{
		x_coord = (stat.ER_boundary_y_min - k) / m;
		y_coord = stat.ER_boundary_y_min;
		if (x_coord >= stat.ER_boundary_x_min && x_coord <= stat.ER_boundary_x_max)
		{
			intersect_ED_determine(ls_dataset[cur_line_index].x_start, ls_dataset[cur_line_index].y_start, x_coord, y_coord, stat.bandwidth, stat);
			return;
		}

		x_coord = stat.ER_boundary_x_min;
		y_coord = m * x_coord + k;
		if (y_coord >= stat.ER_boundary_y_min && y_coord <= stat.ER_boundary_y_max)
		{
			intersect_ED_determine(ls_dataset[cur_line_index].x_start, ls_dataset[cur_line_index].y_start, x_coord, y_coord, stat.bandwidth, stat);
			return;
		}

		stat.intersect_x = -inf - 1;
		stat.intersect_y = -inf - 1;
		return;
	}

	//Region IV
	if ((ls_dataset[cur_line_index].x_start >= stat.ER_boundary_x_min && ls_dataset[cur_line_index].x_start <= stat.ER_boundary_x_max) &&
		ls_dataset[cur_line_index].y_start > stat.ER_boundary_y_max)
	{
		if (m > 0)
		{
			stat.intersect_x = -inf - 1;
			stat.intersect_y = -inf - 1;
			return;
		}

		x_coord = (stat.ER_boundary_y_max - k) / m;
		y_coord = stat.ER_boundary_y_max;
		if (x_coord >= stat.ER_boundary_x_min && x_coord <= stat.ER_boundary_x_max)
		{
			intersect_ED_determine(ls_dataset[cur_line_index].x_start, ls_dataset[cur_line_index].y_start, x_coord, y_coord, stat.bandwidth, stat);
			return;
		}

		stat.intersect_x = -inf - 1;
		stat.intersect_y = -inf - 1;
		return;
	}

	//Region VI
	if ((ls_dataset[cur_line_index].x_start >= stat.ER_boundary_x_min && ls_dataset[cur_line_index].x_start <= stat.ER_boundary_x_max) &&
		ls_dataset[cur_line_index].y_start < stat.ER_boundary_y_min)
	{
		if (m < 0)
		{
			stat.intersect_x = -inf - 1;
			stat.intersect_y = -inf - 1;
			return;
		}

		x_coord = (stat.ER_boundary_y_min - k) / m;
		y_coord = stat.ER_boundary_y_min;
		if (x_coord >= stat.ER_boundary_x_min && x_coord <= stat.ER_boundary_x_max)
		{
			intersect_ED_determine(ls_dataset[cur_line_index].x_start, ls_dataset[cur_line_index].y_start, x_coord, y_coord, stat.bandwidth, stat);
			return;
		}
		stat.intersect_x = -inf - 1;
		stat.intersect_y = -inf - 1;
		return;
	}

	//Region VII
	stat.intersect_x = -inf - 1;
	stat.intersect_y = -inf - 1;
	return;
}

void obtain_line_accumulated_length(statistics& stat)
{
	int cur_line_index = stat.cur_line_index;
	line_segment* ls_dataset = stat.ls_dataset;
	double x_coord, y_coord;
	double ED;
	double m = ls_dataset[cur_line_index].m;
	double k = ls_dataset[cur_line_index].k;
	double threshold;
	int ER_cur_index_x = stat.ER_cur_index_x;
	int ER_cur_index_y = stat.ER_cur_index_y;

	if (ER_cur_index_x >= stat.ER_X || ER_cur_index_y < 0 || ER_cur_index_y >= stat.ER_Y)
		return;

	//Case (1) and (2)
	if (m > 0)
	{
		y_coord = stat.ER_plane[ER_cur_index_x][ER_cur_index_y].y_max;
		x_coord = (y_coord - ls_dataset[cur_line_index].k) / ls_dataset[cur_line_index].m;

		if (x_coord >= stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_min && x_coord <= stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_max)
		{
			threshold = Euclidean_distance(stat.intersect_x, stat.intersect_y, ls_dataset[cur_line_index].x_end, ls_dataset[cur_line_index].y_end);
			ED = intersect_ED_determine(stat.intersect_x, stat.intersect_y, x_coord, y_coord, threshold, stat);

			if (stat.intersect_x < -inf && stat.intersect_y < -inf)
			{
				stat.ER_plane[ER_cur_index_x][ER_cur_index_y].total_length += threshold;
				return;
			}
			else
			{
				stat.ER_plane[ER_cur_index_x][ER_cur_index_y].total_length += ED;
				stat.ER_cur_index_y++;

				//Case (2)
				if (x_coord == stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_max)
					stat.ER_cur_index_x++;

				obtain_line_accumulated_length(stat);
				return;
			}
		}
	}

	//Case (4) and Case (5)
	if (m < 0)
	{
		y_coord = stat.ER_plane[ER_cur_index_x][ER_cur_index_y].y_min;
		x_coord = (y_coord - ls_dataset[cur_line_index].k) / ls_dataset[cur_line_index].m;

		if (x_coord >= stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_min && x_coord <= stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_max)
		{
			threshold = Euclidean_distance(stat.intersect_x, stat.intersect_y, ls_dataset[cur_line_index].x_end, ls_dataset[cur_line_index].y_end);
			ED = intersect_ED_determine(stat.intersect_x, stat.intersect_y, x_coord, y_coord, threshold, stat);

			if (stat.intersect_x < -inf && stat.intersect_y < -inf)
			{
				stat.ER_plane[ER_cur_index_x][ER_cur_index_y].total_length += threshold;
				return;
			}
			else
			{
				stat.ER_plane[ER_cur_index_x][ER_cur_index_y].total_length += ED;
				stat.ER_cur_index_y--;

				//Case (4)
				if (x_coord == stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_max)
					stat.ER_cur_index_x++;

				obtain_line_accumulated_length(stat);
				return;
			}
		}
	}

	//Case (3)
	x_coord = stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_max;
	y_coord = m * stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_max + k;

	if (y_coord >= stat.ER_plane[ER_cur_index_x][ER_cur_index_y].y_min && y_coord <= stat.ER_plane[ER_cur_index_x][ER_cur_index_y].y_max)
	{
		threshold = Euclidean_distance(stat.intersect_x, stat.intersect_y, ls_dataset[cur_line_index].x_end, ls_dataset[cur_line_index].y_end);
		ED = intersect_ED_determine(stat.intersect_x, stat.intersect_y, x_coord, y_coord, threshold, stat);

		if (stat.intersect_x < -inf && stat.intersect_y < -inf)
		{
			stat.ER_plane[ER_cur_index_x][ER_cur_index_y].total_length += threshold;
			return;
		}
		else
		{
			stat.ER_plane[ER_cur_index_x][ER_cur_index_y].total_length += ED;
			stat.ER_cur_index_x++;

			obtain_line_accumulated_length(stat);
			return;
		}
	}
}

//Step 2
void obtain_accumulated_length(statistics& stat)
{
	for (int i = 0;i < stat.n;i++)
	{
		stat.cur_line_index = i;
		find_initial_intersection_point(stat);

		if (stat.intersect_x < -inf && stat.intersect_y < -inf)
			continue;
		else
		{
			stat.ER_cur_index_x = (int)floor((stat.intersect_x - stat.ER_boundary_x_min) / stat.incr_x);
			stat.ER_cur_index_y = (int)floor((stat.intersect_y - stat.ER_boundary_y_min) / stat.incr_y);
			obtain_line_accumulated_length(stat);
		}		
	}
}

//Step 3
void construct_prefix_structure(statistics& stat)
{
	for (int x = 0;x < stat.ER_X;x++)
	{
		for (int y = 0;y < stat.ER_Y;y++)
		{
			if (x == 0 && y == 0)
				stat.ER_prefix_plane[x][y].total_length = stat.ER_plane[x][y].total_length;
			if (x != 0 && y == 0)
				stat.ER_prefix_plane[x][y].total_length = stat.ER_prefix_plane[x - 1][y].total_length + stat.ER_plane[x][y].total_length;
			if (x == 0 && y != 0)
				stat.ER_prefix_plane[x][y].total_length = stat.ER_prefix_plane[x][y - 1].total_length + stat.ER_plane[x][y].total_length;
			if (x != 0 && y != 0)
				stat.ER_prefix_plane[x][y].total_length = stat.ER_prefix_plane[x - 1][y].total_length + stat.ER_prefix_plane[x][y - 1].total_length
					- stat.ER_prefix_plane[x - 1][y - 1].total_length + stat.ER_plane[x][y].total_length;
		}
	}

	for (int y = 0;y < stat.ER_Y;y++)
	{
		for (int x = 0;x < stat.ER_X;x++)
		{
			if (x == 0)
				stat.ER_one_D_prefix_plane[x][y].total_length = stat.ER_plane[x][y].total_length;
			else
				stat.ER_one_D_prefix_plane[x][y].total_length = stat.ER_one_D_prefix_plane[x - 1][y].total_length + stat.ER_plane[x][y].total_length;
		}
	}
}

void construct_LARGE(statistics& stat)
{
	create_extended_region(stat);
	obtain_accumulated_length(stat);
	construct_prefix_structure(stat);
}

void obtain_rectangular_mask(statistics& stat)
{
	double min_length = 0;
	double max_length = 0;
	double new_min_length;
	double new_max_length;

	double diagonal_length = sqrt(stat.incr_x * stat.incr_x + stat.incr_y * stat.incr_y);
	stat.LB_expanded_index = -1;
	stat.UB_expanded_index = -1;

	min_length += min(stat.incr_x, stat.incr_y) / 2.0;
	max_length += diagonal_length / 2.0;
	if (max_length > stat.bandwidth) //LB = 0 in this case.
		return;
	stat.LB_expanded_index++; stat.UB_expanded_index++;

	while (min_length <= stat.bandwidth)
	{
		new_min_length = min_length + min(stat.incr_x, stat.incr_y);
		new_max_length = max_length + diagonal_length;

		stat.UB_expanded_index++;
		if (new_max_length <= stat.bandwidth)
			stat.LB_expanded_index++;
		if (new_min_length > stat.bandwidth)
			break;

		min_length = new_min_length;
		max_length = new_max_length;
	}
}

double bound_rectangle(Pixel& p, statistics& stat, bool is_LB)
{
	int alpha_L, beta_L;
	int alpha_U, beta_U;
	double bound_value = 0;

	if (is_LB == true)
	{
		alpha_L = max(p.ER_region_x_index - stat.LB_expanded_index, 0);
		alpha_U = min(p.ER_region_x_index + stat.LB_expanded_index, stat.ER_X - 1);
		beta_L = max(p.ER_region_y_index - stat.LB_expanded_index, 0);
		beta_U = min(p.ER_region_y_index + stat.LB_expanded_index, stat.ER_Y - 1);
	}
	else
	{
		alpha_L = max(p.ER_region_x_index - stat.UB_expanded_index, 0);
		alpha_U = min(p.ER_region_x_index + stat.UB_expanded_index, stat.ER_X - 1);
		beta_L = max(p.ER_region_y_index - stat.UB_expanded_index, 0);
		beta_U = min(p.ER_region_y_index + stat.UB_expanded_index, stat.ER_Y - 1);
	}

	if (alpha_L == 0 && beta_L == 0)
		bound_value = stat.ER_prefix_plane[alpha_U][beta_U].total_length;
	if (alpha_L != 0 && beta_L == 0)
		bound_value = stat.ER_prefix_plane[alpha_U][beta_U].total_length - stat.ER_prefix_plane[alpha_L - 1][beta_U].total_length;
	if (alpha_L == 0 && beta_L != 0)
		bound_value = stat.ER_prefix_plane[alpha_U][beta_U].total_length - stat.ER_prefix_plane[alpha_U][beta_L - 1].total_length;
	if (alpha_L != 0 && beta_L != 0)
		bound_value = stat.ER_prefix_plane[alpha_U][beta_U].total_length - stat.ER_prefix_plane[alpha_L - 1][beta_U].total_length
		- stat.ER_prefix_plane[alpha_U][beta_L - 1].total_length + stat.ER_prefix_plane[alpha_L - 1][beta_L - 1].total_length;

	return bound_value;
}

int check_grid_status(double center_x, double center_y, statistics& stat)
{
	double ED_dist_upper_left;
	double ED_dist_upper_right;
	double ED_dist_lower_left;
	double ED_dist_lower_right;
	double ED_dist_max;
	double ED_dist_min;

	ED_dist_upper_left = Euclidean_distance(center_x - stat.incr_x / 2.0, center_y + stat.incr_y / 2.0, 0, 0);
	ED_dist_upper_right = Euclidean_distance(center_x + stat.incr_x / 2.0, center_y + stat.incr_y / 2.0, 0, 0);
	ED_dist_lower_left = Euclidean_distance(center_x - stat.incr_x / 2.0, center_y - stat.incr_y / 2.0, 0, 0);
	ED_dist_lower_right= Euclidean_distance(center_x + stat.incr_x / 2.0, center_y - stat.incr_y / 2.0, 0, 0);

	ED_dist_min = min(min(ED_dist_upper_left, ED_dist_upper_right), min(ED_dist_lower_left, ED_dist_lower_right));
	ED_dist_max = max(max(ED_dist_upper_left, ED_dist_upper_right), max(ED_dist_lower_left, ED_dist_lower_right));
	
	if (ED_dist_min > stat.bandwidth)
		return 2; //This grid is not covered by the bandwidth region.
	if (ED_dist_min <= stat.bandwidth && ED_dist_max > stat.bandwidth)
		return 1; //This grid is partly covered by the bandwidth region.
	if (ED_dist_max <= stat.bandwidth)
		return 0; //This grid is fully covered by the bandwidth region.
}

void obtain_arbitrary_mask(statistics& stat)
{
	int x_expanded_index = 0;
	int y_expanded_index = 0;
	int LB_x_expanded_index;
	int UB_x_expanded_index;
	int status;
	bool is_first_start = true;
	bool is_LB = false;

	status = check_grid_status(0, 0, stat);
	if (status != 0)
		return;

	while (true)
	{
		status = check_grid_status(x_expanded_index * stat.incr_x, y_expanded_index * stat.incr_y, stat);
		
		if (status == 0) //This grid is fully covered by the bandwidth region.
		{
			x_expanded_index++;
			LB_x_expanded_index = x_expanded_index;
			UB_x_expanded_index = x_expanded_index;
			is_first_start = false;
			is_LB = true;
		}
			
		if (status == 1) //This grid is partly covered by the bandwidth region.
		{
			x_expanded_index++;
			UB_x_expanded_index = x_expanded_index;
			is_first_start = false;
		}

		if (status == 2) //This grid is not covered by the bandwidth region.
		{
			if (is_first_start == true)
				break;

			if (is_LB == true)
				stat.LB_arbit_expanded_index_vec.push_back(LB_x_expanded_index - 1);
			stat.UB_arbit_expanded_index_vec.push_back(UB_x_expanded_index - 1);
			y_expanded_index++;
			x_expanded_index = 0;
			is_first_start = true;
			is_LB = false;
		}
	}
}

double bound_arbit(Pixel& p, statistics& stat, vector<int>& expanded_index_vec)
{
	int alpha_L, alpha_U;
	double bound_value = 0;

	for (int l = 0;l < expanded_index_vec.size();l++)
	{
		if (l == 0)
		{
			alpha_L = max(p.ER_region_x_index - expanded_index_vec[0], 0);
			alpha_U = min(p.ER_region_x_index + expanded_index_vec[0], stat.ER_X - 1);

			if (alpha_L == 0)
				bound_value += stat.ER_one_D_prefix_plane[alpha_U][p.ER_region_y_index].total_length;
			else
				bound_value += stat.ER_one_D_prefix_plane[alpha_U][p.ER_region_y_index].total_length
				- stat.ER_one_D_prefix_plane[alpha_L - 1][p.ER_region_y_index].total_length;

			continue;
		}

		if (p.ER_region_y_index - l >= 0)
		{
			alpha_L = max(p.ER_region_x_index - expanded_index_vec[l], 0);
			alpha_U = min(p.ER_region_x_index + expanded_index_vec[l], stat.ER_X - 1);

			if (alpha_L == 0)
				bound_value += stat.ER_one_D_prefix_plane[alpha_U][p.ER_region_y_index - l].total_length;
			else
				bound_value += stat.ER_one_D_prefix_plane[alpha_U][p.ER_region_y_index - l].total_length
				- stat.ER_one_D_prefix_plane[alpha_L - 1][p.ER_region_y_index - l].total_length;
		}

		if (p.ER_region_y_index + l < stat.ER_Y)
		{
			alpha_L = max(p.ER_region_x_index - expanded_index_vec[l], 0);
			alpha_U = min(p.ER_region_x_index + expanded_index_vec[l], stat.ER_X - 1);

			if (alpha_L == 0)
				bound_value += stat.ER_one_D_prefix_plane[alpha_U][p.ER_region_y_index + l].total_length;
			else
				bound_value += stat.ER_one_D_prefix_plane[alpha_U][p.ER_region_y_index + l].total_length
				- stat.ER_one_D_prefix_plane[alpha_L - 1][p.ER_region_y_index + l].total_length;
		}
	}

	return bound_value;
}

bool check_condition(statistics& stat, Pixel& p, double LB, double UB)
{
	if (stat.query_type == 0)
	{
		if (LB == 0)
			return false;

		if (UB - LB <= stat.epsilon * (UB + LB))
		{
			p.density_value = (2 * LB * UB) / (LB + UB);
			return true;
		}
		return false;
	}
}

/*void filter_and_refinement(statistics& stat)
{
	double LB_rect, UB_rect;
	double LB_arbit, UB_arbit;
	bool filter_condition;

	obtain_rectangular_mask(stat);
	obtain_arbitrary_mask(stat);
	if (stat.LB_expanded_index == -1) //LB must be 0
	{
		SCAN(stat);
		return;
	}

	for (int x = 0;x < stat.X;x++)
	{
		for (int y = 0;y < stat.Y;y++)
		{
			Pixel& p = stat.plane[x][y];

			//Filter
			LB_rect = bound_rectangle(p, stat, true);
			UB_rect = bound_rectangle(p, stat, false);

			filter_condition = check_condition(stat, p, LB_rect, UB_rect);
			if (filter_condition == true)
				continue;

			LB_arbit = bound_arbit(p, stat, stat.LB_arbit_expanded_index_vec);
			UB_arbit = bound_arbit(p, stat, stat.UB_arbit_expanded_index_vec);

			filter_condition = check_condition(stat, p, LB_arbit, UB_arbit);
			if (filter_condition == true)
				continue;

			//Refinement
			SCAN_one_pixel(stat, x, y);
		}
	}
}*/

void filter_and_refinement(statistics& stat, R_tree& R_tree)
{
	double LB_rect, UB_rect;
	double LB_arbit, UB_arbit;
	bool filter_condition;

	obtain_rectangular_mask(stat);
	obtain_arbitrary_mask(stat);
	if (stat.LB_expanded_index == -1) //LB must be 0
	{
		SCAN(stat);
		return;
	}

	for (int x = 0;x < stat.X;x++)
	{
		for (int y = 0;y < stat.Y;y++)
		{
			Pixel& p = stat.plane[x][y];

			//Filter
			LB_rect = bound_rectangle(p, stat, true);
			UB_rect = bound_rectangle(p, stat, false);

			filter_condition = check_condition(stat, p, LB_rect, UB_rect);
			if (filter_condition == true)
				continue;

			LB_arbit = bound_arbit(p, stat, stat.LB_arbit_expanded_index_vec);
			UB_arbit = bound_arbit(p, stat, stat.UB_arbit_expanded_index_vec);

			filter_condition = check_condition(stat, p, LB_arbit, UB_arbit);
			if (filter_condition == true)
				continue;

			if (stat.method == 1) //LARGE
				SCAN_one_pixel(stat, x, y);

			if (stat.method == 4) //LARGE + R-tree
				p.density_value = R_tree.compute_LDF(p, (R_node*)R_tree.rootNode);
		}
	}
}

//Used for debugging
void output_LARGE(statistics& stat)
{
	//ER_plane
	cout << "The plane of accumulated lengths in the extended region: " << endl;
	for (int x = 0;x < stat.ER_X;x++)
	{
		for (int y = 0;y < stat.ER_Y;y++)
			cout << stat.ER_plane[x][y].total_length << " ";
		cout << endl;
	}
	cout << endl;

	//ER_plane
	cout << "The plane of prefix accumulated lengths in the extended region: " << endl;
	for (int x = 0;x < stat.ER_X;x++)
	{
		for (int y = 0;y < stat.ER_Y;y++)
			cout << stat.ER_prefix_plane[x][y].total_length << " ";
		cout << endl;
	}
}