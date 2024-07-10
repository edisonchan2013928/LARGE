#include "SCAN_line.h"

void init_influence_pixels(statistics& stat)
{
	create_extended_region(stat);

	stat.is_selected_pixel_matrix = new bool* [stat.ER_X];
	for (int x = 0;x < stat.ER_X;x++)
		stat.is_selected_pixel_matrix[x] = new bool[stat.ER_Y];

	for (int x = 0;x < stat.ER_X;x++)
		for (int y = 0;y < stat.ER_Y;y++)
			stat.is_selected_pixel_matrix[x][y] = false;
}

void clear_influence_pixels(statistics& stat)
{
	for (int x = 0;x < stat.ER_X;x++)
		delete[] stat.is_selected_pixel_matrix[x];
	delete[] stat.is_selected_pixel_matrix;
}

void obtain_influence_grid(statistics& stat)
{
	stat.influence_x_grid = (int)ceil(stat.bandwidth / stat.incr_x);
	stat.influence_y_grid = (int)ceil(stat.bandwidth / stat.incr_y);
	pair<int, int> index_pair;

	for (int x_index = stat.ER_cur_index_x - stat.influence_x_grid;x_index < stat.ER_cur_index_x + stat.influence_x_grid;x_index++)
	{
		for (int y_index = stat.ER_cur_index_y - stat.influence_y_grid;y_index < stat.ER_cur_index_y + stat.influence_y_grid;y_index++)
		{
			if (x_index < 0 || y_index < 0 || stat.is_selected_pixel_matrix[x_index][y_index] == true)
				continue;

			index_pair.first = x_index;
			index_pair.second = y_index;
			stat.is_selected_pixel_matrix[x_index][y_index] = true;
			stat.influence_pixel_index_vec.push_back(index_pair);
		}
	}
}

void find_influence_pixels_recur(statistics& stat)
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

	obtain_influence_grid(stat);

	//Case (1) and Case (2)
	if (m > 0)
	{
		y_coord = stat.ER_plane[ER_cur_index_x][ER_cur_index_y].y_max;
		x_coord = (y_coord - ls_dataset[cur_line_index].k) / ls_dataset[cur_line_index].m;

		if (x_coord >= stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_min && x_coord <= stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_max)
		{
			threshold = Euclidean_distance(stat.intersect_x, stat.intersect_y, ls_dataset[cur_line_index].x_end, ls_dataset[cur_line_index].y_end);
			ED = intersect_ED_determine(stat.intersect_x, stat.intersect_y, x_coord, y_coord, threshold, stat);

			if (stat.intersect_x < -inf && stat.intersect_y < -inf)
				return;
			else
			{
				stat.ER_cur_index_y++;

				//Case (2)
				if (x_coord == stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_max)
					stat.ER_cur_index_x++;

				find_influence_pixels_recur(stat);
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
				return;
			else
			{
				stat.ER_cur_index_y--;

				//Case (4)
				if (x_coord == stat.ER_plane[ER_cur_index_x][ER_cur_index_y].x_max)
					stat.ER_cur_index_x++;

				find_influence_pixels_recur(stat);
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
			return;
		else
		{
			stat.ER_cur_index_x++;

			find_influence_pixels_recur(stat);
			return;
		}
	}
}

void find_influence_pixels(statistics& stat)
{
	find_initial_intersection_point(stat);
	if (stat.intersect_x < -inf && stat.intersect_y < -inf)
		return;
	else
	{
		stat.ER_cur_index_x = (int)floor((stat.intersect_x - stat.ER_boundary_x_min) / stat.incr_x);
		stat.ER_cur_index_y = (int)floor((stat.intersect_y - stat.ER_boundary_y_min) / stat.incr_y);

		find_influence_pixels_recur(stat);
	}
}

void SCAN_single_line(statistics& stat)
{
	int x_index, y_index;
	int ori_x_index, ori_y_index;
	line_segment& l = stat.ls_dataset[stat.cur_line_index];

	find_influence_pixels(stat);
	for (int p_index = 0;p_index < stat.influence_pixel_index_vec.size();p_index++)
	{
		x_index = stat.influence_pixel_index_vec[p_index].first;
		y_index = stat.influence_pixel_index_vec[p_index].second;
		
		ori_x_index = x_index - stat.os_appended_x_grid_num;
		ori_y_index = y_index - stat.os_appended_y_grid_num;
		Pixel& p = stat.plane[ori_x_index][ori_y_index];

		if ((ori_x_index >= 0 && ori_x_index < stat.X) && (ori_y_index >= 0 && ori_y_index < stat.Y))
			p.density_value += compute_length(p, l, stat.bandwidth);

		stat.is_selected_pixel_matrix[x_index][y_index] = false;
	}

	stat.influence_pixel_index_vec.clear();
}

void SCAN_line(statistics& stat)
{
	init_influence_pixels(stat);

	for (int l_index = 0;l_index < stat.n;l_index++)
	{
		stat.cur_line_index = l_index;
		SCAN_single_line(stat);
	}

	//clear_influence_pixels(stat);
}