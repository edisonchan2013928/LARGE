#include "PMR_QUAD_Tree.h"

pmr_quad_node::pmr_quad_node()
{
	boundary = new double* [2];
	for (int b = 0;b < 2;b++)
		boundary[b] = new double[2];
}

int pmr_quad_node::search(double* q, double bandwidth)
{
	double dist_a, dist_b, dist_c, dist_d;
	double min_dist = 0;
	double max_dist = 0;

	dist_a = Euclidean_distance(q[0], q[1], boundary[0][0], boundary[1][0]);
	dist_b = Euclidean_distance(q[0], q[1], boundary[0][0], boundary[1][1]);
	dist_c = Euclidean_distance(q[0], q[1], boundary[0][1], boundary[1][1]);
	dist_d = Euclidean_distance(q[0], q[1], boundary[0][1], boundary[1][0]);

	//Region V
	if ((q[0] >= boundary[0][0] && q[0] <= boundary[0][1]) && (q[1] >= boundary[1][0] && q[1] <= boundary[1][1]))
	{
		min_dist = min(min(dist_a, dist_b), min(dist_c, dist_d));
		max_dist = max(max(dist_a, dist_b), max(dist_c, dist_d));

		if (bandwidth > max_dist)
			return 0;
		else
			return 2;
	}

	//Region I
	if (q[0] < boundary[0][0] && q[1] > boundary[1][1])
	{
		min_dist = dist_b;
		max_dist = dist_d;
	}

	//Region II
	if (q[0] < boundary[0][0] && (q[1] >= boundary[1][0] && q[1] <= boundary[1][1]))
	{
		min_dist = boundary[0][0] - q[0];
		max_dist = max(dist_c, dist_d);
	}

	//Region III
	if (q[0] < boundary[0][0] && q[1] < boundary[1][0])
	{
		min_dist = dist_a;
		max_dist = dist_c;
	}

	//Region IV
	if ((q[0] >= boundary[0][0] && q[0] <= boundary[0][1]) && q[1] > boundary[1][1])
	{
		min_dist = q[1] - boundary[1][1];
		max_dist = max(dist_a, dist_d);
	}

	//Region VI
	if ((q[0] >= boundary[0][0] && q[0] <= boundary[0][1]) && q[1] < boundary[1][0])
	{
		min_dist = boundary[1][0] - q[1];
		max_dist = max(dist_b, dist_c);
	}

	//Region VII
	if (q[0] > boundary[0][1] && q[1] > boundary[1][1])
	{
		min_dist = dist_c;
		max_dist = dist_a;
	}

	//Region VIII
	if (q[0] > boundary[0][1] && (q[1] >= boundary[1][0] && q[1] <= boundary[1][1]))
	{
		min_dist = q[0] - boundary[0][1];
		max_dist = max(dist_a, dist_b);
	}

	//Region IX
	if (q[0] > boundary[0][1] && q[1] < boundary[1][0])
	{
		min_dist = dist_d;
		max_dist = dist_b;
	}

	if (max_dist < bandwidth) //The search region fully covers the grid region.
		return 0;
	if (min_dist > bandwidth) //The search region never covers the grid region.
		return 1;

	return 2; //The search region partially covers the grid region.

	/*double dist_a, dist_b, dist_c, dist_d;
	double min_dist, max_dist;

	dist_a = Euclidean_distance(q[0], q[1], boundary[0][0], boundary[1][0]);
	dist_b = Euclidean_distance(q[0], q[1], boundary[0][0], boundary[1][1]);
	dist_c = Euclidean_distance(q[0], q[1], boundary[0][1], boundary[1][1]);
	dist_d = Euclidean_distance(q[0], q[1], boundary[0][1], boundary[1][0]);

	min_dist = min(min(dist_a, dist_b), min(dist_c, dist_d));
	max_dist = max(max(dist_a, dist_b), max(dist_c, dist_d));

	if (max_dist <= bandwidth) //The search region fully covers the grid region.
		return 0;

	if (min_dist > bandwidth) //The search region never covers the grid region.
		return 1;

	return 2; //The search region partially covers the grid region.*/
}

bool intersect_ED_determine(double line_x_start, double line_y_start, double x_coord, double y_coord, double threshold)
{
	double ED = Euclidean_distance(line_x_start, line_y_start, x_coord, y_coord);
	if (ED > threshold) //Cannot intersect the boundary
		return false;
	else
		return true;
}

bool pmr_quad_node::check_line_intersection(line_segment* ls_dataset, int lid)
{
	double m = ls_dataset[lid].m;
	double k = ls_dataset[lid].k;
	double x_coord, y_coord;
	double x_min, y_min;
	double x_max, y_max;
	double length;

	x_min = boundary[0][0]; x_max = boundary[0][1];
	y_min = boundary[1][0]; y_max = boundary[1][1];
	length = Euclidean_distance(ls_dataset[lid].x_start, ls_dataset[lid].y_start, ls_dataset[lid].x_end, ls_dataset[lid].y_end);

	//Region I
	if (ls_dataset[lid].x_start < x_min && ls_dataset[lid].y_start > y_max)
	{
		x_coord = (y_max - k) / m;
		y_coord = y_max;

		if (x_coord >= x_min && x_coord <= x_max)
			return intersect_ED_determine(ls_dataset[lid].x_start, ls_dataset[lid].y_start, x_coord, y_coord, length);

		x_coord = x_min;
		y_coord = m * x_coord + k;
		if (y_coord >= y_min && y_coord <= y_max)
			return intersect_ED_determine(ls_dataset[lid].x_start, ls_dataset[lid].y_start, x_coord, y_coord, length);

		return false;
	}

	//Region II
	if (ls_dataset[lid].x_start < x_min &&
		(ls_dataset[lid].y_start >= y_min && ls_dataset[lid].y_start <= y_max))
	{
		x_coord = x_min;
		y_coord = m * x_coord + k;
		if (y_coord >= y_min && y_coord <= y_max)
			return intersect_ED_determine(ls_dataset[lid].x_start, ls_dataset[lid].y_start, x_coord, y_coord, length);

		return false;
	}

	//Region III
	if (ls_dataset[lid].x_start < x_min && ls_dataset[lid].y_start < y_min)
	{
		x_coord = (y_min - k) / m;
		y_coord = y_min;
		if (x_coord >= x_min && x_coord <= x_max)
			return intersect_ED_determine(ls_dataset[lid].x_start, ls_dataset[lid].y_start, x_coord, y_coord, length);

		x_coord = x_min;
		y_coord = m * x_coord + k;
		if (y_coord >= y_min && y_coord <= y_max)
			return intersect_ED_determine(ls_dataset[lid].x_start, ls_dataset[lid].y_start, x_coord, y_coord, length);

		return false;
	}

	//Region IV
	if ((ls_dataset[lid].x_start >= x_min && ls_dataset[lid].x_start <= x_max) &&
		ls_dataset[lid].y_start > y_max)
	{
		if (m > 0)
			return false;

		x_coord = (y_max - k) / m;
		y_coord = y_max;
		if (x_coord >= x_min && x_coord <= x_max)
			return intersect_ED_determine(ls_dataset[lid].x_start, ls_dataset[lid].y_start, x_coord, y_coord, length);

		return false;
	}

	//Region V
	if ((ls_dataset[lid].x_start >= x_min && ls_dataset[lid].x_start <= x_max) &&
		(ls_dataset[lid].y_start >= y_min && ls_dataset[lid].y_start <= y_max))
		return true;

	//Region VI
	if ((ls_dataset[lid].x_start >= x_min && ls_dataset[lid].x_start <= x_max) &&
		ls_dataset[lid].y_start < y_min)
	{
		if (m < 0)
			return false;

		x_coord = (y_min - k) / m;
		y_coord = y_min;
		if (x_coord >= x_min && x_coord <= x_max)
			return intersect_ED_determine(ls_dataset[lid].x_start, ls_dataset[lid].y_start, x_coord, y_coord, length);
		
		return false;
	}

	//Region VII
	return false;
}

pmr_quad_node* pmr_quad_node::createNode()
{
	return new pmr_quad_node();
}

pmr_quad_tree::pmr_quad_tree(line_segment* ls_dataset, int n, int ls_num_split)
{
	this->ls_dataset = ls_dataset;
	this->n = n;
	this->ls_num_split = ls_num_split;

	this->sel_id_vec = new bool[n];

	for (int i = 0;i < n;i++)
		this->sel_id_vec[i] = false;
}

void pmr_quad_tree::insert_one_line_segment(pmr_quad_node* node, int lid)
{
	if (node->childVector.size() != 0)
	{
		for (int c = 0;c < node->childVector.size();c++)
		{
			if (((pmr_quad_node*)node->childVector[c])->check_line_intersection(ls_dataset, lid) == true)
				insert_one_line_segment((pmr_quad_node*)node->childVector[c], lid);
		}
		return;
	}

	node->idList.push_back(lid);

	if (node->idList.size() <= ls_num_split)
		return;
	else //split this node into four child nodes
	{
		//(1) insert four child nodes and (2) remove the idList for this node.
		pmr_quad_node* node_left_low = new pmr_quad_node();
		pmr_quad_node* node_left_up = new pmr_quad_node();
		pmr_quad_node* node_right_low = new pmr_quad_node();
		pmr_quad_node* node_right_up = new pmr_quad_node();
		
		node_left_low->boundary[0][0] = node->boundary[0][0];
		node_left_low->boundary[0][1] = (node->boundary[0][0] + node->boundary[0][1]) / 2.0;
		node_left_low->boundary[1][0] = node->boundary[1][0];
		node_left_low->boundary[1][1] = (node->boundary[1][0] + node->boundary[1][1]) / 2.0;

		node_left_up->boundary[0][0] = node->boundary[0][0];
		node_left_up->boundary[0][1] = (node->boundary[0][0] + node->boundary[0][1]) / 2.0;
		node_left_up->boundary[1][0] = (node->boundary[1][0] + node->boundary[1][1]) / 2.0;
		node_left_up->boundary[1][1] = node->boundary[1][1];

		node_right_low->boundary[0][0] = (node->boundary[0][0] + node->boundary[0][1]) / 2.0;
		node_right_low->boundary[0][1] = node->boundary[0][1];
		node_right_low->boundary[1][0] = node->boundary[1][0];
		node_right_low->boundary[1][1] = (node->boundary[1][0] + node->boundary[1][1]) / 2.0;

		node_right_up->boundary[0][0] = (node->boundary[0][0] + node->boundary[0][1]) / 2.0;
		node_right_up->boundary[0][1] = node->boundary[0][1];
		node_right_up->boundary[1][0] = (node->boundary[1][0] + node->boundary[1][1]) / 2.0;
		node_right_up->boundary[1][1] = node->boundary[1][1];

		for (int id = 0;id < node->idList.size();id++)
		{
			if (node_left_low->check_line_intersection(ls_dataset, node->idList[id]) == true)
				node_left_low->idList.push_back(node->idList[id]);
			if (node_left_up->check_line_intersection(ls_dataset, node->idList[id]) == true)
				node_left_up->idList.push_back(node->idList[id]);
			if (node_right_low->check_line_intersection(ls_dataset, node->idList[id]) == true)
				node_right_low->idList.push_back(node->idList[id]);
			if (node_right_up->check_line_intersection(ls_dataset, node->idList[id]) == true)
				node_right_up->idList.push_back(node->idList[id]);
		}

		node->idList.clear();
		node->childVector.push_back(node_left_low);
		node->childVector.push_back(node_left_up);
		node->childVector.push_back(node_right_low);
		node->childVector.push_back(node_right_up);
	}
}

void pmr_quad_tree::construct_pmr_quad_tree()
{
	//Obtain the boundary of the root node
	rootNode = new pmr_quad_node();

	/*rootNode->boundary = new double* [2];
	for (int b = 0;b < 2;b++)
		rootNode->boundary[b] = new double[2];*/

	rootNode->boundary[0][0] = inf; rootNode->boundary[0][1] = -inf;
	rootNode->boundary[1][0] = inf; rootNode->boundary[1][1] = -inf;

	for (int lid = 0;lid < n;lid++)
	{
		if (ls_dataset[lid].x_start < rootNode->boundary[0][0])
			rootNode->boundary[0][0] = ls_dataset[lid].x_start;
		if (ls_dataset[lid].x_end > rootNode->boundary[0][1])
			rootNode->boundary[0][1] = ls_dataset[lid].x_end;
		if (min(ls_dataset[lid].y_start, ls_dataset[lid].y_end) < rootNode->boundary[1][0])
			rootNode->boundary[1][0] = min(ls_dataset[lid].y_start, ls_dataset[lid].y_end);
		if (max(ls_dataset[lid].y_start, ls_dataset[lid].y_end) > rootNode->boundary[1][1])
			rootNode->boundary[1][1] = max(ls_dataset[lid].y_start, ls_dataset[lid].y_end);
	}

	//construct the PMR Quad-tree structure
	for (int lid = 0;lid < n;lid++)
		insert_one_line_segment((pmr_quad_node*)rootNode, lid);
}

void pmr_quad_tree::obtain_idList(double* q, double bandwidth, pmr_quad_node* node)
{
	int intersect_property;

	intersect_property = node->search(q, bandwidth); //0: fully cover, 1: never cover, 2: partially cover

	/*if (intersect_property == 0) //fully cover
	{
		for (int id = 0;id < node->idList.size();id++)
			result_idList.push_back(node->idList[id]);
		return;
	}*/

	if (intersect_property == 1) //never cover
		return;

	if (intersect_property == 0 || intersect_property == 2) //fully cover or partially cover
	{
		if (node->childVector.size() != 0) //non-leaf node
		{
			for (int c = 0;c < node->childVector.size();c++)
				obtain_idList(q, bandwidth, (pmr_quad_node*)node->childVector[c]);

			return;
		}

		//leaf node
		for (int id = 0;id < node->idList.size();id++)
		{
			if (sel_id_vec[node->idList[id]] == false)
				result_idList.push_back(node->idList[id]);
			sel_id_vec[node->idList[id]] = true;
		}
	}
}

double pmr_quad_tree::compute_LDF(double* q, Pixel& p, double bandwidth)
{
	double LDF_value = 0;
	int lid;
	obtain_idList(q, bandwidth, (pmr_quad_node*)rootNode);

	//compute LDF_value
	for (int id = 0;id < result_idList.size();id++)
	{
		lid = result_idList[id];
		line_segment& l = ls_dataset[lid];
		LDF_value += compute_length(p, l, bandwidth);

		/*if (sel_id_vec[id] == true)
			continue;
		else
		{
			line_segment& l = ls_dataset[lid];
			p.density_value += compute_length(p, l, bandwidth);
			sel_id_vec[lid] = true;
		}*/
	}

	for (int id = 0;id < result_idList.size();id++)
		sel_id_vec[result_idList[id]] = false;
	result_idList.clear();

	return LDF_value;
}

void pmr_quad_tree::generate_LDV(statistics& stat)
{
	double* q = new double[2];

	for (int x_index = 0;x_index < stat.X;x_index++)
	{
		for (int y_index = 0;y_index < stat.Y;y_index++)
		{
			Pixel& p = stat.plane[x_index][y_index];
			q[0] = p.x_center; q[1] = p.y_center;
			p.density_value = compute_LDF(q, p, stat.bandwidth);
		}
	}
}