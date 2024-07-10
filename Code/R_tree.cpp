#include "R_tree.h"

R_node* R_node::createNode()
{
	return new R_node();
}

R_tree::R_tree(int n, line_segment* ls_dataset, int leafCapacity, double bandwidth)
{
	this->n = n;
	this->ls_dataset = ls_dataset;
	this->leafCapacity = leafCapacity;
	this->bandwidth = bandwidth;
}

void R_tree::construct_ls_rect_data()
{
	R_node* data_node;

	for (int l = 0;l < n;l++)
	{
		data_node = new R_node();
		data_node->boundary = new double* [2];
		data_node->boundary[0] = new double[2];	data_node->boundary[1] = new double[2];
		data_node->center = new double[2];

		data_node->boundary[0][0] = min(ls_dataset[l].x_start, ls_dataset[l].x_end);
		data_node->boundary[0][1] = max(ls_dataset[l].x_start, ls_dataset[l].x_end);
		data_node->boundary[1][0] = min(ls_dataset[l].y_start, ls_dataset[l].y_end);
		data_node->boundary[1][1] = max(ls_dataset[l].y_start, ls_dataset[l].y_end);

		data_node->center[0] = (data_node->boundary[0][0] + data_node->boundary[0][1]) / 2.0;
		data_node->center[1] = (data_node->boundary[1][0] + data_node->boundary[1][1]) / 2.0;
		//data_node->id = l;

		ls_rect_dataset.push_back(data_node);
	}
}

void R_tree::build_RTree()
{
	rootNode = new R_node();
	for (int i = 0; i < n; i++)
		rootNode->idList.push_back(i);

	update_R_node_info((R_node*)rootNode);
	R_Tree_Recur((R_node*)rootNode, 0);
}

double R_tree::refinement(Pixel& p, R_node* r_node)
{
	int id;
	double density_value = 0;
	for (int i = 0;i < r_node->idList.size();i++)
	{
		id = r_node->idList[i];
		line_segment& l = ls_dataset[id];
		density_value += compute_length(p, l, bandwidth);
	}

	return density_value;
}

double R_tree::compute_LDF(Pixel& p, R_node* r_node)
{
	double x_dist, y_dist;
	double min_dist;
	double density_value;

	if (r_node->childVector.size() == 0) //leaf node
		return refinement(p, r_node);

	if (p.x_center >= r_node->boundary[0][0] && p.x_center <= r_node->boundary[0][1])
		x_dist = 0;
	else
		x_dist = min(fabs(p.x_center - r_node->boundary[0][0]), fabs(p.x_center - r_node->boundary[0][1]));

	if (p.y_center >= r_node->boundary[1][0] && p.y_center <= r_node->boundary[1][1])
		y_dist = 0;
	else
		y_dist = min(fabs(p.y_center - r_node->boundary[1][0]), fabs(p.y_center - r_node->boundary[1][1]));

	min_dist = sqrt(x_dist * x_dist + y_dist * y_dist);

	if (min_dist > bandwidth)
		return 0;
	else
	{
		density_value = 0;
		for (int c = 0;c < r_node->childVector.size();c++)
			density_value += compute_LDF(p, (R_node*)r_node->childVector[c]);

		return density_value;
	}
}

void R_tree::generate_LDV(statistics& stat)
{
	for (int x_index = 0;x_index < stat.X;x_index++)
	{
		for (int y_index = 0;y_index < stat.Y;y_index++)
		{
			Pixel& p = stat.plane[x_index][y_index];
			p.density_value = compute_LDF(p, (R_node*)rootNode);
		}
	}
}

void R_tree::R_Tree_Recur(R_node* r_node, int split_Dim)
{
	int id;
	int counter;
	
	//base case
	if ((int)r_node->idList.size() <= leafCapacity)
		return;

	double splitValue = obtain_SplitValue(r_node, split_Dim);

	//create two children
	R_node* leftNode;
	R_node* rightNode;

	leftNode = r_node->createNode();
	rightNode = r_node->createNode();

	counter = 0;
	int halfSize = ((int)r_node->idList.size()) / 2;
	for (int i = 0; i < (int)r_node->idList.size(); i++)
	{
		id = r_node->idList[i];
		
		if (ls_rect_dataset[id]->center[split_Dim] <= splitValue && counter <= halfSize)
		{
			leftNode->idList.push_back(id);
			counter++;
		}
		else
			rightNode->idList.push_back(id);
	}

	update_R_node_info(leftNode);
	update_R_node_info(rightNode);

	R_Tree_Recur(leftNode, (split_Dim + 1) % 2);
	R_Tree_Recur(rightNode, (split_Dim + 1) % 2);

	r_node->childVector.push_back(leftNode);
	r_node->childVector.push_back(rightNode);
}

void R_tree::update_R_node_info(R_node* r_node)
{
	int id;
	r_node->boundary = new double* [2];
	r_node->boundary[0] = new double[2]; r_node->boundary[1] = new double[2];
	r_node->center = new double[2];

	r_node->boundary[0][0] = inf; r_node->boundary[0][1] = -inf;
	r_node->boundary[1][0] = inf; r_node->boundary[1][1] = -inf;

	for (int i = 0;i < r_node->idList.size();i++)
	{
		id = r_node->idList[i];

		if (ls_rect_dataset[id]->boundary[0][0] < r_node->boundary[0][0])
			r_node->boundary[0][0] = ls_rect_dataset[id]->boundary[0][0];
		if (ls_rect_dataset[id]->boundary[0][1] > r_node->boundary[0][1])
			r_node->boundary[0][1] = ls_rect_dataset[id]->boundary[0][1];
		if (ls_rect_dataset[id]->boundary[1][0] < r_node->boundary[1][0])
			r_node->boundary[1][0] = ls_rect_dataset[id]->boundary[1][0];
		if (ls_rect_dataset[id]->boundary[1][1] > r_node->boundary[1][1])
			r_node->boundary[1][1] = ls_rect_dataset[id]->boundary[1][1];
	}

	r_node->center[0] = (r_node->boundary[0][0] + r_node->boundary[0][1]) / 2.0;
	r_node->center[1] = (r_node->boundary[1][0] + r_node->boundary[1][1]) / 2.0;
}

double R_tree::obtain_SplitValue(R_node* r_node, int split_Dim)
{
	vector<double> tempVector;
	int id;
	int middle_left, middle_right, middle;
	
	for (int i = 0; i < (int)r_node->idList.size(); i++)
	{
		id = r_node->idList[i];
		tempVector.push_back(ls_rect_dataset[id]->center[split_Dim]);
	}

	sort(tempVector.begin(), tempVector.end());

	if ((int)tempVector.size() % 2 == 0)//even number
	{
		middle_right = (int)tempVector.size() / 2;
		middle_left = middle_right - 1;

		return ((tempVector[middle_left] + tempVector[middle_right]) / 2.0);
	}
	else
	{
		middle = ((int)tempVector.size() - 1) / 2;
		return tempVector[middle];
	}

	tempVector.clear();
}