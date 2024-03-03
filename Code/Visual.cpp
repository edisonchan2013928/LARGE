#include "Visual.h"

void load_data(statistics& stat)
{
	double x_first, y_first;
	double x_second, y_second;
	fstream input_data_file;
	input_data_file.open(stat.input_data_fileName);

	if (input_data_file.is_open() == false)
	{
		cout << "Cannnot open input data file!" << endl;
		exit(0);
	}
	
	input_data_file >> stat.n;
	stat.ls_dataset = new line_segment[stat.n];
	for (int i = 0;i < stat.n;i++)
	{
		input_data_file >> x_first; input_data_file >> y_first;
		input_data_file >> x_second; input_data_file >> y_second;

		if (x_first < x_second)
		{
			stat.ls_dataset[i].x_start = x_first;
			stat.ls_dataset[i].y_start = y_first;
			stat.ls_dataset[i].x_end = x_second;
			stat.ls_dataset[i].y_end = y_second;
		}
		else
		{
			stat.ls_dataset[i].x_start = x_second;
			stat.ls_dataset[i].y_start = y_second;
			stat.ls_dataset[i].x_end = x_first;
			stat.ls_dataset[i].y_end = y_first;
		}

		if (fabs(stat.ls_dataset[i].x_end - stat.ls_dataset[i].x_start) < eps)
			stat.ls_dataset[i].m = inf;
		else
		{
			stat.ls_dataset[i].m = (stat.ls_dataset[i].y_end - stat.ls_dataset[i].y_start) /
				(stat.ls_dataset[i].x_end - stat.ls_dataset[i].x_start);
		}

		stat.ls_dataset[i].k = stat.ls_dataset[i].y_end - stat.ls_dataset[i].m * stat.ls_dataset[i].x_end;
	}

	input_data_file.close();
}

void init_pixel(statistics& stat)
{
	double cur_x;
	double cur_y;

	stat.plane = new Pixel*[stat.X];
	for (int x = 0;x < stat.X;x++)
		stat.plane[x] = new Pixel[stat.Y];

	stat.incr_x = (stat.boundary_x_max - stat.boundary_x_min) / stat.X;
	stat.incr_y = (stat.boundary_y_max - stat.boundary_y_min) / stat.Y;

	for (int x = 0;x < stat.X;x++)
	{
		cur_x = stat.boundary_x_min + x * stat.incr_x;
		for (int y = 0;y < stat.Y;y++)
		{
			cur_y = stat.boundary_y_min + y * stat.incr_y;
			
			stat.plane[x][y].x_min = cur_x;
			stat.plane[x][y].x_max = cur_x + stat.incr_x;
			stat.plane[x][y].y_min = cur_y;
			stat.plane[x][y].y_max = cur_y + stat.incr_y;
			stat.plane[x][y].x_center = cur_x + stat.incr_x / 2.0;
			stat.plane[x][y].y_center = cur_y + stat.incr_y / 2.0;
			stat.plane[x][y].density_value = 0;
		}
	}
}

void obtain_boundary(statistics& stat)
{
	stat.boundary_x_min = inf; stat.boundary_y_min = inf;
	stat.boundary_x_max = -inf; stat.boundary_y_max = -inf;

	for (int i = 0;i < stat.n;i++)
	{
		if (stat.ls_dataset[i].x_start < stat.boundary_x_min)
			stat.boundary_x_min = stat.ls_dataset[i].x_start;
		if (stat.ls_dataset[i].x_end > stat.boundary_x_max)
			stat.boundary_x_max = stat.ls_dataset[i].x_end;
		if (min(stat.ls_dataset[i].y_start, stat.ls_dataset[i].y_end) < stat.boundary_y_min)
			stat.boundary_y_min = min(stat.ls_dataset[i].y_start, stat.ls_dataset[i].y_end);
		if (max(stat.ls_dataset[i].y_start, stat.ls_dataset[i].y_end) > stat.boundary_y_max)
			stat.boundary_y_max = max(stat.ls_dataset[i].y_start, stat.ls_dataset[i].y_end);
	}
}

/*void init_tree(statistics& stat)
{
	stat.q = new double[2];
	stat.ls_set_select_vec = new bool[stat.n];
	for (int i = 0;i < stat.n;i++)
		stat.ls_set_select_vec[i] = false;
	stat.point_dataset = new Point[2 * stat.n];

	//Used for kd-tree
	for (int i = 0;i < stat.n;i++)
	{
		stat.point_dataset[2 * i].p = new double[2];
		stat.point_dataset[2 * i].line_index = i;
		stat.point_dataset[2 * i].p[0] = stat.ls_dataset[i].x_start;
		stat.point_dataset[2 * i].p[1] = stat.ls_dataset[i].y_start;
		stat.point_dataset[2 * i + 1].p = new double[2];
		stat.point_dataset[2 * i + 1].line_index = i;
		stat.point_dataset[2 * i + 1].p[0] = stat.ls_dataset[i].x_end;
		stat.point_dataset[2 * i + 1].p[1] = stat.ls_dataset[i].y_end;
	}
}*/

void init(statistics& stat, int argc, char** argv)
{
	stat.input_data_fileName = argv[1];
	stat.output_visual_fileName = argv[2];
	stat.query_type = atoi(argv[3]);
	stat.method = atoi(argv[4]);
	stat.X = atoi(argv[5]);
	stat.Y = atoi(argv[6]);
	stat.bandwidth = atof(argv[7]);
	if (stat.query_type == 0)
		stat.epsilon = atof(argv[8]);
	else
	{
		stat.num_tau = atoi(argv[8]);
		for (int t = 0;t < stat.num_tau;t++)
			stat.tau_vector.push_back(atof(argv[9 + t]));
	}

	//Debug
	/*stat.input_data_fileName = (char*)"../../../Datasets/San_Francisco_taxi/San_Francisco_taxi_S5";
	stat.output_visual_fileName = (char*)"./Results/San_Francisco_taxi_S5_m3_q0_X16_Y16_b1000_e0.2";
	stat.query_type = 0;
	stat.method = 3;
	stat.X = 16;
	stat.Y = 16;
	stat.bandwidth = 1000;
	stat.epsilon = 0.2;*/

	load_data(stat);
	obtain_boundary(stat);
	init_pixel(stat);

	//if (stat.method == 2)
	//	init_tree(stat);

}

void visual_algorithm(statistics& stat)
{
	double run_time;
	//kd_tree kdTree(stat.ls_dataset, stat.point_dataset, leafCapacity); //Remember to construct the point dataset
	R_tree R_tree(stat.n, stat.ls_dataset, leafCapacity, stat.bandwidth);
	pmr_quad_tree pmr_quadTree(stat.ls_dataset, stat.n, (int)sqrt(stat.n));

	auto start_s = chrono::high_resolution_clock::now();
	if (stat.method == 0) //SCAN
		SCAN(stat);
	if (stat.method == 1) //LARGE + SCAN (as a refinement method)
	{
		construct_LARGE(stat);
		filter_and_refinement(stat, R_tree);
	}
	if (stat.method == 2) //pmr_quad_tree
	{
		pmr_quadTree.construct_pmr_quad_tree();
		pmr_quadTree.generate_LDV(stat);
	}
	if (stat.method == 3) //R_tree
	{
		R_tree.construct_ls_rect_data();
		R_tree.build_RTree();
		R_tree.generate_LDV(stat);
	}
	if (stat.method == 4) //LARGE + R_tree (as a refinement method)
	{
		construct_LARGE(stat);
		R_tree.construct_ls_rect_data();
		R_tree.build_RTree();
		filter_and_refinement(stat, R_tree);
	}

	auto end_s = chrono::high_resolution_clock::now();

	run_time = (chrono::duration_cast<chrono::nanoseconds>(end_s - start_s).count()) / 1000000000.0;

	cout << "method " << stat.method << ":" << run_time << endl;
}

void output_visual(statistics& stat)
{
	fstream output_visual_file;
	output_visual_file.open(stat.output_visual_fileName, ios::in | ios::out | ios::trunc);
	if (output_visual_file.is_open() == false)
	{
		cout << "Cannot open output file!" << endl;
		exit(0);
	}

	for (int x = 0;x < stat.X;x++)
	{
		for (int y = 0;y < stat.Y;y++)
		{
			output_visual_file << stat.plane[x][y].x_center << " " << stat.plane[x][y].y_center
				<< " " << stat.plane[x][y].density_value << endl;
		}
	}
}