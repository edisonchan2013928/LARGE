g++ -c SCAN.cpp -o SCAN.o
g++ -c R_tree.cpp -w R_tree.o
g++ -c PMR_QUAD_Tree.cpp -w PMR_QUAD_Tree.o
g++ -c LARGE.cpp -w LARGE.o
g++ -c Visual.cpp -w Visual.o
g++ main.cpp -O3 -o main SCAN.o R_tree.o PMR_QUAD_Tree.o LARGE.o Visual.o
exit #Users can remove this line if they want to call our code.

#stat.input_data_fileName = argv[1]; //The name of the input dataset
#stat.output_visual_fileName = argv[2]; //The name of the output file
#stat.method = atoi(argv[3]); //The method name (method = 0: SCAN, method = 2: PMR quadtree, method = 3: R-tree, method = 4: LARGE) 
#stat.X = atoi(argv[4]); //The resolution size of the x-axis (e.g., 320)
#stat.Y = atoi(argv[5]); //The resolution size of the y-axis (e.g., 240)
#stat.bandwidth = atof(argv[6]); //The bandwidth parameter (in terms of meters)
#stat.epsilon = atof(argv[7]); //The relative error

#An example for calling our C++ code.
dataset="San_Francisco_taxi"
method=4
X=320
Y=240
bandwidth=1000
epsilon=0.1

./main $dataset "output_result_"$dataset $method $X $Y $bandwidth $epsilon