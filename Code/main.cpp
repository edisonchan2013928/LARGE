#include "Visual.h"

int main(int argc, char**argv)
{
	statistics stat;
	init(stat, argc, argv);
	visual_algorithm(stat);
	output_visual(stat);
}