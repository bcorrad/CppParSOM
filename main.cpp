
#include "globals.h"
#include "somClasses.cpp"

#include <omp.h>
#include <iostream>
#include <string>
#include "chrono"

// Self Organizing Maps algorithm
// 1. Each node's weights are initialized.
// 2. A std::vector is chosen at random from the set of training data and presented to the lattice.
// 3. Every node is examined to calculate which one's weights are most like the input std::vector. The winning
//      node is commonly known as the Best Matching Unit (BMU).
// 4. The radius of the neighbourhood of the BMU is now calculated. This is a value that starts large,
//      typically set to the 'radius' of the lattice,  but diminishes each time-step. Any nodes found
//      within this radius are deemed to be inside the BMU's neighbourhood.
// 5. Each neighbouring node's (the nodes found in step 4) weights are adjusted to make them more like
//      the input std::vector. The closer a node is to the BMU, the more its weights get altered.
// Repeat step 2 for N iterations.
//
//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "openmp-use-default-none"

int main() {
	
	std::srand((unsigned int)time(NULL));
	// Initializing timer
	somTimer timer = somTimer();
	inputNodes nodesInput = inputNodes(N_INPUTS);
	somGrid nodesGrid = somGrid(ROWS, COLS);
//	omp_set_dynamic(0);
	omp_set_num_threads(N_THREADS);
	std::cout << "NUM_THREADS = " << omp_get_num_threads() << std::endl;
	std::cout << "MAX_THREADS = " << omp_get_max_threads() << std::endl;
	std::cout << "NODESGRID AT t=O" << std::endl;
//	nodesGrid.printGrid();
//	system("PAUSE");
	timer.tic();
	nodesGrid.somTrain(nodesInput);
	timer.toc("", true);
	std::cout << "NODESGRID AT t=100" << std::endl;
//	nodesGrid.printGrid();
	
	return 0;
}
//#pragma clang diagnostic pop
