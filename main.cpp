
#include "globals.h"
#include "metrics.h"
#include "somClasses.cpp"
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

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

int main() {
	
	std::srand((unsigned int)time(NULL));
	// Initializing timer
	somTimer timer = somTimer();
	inputNodes nodesInput = inputNodes(N_INPUTS);
	somGrid nodesGrid = somGrid(ROWS, COLS);
//	omp_set_dynamic(0);

//	std::cout << "NUM_THREADS = " << omp_get_num_threads() << std::endl;
//	std::cout << "MAX_THREADS = " << omp_get_max_threads() << std::endl;
//	std::cout << "NODESGRID AT t=O" << std::endl;
//	nodesGrid.printGrid();
//	system("PAUSE");
//	for(int i = 1; i <= N_THREADS; i++) {
	int i = N_THREADS;
	std::cout << i << " THREADS" << std::endl;
	omp_set_dynamic(0);
	omp_set_num_threads(i);
	timer.tic();
	nodesGrid.somTrain(nodesInput);
	timer.toc();
	std::cout << "END " << i << " THREADS" << std::endl;
		std::cout << "==================================================" << std::endl;
//	}
//	std::cout << "NODESGRID AT t=100" << std::endl;
	timer.printDeltaT();
	long double su = speedUp(nodesGrid.getParalTime(), timer.getDeltaT(), N_THREADS);
	std::cout << "SPEED UP = " << su << std::endl;
	long double effic = efficiency(nodesGrid.getParalTime(), timer.getDeltaT(), N_THREADS);
	std::cout << "EFFICIENCY = " << effic << std::endl;
	long double costo = cost(nodesGrid.getParalTime(), N_THREADS);
	std::cout << "COST = " << costo << std::endl;
	return 0;
	
}
#pragma clang diagnostic pop
