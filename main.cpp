
#include "globals.h"
#include "metrics.h"
#include "somClasses.cpp"
#include <iostream>
#include <fstream>
#include "iomanip"
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
	somGrid nodesGrid = somGrid(GRID_ROWS, GRID_COLS);
	
	std::cout << N_THREADS << " THREADS" << std::endl;
	omp_set_num_threads(N_THREADS);
	timer.tic();
	nodesGrid.somTrain(nodesInput);
	timer.toc();
	
	std::cout << "*************************************************" << std::endl;
	std::cout << "***********************END***********************" << std::endl;
	std::cout << "*************************************************" << std::endl;
	std::cout << "STATS FOR " << N_THREADS << " THREADS" << std::endl;
	std::cout << "WEIGHT_SIZE = " << WEIGHT_SIZE << std::endl;
	std::cout << "EPOCHS = " << EPOCHS << std::endl;
	std::cout << "GRID_GRID_ROWS = " << GRID_ROWS << std::endl;
	std::cout << "GRID_GRID_COLS = " << GRID_COLS << std::endl;
	std::cout << "N_INPUTS = " << N_INPUTS << std::endl;
	
	long double totExecTime = timer.getDeltaT();
	std::cout << "TOTAL EXECUTION TIME [microseconds] = " << std::setprecision(15) << totExecTime << std::endl;
	long double paralDeltaT = nodesGrid.getParalTime();
	std::cout << "PARALLEL EXECUTION TIME [microseconds] = " << nodesGrid.getParalTime() << std::endl;
	long double su = speedUp(nodesGrid.getParalTime(), timer.getDeltaT());
	std::cout << "SPEED UP = " << su << std::endl;
	long double effic = efficiency(nodesGrid.getParalTime(), timer.getDeltaT(), N_THREADS);
	std::cout << "EFFICIENCY = " << effic << std::endl;
	long double costo = cost(nodesGrid.getParalTime(), N_THREADS);
	std::cout << "COST = " << costo << std::endl;
	
	std::ofstream my_file;
	std::string reportPath = "C:\\Users\\barba\\CLionProjects\\PPF_SOM\\report.txt";
	my_file.open(reportPath, std::ios_base::app); // append instead of overwrite
	my_file << N_THREADS << ","
	        << WEIGHT_SIZE << ","
	        << EPOCHS << ","
	        << GRID_ROWS << ","
	        << GRID_COLS << ","
	        << N_INPUTS << ","
	        << std::setprecision(15) << totExecTime << ","
	        << std::setprecision(15) << paralDeltaT << ","
	        << std::setprecision(15) << su << ","
	        << effic << ","
	        << costo << ",\n";
	my_file.close();
	
	return 0;
}
#pragma clang diagnostic pop
