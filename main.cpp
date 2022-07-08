
#include "metrics.h"
#include "somClasses.cpp"
#include "somClassesAdv.cpp"
#include <iostream>
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

constexpr int MAX_THREADS = 4;
constexpr const int OMP_ADV = true;
constexpr const int EPOCHS = 10;
int N_THREADS;
int WEIGHT_SIZE;
int GRID_ROWS;
int GRID_COLS;
int N_INPUTS;
std::string LANG;

int main() {
	
	// setting parameters for experiments
	std::vector<int> weight_size{10, 100};
	std::vector<int> grid_rows{10, 100};
	std::vector<int> grid_cols{10, 100};
	std::vector<int> n_inputs{10, 100};
	
	for(auto & w:weight_size){
		for(auto & col : grid_cols){
			for(auto & row : grid_rows){
				for(auto & inp : n_inputs){
					std::vector<long double> T_;
					for(int i = 1; i <= MAX_THREADS; i++){
						WEIGHT_SIZE = w; // size of weight vector of each node
						GRID_ROWS = row;
						GRID_COLS = col;
						N_INPUTS = inp;
						N_THREADS = i;
						std::cout << N_THREADS << " THREADS" << std::endl;
						omp_set_num_threads(N_THREADS);
						std::srand((unsigned int) time(NULL));
						
						// initializing timer
						long double totExecTime;
						somTimer timer = somTimer();
						
						if(OMP_ADV == false){
							LANG = "cpp";
							
							inputNodes nodesInput = inputNodes(N_THREADS, N_INPUTS, WEIGHT_SIZE);
							somGrid nodesGrid = somGrid(row, col, N_THREADS, w);
							
							timer.tic();
							nodesGrid.somTrain(nodesInput, EPOCHS);
							timer.toc();
							
							totExecTime = timer.getDeltaT();
							T_.push_back(totExecTime);
						}
						else if(OMP_ADV == true){
							LANG = "cpp_adv";
							
							inputNodesAdv nodesInput = inputNodesAdv(N_THREADS, N_INPUTS, WEIGHT_SIZE);
							somGridAdv nodesGridAdv = somGridAdv(row, col, N_THREADS, w);
							
							timer.tic();
							nodesGridAdv.somTrain(nodesInput, EPOCHS);
							timer.toc();
							
							totExecTime = timer.getDeltaT();
							T_.push_back(totExecTime);
						}
						std::cout << "*************************************************" << std::endl;
						std::cout << "***********************END***********************" << std::endl;
						std::cout << "*************************************************" << std::endl;
						std::cout << "STATS FOR " << N_THREADS << " THREADS" << std::endl;
						std::cout << "WEIGHT_SIZE = " << WEIGHT_SIZE << std::endl;
						std::cout << "EPOCHS = " << EPOCHS << std::endl;
						std::cout << "GRID_GRID_ROWS = " << GRID_ROWS << std::endl;
						std::cout << "GRID_GRID_COLS = " << GRID_COLS << std::endl;
						std::cout << "N_INPUTS = " << N_INPUTS << std::endl;
						std::cout << "TOTAL EXECUTION TIME [microseconds] = " << std::setprecision(15) << totExecTime << std::endl;
					}
					report(T_, LANG, WEIGHT_SIZE, GRID_ROWS, GRID_COLS, N_INPUTS, EPOCHS);
				}
			}
		}
	}
	return 0;
}

#pragma clang diagnostic pop
