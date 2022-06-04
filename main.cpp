#include <omp.h>
#include <iostream>
#include <string>
#include "somClasses.cpp"

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

constexpr int ROWS = 3;
constexpr int COLS = 3;
constexpr int NTHREADS = 4;
// number of weights each node must contain. One for each element of
// the input std::vector. In this example it is 3 because a color is
// represented by its red, green and blue components. (RGB)
constexpr int INPUTDIM = 3;
//the number of epochs desired for the training
//constexpr int EPOCHS = 10;
////the value of the learning rate at the start of training
//constexpr double START_LR = 0.1;

class somTimer {
	private:
		double itime_d;
		double ftime_d;
		double exectime_d;
	public:
		somTimer(): itime_d((double) 0), ftime_d((double) 0), exectime_d((double) 0) {};
		
		double tic() {
			this->itime_d = omp_get_wtime();
			return itime_d;
		};
		
		double toc() {
			this->ftime_d = omp_get_wtime();
			return this->ftime_d;
		};
		
		double execTime() {
			this->exectime_d = this->ftime_d - this->itime_d;
			printf("\n\nTime taken is %f", this->exectime_d);
			return this->exectime_d;
		};
};


#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
int main() {
	
	std::srand((unsigned int)time(NULL));
	// Initializing timer
	somTimer timer = somTimer();
	inputNodes nodesInput = inputNodes(INPUTDIM);
	somGrid nodesGrid = somGrid(ROWS, COLS);
	
	timer.tic();
	omp_set_num_threads(NTHREADS);
	nodesGrid.somTrain(nodesInput);
	timer.toc();
	timer.execTime();
	
	return 0;
}
#pragma clang diagnostic pop
