#include <omp.h>
#include "somNode.cpp"

// Self Organizing Maps algorithm
// 1. Each node's weights are initialized.
// 2. A vector is chosen at random from the set of training data and presented to the lattice.
// 3. Every node is examined to calculate which one's weights are most like the input vector. The winning
//      node is commonly known as the Best Matching Unit (BMU).
// 4. The radius of the neighbourhood of the BMU is now calculated. This is a value that starts large,
//      typically set to the 'radius' of the lattice,  but diminishes each time-step. Any nodes found
//      within this radius are deemed to be inside the BMU's neighbourhood.
// 5. Each neighbouring node's (the nodes found in step 4) weights are adjusted to make them more like
//      the input vector. The closer a node is to the BMU, the more its weights get altered.
// Repeat step 2 for N iterations.

constexpr int ROWS = 100;
constexpr int COLS = 100;
constexpr int NTHREADS = 4;

class somTimer {
	private:
		double itime;
		double ftime;
		double exectime;
	public:
		somTimer() {};
		double tic() {
			this->itime = omp_get_wtime();
			return itime;
		};
		double toc() {
			this->ftime = omp_get_wtime();
			this->exectime = this->ftime - this->itime;
			printf("\n\nTime taken is %f", this->exectime);
			return this->exectime;
		};
};



#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
int main() {
	somTimer timer = somTimer();
	
	timer.tic();
	omp_set_num_threads(NTHREADS);
	// Initializing The Weights
	#pragma omp parallel for collapse(2)
	for(int row = 0; row < ROWS; row++){
		for(int col = 0; col < COLS; col++){
			// std::cout << omp_get_thread_num() << std::endl;
			somNode n = somNode(2, 2, 2, 2, 10);
		}
	}
	timer.toc();
	// Calculating the Best Matching Unit
	return 0;
	
}
#pragma clang diagnostic pop
