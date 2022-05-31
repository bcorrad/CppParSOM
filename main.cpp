#include <iostream>
#include <cstdlib>     /* srand, rand */
#include <ctime>
#include "vector"
#include "omp.h"

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

constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 1;

class somNode{
private:
	std::vector<double> weights_d;
	double posX, posY;
	int edgeL_i, edgeR_i, edgeT_i, edgeB_i; // edges of the node
public:
#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
	somNode(int left, int right, int top, int bottom, int weightsSize) : edgeL_i(left),
                                                                         edgeR_i(right),
                                                                         edgeT_i(top),
                                                                         edgeB_i(bottom)
	{
		std::srand(std::time(nullptr));
		// weights initialization
		#pragma omp parallel for num_threads(4)
			for (int i = 0; i < weightsSize; i++){
				float randNum = FLOAT_MIN + (float)(rand()) / ((float)(RAND_MAX/(FLOAT_MAX - FLOAT_MIN)));
				// std::cout << randNum << std::endl;
				weights_d.push_back(randNum);
			}
		//calculate the node's center
		posX = edgeL_i + (double)(edgeR_i - edgeL_i)/2;
		posY = edgeT_i  + (double)(edgeB_i - edgeT_i)/2;
	}
	#pragma clang diagnostic pop
};

int main() {
	somNode n = somNode(2, 3, 4, 5, 6);
	return 0;
}
