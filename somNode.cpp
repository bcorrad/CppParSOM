#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include <iostream>
#include <cstdlib>     /* srand, rand */
#include <ctime>
#include "math.h"
#include "vector"

constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 1;

class somNode {
	private:
	
		std::vector<double> weights_vd;
		double posX, posY;
		int edgeL_i, edgeR_i, edgeT_i, edgeB_i; // edges of the node
		
	public:
	
		somNode(int left, int right, int top, int bottom, int weightsSize) : edgeL_i(left), edgeR_i(right), edgeT_i(top), edgeB_i(bottom)
		{
			std::srand(std::time(nullptr));
			// weights initialization
			#pragma omp parallel for
			for(int i = 0; i < weightsSize; i++) {
				float randNum = FLOAT_MIN + (float)(rand()) / ((float)(RAND_MAX/(FLOAT_MAX - FLOAT_MIN)));
				// std::cout << randNum << std::endl;
				weights_vd.push_back(randNum);
			}
			//calculate the node's center
			posX = edgeL_i + (double)(edgeR_i - edgeL_i)/2;
			posY = edgeT_i  + (double)(edgeB_i - edgeT_i)/2;
		};
		
		double euclDistance(const std::vector<double> &inputWeights_vd) {
			double distance = 0;
			#pragma omp parallel for
			for(int i = 0; i < weights_vd.size(); ++i) {
				distance += (inputWeights_vd[i] - weights_vd[i]) * (inputWeights_vd[i] - weights_vd[i]);
			}
			return sqrt(distance);
		};
};
#pragma clang diagnostic pop