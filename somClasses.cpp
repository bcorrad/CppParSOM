
#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include "somTimer.cpp"
#include "operators.h"

#include <iostream>
#include <cstdlib>     /* srand, rand */
#include <cmath>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include "algorithm"
#include "iterator"
#include "omp.h"

constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 1;
constexpr double START_LR = 0.001;

class somNode {
	private:
	
		int N_THREADS_I;
		int posX_i, posY_i;
		double dist_d;
		bool isNeigh_b;
	public:
	
		std::vector<double> weights_vd;
		
		// operator overloading for vectors
		template <typename T>
			void operator =(const std::vector<T>& v2)
			{
				for(int i = 0; i < getWeights().size(); i++) {
					getWeights()[i] = v2[i];
				}
			}

		somNode()=default;
		
		somNode(int N_THREADS, int nodeX, int nodeY, int weightsSize) :
				N_THREADS_I(N_THREADS), posX_i(nodeX), posY_i(nodeY), dist_d(0), isNeigh_b(false) {
			initWeights(weightsSize);
		};
		
		void initWeights(int weightsSize) {
			/**
			 * Weights initialization function for a node.
			 * A vector of length weightSize with random numbers.
			 */
			float randNum;
			for(int i = 0; i < weightsSize; i++) {
				randNum = (float) FLOAT_MIN + (float)(rand()) / ((float)(RAND_MAX/(FLOAT_MAX - FLOAT_MIN)));
				weights_vd.push_back(randNum);
			}
		};
	
		double euclDistance(std::vector<double> inputWeights_vd) {
			/**
			 * Function which calculates Euclidean distance
			 * between the current node and input node weights.
			 */
			auto distance = (double) 0;
			int i;
			omp_set_num_threads(N_THREADS_I);
			#pragma omp parallel for shared(distance) private(i)
			for(i = 0; i < weights_vd.size(); ++i) {
				 distance += (inputWeights_vd[i] - weights_vd[i]) * (inputWeights_vd[i] - weights_vd[i]);
			}
			return sqrt(distance);
		};
		
		std::vector<double> getWeights() const {
			return weights_vd;
		};
		
		int getPosX() const {
			return posX_i;
		};
		
		int getPosY() const {
			return posY_i;
		};
		
		double getDist() const {
			return dist_d;
		};
		
		bool getIsNb() const {
			return isNeigh_b;
		};
		
		void setDist(double dist) {
			dist_d = dist;
		};
		
		void setIsNb(bool isN) {
			isNeigh_b = isN;
		};
	
		~somNode()= default;
};


class inputNodes {
	private:
	
		std::vector<somNode> nodes_vs;
		int N_THREADS_I;
	public:
	
		explicit inputNodes(int N_THREADS, int nInput, int WEIGHT_SIZE, bool fromFile = false){
			N_THREADS_I = N_THREADS;
			// variables for text file parsing
			std::vector<double> content;
			std::vector<double> txtLine;
			std::string line, word;
			auto txtLineMax = (double) 0;
			int nLine = 0;
			if(fromFile) {
				std::string fname = "C:\\Users\\barba\\CLionProjects\\SOM_OMP\\iris.csv";
				std::fstream file (fname, std::ios::in);
				if(file.is_open())
				{
					while(getline(file, line) && nLine < nInput)
					{
						txtLine.clear();
						std::stringstream str(line);
						while(getline(str, word, ',')){
							if(std::stod(word) > txtLineMax) txtLineMax = std::stod(word);
							txtLine.push_back(std::stod(word));
						}
						content = (txtLine/txtLineMax);
						
						somNode n = somNode(N_THREADS_I,0,0,WEIGHT_SIZE);
						n.weights_vd = content;
						addNode(n);
						nLine++;
					}
				}
				else
					std::cout << "Could not open the file" << std::endl;
			}
			else {
				for(int i = 0; i < nInput; i++) {
					somNode n = somNode(N_THREADS_I,0,0,WEIGHT_SIZE);
					addNode(n);
				}
			}
		};
		
		void addNode(somNode &node) {
			nodes_vs.push_back(node);
		};
		
		std::vector<somNode> getNodes() {
			return nodes_vs;
		};
	
		~inputNodes()= default;
};


class somGrid {
	private:
	
		int N_THREADS_I;
		int somGRID_ROWS, somGRID_COLS;     // Rows and columns of SOM grid
		std::vector<somNode> nodes_vs;
		somTimer t_neighbExplorer, t_inputVsGrid, t_adjustWeights, t_resetGrid;
		int winNodeX_i, winNodeY_i; // BMU coordinates
		double winNodeDist_d;       // BMU distance from input node
		double somInitRadius;   // Initial neighborhood radius
		double neighbRadius = (double) 1;   // Updated neighborhood radius calculated as
											// somInitRadius * exp(-(double)epoch/timeConst); (aka sigma0)
		double learnRate;       // Initial learning rate
	public:
	
		somGrid()= default;
		
		somGrid(int GRID_ROWS, int GRID_COLS, int N_THREADS, int WEIGHT_SIZE) {
			N_THREADS_I = N_THREADS;
			t_inputVsGrid.setFunc("inputVsGrid");
			t_neighbExplorer.setFunc("neighbExplorer");
			t_adjustWeights.setFunc("adjustWeights");
			t_resetGrid.setFunc("resetGrid");
			somGRID_ROWS = GRID_ROWS;
			somGRID_COLS = GRID_COLS;
			for(int r = 0; r < GRID_ROWS; r++) {
				for(int c = 0; c < GRID_COLS; c++) {
					somNode n = somNode(N_THREADS_I, r, c, WEIGHT_SIZE);
					nodes_vs.push_back(n);
				}
			}
		}
		
		void inputVsGrid(somNode inputNode) {
			/**
			 * Function which compares the input node with each
			 * node of the grid, calculates the euclidean distance
			 * and finds BMU.
			 */
			t_inputVsGrid.tic();
			omp_set_num_threads(N_THREADS_I);
			int i;
			#pragma omp parallel for private(i)
			for(int i = 0; i < nodes_vs.size(); i++){
				double dist = nodes_vs[i].euclDistance(inputNode.getWeights());
				nodes_vs[i].setDist(dist);
			}
			t_inputVsGrid.toc();
			BMU();
		}
		
		void BMU() {
			/**
			 * Function which finds the BMU as the minimum-distance node
			 * from the input node.
			 */
			auto minDist = (double) 100000;
			for(int i = 0; i < nodes_vs.size(); i++) {
				if(nodes_vs[i].getDist() < minDist) {
					winNodeX_i = nodes_vs[i].getPosX();
					winNodeY_i = nodes_vs[i].getPosY();
					winNodeDist_d = nodes_vs[i].getDist();
					minDist = nodes_vs[i].getDist();
				}
			}
		}
	
		double bmuEuclideanDist(somNode& node) {
			double dist = pow((double) node.getPosX() - (double) winNodeX_i, 2) +
					pow((double) node.getPosY() - (double) winNodeY_i, 2);
			return dist;
		}
		
		void neighbExplorer(double radius) {
			/**
			 * Function which finds the neighborhood of the BMU
			 * based on the current radius.
			 */
			double D;
			int i;
			omp_set_num_threads(N_THREADS_I);
			t_neighbExplorer.tic();
			#pragma omp parallel for private(i, D)
			for(i = 0; i < nodes_vs.size(); i++) {
				D = bmuEuclideanDist(nodes_vs[i]);
				if(D < radius)
					nodes_vs[i].setIsNb(true);
			}
			t_neighbExplorer.toc();
		}
		
		void adjustWeights(somNode inputnode, double lr) {
			int i;
			double adjustRate;
			omp_set_num_threads(N_THREADS_I);
			t_adjustWeights.tic();
			#pragma omp parallel for private(i, winNodeDist_d)
			for(i = 0; i < nodes_vs.size(); i++) {
				if(nodes_vs[i].getIsNb()) {
					// Calculate by how much its weights are adjusted
					adjustRate = exp(-(winNodeDist_d) / (2*(neighbRadius*neighbRadius)));
					auto miao = nodes_vs[i].getWeights() + lr * adjustRate * (inputnode.getWeights() - nodes_vs[i].getWeights());
					nodes_vs[i].weights_vd = miao;
				}
			}
			t_adjustWeights.toc();
		}
	
		void resetGrid() {
			/**
			 * Function which resets the neighbours of the grid.
			 */
			int i;
			t_resetGrid.tic();
			omp_set_num_threads(N_THREADS_I);
			#pragma omp parallel for private(i)
			for(i = 0; i < nodes_vs.size(); i++){
				nodes_vs[i].setIsNb(false);
			}
			t_resetGrid.toc();
		}
		
		void printGridStats() {
			t_inputVsGrid.printDeltaT();
			t_neighbExplorer.printDeltaT();
			t_adjustWeights.printDeltaT();
			t_resetGrid.printDeltaT();
		}
		
		bool somTrain(inputNodes inputs, int epochs) {
			/*
			1. create n x n map with random node std::vector values
			2. loop while s < StepsMax times
				compute what a "close" node means, based on s
				compute a learn rate, based on s
				pick a random data item
				determine the map node closest to data item (BMU)
				for-each node close to the BMU
				adjust node std::vector values towards data item
				end-loop
			 */
			
			omp_set_num_threads(N_THREADS_I);

			// Calculate the biggest possible initial radius (half of width either height (delta_0))
			somInitRadius = (double) std::max(somGRID_ROWS, somGRID_COLS) / 2;
			// Time constant (lambda) used in the calculation of the neighbourhood width
			double timeConst = epochs/log(somInitRadius);
			
			for(int epoch = 0; epoch < epochs; epoch++) {
				std::cout << "========= EPOCH " << epoch << "/" << epochs << " =========" << std::endl;
				
				// Calculate the width of the neighbourhood for this timestep
				neighbRadius = somInitRadius * exp(-(double)epoch/timeConst);
				// Reduce the learning rate
				learnRate = START_LR * exp(-(double)epoch/epochs);
				
				for(auto & inputnode : inputs.getNodes()){
					// Calculating the Best Matching Unit
					inputVsGrid(inputnode);
					// BMU neighborhood exploration
					neighbExplorer(neighbRadius);
					// Adjust weights
					adjustWeights(inputnode, learnRate);
					// Reset grid
					resetGrid();
				}
				printGridStats();
			}
			return true;
		}
		
		~somGrid()= default;
};

#pragma clang diagnostic pop