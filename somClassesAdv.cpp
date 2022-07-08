
#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include "operators.h"

#include <iostream>
#include <cstdlib>     /* srand, rand */
#include <cmath>
#include <vector>
#include <string>
#include "iterator"
#include "omp.h"

class somNodeAdv {
	private:
		int N_THREADS_I;
		int nodeId_i;
		int posX_i, posY_i;
		double dist_d;
		bool isNeigh_b;
	
	public:
		std::vector<double> weights_vd;
		
		template <typename T>
			void operator =(const std::vector<T>& v2)
			{
				for(int i = 0; i < getWeights().size(); i++) {
					getWeights()[i] = v2[i];
				}
			}

		somNodeAdv()=default;
		
		somNodeAdv(int N_THREADS, int nodeId, int nodeX, int nodeY, int weightsSize) :
				N_THREADS_I(N_THREADS), nodeId_i(nodeId), posX_i(nodeX), posY_i(nodeY), dist_d(0),
				isNeigh_b(false) {
			initWeights(weightsSize);
		};
		
		void initWeights(int weightsSize) {
			float randNum;
			for(int i = 0; i < weightsSize; i++) {
				randNum = (float) FLOAT_MIN + (float)(rand()) / ((float)(RAND_MAX/(FLOAT_MAX - FLOAT_MIN)));
				weights_vd.push_back(randNum);
			}
		};
	
		double euclDistance(std::vector<double> inputWeights_vd) {
			auto distance = (double) 0;
			int i;
			int n_per_thread = inputWeights_vd.size()/N_THREADS_I;
			omp_set_num_threads(N_THREADS_I);
			#pragma omp parallel for reduction(+ : distance) private(i) schedule(static, n_per_thread)
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
		
		int getId() const {
			return nodeId_i;
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
	
		~somNodeAdv()= default;
};


class inputNodesAdv {
private:
	std::vector<somNodeAdv> nodes_vs;
	int N_THREADS_I;
public:
	
	explicit inputNodesAdv(int N_THREADS, int nInput, int WEIGHT_SIZE, bool fromFile = false){
		
		N_THREADS_I = N_THREADS;
		std::vector<double> content;
		std::vector<double> row;
		std::string line, word;
		auto rowMax = (double) 0;
		int id = 0;
		if(fromFile) {
			std::string fname = "C:\\Users\\barba\\CLionProjects\\SOM_OMP\\iris.csv";
			std::fstream file (fname, std::ios::in);
			if(file.is_open())
			{
				while(getline(file, line) && id < nInput)
				{
					row.clear();
					std::stringstream str(line);
					while(getline(str, word, ',')){
						if(std::stod(word) > rowMax) rowMax = std::stod(word);
						row.push_back(std::stod(word));
					}
					content = (row/rowMax);
					
					somNodeAdv n = somNodeAdv(N_THREADS_I,id,0,0,WEIGHT_SIZE);
					n.weights_vd = content;
					addNode(n);
					id++;
				}
			}
			else
				std::cout << "Could not open the file" << std::endl;
		}
		else {
			for(int i = 0; i < nInput; i++) {
				somNodeAdv n = somNodeAdv(N_THREADS_I,
				                    id,
				                    0,
				                    0,
				                    WEIGHT_SIZE);
				addNode(n);
				id++;
			}
		}
	};
	
	void addNode(somNodeAdv &node) {
		nodes_vs.push_back(node);
	};
	
	std::vector<somNodeAdv> getNodes(bool verbose = false) {
		if(verbose) {
			for(auto & nodes_v : nodes_vs) {
				std::cout << "INPUT NODE ("
				          << nodes_v.getId()
				          << ") HAS WEIGHTS -> "
				          << nodes_v.getWeights()
				          << std::endl;
			}
		}
		
		return nodes_vs;
	};
	
	~inputNodesAdv()= default;
};

class somGridAdv {
	private:
	
		int N_THREADS_I;
		int somGRID_ROWS, somGRID_COLS;
		std::vector<somNodeAdv> nodes_vs;
		somTimer t_neighbExplorer, t_inputVsGrid;
		int winNodeX_i, winNodeY_i;
		double winNodeDist_d;
		double somInitRadius;
		double neighbRadius = (double) 1; // sigma0
		double learnRate;
	public:
	
		somGridAdv()= default;
		
		somGridAdv(int GRID_ROWS, int GRID_COLS, int N_THREADS, int WEIGHT_SIZE) {
			N_THREADS_I = N_THREADS;
			t_inputVsGrid = somTimer("inputVsGrid");
			t_neighbExplorer = somTimer("neighbExplorer");
			somGRID_ROWS = GRID_ROWS;
			somGRID_COLS = GRID_COLS;
			int id = 0;
			for(int r = 0; r < GRID_ROWS; r++) {
				for(int c = 0; c < GRID_COLS; c++) {
					somNodeAdv n = somNodeAdv(N_THREADS_I,id,r,c, WEIGHT_SIZE);
					nodes_vs.push_back(n);
					id++;
				}
			}
		};
		
		std::vector<somNodeAdv> getNodes(bool verbose = false) {
			if(verbose) {
				for(auto & nodes_v : nodes_vs) {
					std::cout << "NODE ("
					          << nodes_v.getPosX() << ","
					          << nodes_v.getPosY()
					          << ") ID " << nodes_v.getId()
					          <<" HAS WEIGHTS -> "
					          << nodes_v.getWeights()
					          << std::endl;
				}
			}
			return nodes_vs;
		};
		
		void inputVsGrid(somNodeAdv inputNode) {
			t_inputVsGrid.tic();
			omp_set_num_threads(N_THREADS_I);
			int i;
			double dist;
			int n_per_thread = nodes_vs.size()/N_THREADS_I;
			#pragma omp parallel for private(i, dist) schedule(static, n_per_thread)
			for(int i = 0; i < nodes_vs.size(); i++){
				dist = nodes_vs[i].euclDistance(inputNode.getWeights());
				nodes_vs[i].setDist(dist);
			}
			t_inputVsGrid.toc();
			BMU();
		};
		
		void BMU() {
			auto minDist = (double) 10000;
			for(int i = 0; i < nodes_vs.size(); i++) {
				if(nodes_vs[i].getDist() < minDist) {
//					winNodeId_i = nodes_vs[i].getId();
					winNodeX_i = nodes_vs[i].getPosX();
					winNodeY_i = nodes_vs[i].getPosY();
					winNodeDist_d = nodes_vs[i].getDist();
					minDist = nodes_vs[i].getDist();
				}
			}
		};
		
		double bmuEuclideanDist(somNodeAdv& node) {
			double dist = pow((double) node.getPosX() - (double) winNodeX_i, 2) +
			              pow((double) node.getPosY() - (double) winNodeY_i, 2);
			return dist;
		};
		
		void neighbExplorer(double radius) {
			double D;
			int i;
			int n_per_thread = nodes_vs.size()/N_THREADS_I;
			omp_set_num_threads(N_THREADS_I);
			t_neighbExplorer.tic();
			#pragma omp parallel for private(i, D) schedule(static, n_per_thread)
			for(i = 0; i < nodes_vs.size(); i++) {
				D = bmuEuclideanDist(nodes_vs[i]);
				if(D < radius)
					nodes_vs[i].setIsNb(true);
			}
			t_neighbExplorer.toc();
		};
		
		void adjustWeights(somNodeAdv inputnode, double lr) {
			int i;
			double adjustRate;
			int n_per_thread = nodes_vs.size()/N_THREADS_I;
			omp_set_num_threads(N_THREADS_I);
			#pragma omp parallel for private(i, winNodeDist_d) schedule(static, n_per_thread)
			for(i = 0; i < nodes_vs.size(); i++) {
				if(nodes_vs[i].getIsNb()) {
					// Calculate by how much its weights are adjusted
					adjustRate = exp(-(winNodeDist_d) / (2*(neighbRadius*neighbRadius)));
					auto miao = nodes_vs[i].getWeights() + lr * adjustRate * (inputnode.getWeights() - nodes_vs[i].getWeights());
					nodes_vs[i].weights_vd = miao;
				}
			}
		};
		
		void resetGrid() {
			int i;
			int n_per_thread = nodes_vs.size()/N_THREADS_I;
			omp_set_num_threads(N_THREADS_I);
			#pragma omp parallel for private(i) schedule(static, n_per_thread)
			for(i = 0; i < nodes_vs.size(); i++){
				nodes_vs[i].setIsNb(false);
			}
		};
		
		void printGridStats() {
			t_inputVsGrid.printDeltaT();
			t_neighbExplorer.printDeltaT();
		};
		
		bool somTrain(inputNodesAdv inputs, int epochs) {
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
				
				for(auto & inputnode : inputs.getNodes(false)){
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
		};
		
		~somGridAdv()= default;
		
};

#pragma clang diagnostic pop