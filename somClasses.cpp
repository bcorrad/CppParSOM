
#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include "somTimer.cpp"

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

template <typename T>
	std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
	{
		for (int i = 0; i < v.size(); ++i) {
			os << v[i];
			if (i != v.size() - 1)
				os << ", ";
		}
		return os;
	}

template <typename T>
	std::vector<T> operator -(const std::vector<T>& v1, const std::vector<T>& v2)
	{
		int i;
		std::vector<T> result_v(v1.size());
		for(i = 0; i < v1.size(); i++) {
			result_v[i] = v1[i]-v2[i];
		}
		return result_v;
	}

template <typename T>
	std::vector<T> operator *(const T factor, const std::vector<T>& v)
	{
		int i;
		std::vector<T> result_v(v.size());
		for(i = 0; i < v.size(); i++) {
			result_v[i] = v[i]*factor;
		}
		return result_v;
	}


template <typename T>
	std::vector<T> operator /(const std::vector<T>& v, const T factor)
	{
		int i;
		std::vector<T> result_v(v.size());
		for(i = 0; i < v.size(); i++) {
			result_v[i] = v[i]/factor;
		}
		return result_v;
	}


template <typename T>
	std::vector<T> operator +(const std::vector<T>& v1, const std::vector<T>& v2)
	{
		std::vector<T> result_v;
		for(int i = 0; i < v1.size(); i++) {
			result_v.push_back(v1[i]+v2[i]);
		}
		return result_v;
	}


	
class somNode {
	private:
		
		int N_THREADS_I;
		int nodeId_i;
		int nodeCreatorId_i;
		
		int weightsSize_i;
		int posX_i, posY_i;
		double dist_d;
		bool isNeigh_b;
		bool verbose = false;
	
	public:
		std::vector<double> weights_vd;
		
		template <typename T>
			void operator =(const std::vector<T>& v2)
			{
				for(int i = 0; i < getWeights().size(); i++) {
					getWeights()[i] = v2[i];
				}
			}

		somNode()=default;
		
		somNode(int N_THREADS, int nodeId, int nodeX, int nodeY, int threadId, int weightsSize) :
				N_THREADS_I(N_THREADS), nodeId_i(nodeId), posX_i(nodeX), posY_i(nodeY), dist_d(0),
				nodeCreatorId_i(threadId), weightsSize_i(weightsSize), isNeigh_b(false) {
			if(verbose) std::cout << "Node " << nodeId << " initialized by thread " << threadId << std::endl;
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
			#pragma omp parallel for shared(distance) private(i) schedule(static, n_per_thread)
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
		
		void setId(int id) {
			nodeId_i = id;
		};
		
		void setPosX (int p) {
			posX_i = p;
		};
	
		void setPosY (int p) {
			posY_i = p;
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
					
					somNode n = somNode(N_THREADS_I,
										id,
					                    0,
					                    0,
										1,
					                    WEIGHT_SIZE);
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
				somNode n = somNode(N_THREADS_I,
									id,
				                    0,
				                    0,
				                    1,
				                    WEIGHT_SIZE);
				addNode(n);
				id++;
			}
		}
	};
	
	void addNode(somNode &node) {
		nodes_vs.push_back(node);
	};
	
	std::vector<somNode> getNodes(bool verbose = false) {
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
	
	~inputNodes()= default;
	
};


class somGrid {
	private:
		int N_THREADS_I;
		int somGRID_ROWS, somGRID_COLS;
		std::vector<somNode> nodes_vs;
		somTimer t_neighbExplorer, t_inputVsGrid, t_adjustWeights;
		long long paralTime;
		
		int winNodeId_i, winNodeX_i, winNodeY_i;
		double winNodeDist_d;
		
		double somInitRadius;
		double neighbRadius = (double) 1; // sigma0
		double learnRate;
	public:
		somGrid()= default;
		
		somGrid(int GRID_ROWS, int GRID_COLS, int N_THREADS, int WEIGHT_SIZE) {
			N_THREADS_I = N_THREADS;
			paralTime = (double) 0;
			t_inputVsGrid = somTimer("inputVsGrid");
			t_neighbExplorer = somTimer("neighbExplorer");
			this->somGRID_ROWS = GRID_ROWS;
			this->somGRID_COLS = GRID_COLS;
			int id = 0;
			for(int r = 0; r < GRID_ROWS; r++) {
				for(int c = 0; c < GRID_COLS; c++) {
					somNode n = somNode(N_THREADS_I,
										id,
										r,
										c,
										1, 
										WEIGHT_SIZE);
					nodes_vs.push_back(n);
					id++;
				}
			}
		};
		
		std::vector<somNode> getNodes(bool verbose = false) {
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
		
		void inputVsGrid(somNode inputNode, bool verbose = false) {
			t_inputVsGrid.tic();
			omp_set_num_threads(N_THREADS_I);
			int i;
			#pragma omp parallel for private(i)
			for(int i = 0; i < nodes_vs.size(); i++){
//			for(auto & nodes_v : nodes_vs) {
				double dist = nodes_vs[i].euclDistance(inputNode.getWeights());
				nodes_vs[i].setDist(dist);
				if(verbose) {
					std::cout << "INPUT NODE "
					          << inputNode.getId()
					          << " HAS DISTANCE -> "
					          << nodes_vs[i].getDist()
					          << " WRT NODE ("
					          << nodes_vs[i].getPosX() << ","
					          << nodes_vs[i].getPosY() << ")"
					          << std::endl;
				}
			}
			t_inputVsGrid.toc();
			BMU();
		};
		
		void BMU(bool verbose = false) {
			auto minDist = (double) 100000;
			for(int i = 0; i < nodes_vs.size(); i++) {
				if(nodes_vs[i].getDist() < minDist) {
					winNodeId_i = nodes_vs[i].getId();
					winNodeX_i = nodes_vs[i].getPosX();
					winNodeY_i = nodes_vs[i].getPosY();
					winNodeDist_d = nodes_vs[i].getDist();
					minDist = nodes_vs[i].getDist();
				}
			}
			if(verbose) {
				std::cout << "************************" << std::endl;
				std::cout << "MIN DIST IS "
				          << winNodeDist_d
				          << " AT NODE ("
				          << winNodeX_i
				          << "," << winNodeY_i
				          << ") ID "
				          << winNodeId_i
				          << std::endl;
			}
		};
		
		double getWinNodeDist() const {
			return winNodeDist_d;
		};
	
		double bmuEuclideanDist(somNode& node) {
			double dist = pow((double) node.getPosX() - (double) winNodeX_i, 2) +
					pow((double) node.getPosY() - (double) winNodeY_i, 2);
			return dist;
		};
		
		void neighbExplorer(double radius) {
			double D;
			int i;
//			int n_per_thread = nodes_vs.size()/N_THREADS_I;
			omp_set_num_threads(N_THREADS_I);
			t_neighbExplorer.tic();
			#pragma omp parallel for private(i, D)
			for(i = 0; i < nodes_vs.size(); i++) {
				D = bmuEuclideanDist(nodes_vs[i]);
				if(D < radius)
					nodes_vs[i].setIsNb(true);
			}
			t_neighbExplorer.toc();
		};
		
		void adjustWeights(somNode inputnode, double lr, bool verbose = false) {
			int i;
			double adjustRate;
//			int n_per_thread = nodes_vs.size()/N_THREADS_I;
			omp_set_num_threads(N_THREADS_I);
			t_adjustWeights.tic();
			#pragma omp parallel for private(i, winNodeDist_d)
			for(i = 0; i < nodes_vs.size(); i++) {
				if(nodes_vs[i].getIsNb()) {
					if(verbose){
						std::cout << "WEIGHTS OF NODE ID " << nodes_vs[i].getId() << " BEFORE ADJUSTMENT" << std::endl;
						std::cout << nodes_vs[i].getWeights() << std::endl;
					}
					// Calculate by how much its weights are adjusted
					adjustRate = exp(-(winNodeDist_d) / (2*(neighbRadius*neighbRadius)));
					auto miao = nodes_vs[i].getWeights() + lr * adjustRate * (inputnode.getWeights() - nodes_vs[i].getWeights());
					nodes_vs[i].weights_vd = miao;
					if(verbose){
						std::cout << "WEIGHTS AFTER ADJUSTMENT" << std::endl;
						std::cout << nodes_vs[i].getWeights() << std::endl;
						std::cout << "************************" << std::endl;
					}
				}
			}
			t_adjustWeights.toc();
		};
	
		void resetGrid() {
			int i;
//			int n_per_thread = nodes_vs.size()/N_THREADS_I;
			omp_set_num_threads(N_THREADS_I);
			#pragma omp parallel for private(i)
			for(i = 0; i < nodes_vs.size(); i++){
				nodes_vs[i].setIsNb(false);
			}
		};
		
		void printGridStats() {
			t_inputVsGrid.printDeltaT();
			t_neighbExplorer.printDeltaT();
			t_adjustWeights.printDeltaT();
		};
		
		long long getParalTime() {
			return paralTime;
		};
		
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
			
			learnRate = START_LR;
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
//				printGridStats();
//				paralTime += t_neighbExplorer.getDeltaT() + t_inputVsGrid.getDeltaT();
			}
			return true;
		};
		
		~somGrid()= default;
		
};


class somGridAdv {
private:
	int N_THREADS_I;
	int somGRID_ROWS, somGRID_COLS;
	std::vector<somNode> nodes_vs;
	somTimer t_neighbExplorer, t_inputVsGrid;
	long long paralTime;
	
	int winNodeId_i, winNodeX_i, winNodeY_i;
	double winNodeDist_d;
	
	double somInitRadius;
	double neighbRadius = (double) 1; // sigma0
	double learnRate;
public:
	somGridAdv()= default;
	
	somGridAdv(int GRID_ROWS, int GRID_COLS, int N_THREADS, int WEIGHT_SIZE) {
		N_THREADS_I = N_THREADS;
		paralTime = (double) 0;
		t_inputVsGrid = somTimer("inputVsGrid");
		t_neighbExplorer = somTimer("neighbExplorer");
		this->somGRID_ROWS = GRID_ROWS;
		this->somGRID_COLS = GRID_COLS;
		int id = 0;
		for(int r = 0; r < GRID_ROWS; r++) {
			for(int c = 0; c < GRID_COLS; c++) {
				somNode n = somNode(N_THREADS_I,
				                    id,
				                    r,
				                    c,
				                    1,
				                    WEIGHT_SIZE);
				nodes_vs.push_back(n);
				id++;
			}
		}
	};
	
	std::vector<somNode> getNodes(bool verbose = false) {
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
	
	void inputVsGrid(somNode inputNode, bool verbose = false) {
		t_inputVsGrid.tic();
		omp_set_num_threads(N_THREADS_I);
		int i;
		#pragma omp parallel for private(i)
		for(int i = 0; i < nodes_vs.size(); i++){
//			for(auto & nodes_v : nodes_vs) {
			double dist = nodes_vs[i].euclDistance(inputNode.getWeights());
			nodes_vs[i].setDist(dist);
			if(verbose) {
				std::cout << "INPUT NODE "
				          << inputNode.getId()
				          << " HAS DISTANCE -> "
				          << nodes_vs[i].getDist()
				          << " WRT NODE ("
				          << nodes_vs[i].getPosX() << ","
				          << nodes_vs[i].getPosY() << ")"
				          << std::endl;
			}
		}
		t_inputVsGrid.toc();
		BMU();
	};
	
	void BMU(bool verbose = false) {
		auto minDist = (double) 10000;
		for(int i = 0; i < nodes_vs.size(); i++) {
			if(nodes_vs[i].getDist() < minDist) {
				winNodeId_i = nodes_vs[i].getId();
				winNodeX_i = nodes_vs[i].getPosX();
				winNodeY_i = nodes_vs[i].getPosY();
				winNodeDist_d = nodes_vs[i].getDist();
				minDist = nodes_vs[i].getDist();
			}
		}
		if(verbose) {
			std::cout << "************************" << std::endl;
			std::cout << "MIN DIST IS "
			          << winNodeDist_d
			          << " AT NODE ("
			          << winNodeX_i
			          << "," << winNodeY_i
			          << ") ID "
			          << winNodeId_i
			          << std::endl;
		}
	};
	
	double getWinNodeDist() const {
		return winNodeDist_d;
	};
	
	double bmuEuclideanDist(somNode& node) {
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
	
	void adjustWeights(somNode inputnode, double lr, bool verbose = false) {
		int i;
		double adjustRate;
		int n_per_thread = nodes_vs.size()/N_THREADS_I;
		omp_set_num_threads(N_THREADS_I);
		#pragma omp parallel for private(i, winNodeDist_d) schedule(static, n_per_thread)
		for(i = 0; i < nodes_vs.size(); i++) {
			if(nodes_vs[i].getIsNb()) {
				if(verbose){
					std::cout << "WEIGHTS OF NODE ID " << nodes_vs[i].getId() << " BEFORE ADJUSTMENT" << std::endl;
					std::cout << nodes_vs[i].getWeights() << std::endl;
				}
				// Calculate by how much its weights are adjusted
				adjustRate = exp(-(winNodeDist_d) / (2*(neighbRadius*neighbRadius)));
				auto miao = nodes_vs[i].getWeights() + lr * adjustRate * (inputnode.getWeights() - nodes_vs[i].getWeights());
				nodes_vs[i].weights_vd = miao;
				if(verbose){
					std::cout << "WEIGHTS AFTER ADJUSTMENT" << std::endl;
					std::cout << nodes_vs[i].getWeights() << std::endl;
					std::cout << "************************" << std::endl;
				}
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
	
	long long getParalTime() {
		return paralTime;
	};
	
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
		
		learnRate = START_LR;
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
//			paralTime += t_neighbExplorer.getDeltaT() + t_inputVsGrid.getDeltaT();
		}
		return true;
	};
	
	~somGridAdv()= default;
	
};

#pragma clang diagnostic pop