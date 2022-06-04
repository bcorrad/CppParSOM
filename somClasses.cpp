#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include <iostream>
#include <cstdlib>     /* srand, rand */
#include <cmath>
#include <utility>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "algorithm"
#include "iterator"
#include "omp.h"

constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 1;
constexpr int WEIGHT_SIZE = 4;

//the number of epochs desired for the training
constexpr int EPOCHS = 10000;
//the value of the learning rate at the start of training
constexpr double START_LR = 0.1;

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
		std::vector<T> result_v;
		for(int i = 0; i < v1.size(); i++) {
			result_v.push_back(v1[i]-v2[i]);
		}
		return result_v;
	}

template <typename T>
	std::vector<T> operator *(const T factor, const std::vector<T>& v)
	{
		std::vector<T> result_v;
		for(int i = 0; i < v.size(); i++) {
			result_v.push_back(v[i]*factor);
		}
		return result_v;
	}


template <typename T>
	std::vector<T> operator /(const std::vector<T>& v, const T factor)
	{
		std::vector<T> result_v;
		for(int i = 0; i < v.size(); i++) {
			result_v.push_back(v[i]/factor);
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
		
		somNode(int nodeId, int nodeX, int nodeY, int threadId, int weightsSize) :
				nodeId_i(nodeId), posX_i(nodeX), posY_i(nodeY), dist_d(0),
				nodeCreatorId_i(threadId), weightsSize_i(weightsSize), isNeigh_b(false) {
			if(verbose) std::cout << "Node " << nodeId << " initialized by thread " << threadId << std::endl;
			initWeights(weightsSize);
		};
		
		void initWeights(int weightsSize) {
//			#pragma omp parallel for
				for(int i = 0; i < weightsSize; i++) {
					float randNum = FLOAT_MIN + (float)(rand()) / ((float)(RAND_MAX/(FLOAT_MAX - FLOAT_MIN)));
					weights_vd.push_back(randNum);
				}
		};
		
		double euclDistance(std::vector<double> inputWeights_vd) {
			auto distance = (double) 0;
//			#pragma omp parallel for
				for(int i = 0; i < weights_vd.size(); ++i) {
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
public:
	
	explicit inputNodes(int nInput){
		
		std::string fname = "C:\\Users\\barba\\CLionProjects\\PPF_SOM\\iris.csv";
		std::vector<double> content;
		std::vector<double> row;
		std::string line, word;
		auto rowMax = (double) 0;
		std::fstream file (fname, std::ios::in);
		int id = 0;
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
				
				somNode n = somNode(id,
				                    0,
				                    0,
				                    omp_get_thread_num(),
				                    WEIGHT_SIZE);
				n.weights_vd = content;
				addNode(n);
			}
		}
		else
			std::cout << "Could not open the file" << std::endl;
		
		/*for(int i=0;i<content.size();i++)
		{
			for(int j=0;j<content[i].size();j++)
			{
				std::cout<<content[i][j]<<" ";
			}
			std::cout<<"\n";
		}*/
		
//			for(int d = 0; d < nInput; d++) {
//				somNode n = somNode(d,
//				                    0,
//				                    0,
//				                    omp_get_thread_num(),
//				                    WEIGHT_SIZE);
//				n.weights_vd = content;
//				addNode(n);
//			}
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
		int somRows, somCols;
		std::vector<somNode> nodes_vs;
		
		int winNodeId_i, winNodeX_i, winNodeY_i;
		double winNodeDist_d;
		
		double somInitRadius;
		double neighbRadius; // sigma0
		double learnRate;
	public:
		somGrid()= default;
		
		somGrid(int rows, int cols) : somRows(rows), somCols(cols) {
			int id = 0;
			for(int r = 0; r < rows; r++) {
				for(int c = 0; c < cols; c++) {
					somNode n = somNode(id,
										r,
										c,
										omp_get_thread_num(),
										WEIGHT_SIZE);
					addNode(n);
					id++;
				}
			}
		};
		
		double getNeighbRadius() const {
			return neighbRadius;
		};
		
		void addNode(somNode &node) {
			nodes_vs.push_back(node);
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
			for(auto & nodes_v : nodes_vs) {
				double dist = nodes_v.euclDistance(inputNode.getWeights());
				nodes_v.setDist(dist);
				if(verbose) {
					std::cout << "INPUT NODE "
					          << inputNode.getId()
					          << " HAS DISTANCE -> "
					          << nodes_v.getDist()
					          << " WRT NODE ("
					          << nodes_v.getPosX() << ","
					          << nodes_v.getPosY() << ")"
					          << std::endl;
				}
			}
			BMU();
		};
		
		void BMU(bool verbose = false) {
			auto minDist = (double) 10000;
			for(auto & nodes_v : nodes_vs) {
				if(nodes_v.getDist() < minDist) {
					winNodeId_i = nodes_v.getId();
					winNodeX_i = nodes_v.getPosX();
					winNodeY_i = nodes_v.getPosY();
					winNodeDist_d = nodes_v.getDist();
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
					pow((double) node.getPosY() -(double) winNodeY_i, 2);
			return dist;
		};
		
		void neighbExplorer(double radius) {
			long double start = omp_get_wtime();
			#pragma omp parallel for
				for(auto & node_v : nodes_vs) {
					double D = bmuEuclideanDist(node_v);
					if(D < radius)
						node_v.setIsNb(true);
				}
			long double end = omp_get_wtime();
			printf("Work took %ld seconds\n", end - start);
		};
		
		void adjustWeights(somNode inputnode, double lr, bool verbose = false) {
			#pragma omp parallel for
				for(auto & node_v : nodes_vs) {
					if(node_v.getIsNb()) {
						if(verbose){
							std::cout << "WEIGHTS OF NODE ID " << node_v.getId() << " BEFORE ADJUSTMENT" << std::endl;
							std::cout << node_v.getWeights() << std::endl;
						}
						// Calculate by how much its weights are adjusted
						double adjustRate = exp(-(getWinNodeDist()) / (2*(neighbRadius*neighbRadius)));
						auto miao = node_v.getWeights() + lr * adjustRate * (inputnode.getWeights() - node_v.getWeights());
						node_v.weights_vd = miao;
						if(verbose){
							std::cout << "WEIGHTS AFTER ADJUSTMENT" << std::endl;
							std::cout << node_v.getWeights() << std::endl;
							std::cout << "************************" << std::endl;
						}
					}
				}
		};
	
		void resetGrid() {
			#pragma omp parallel for
				for(auto & node_v : nodes_vs){
					node_v.setIsNb(false);
				}
		};
		
		void printGrid() {
			for(auto & node : nodes_vs) {
				std::cout << node.getWeights() << std::endl;
			}
		};
		
		bool somTrain(inputNodes inputs) {
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
			// Time constant (lambda) used in the calculation of the neighbourhood width
			double timeConst = EPOCHS/log(neighbRadius);
			// Calculate the biggest possible initial radius (half of width either height (delta_0))
			somInitRadius = (double) std::max(somRows, somCols) / 2;
			
			for(int epoch = 0; epoch < EPOCHS; epoch++) {
				if(epoch % (EPOCHS/10) == 0) std::cout << "========= EPOCH "
				                                       << epoch << "/" << EPOCHS
				                                       << " =========" << std::endl;
				
				// Calculate the width of the neighbourhood for this timestep
				neighbRadius = somInitRadius * exp(-(double)epoch/timeConst);
				// Reduce the learning rate
				learnRate = START_LR * exp(-(double)epoch/EPOCHS);
				
				for(auto & inputnode : inputs.getNodes(0)){
					// Calculating the Best Matching Unit
					inputVsGrid(inputnode);
					// BMU neighborhood exploration
					neighbExplorer(neighbRadius);
					// Adjust weights
					adjustWeights(inputnode, learnRate);
				}
				std::cout << "=================================" << std::endl;
				printGrid();
				resetGrid();
			}
			return true;
		};
		
		~somGrid()= default;
	
	
};

//#pragma clang diagnostic pop