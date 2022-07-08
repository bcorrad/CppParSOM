#ifndef UNTITLED_METRICS_H
#define UNTITLED_METRICS_H

#include <iomanip>
#include <fstream>
#include "vector"
#include "iostream"
#include "string"

long double speedUp(long double paralTime_lld, long double totalTime_lld) {
	/// speedup(P) = T_1 / T_P
	/// Where:
	/// - speedup(P) is the max of how many times that parallel computing
	/// could be faster than single-thread performance.
	/// - paralTime_lld  is the execution time when using more than one thread.
	/// - totalTime_lld is the execution time when using only one thread.
	long double speedup = (long double) totalTime_lld / (long double) paralTime_lld;
	return speedup;
}

long double efficiency(long double paralTime_lld, long double totalTime_lld, int nThr) {
	/// E(P) = T_1 / (T_P * P)
	/// Where:
	/// - E(P) is the efficiency depending on number of threads P.
	/// - totalTime_lld or T_1 is the serial time (P = 1).
	/// - paralTime_lld or T_P is the parallel time (when using threads).
	int nWorkers = nThr;
	auto efficiency = (long double)  (totalTime_lld) / (long double) (paralTime_lld * nWorkers);
	return efficiency;
}

long double cost(long double paralTime_lld, int nThr) {
	/// C(P) = T_P * P
	return (long double) paralTime_lld * nThr;
}

void report(std::vector<long double> T_, std::string lang, int weightsize, int gridrows, int gridcols, int n_inputs, int epochs){
	for(int P = 0; P < T_.size(); P++) {
		int N = P+1;
		long double T_1 = T_[0];
		long double T_P = T_[P];
		long double su = speedUp(T_P, T_1);
		long double effic = efficiency(T_P, T_1, N);
		long double costo = cost(T_P, N);
		
		std::cout <<	"  /\\_/\\ \n"
		                " ( o.o ) \n"
		                "  > ^ <  \n" << std::endl;
		std::cout << N << " THREADS" << std::endl;
		std::cout << "EXECUTION TIME [microseconds] = " << T_P << std::endl;
		std::cout << "SPEED UP = " << su << std::endl;
		std::cout << "EFFICIENCY = " << effic << std::endl;
		std::cout << "COST = " << costo << std::endl;
		std::ofstream my_file;
		std::string reportPath = "C:\\Users\\barba\\CLionProjects\\SOM_OMP\\report_"+lang+".csv";
		my_file.open(reportPath, std::ios_base::app); // append instead of overwrite
		my_file << N << ","
		        << weightsize << ","
		        << epochs << ","
		        << gridrows << ","
		        << gridcols << ","
		        << n_inputs << ","
		        << std::setprecision(15) << T_1 << ","
		        << std::setprecision(15) << T_P << ","
		        << std::setprecision(15) << su << ","
		        << std::setprecision(15) <<  effic << ","
		        << std::setprecision(15) << costo << ","
		        << lang << std::endl;
		my_file.close();
	}
}

#endif //UNTITLED_METRICS_H
