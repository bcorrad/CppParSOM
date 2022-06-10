#ifndef UNTITLED_METRICS_H
#define UNTITLED_METRICS_H

//long double speedUp(long long paralTime_lld, long long totalTime_lld, int nThr) {
//	/// speedup = 1/((1-f)+(f/N_workers))
//	/// Where:
//	/// - speedup is the max of how many times that parallel computing
//	/// could be faster than single-thread performance.
//	/// - f is the fraction of the algorithm that can be parallelized.
//	/// - N_workers is the number of threads currently used (a.k.a. processors).
//	long double f = (long double) paralTime_lld/(long double) totalTime_lld;
//	auto nWorkers = (double) nThr;
//	auto speedup = 1/((1-f)+(f/(long double) nWorkers));
//	return speedup;
//}

long double speedUp(long long paralTime_lld, long long totalTime_lld) {
	/// speedup(P) = T_1/T_P
	/// Where:
	/// - speedup(P) is the max of how many times that parallel computing
	/// could be faster than single-thread performance.
	/// - paralTime_lld  is the execution time when using more than one thread.
	/// - totalTime_lld is the execution time when using only one thread.
	long double speedup = (long double) totalTime_lld/(long double) paralTime_lld;
	return speedup;
}

long double efficiency(long long paralTime_lld, long long totalTime_lld, int nThr) {
	/// E(P) = T_1 / (T_P * P)
	/// Where:
	/// - E(P) is the efficiency depending on number of threads P.
	/// - totalTime_lld or T_1 is the serial time (P = 1).
	/// - paralTime_lld or T_P is the parallel time (when using threads).
	int nWorkers = nThr;
	auto efficiency = (totalTime_lld-paralTime_lld)/(paralTime_lld*nWorkers);
	return efficiency;
}

long double cost(long long paralTime_lld, int nThr) {
	return paralTime_lld * nThr;
}

#endif //UNTITLED_METRICS_H
