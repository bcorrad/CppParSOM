#ifndef UNTITLED_METRICS_H
#define UNTITLED_METRICS_H

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

#endif //UNTITLED_METRICS_H
