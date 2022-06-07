#include <iostream>
#include "chrono"

class somTimer {
private:
	std::chrono::steady_clock::time_point itime_tp;
	std::chrono::steady_clock::time_point ftime_d;
//	long exectime_ld;
public:
	somTimer()= default;
	
	void tic() {
		this->itime_tp = std::chrono::steady_clock::now();
	};
	
	void toc(std::string optional = "", bool verbose = false) {
		this->ftime_d = std::chrono::steady_clock::now();
//		this->exectime_ld = std::chrono::duration_cast<std::chrono::microseconds>(ftime_d - itime_tp).count();
		if(verbose) std::cout << "EXECUTION TIME " << optional << " = " <<
		                      std::chrono::duration_cast<std::chrono::milliseconds>(ftime_d - itime_tp).count() << " [ms]" << std::endl;
	};
	
};
