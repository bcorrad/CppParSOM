#include <iostream>
#include "chrono"

class somTimer {
	private:
	
		std::chrono::steady_clock::time_point itime_tp;
		std::chrono::steady_clock::time_point ftime_d;
		std::string funcName_s;
	public:
	
		somTimer()= default;
		
		explicit somTimer(std::string funcName) {
			funcName_s = funcName;
		}
		
		void setFunc(std::string funcName) {
			funcName_s = funcName;
		}
		
		void tic() {
			this->itime_tp = std::chrono::steady_clock::now();
		}
		
		void toc() {
			this->ftime_d = std::chrono::steady_clock::now();
		}
		
		void printDeltaT() {
			std::cout << "EXECUTION TIME " << this->funcName_s << " = " <<
			                      std::chrono::duration_cast<std::chrono::microseconds>(ftime_d - itime_tp).count()
								  << " [micros]" << std::endl;
		}
		
		long long getDeltaT() {
			auto delta = std::chrono::duration_cast<std::chrono::microseconds>(ftime_d - itime_tp);
			return delta.count();
		}
	
		~somTimer()= default;
};
