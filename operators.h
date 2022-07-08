#ifndef UNTITLED_OPERATORS_H
#define UNTITLED_OPERATORS_H

#include "iostream"
#include "vector"

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
#endif //UNTITLED_OPERATORS_H
