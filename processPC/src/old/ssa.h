#ifndef CODES_POSTPROCESSING_PROCESSPCG_SRC_SSA_H_
#define CODES_POSTPROCESSING_PROCESSPCG_SRC_SSA_H_

#include <vector>

#include "datafileinput.h"

bool SortType(double i, double j);

inline double covariance(int i, int j, std::vector<double> X, int N);
inline double mean(std::vector<double> X);
inline double variance(std::vector<double> X);
inline double sd(std::vector<double> X);
int ssa(DataFileInput &cInput,
		std::vector<double> &time,
		std::vector<double> &pressure,
		int M,
		int kmax,
		int colmn);

#endif
