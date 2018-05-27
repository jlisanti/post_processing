#ifndef CODES_POSTPROCESSING_PROCESSPC_SRC_FITPOLY_H_
#define CODES_POSTPROCESSING_PROCESSPC_SRC_FITPOLY_H_

#include <vector>

#include "datafileinput.h"

void curve_fit(std::vector<double> x_col,
		       std::vector<double> y_col, 
			   std::vector<double> &coefficients);

#endif
