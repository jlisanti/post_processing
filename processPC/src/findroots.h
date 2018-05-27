#include <vector>

/*
void find_max(std::vector<double> coefficients, double &max, double &phase_max,
		                                        double &min, double &phase_min);
												*/
void find_min(std::vector<double> x_col,
		      std::vector<double> y_col,
			  double &min,
			  double &phase_min);

void find_max(std::vector<double> x_col,
		      std::vector<double> y_col,
			  double &max,
			  double &phase_max);
