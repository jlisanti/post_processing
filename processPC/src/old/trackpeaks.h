#ifndef CODES_POSTPROCESSING_PROCESSPC_SRC_TRACKPEAKS_H_ 
#define CODES_POSTPROCESSING_PROCESSPC_SRC_TRACKPEAKS_H_

#include <vector>

void track_peaks(DataFileInput &cInput,
		         std::vector<double> &p_vector,
				 std::vector<double> &e_vector,
				 std::vector<double> &i_vector,
		         int e_colmn,
				 int p_colmn,
				 int i_colmn,
				 int t_colmn,
				 std::vector<double> &e_max_p,
				 std::vector<double> &e_min_p,
				 std::vector<double> &e_max_i,
				 std::vector<double> &e_min_i,
				 std::vector<double> &t_max_p,
				 std::vector<double> &t_min_p,
				 std::vector<double> &t_max_i,
				 std::vector<double> &t_min_i,
				 std::vector<double> &p_max,
				 std::vector<double> &p_min,
				 std::vector<double> &i_max,
				 std::vector<double> &i_min);

#endif
