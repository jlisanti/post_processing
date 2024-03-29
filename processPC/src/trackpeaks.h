#ifndef CODES_POSTPROCESSING_PROCESSPC_SRC_TRACKPEAKS_H_ 
#define CODES_POSTPROCESSING_PROCESSPC_SRC_TRACKPEAKS_H_

#include <vector>

void track_peaks(DataFileInput &cInput,
		         std::vector<double> &p_vector,
				 std::vector<double> &i_vector,
				 int p_colmn,
				 int i_colmn,
				 int t_colmn,
				 std::vector<double> &t_max_p,
				 std::vector<double> &t_min_p,
				 std::vector<double> &t_max_i,
				 std::vector<double> &t_min_i,
				 std::vector<double> &p_max,
				 std::vector<double> &p_min,
				 std::vector<double> &i_max,
				 std::vector<double> &i_min,
				 std::vector<double> &downCrossOver,
				 std::vector<double> &upCrossOver);
void track_peaks_active(DataFileInput &cInput,
                        std::vector<double> &p_vector,
                        std::vector<double> &i_vector,
                        int p_colmn,
                        int i_colmn,
                        int t_colmn,
                        std::vector<double> &vt_max_p,
                        std::vector<double> &vt_min_p,
                        std::vector<double> &vt_max_i,
                        std::vector<double> &vt_min_i,
                        std::vector<double> &vp_max,
                        std::vector<double> &vp_min,
                        std::vector<double> &vi_max,
                        std::vector<double> &vi_min,
                        std::vector<double> &downCrossOver,
                        std::vector<double> &upCrossOver);

#endif
