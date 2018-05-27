#ifndef CODES_POSTPROCESSING_PROCESSPC_SRC_AVERAGECYCLE_H_
#define CODES_POSTPROCESSING_PROCESSPC_SRC_AVERAGECYCLE_H_

void average_cycle(DataFileInput &cInput,
		           std::vector<double> &time_ssa,
				   std::vector<double> &pressure_ssa,
				   std::vector<double> &ion_ssa,
		           int encoder_colmn,
				   int pressur_colmn,
				   double window,
				   std::vector<double> &e_vector,
				   std::vector<double> &p_vector,
				   std::vector<double> &i_vector,
				   std::vector<double> &pressure_scatter,
				   std::vector<double> &encoder_scatter,
				   std::vector<double> &ion_scatter);

#endif
