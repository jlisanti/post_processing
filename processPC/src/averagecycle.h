#ifndef CODES_POSTPROCESSING_PROCESSPC_SRC_AVERAGECYCLE_H_
#define CODES_POSTPROCESSING_PROCESSPC_SRC_AVERAGECYCLE_H_

void average_cycle_active(DataFileInput &cInput,
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

void average_cycle_passive(DataFileInput &cInput,
		                   int pressur_colmn,
						   int ion_colmn,
						   double window,
                           std::vector<double> &position_average,
						   std::vector<double> &pressure_average,
						   std::vector<double> &ion_average,
						   std::vector<double> &position_scatter,
						   std::vector<double> &pressure_scatter,
						   std::vector<double> &ion_scatter,
                           double static_pressure);

#endif
