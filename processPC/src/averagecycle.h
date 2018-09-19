#ifndef CODES_POSTPROCESSING_PROCESSPC_SRC_AVERAGECYCLE_H_
#define CODES_POSTPROCESSING_PROCESSPC_SRC_AVERAGECYCLE_H_

void average_cycle_active(DataFileInput &cInput,
		                  double window,
						  double stepSize,
						  std::vector<double> &time_ssa,
						  std::vector<double> &pressure_ssa,
						  std::vector<double> &ion_ssa,
						  int encoder_colmn,
						  int pressur_colmn,
						  std::vector<double> &e_vector,
						  std::vector<double> &p_vector,
						  std::vector<double> &i_vector,
						  std::vector<double> &pressure_scatter,
						  std::vector<double> &encoder_scatter,
						  std::vector<double> &ion_scatter,
						  std::vector<double> &up_cross_over,
						  std::vector<double> &down_cross_over,
						  std::vector<double> &areaFunction,
						  std::string printFromStart);

void average_cycle_passive(DataFileInput &cInput,
		                   int pressur_colmn,
						   int ion_colmn,
						   double window,
						   double step,
                           std::vector<double> &phase_average,
						   std::vector<double> &time_average, 
						   std::vector<double> &pressure_average,
						   std::vector<double> &ion_average,
						   std::vector<double> &phase_scatter,
						   std::vector<double> &time_scatter,
						   std::vector<double> &pressure_scatter,
						   std::vector<double> &ion_scatter,
						   std::vector<double> &down_cross_over,
						   std::vector<double> &up_cross_over,
						   double frequency,
                           double static_pressure,
						   std::string printFromStart);

void average_files(std::vector<double> time_scatter,
		           std::vector<double> pressure_scatter,
				   std::vector<double> ion_scatter,
				   std::vector<double> &time_average,
				   std::vector<double> &pressure_average,
				   std::vector<double> &ion_average,
				   std::vector<double> &phase_average);

void average_files_active(std::vector<double> encoder_scatter,
		                  std::vector<double> pressure_scatter,
				          std::vector<double> ion_scatter,
				          std::vector<double> &time_average,
				          std::vector<double> &pressure_average,
				          std::vector<double> &ion_average,
				          std::string ion_probe);

#endif
