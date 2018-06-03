#ifndef CODES_POSTPROCESSING_PROCESSPCG_SRC_FFT_H_
#define CODES_POSTPROCESSING_PROCESSPCG_SRC_FFT_H_

#include "datafileinput.h"

void compute_frequency_spectrum(DataFileInput &cInput,
		                        std::vector<double> &magnitude,
								std::vector<double> &frequency,
								int col,
								double &peak_frequency);


#endif
