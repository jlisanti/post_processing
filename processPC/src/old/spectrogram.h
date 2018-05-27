#ifndef CODES_POSTPROCESSING_PROCESSPCG_SRC_SPECTROGRAM_H_
#define CODES_POSTPROCESSING_PROCESSPCG_SRC_SPECTROGRAM_H_

#include <vector>

#include "datafileinput.h"

void build_spectrogram(DataFileInput &cInput,
					   double &freq,
					   std::vector<std::vector<double>> &spectrogram_f,
					   std::vector<std::vector<double>> &spectrogram_m,
					   int &n,
					   double &dt,
					   std::vector<double> &time);


#endif
