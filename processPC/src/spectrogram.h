#ifndef CODES_POSTPROCESSING_PROCESSPCG_SRC_SPECTROGRAM_H_
#define CODES_POSTPROCESSING_PROCESSPCG_SRC_SPECTROGRAM_H_

#include <vector>

#include "datafileinput.h"

void build_spectrogram(DataFileInput &cInput,
					   double &freq,
					   std::vector<std::vector<double>> &spectrogram_f,
					   std::vector<std::vector<double>> &spectrogram_m,
					   double numberCyclesSpace,
					   double numberCyclesWindow,
					   int &n,
					   int col,
					   std::vector<double> &time);


#endif
