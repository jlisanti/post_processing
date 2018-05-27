#ifndef CODES_POSTPROCESSING_SMOOTHENCODER_H_
#define CODES_POSTPROCESSING_SMOOTHENCODER_H_

#include <vector>
#include <string>

#include "datafileinput.h"

std::vector<double> smooth_encoder (DataFileInput &cInput,
		                            int enco_colmn);
#endif
