#ifndef CODES_POSTPROCESSING_PROCESSPC_SRC_COMPUTELIQUIDFUELFLOWRATE_H_
#define CODES_POSTPROCESSING_PROCESSPC_SRC_COMPUTELIQUIDFUELFLOWRATE_H_

#include <vector>

double computeLiquidFuelFlowRate(double frequency,
		                             double fuelPressure,
						             std::vector<double> pressure);

double qm(double dP);

#endif
