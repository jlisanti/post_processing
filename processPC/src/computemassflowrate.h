#ifndef CODES_POSTPROCESSING_PROCESSPC_SRC_COMPUTEMASSFLOWRATE_H_
#define CODES_POSTPROCESSING_PROCESSPC_SRC_COMPUTEMASSFLOWRATE_H_

void computeMassFlowRate(double &inletTemperature,
		                 double &staticPressure,
						 double &gasConstant,
						 std::vector<double> areaFunction,
						 double frequency,
						 double &inletLength,
						 double &C,
						 double &eta,
						 double &averageMassFlowRate,
						 std::vector<double> &massFlow,
						 std::vector<double> &massFlowTime,
						 std::vector<double> pressure,
						 std::vector<double> time);

double massFlowRate(double dt,
		            double cfi,
				    double deltal,
			        double mdoti,
				    double rho,
				    double area,
				    double deltap,
					double eta);

double Cf(double deltaP,
		  double area,
		  double rho,
		  double massFlowRate,
		  double C);

#endif
