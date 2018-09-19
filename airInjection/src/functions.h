#ifndef CODES_POSTPROCESSING_AIRINJECTION_SRC_FUNCTIONS_H_
#define CODES_POSTPROCESSING_AIRINJECTION_SRC_FUNCTIONS_H_

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
