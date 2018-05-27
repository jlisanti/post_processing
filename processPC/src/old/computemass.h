#ifndef CODES_PROCESSPC_SRC_COMPUTEMASS_H_
#define CODES_PROCESSPC_SRC_COMPUTEMASS_H_

#include <vector>

void compute_mass_flow_rate(std::vector<double> p,
		                    std::vector<double> e,
							double mdot_fuel_1,
							double mdot_fuel_2,
							double r_b,
							double r_i,
							double l,
							double freq,
							double &mdot_air,
							double &phi,
							double &AF_ratio,
							std::vector<double> &mdot_inlet,
							std::vector<double> &area_inlet,
							std::vector<double> &area_inlet_2);
double area_func(double l,
		         double sr,
				 double lr,
				 double alpha);
double d_func(double l,
		      double r_b,
			  double r_i,
			  double alpha);
#endif
