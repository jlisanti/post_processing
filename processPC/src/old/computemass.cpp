#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

#include "computemass.h"

void compute_mass_flow_rate(std::vector<double> p, 
		                    std::vector<double> enc,
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
							std::vector<double> &area_inlet_2)
{
	double p_mean = 100000.0*accumulate( p.begin(), p.end(), 0.0)/p.size(); 
	double air_mass = 0.0;
	double rho = p_mean/(287.15*300.0);
	for (int i = 0; i < p.size()-1; i++)
	{
		double p_c = 100000.0*p[i];
		double alpha = (enc[i]/2.5)*2.0*M_PI;
		double rate = 0.0; 
		if ((p_mean - p_c) > 0.0) 
			rate = rho*sqrt((2.0/rho)*((p_mean-p_c)))
			          *area_func(l,r_b,r_i,alpha);
		else 
			rate = -rho*sqrt((2.0/rho)*((p_c-p_mean)))
			          *area_func(l,r_b,r_i,alpha);
		double dt = (((enc[i+1]/2.5)*2.0*M_PI) - ((enc[i]/2.5)*2.0*M_PI))/freq;
		mdot_inlet.push_back(rate);

		area_inlet.push_back(area_func(l,r_b,r_i,alpha));
		//area_inlet_2.push_back(area_func(60.0/1000.0,35.0/2000.0,35.0/2000.0,alpha));

	    air_mass = air_mass + dt*rate; 	
	}
	mdot_air = air_mass/(1.0/freq);
	double M_H2 = 2.016;
	double M_CH4 = 16.04;
	double M_air = 28.97;
	double H2_CH4_ratio = (mdot_fuel_2/M_H2)/(mdot_fuel_1/M_CH4);
	double a = H2_CH4_ratio/(1.0+H2_CH4_ratio);
	double b = 1.0 - a;
	double d = b;
	double e = a + 2.0*b;
	double c = 2.0*d + e;
	double AF_st = (c*M_air)/((b*M_CH4) + (a*M_H2)); 
    double AF = (mdot_air)/(mdot_fuel_1/1000.0 + mdot_fuel_2/1000.0);
	phi = AF/AF_st;
}

double area_func(double l,
		         double sr,
				 double lr,
				 double alpha)
{
	if (sr > lr)
	{
		double tmp = lr;
		lr = sr;
		sr = tmp;
	}
	double d = d_func(l,sr,lr,alpha);
	if (d > (sr+lr))
		return 0.0;
	else if (d == 0.0)
		return M_PI*pow(sr,2);
	else if (d < (lr - sr))
	{
		return M_PI*pow(sr,2);
	}
	else
		return pow(sr,2)*acos((pow(d,2)+pow(sr,2)-pow(lr,2))
				/(2.0*d*sr)) 
							+ pow(lr,2)*acos((pow(d,2)+pow(lr,2)-pow(sr,2))
									/(2.0*d*lr)) 
							- 0.5*sqrt((-d+sr+lr)*(d+sr-lr)*(d-sr+lr)*(d+sr+lr))
							;
}

double d_func(double l,
		      double r_b,
			  double r_i,
			  double alpha)
{
	double h = sqrt(pow(l,2) + pow(2.0*r_i,2));
	double alpha_close = 4.0*asin((2.0*r_i)/h);
	if (alpha < alpha_close)
	{
		return ((r_i+r_b)/alpha_close)*alpha;
	}
	else if (alpha > ((360.0/180.0)*M_PI - alpha_close)) 
	{
		return (r_i+r_b) + (-((r_i+r_b)/alpha_close)
				*(alpha-(360.0/180.0*M_PI)) - (r_i+r_b));
	}
	else
		return 10000.0*r_i;
}
