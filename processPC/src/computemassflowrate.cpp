#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include "menu.h"
#include "datafileinput.h"
#include "computeballvalvearea.h"
#include "computemassflowrate.h"

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
						 std::vector<double> encoder)
{
	/* set parameters */
	/*
	double Tin = 300.0;
	double Pstatic = 1.00944683*1.0e5;
	double R = 287.15;
	//double area = pow((9.0/2000.0),2)*M_PI*10.0;
	double area = pow(12.0/1000.0,2.0)*M_PI;
	double deltaL = 46.0/1000.0;

	double C = 10.0;
	double eta = 0.15; 
	*/

	/* compute time vector */
	double cycleLength = 1.0/frequency;
	double dt = cycleLength/double(pressure.size());

	double thetaStep = M_PI/double(pressure.size());

 	double mdoti = 0.0; 
	double rho = staticPressure/(inletTemperature*gasConstant);

	std::vector<double> mdot;
	std::vector<double> t;
	std::vector<double> dt_vec;

	/* loop through pressure trace to compute mdot */
	for (int i = 0; i < pressure.size(); i++)
	{
		//double dt = time[i+1] - time[i];

		double deltaP = - staticPressure + pressure[i]*1.0e5;
		double inletArea = ballValveArea(double(i)*thetaStep,
				                         areaFunction);
				  

		/*
		std::cout << std::endl;
		std::cout << "dp: "   << deltaP << '\n'
			      << "area: " << inletArea << std::endl; 
		std::cout << std::endl;
		*/
		if((dt > 0.0))
		{
			if(inletArea > 0.0)
			{
		        double Cfi = Cf(deltaP,
		        			    inletArea,
		        				rho,
		   	    				mdoti,
		   	    				C);
			    //Cfi = 10.0;
			    double mdotii = massFlowRate(dt,
			    		                     Cfi,
			    			    			 inletLength,
			    			    			 mdoti,
			    			    			 rho,
			    			    			 inletArea,
			    			    			 deltaP,
			    							 eta);
				/*
			    if(mdotii < 0.0)
			    	mdotii = 0.0;
					*/

			    mdot.push_back(mdotii);
			    t.push_back(double(i)*dt);
			    dt_vec.push_back(dt);
			    mdoti = mdotii;
			}
			else
			{
			    double mdotii = 0.0;

			    mdot.push_back(mdotii);
			    t.push_back(double(i)*dt);
			    dt_vec.push_back(dt);
			    mdoti = mdotii;
			}
		}
	}

	double totalMass = 0.0;

	for (int i = 0; i < mdot.size(); i++)
		totalMass = totalMass + mdot[i]*dt_vec[i];

	massFlow = mdot;
	massFlowTime = t;

	averageMassFlowRate = totalMass/(t[t.size()-1] - t[0]);

}

double massFlowRate(double dt,
		            double Cfi,
				    double deltaL,
			        double mdoti,
				    double rho,
				    double area,
				    double deltaP,
					double eta)
{
	return mdoti - (eta*(((dt*Cfi)/(2.0*deltaL))
		            *(mdoti/(rho*area))
				    *std::abs(mdoti)
				 + ((dt*deltaP*area)/deltaL)));
}

double Cf(double deltaP,
		  double area,
		  double rho,
		  double massFlowRate,
		  double C)
{
	if(massFlowRate == 0.0)
		return 0.0;
	else
		return C*(2.0*rho*pow(area,2)*(deltaP/pow(massFlowRate,2)));
}
