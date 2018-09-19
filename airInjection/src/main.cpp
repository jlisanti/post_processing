#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include "menu.h"
#include "datafileinput.h"
#include "functions.h"

int main(int argc, char*argv[])
{
	Menu cMenu("->");

	/* read processed input file */
	if (argc <= 1)
	{
	    cMenu.display_message("Error : no input file");	
		exit(1);
	}
	else if(argc == 2)
	{
		int skip_lines = 0;
		DataFileInput cInput(argv[1],skip_lines);

		/* set parameters */
		double Tin = 300.0;
		double Pstatic = 1.00944683*1.0e5;
		double R = 287.15;
		double rho = Pstatic/(Tin*R);
		//double area = pow((9.0/2000.0),2)*M_PI*10.0;
		double area = pow(12.0/1000.0,2.0)*M_PI;
		double deltaL = 46.0/1000.0;

 		double mdoti = 0.0; 
		double C = 10.0;
		double eta = 0.15; 

		std::vector<double> mdot;
		std::vector<double> t;
		std::vector<double> dt_vec;

		/* loop through pressure trace to compute mdot */
		for (int i = 0; i < cInput.file_length()-1; i++)
		{
			double dt = cInput.table_value(i+1,0)
				        - cInput.table_value(i,0);

			double deltaP = -Pstatic 
				             + cInput.table_value(i,1)*1.0e5;

			if((dt > 0.0))
			{
		    	double Cfi = Cf(deltaP,
		    			        area,
		    					rho,
		   					    mdoti,
		   					    C);
				//Cfi = 5.0;
			    double mdotii = massFlowRate(dt,
					                         Cfi,
						    				 deltaL,
						    				 mdoti,
						    				 rho,
						    				 area,
						    				 deltaP,
											 eta);
				if(mdotii < 0.0)
					mdotii = 0.0;
			    std::cout << dt     << '\t'
				          << deltaP << '\t'
					      << Cfi    << '\t'
					      << mdotii << std::endl;
			    mdot.push_back(mdotii);
			    t.push_back(cInput.table_value(i,0));
				dt_vec.push_back(dt);
			    mdoti = mdotii;
			}
		}
		std::ofstream fout("massFlowRate.dat");

		double totalMass = 0.0;

		for (int i = 0; i < mdot.size(); i++)
		{
			fout << t[i] << '\t' << mdot[i] << std::endl; 
			totalMass = totalMass + mdot[i]*dt_vec[i];
		}
		std::cout << "Mean mass flow rate : " 
			      << totalMass/(t[t.size()-1]-t[0]) 
				  << std::endl;
	}
	else
	{
		cMenu.display_message("Error!");
	}
	return 0;
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

	/*
	return  - ((dt*Cfi)/(2.0*deltaL))
		           *(mdoti/(rho*area))
				   *std::abs(mdoti)
				 - ((dt*deltaP*area)/deltaL);
				 */
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
