#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include "menu.h"
#include "datafileinput.h"
#include "computeliquidfuelflowrate.h"

double computeLiquidFuelFlowRate(double frequency,
                                 double fuelPressure,
						         std::vector<double> pressure)
{

	std::cout << " Computing mass flow rate" << std::endl;

	double cycleLength = 1.0/frequency;
    double totalMass = 0.0;

	for (int i = 0; i < pressure.size(); i++)
    {
        double dP = fuelPressure - pressure[i]*100000.0;
		totalMass = totalMass + qm(dP);
    }

	return totalMass/cycleLength;
}

double qm(double dP)
{
	double Cd = 0.7;
    double rho = 748.9; // kg/m3
	double d = 0.0007;  // m
    double D = 0.004;   // m
    double beta = d/D;
	return (Cd/sqrt(1.0-pow(beta,4)))*(M_PI/4.0)*pow(d,2)*sqrt(2.0*dP*rho);
}
