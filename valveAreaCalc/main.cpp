#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "computeballvalvearea.h"

double theta(double freq,
             double time)
{
    double tCycle = 1.0/freq;
    int nCycles   = time/tCycle;
    double tIntoCycle = time - (tCycle * double(nCycles));
    
    return (tIntoCycle/tCycle)*M_PI;
}

int main(int argc, char*argv[])
{
    std::string freqString = argv[1];
    double freq = std::stod(freqString);

    std::vector<double> areaFunction;

    valveProperties valveGeom;
	valveGeom.spanCount = 1000;
	valveGeom.inletRadius = 0.0098;
	valveGeom.ballInternalRadius = 0.0098;
	valveGeom.h = 0.0425;

	std::cout << "span "  << valveGeom.spanCount << std::endl;
	std::cout << "inlet " << valveGeom.inletRadius << std::endl;
	std::cout << "ball "  << valveGeom.ballInternalRadius << std::endl;
	std::cout << "h "     << valveGeom.h << std::endl;

	computeBallValveAreaFunction(valveGeom,
				                 1000,
								 areaFunction);

    double tFinal = 1.0;
    int steps = 1000000;
    double dt = tFinal/double(steps);

    double areaSum = 0.0;
    double time = 0.0;
   
    std::cout << "freq = " << freq << std::endl;

    std::string filename;

    filename = "area_" + std::to_string(freq) + "Hz.dat";
    double areai = 0.0;
    double areaPrev = 0.0;

    std::ofstream fout1(filename);
    for (int i = 0; i < steps; i++)
    {
        areaPrev = areai;
        areai = ballValveArea(theta(freq,time), areaFunction);
        if(i != 0)
            areaSum += dt*((areai - areaPrev)/2.0);

        fout1 << time << '\t' << areai << '\t' << theta(freq,time) << std::endl;
        time += dt;
    }

    std::cout << "freq = " << freq << " area per second = " << areaSum << std::endl;
    return 0;
}
