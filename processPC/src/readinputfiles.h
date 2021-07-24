#ifndef CODES_POSTPROCESSING_PROCESSPCG_SRC_READINPUTFILES_H_
#define CODES_POSTPROCESSING_PROCESSPCG_SRC_READINPUTFILES_H_

#include <string>
#include <vector>
#include <iostream>

struct main_menu {
	std::string experiment;
	std::string combustorType;
	std::string fuelType;
	std::string ionProbe;
	int numberFuels;
	std::string fuelA;
	std::string fuelB;
	std::vector<double> fuelPercent;
	double gasoline;
	double heptane;
	double ethanol;
	double diesel;
};


struct data_conditioning_menu {
	std::string SSA;
	int vectorDimensionM;
	int kMax;
	double windowSize;
	double stepSize;
	double pressureOffset;
	std::string shiftCurve;
	std::string shiftFixedCurve;
	double shiftFixed;
};

struct data_analysis_menu {
	std::string findPeaks;
	std::string averageCycle;
	std::string freqSpectrum;
	std::string spectrogram;
	double spectNcyclesWindow;
	double spectNcyclesSpace;
	std::string trackPeaks;
	std::string bodePlots;
	std::string computeRMS;
	std::string computeMassFlowRate;
	double inletTemperature;
	double gasConstant;
	double inletArea;
	double inletLength;
	double C;
	double eta;
	int spanCount;
	double inletRadius;
	double ballInternalRadius;
	double lengthBall;
	double staticPressure;
	std::string setMeanP;
	double meanP;
};

struct output_menu {
	std::string printFromOpen;
	std::string shiftCurve;
};

struct options {
	main_menu mainMenu;
	data_conditioning_menu dataCondMenu;
	data_analysis_menu dataAnalysisMenu;
	output_menu outputMenu;
};


void read_input_files(std::string file, std::vector<std::string> &files_names);
void read_settings_file(std::string file,
		                options &optionsMenu);

#endif
