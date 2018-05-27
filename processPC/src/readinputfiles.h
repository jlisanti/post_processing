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
};


struct data_conditioning_menu {
	std::string SSA;
	int vectorDimensionM;
	int kMax;
	double windowSize;
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
};

struct output_menu {
	std::string printFromOpen;
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
