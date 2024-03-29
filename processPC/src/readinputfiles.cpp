#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>

#include <boost/algorithm/string.hpp>

#include "readinputfiles.h"

void read_input_files(std::string file, std::vector<std::string> &file_name)
{
	std::cout << "reading input files" << std::endl;
	std::ifstream input_file_stream(file);

	while (true)
	{
		std::string line;
		std::getline (input_file_stream, line);

		if (!input_file_stream)
			break;

		if (strncmp(line.c_str(),"files {",7)==0)
		{
			while (!(strncmp(line.c_str(),"}",1)==0))
			{
				std::getline (input_file_stream, line);
				if (!(strncmp(line.c_str(),"}",1)==0))
					file_name.push_back(line);
			}
		}

	}
}

void read_settings_file(std::string file,
		                options &optionsMenu)
{
	std::ifstream input_file_stream(file);
	while (true)
	{
		std::string line;
		std::getline (input_file_stream, line, '\t');
		boost::trim(line);

		if (!input_file_stream)
			break;
		if (strncmp(line.c_str(),"dataType {",10)==0)
		{

			/* Read primary settings */
			while (!(strncmp(line.c_str(),"}",1)==0))
			{
				std::string setting;
				std::string value;
				std::getline (input_file_stream, line, ';');
				std::istringstream ss(line);
				boost::trim(line);
				if(strncmp(line.c_str(),"}",1)==0)
					break;
				ss >> setting >> value;
				boost::trim(setting);
				boost::trim(value);
				if(setting=="experiment")
					optionsMenu.mainMenu.experiment=value;
				if(setting=="combustorType")
					optionsMenu.mainMenu.combustorType=value;
				if(setting=="fuelType")
					optionsMenu.mainMenu.fuelType=value;
				if(setting=="ionProbe")
					optionsMenu.mainMenu.ionProbe=value;
				if(setting=="numberFuels")
					optionsMenu.mainMenu.numberFuels=atoi(value.c_str());
				if(setting=="fuelA")
					optionsMenu.mainMenu.fuelA=value;
				if(setting=="fuelB")
					optionsMenu.mainMenu.fuelB=value;
				if(setting==optionsMenu.mainMenu.fuelA)
					optionsMenu.mainMenu.fuelPercent.push_back(atof(value.c_str()));
				if(setting==optionsMenu.mainMenu.fuelB)
					optionsMenu.mainMenu.fuelPercent.push_back(atof(value.c_str()));
				if(setting=="gasoline")
					optionsMenu.mainMenu.gasoline=atof(value.c_str());
				if(setting=="heptane")
					optionsMenu.mainMenu.heptane=atof(value.c_str());
				if(setting=="ethanol")
					optionsMenu.mainMenu.ethanol=atof(value.c_str());
				if(setting=="diesel")
					optionsMenu.mainMenu.diesel=atof(value.c_str());
			}
		}

		if (strncmp(line.c_str(),"dataConditioning {",17)==0)
		{
			/* read data conditioning settings */
			while (!(strncmp(line.c_str(),"}",1)==0))
			{
				std::string setting;
				std::string value;
				std::getline (input_file_stream, line, ';');
				std::istringstream ss(line);
				boost::trim(line);
				if(strncmp(line.c_str(),"}",1)==0)
					break;
				ss >> setting >> value;
				boost::trim(setting);
				boost::trim(value);
				if(setting=="SSA")
					optionsMenu.dataCondMenu.SSA=value;
				if(setting=="vectorDimensionM")
					optionsMenu.dataCondMenu.vectorDimensionM=atof(value.c_str());
				if(setting=="kMax")
					optionsMenu.dataCondMenu.kMax=atof(value.c_str());
				if(setting=="windowSize")
					optionsMenu.dataCondMenu.windowSize=atof(value.c_str());
				if(setting=="stepSize")
					optionsMenu.dataCondMenu.stepSize=atof(value.c_str());
				if(setting=="pressureOffset")
					optionsMenu.dataCondMenu.pressureOffset=atof(value.c_str());
				if(setting=="shiftCurve")
					optionsMenu.dataCondMenu.shiftCurve=value;
				if(setting=="shiftFixedCurve")
					optionsMenu.dataCondMenu.shiftFixedCurve=value;
				if(setting=="shiftFixed")
					optionsMenu.dataCondMenu.shiftFixed=atof(value.c_str());
			}
		}

		if (strncmp(line.c_str(),"dataAnalysis {",14)==0)
		{
			/* Read data analysis settings */
			while (!(strncmp(line.c_str(),"}",1)==0))
			{
				std::string setting;
				std::string value;
				std::getline (input_file_stream, line, ';');
				std::istringstream ss(line);
				boost::trim(line);
				if(strncmp(line.c_str(),"}",1)==0)
					break;
				ss >> setting >> value;
				boost::trim(setting);
				boost::trim(value);
				if(setting=="findPeaks")
					optionsMenu.dataAnalysisMenu.findPeaks=value;
				if(setting=="averageCycle")
					optionsMenu.dataAnalysisMenu.averageCycle=value;
				if(setting=="freqSpectrum")
					optionsMenu.dataAnalysisMenu.freqSpectrum=value;
				if(setting=="spectrogram")
					optionsMenu.dataAnalysisMenu.spectrogram=value;
				if(setting=="spectNcyclesWindow")
					optionsMenu.dataAnalysisMenu.spectNcyclesWindow=atof(value.c_str());
				if(setting=="spectNcyclesSpace")
					optionsMenu.dataAnalysisMenu.spectNcyclesSpace=atof(value.c_str());
				if(setting=="trackPeaks")
					optionsMenu.dataAnalysisMenu.trackPeaks=value;
				if(setting=="bodePlot")
					optionsMenu.dataAnalysisMenu.bodePlots=value;
				if(setting=="computeRMS")
					optionsMenu.dataAnalysisMenu.computeRMS=value;
				if(setting=="computeMassFlowRate")
					optionsMenu.dataAnalysisMenu.computeMassFlowRate=value;
				if(setting=="inletTemperature")
					optionsMenu.dataAnalysisMenu.inletTemperature=atof(value.c_str());
				if(setting=="gasConstant")
					optionsMenu.dataAnalysisMenu.gasConstant=atof(value.c_str());
				if(setting=="inletArea")
					optionsMenu.dataAnalysisMenu.inletArea=atof(value.c_str());
				if(setting=="inletLength")
					optionsMenu.dataAnalysisMenu.inletLength=atof(value.c_str());
				if(setting=="C")
					optionsMenu.dataAnalysisMenu.C=atof(value.c_str());
				if(setting=="eta")
					optionsMenu.dataAnalysisMenu.eta=atof(value.c_str());
				if(setting=="spanCount")
					optionsMenu.dataAnalysisMenu.spanCount=atoi(value.c_str());
				if(setting=="inletRadius")
					optionsMenu.dataAnalysisMenu.inletRadius=atof(value.c_str());
				if(setting=="ballInternalRadius")
					optionsMenu.dataAnalysisMenu.ballInternalRadius=atof(value.c_str());
				if(setting=="lengthBall")
					optionsMenu.dataAnalysisMenu.lengthBall=atof(value.c_str());
				if(setting=="inletTemperature")
					optionsMenu.dataAnalysisMenu.inletTemperature=atof(value.c_str());
				if(setting=="staticPressure")
					optionsMenu.dataAnalysisMenu.staticPressure=atof(value.c_str());
				if(setting=="gasConstant")
					optionsMenu.dataAnalysisMenu.gasConstant=atof(value.c_str());
				if(setting=="inletLength")
					optionsMenu.dataAnalysisMenu.inletLength=atof(value.c_str());
				if(setting=="setMeanP")
					optionsMenu.dataAnalysisMenu.setMeanP=value;
				if(setting=="meanP")
					optionsMenu.dataAnalysisMenu.meanP=atof(value.c_str());
			}
		}

		if (strncmp(line.c_str(),"output {",8)==0)
		{
			/* Read output settings */
			while (!(strncmp(line.c_str(),"}",1)==0))
			{
				std::string setting;
				std::string value;
				std::getline (input_file_stream, line, ';');
				std::istringstream ss(line);
				boost::trim(line);
				if(strncmp(line.c_str(),"}",1)==0)
					break;
				ss >> setting >> value;
				boost::trim(setting);
				boost::trim(value);
				if(setting=="printFromOpen")
					optionsMenu.outputMenu.printFromOpen=value;
				if(setting=="shiftCurve")
					optionsMenu.outputMenu.shiftCurve=value;
			}
		}
	}
}
