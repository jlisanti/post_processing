#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <locale>

#include <gsl/gsl_statistics.h>

#include "datafileinput.h"
#include "readinputfiles.h"
#include "processfile.h"

#include "ssa.h"
#include "computerms.h"
#include "fft.h"
#include "averagecycle.h"
#include "scatterplot.h"
#include "smoothencoder.h"
#include "curvefit.h"
#include "findroots.h"
#include "spectrogram.h"
#include "trackpeaks.h"
#include "computemassflowrate.h"
#include "computeballvalvearea.h"


void process_file(std::string file, output_data &output, options &optionsMenu)
{

	output.file_name = file;
	int skipLines = 0;
	std::cout << "reading header file" << std::endl;
	read_header_file_data(file,output,skipLines);
	std::cout << "Read data for processing" << std::endl;
	DataFileInput cInput(file,skipLines);

	if(optionsMenu.dataCondMenu.SSA=="true")
		ssa(cInput,
			output.time_ssa,
			output.pressure_ssa,
			output.M,
			output.kmax,
			1);

	if((optionsMenu.dataCondMenu.SSA=="true") &&
	   (optionsMenu.mainMenu.ionProbe=="true"))
		ssa(cInput,
			output.time_ssa,
			output.ion_ssa,
			output.M,
			output.kmax,
			3);

	/* process passive valve data */
	if(optionsMenu.mainMenu.combustorType=="passive")
	{
		output.p_rms = compute_prms(cInput,output.p_static);

		compute_frequency_spectrum(cInput,
								   output.spectrum_magnitude,
                                   output.spectrum_frequency,
								   1,
                                   output.combustor_frequency);

		compute_frequency_spectrum(cInput,
								   output.ion_spectrum_magnitude,
                                   output.ion_spectrum_frequency,
								   2,
                                   output.ion_frequency);

		average_cycle_passive(cInput,
							  1,
							  2,
							  optionsMenu.dataCondMenu.windowSize,
							  optionsMenu.dataCondMenu.stepSize,
							  output.phase_vector_smooth,
							  output.time_vector_smooth,
							  output.pressure_vector_smooth,
							  output.ion_vector_smooth,
							  output.phase_scatter,
							  output.time_scatter,
							  output.pressure_scatter,
							  output.ion_scatter,
							  output.downCrossOver,
							  output.upCrossOver,
							  output.combustor_frequency,
							  output.static_pressure,
							  optionsMenu.outputMenu.printFromOpen);

		std::cout << "successfully averaged file" << std::endl;

		find_min(output.time_vector_smooth,
				 output.pressure_vector_smooth,
				 output.pressure_min,
				 output.phase_min);

		find_max(output.time_vector_smooth,
				 output.pressure_vector_smooth,
				 output.pressure_max,
				 output.phase_max);

		if(optionsMenu.mainMenu.ionProbe=="true")
			find_min(output.time_vector_smooth,
					 output.ion_vector_smooth,
					 output.ion_min,
					 output.ion_phase);

		std::cout << "successfully found min and max" << std::endl;

		if(optionsMenu.dataAnalysisMenu.spectrogram=="true")
		{
			build_spectrogram(cInput,
							  output.combustor_frequency,
							  output.spectrogram_f,
						      output.spectrogram_m,
						      optionsMenu.dataAnalysisMenu.spectNcyclesSpace,
						      optionsMenu.dataAnalysisMenu.spectNcyclesWindow,
						      output.spectrogram_n,
							  1,
						      output.time);

			build_spectrogram(cInput,
							  output.combustor_frequency,
							  output.ion_spectrogram_f,
						      output.ion_spectrogram_m,
						      optionsMenu.dataAnalysisMenu.spectNcyclesSpace,
						      optionsMenu.dataAnalysisMenu.spectNcyclesWindow,
						      output.spectrogram_n,
							  2,
						      output.ion_time);
		}
		if(optionsMenu.dataAnalysisMenu.trackPeaks=="true")
			track_peaks(cInput,
						output.pressure_vector,
						output.ion_vector,
						1,
						2,
						0,
						output.t_max_p,
						output.t_min_p,
						output.t_max_i,
						output.t_min_i,
						output.p_max,
					    output.p_min,
						output.i_max,
						output.i_min,
						output.downCrossOver,
						output.upCrossOver);

	}
	/* process active valve data */
	else if(optionsMenu.mainMenu.combustorType=="active")
	{
		std::cout << " Processing active valve data"  << std::endl;

		valveProperties valveGeom; 
		valveGeom.spanCount = optionsMenu.dataAnalysisMenu.spanCount;
		valveGeom.inletRadius = optionsMenu.dataAnalysisMenu.inletRadius;
		valveGeom.ballInternalRadius = optionsMenu.dataAnalysisMenu.ballInternalRadius;
		valveGeom.h = optionsMenu.dataAnalysisMenu.lengthBall;

		//std::vector<double> areaFunction;
		std::cout << "computing valve area" << std::endl;
		computeBallValveAreaFunction(valveGeom,
				                     1000,
									 output.areaFunction);
		std::cout << "computed valve area" << std::endl;

		output.p_rms = compute_prms(cInput,output.p_static);

		compute_frequency_spectrum(cInput,
								   output.spectrum_magnitude,
                                   output.spectrum_frequency,
								   output.pressure_colmn,
                                   output.combustor_frequency);

		/*
		compute_frequency_spectrum(cInput,
								   output.ion_spectrum_magnitude,
                                   output.ion_spectrum_frequency,
								   2,
                                   output.ion_frequency);
								   */
		/* need to determine pressure and econder coloumns dynamically */

		std::vector<double> tmp1 = cInput.table_column(0);
		std::vector<double> tmp2 = cInput.table_column(output.pressure_colmn);
		std::vector<double> tmp3 = cInput.table_column(output.encoder_colmn);

		std::cout << "Encoder colmn: " << output.encoder_colmn << '\t'
			      << "Pressur colmn: " << output.pressure_colmn << std::endl;

		average_cycle_active(cInput,
							  optionsMenu.dataCondMenu.windowSize,
							  optionsMenu.dataCondMenu.stepSize,
							  tmp1,
							  tmp2,
							  tmp3,
							  output.encoder_colmn,
							  output.pressure_colmn,
							  output.encoder_vector_smooth,
							  output.pressure_vector_smooth,
							  output.ion_vector_smooth,
							  output.pressure_scatter,
							  output.encoder_scatter,
							  output.ion_scatter,
							  output.upCrossOver,
							  output.downCrossOver,
							  output.areaFunction,
							  optionsMenu.outputMenu.printFromOpen);

		std::cout << std::endl;
		std::cout << "finished averaging cycle" << std::endl;
		std::cout << std::endl;


		find_min(output.encoder_vector_smooth,
				 output.pressure_vector_smooth,
				 output.pressure_min,
				 output.phase_min);

		find_max(output.encoder_vector_smooth,
				 output.pressure_vector_smooth,
				 output.pressure_max,
				 output.phase_max);

		if(optionsMenu.mainMenu.ionProbe=="true")
			find_min(output.time_vector_smooth,
					 output.ion_vector_smooth,
					 output.ion_min,
					 output.ion_phase);

		if(optionsMenu.dataAnalysisMenu.spectrogram=="true")
		{
			build_spectrogram(cInput,
							  output.combustor_frequency,
							  output.spectrogram_f,
						      output.spectrogram_m,
						      optionsMenu.dataAnalysisMenu.spectNcyclesSpace,
						      optionsMenu.dataAnalysisMenu.spectNcyclesWindow,
						      output.spectrogram_n,
							  1,
						      output.time);

			build_spectrogram(cInput,
							  output.combustor_frequency,
							  output.ion_spectrogram_f,
						      output.ion_spectrogram_m,
						      optionsMenu.dataAnalysisMenu.spectNcyclesSpace,
						      optionsMenu.dataAnalysisMenu.spectNcyclesWindow,
						      output.spectrogram_n,
							  2,
						      output.ion_time);
		}
		if(optionsMenu.dataAnalysisMenu.trackPeaks=="true")
			track_peaks(cInput,
						output.pressure_vector,
						output.ion_vector,
						1,
						2,
						0,
						output.t_max_p,
						output.t_min_p,
						output.t_max_i,
						output.t_min_i,
						output.p_max,
					    output.p_min,
						output.i_max,
						output.i_min,
						output.downCrossOver,
						output.upCrossOver);

        double Tin = 300.0;
        double Pstatic = 1.00944683*1.0e5;
        double R = 287.15;
        double deltaL = 66.6/1000.0; //66.0/1000.0

        //double C = 0.9;
        //double eta = 0.079; 

		double massFlowAverage = 0;
		std::vector<double> massFlow;
		std::vector<double> massFlowTime;

		computeMassFlowRate(optionsMenu.dataAnalysisMenu.inletTemperature,
				            optionsMenu.dataAnalysisMenu.staticPressure,
							optionsMenu.dataAnalysisMenu.gasConstant,
							output.areaFunction,
							output.combustor_frequency,
							optionsMenu.dataAnalysisMenu.inletLength,
							optionsMenu.dataAnalysisMenu.C,
							optionsMenu.dataAnalysisMenu.eta,
							output.massFlowRateAir,
							massFlow,
							massFlowTime,
							output.pressure_vector_smooth,
							output.encoder_vector_smooth);

        if (output.mass_flow_rate_fuel_1 > 10.0)
            output.mass_flow_rate_fuel_1 = output.mass_flow_rate_fuel_1/1000.0;

        if (output.mass_flow_rate_fuel_1 > 10.0)
            output.mass_flow_rate_fuel_1 = output.mass_flow_rate_fuel_1/1000.0;
	
		output.equivalenceRatio = (output.mass_flow_rate_fuel_1
				                   /output.massFlowRateAir)/(1.0/3.0984);

		std::cout << "air mass flow average: " << output.massFlowRateAir << '\t'
			      << "fuel mass flow rate: "   << output.mass_flow_rate_fuel_1 << '\t'
				  << "equilance ratio: "       << output.equivalenceRatio << std::endl;
		std::ofstream foutTest("testMFR.dat");
		for (int i = 0; i < output.pressure_vector_smooth.size(); i++)
		{
		    foutTest << output.encoder_vector_smooth[i]  << '\t'
				     << output.pressure_vector_smooth[i] << '\t'
					 << massFlow[i] << std::endl;
		}

	}
	else
	{
	}

	std::cout << "Exiting proccessfile..." << std::endl;
}

static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

void read_header_file_data(std::string file, output_data &output, int &skip_lines)
{

	int encoder_colmn = 0;

	std::ifstream input_file_stream(file);

	while (true)
	{
		std::string line;
		std::getline (input_file_stream, line);

		bool exit_loop = false;

		if (!input_file_stream)
			break;

		if ((strncmp(line.c_str(),"Mass Flow Rate Air",18)==0 ||
            (strncmp(line.c_str(),"Air mass flow rate",18)==0)))
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+2);
			output.mass_flow_rate_air = std::stod(tmp);
		}
		else if ((strncmp(line.c_str(),"Phase",5)==0))
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+2);
			output.fuel_injection_phase = std::stod(tmp);
		}
		else if ((strncmp(line.c_str(),"Pulse Width",11)==0))
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+2);
			output.fuel_injection_pulse_width = std::stod(tmp);
		}
		else if ((strncmp(line.c_str(),"Mass Flow Rate Fuel",19)==0) ||
				 (strncmp(line.c_str(),"Fuel mass flow rate",19)==0) ||
				 (strncmp(line.c_str(),"CH4  mass flow rate", 18)==0) ||
				 (strncmp(line.c_str(),"Fuel mass flow rate   [g/s]:",28)==0) ||
				 (strncmp(line.c_str(),"C2H4 flow rate:", 15)==0)) 
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+2);
			output.mass_flow_rate_fuel_1 = std::stod(tmp);
		}
		else if ((strncmp(line.c_str(),"H2 mass flow rate     [g/s]:", 18)==0))
		{

			size_t position = line.find(":");
			std::string tmp = line.substr(position+2);
			output.mass_flow_rate_fuel_2 = std::stod(tmp);
		}
		else if ((strncmp(line.c_str(),"Valve set point",15)==0) ||
		         (strncmp(line.c_str(),"Valve frequency:",15)==0) ||
		         (strncmp(line.c_str(),"Frequency:",10)==0))
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+2);
			output.combustor_frequency = std::stod(tmp);
		}
		else if (strncmp(line.c_str(),"Fuel pressure:",14)==0)
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+2);
			output.fuel_pressure = std::stod(tmp);
		}
		else if (strncmp(line.c_str(),"Static pressure:",16)==0)
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+2);
			output.static_pressure = std::stod(tmp);
		}
		else if (strncmp(line.c_str(),"Total pressure:",14)==0)
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+2);
			output.total_pressure = std::stod(tmp);
		}
		else if (strncmp(line.c_str(),"Total temperature:",16)==0)
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+2);
			output.total_temperature = std::stod(tmp);
		}
		else if (strncmp(line.c_str(),"UHC:",4)==0)
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+3);
			output.UHC = std::stod(tmp);
		}
		else if (strncmp(line.c_str(),"NOx:",4)==0)
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+3);
			output.NOx = std::stod(tmp);
		}
		else if (strncmp(line.c_str(),"CO:",3)==0)
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+3);
			output.CO = std::stod(tmp);
		}
		else if (strncmp(line.c_str(),"CO2:",4)==0)
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+3);
			output.CO2 = std::stod(tmp);
		}
		else if (strncmp(line.c_str(),"O2:",3)==0)
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+3);
			output.O2 = std::stod(tmp);
		}
        else if (strncmp(line.c_str(),"Data collected on",17)==0)
        {
		    std::string dateInfo = line;
			std::string timeInfo (" at ");
			rtrim(dateInfo);
		    if(dateInfo.find(timeInfo) != std::string::npos)
			{
                size_t position = line.find(":");
                std::string time = line.substr(position-2);
                std::istringstream ss(time);
                std::string token;
                std::vector<double> t;

                while(std::getline(ss, token, ':')) {
                    t.push_back(std::stod(token));
                }
 
                double tmp = t[0]*60*60 + t[1]*60.0 + t[2];
                output.acquisition_time = tmp;
			}

        }
		else if ((strncmp(line.c_str(),
				"Time [s]            Pressure [bar]",34)==0) ||
				 (strncmp(line.c_str(),
				 "Time [s]                Pressure [bar]",36)==0) ||
				 (strncmp(line.c_str(),
				 "Time [s]        Pressure [bar]  Position [volts]",45)==0) ||
				 (strncmp(line.c_str(),
				 "Time [s]        Pressure [bar]  Ion Probe [volts])",45)==0))
		{

			if ((strncmp(line.c_str(),
					"Time [s]            Pressure [bar]      Injector [volts]",56)==0))
			{
				encoder_colmn = 3;
				output.pressure_colmn = 1;
			}
			if ((strncmp(line.c_str(),
			    "Time [s]        Pressure [bar]  Position [volts]",48)==0))
			{
				encoder_colmn = 3;
				output.pressure_colmn = 1;
			}
			else
			{
				encoder_colmn = 2;
				output.pressure_colmn = 3;
			}
			if ((strncmp(line.c_str(), "Time [s]                Pressure [bar]          Position [volts]",
							60)==0))
			{
				encoder_colmn = 3;
				output.pressure_colmn = 1;
			}
			output.encoder_colmn = encoder_colmn-1;
			exit_loop = true;
		}

		skip_lines++;

		if (exit_loop)
			break;
	}
}
