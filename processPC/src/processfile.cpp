#include <string>
#include <iostream>
#include <fstream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <locale>
#include <vector>

#include <gsl/gsl_statistics.h>

#include "datafileinput.h"
#include "readinputfiles.h"
#include "processfile.h"
#include "menu.h"

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
#include "bodeplot.h"
#include "computeliquidfuelflowrate.h"


void process_file(std::string file, output_data &output, options &optionsMenu, Menu &cMenu)
{

	output.file_name = file;
	int skipLines = 0;
	if(optionsMenu.mainMenu.combustorType!="average")
	{
	    std::cout << "reading header file" << std::endl;
	    read_header_file_data(file,output,skipLines);
	}
	std::cout << "Read data for processing" << std::endl;
	DataFileInput  cInput(file,skipLines);

	/*
	std::ofstream fout1("orig.dat");
	for (int i = 0; i < cInput.size(); i++)
		fout1 << cInput.table_value(i,0) << '\t'
			  << cInput.table_value(i,1) << '\t'
			  << cInput.table_value(i,2) << std::endl;
			  */


	if(optionsMenu.dataCondMenu.SSA=="true")
	{
		ssa(cInput,
			output.time_ssa,
			output.pressure_ssa,
			optionsMenu.dataCondMenu.vectorDimensionM,
			optionsMenu.dataCondMenu.kMax,
			output.pressure_colmn);

	    /* create temporary file and replace cInput */
	    std::ofstream fout("tmp");
		if(optionsMenu.mainMenu.ionProbe=="true")
		{
		    ssa(cInput,
			    output.time_ssa,
			    output.ion_ssa,
			    output.M,
			    output.kmax,
			    output.ion_probe_colmn);
	        for (int i = 0; i < output.time_ssa.size(); i++)
	            fout << output.time_ssa[i]                         << '\t'
			         << output.pressure_ssa[i]                     << '\t'
			         << cInput.table_value(i,output.encoder_colmn) << '\t' 
					 << output.ion_ssa[i]                          << std::endl;
		}
		else
		{

	        for (int i = 0; i < output.time_ssa.size(); i++)
	            fout << output.time_ssa[i]                         << '\t'
			         << output.pressure_ssa[i]                     << '\t'
			         << cInput.table_value(i,output.encoder_colmn) << std::endl;

		}
		cInput.~DataFileInput();
		new(&cInput) DataFileInput("tmp",0);
		std::remove ("tmp");
	}

	/*
	std::ofstream fout2("smooth.dat");

	for (int i = 0; i < cInput.size(); i++)
		fout2 << cInput.table_value(i,0) << '\t'
			  << cInput.table_value(i,1) << '\t' 
			  << cInput.table_value(i,2) << std::endl;
			  */



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
		cMenu.display_message("Processing active valve data");
		cMenu.display_message("Setting valve geometry params...");
		valveProperties valveGeom; 
		valveGeom.spanCount = optionsMenu.dataAnalysisMenu.spanCount;
		valveGeom.inletRadius = optionsMenu.dataAnalysisMenu.inletRadius;
		valveGeom.ballInternalRadius = optionsMenu.dataAnalysisMenu.ballInternalRadius;
		valveGeom.h = optionsMenu.dataAnalysisMenu.lengthBall;

		std::cout << "span " << valveGeom.spanCount << std::endl;
		std::cout << "inlet " << valveGeom.inletRadius << std::endl;
		std::cout << "ball "  << valveGeom.ballInternalRadius << std::endl;
		std::cout << "h "     << valveGeom.h << std::endl;

		cMenu.display_message("Computing valve area function");
		computeBallValveAreaFunction(valveGeom,
				                     1000,
									 output.areaFunction);


	    std::ofstream fout1("area_function.dat");
        for (int i = 0; i < 1000; i++)
        {
		fout1 << i << '\t'
			  << ballValveArea((M_PI/1000.0)*i,output.areaFunction) << std::endl;
        }

		output.ion_probe_colmn = 2;
		
		std::cout << "file length before: " << cInput.file_length() << std::endl;

		int tmpIter = cInput.file_length();

        /* Clear bad data */
		for (int i = 0; i < tmpIter; i++)
		{
			if((cInput.table_value(i,output.pressure_colmn)==0) || 
               (cInput.table_value(i,output.encoder_colmn)==0))
				cInput.delete_row(i);
		}

		std::cout << "file length after: " << cInput.file_length() << std::endl;

		/*
		for (int i = 0; i < output.areaFunction.size(); i++)
		    std::cout << output.areaFunction[i] << std::endl;
			*/



		if(optionsMenu.dataAnalysisMenu.computeRMS=="true")
		{
		     cMenu.display_message("Computing pressure RMS");
		     output.p_rms = compute_prms(cInput,output.p_static);
		}
/*
		if(optionsMenu.dataAnalysisMenu.setMeanP=="true")
		{
		     cMenu.display_message("Computing pressure RMS");
		     output.p_rms = compute_prms(cInput,output.p_static);
		}
*/

		double tmp = 0.0;

		if(optionsMenu.dataAnalysisMenu.freqSpectrum=="true")
		{
		    cMenu.display_message("Computing pressure frequency spectrum");
		    compute_frequency_spectrum(cInput,
								       output.spectrum_magnitude,
                                       output.spectrum_frequency,
								       output.pressure_colmn,
                                       tmp);
		}


		if((optionsMenu.dataAnalysisMenu.freqSpectrum=="true") &&
				(optionsMenu.mainMenu.ionProbe=="true"))
		{
		    cMenu.display_message("Computing ion probe frequency spectrum");
		    compute_frequency_spectrum(cInput,
								       output.ion_spectrum_magnitude,
                                       output.ion_spectrum_frequency,
								       output.ion_probe_colmn,
                                       output.ion_frequency);
		}

		std::vector<double> tmp1 = cInput.table_column(0);
		std::vector<double> tmp2 = cInput.table_column(output.pressure_colmn);
		std::vector<double> tmp3 = cInput.table_column(output.ion_probe_colmn);

        double meanP = compute_pmean(cInput,output.p_static);
        double meanPshift = meanP-optionsMenu.dataAnalysisMenu.meanP;

		meanPshift = 0.0;

        if(optionsMenu.dataAnalysisMenu.averageCycle=="true")
		{
			cMenu.display_message("Averaging cycle...");
			average_cycle_active(cInput,
							  optionsMenu.dataCondMenu.windowSize,
							  optionsMenu.dataCondMenu.stepSize,
							  optionsMenu.dataCondMenu.pressureOffset,
							  meanPshift,
							  tmp1,
							  tmp2,
							  tmp3,
							  output.encoder_colmn,
							  output.pressure_colmn,
							  output.ion_probe_colmn,
							  output.encoder_vector_smooth,
							  output.pressure_vector_smooth,
							  output.ion_vector_smooth,
							  output.pressure_scatter,
							  output.encoder_scatter,
							  output.ion_scatter,
							  output.upCrossOver,
							  output.downCrossOver,
							  output.areaFunction,
							  optionsMenu.dataCondMenu.shiftCurve,
							  optionsMenu.dataCondMenu.shiftFixedCurve,
							  optionsMenu.dataCondMenu.shiftFixed,
							  optionsMenu.outputMenu.printFromOpen);
		}
		else
        {
			output.encoder_vector_smooth.assign(250,0.0);
            output.pressure_vector_smooth.assign(250,0.0);
            output.airMassFlow.assign(250,0.0);
        }
		/*

		for (int i =0; i < output.ion_vector_smooth.size(); i++)
			std::cout << i << '\t' << output.ion_scatter[i] << std::endl;
			*/



		if(optionsMenu.dataAnalysisMenu.findPeaks=="true")
		{
		    cMenu.display_message("Finding pressure peaks");
		    find_min(output.encoder_vector_smooth,
				     output.pressure_vector_smooth,
				     output.pressure_min,
				     output.phase_min);

		    find_max(output.encoder_vector_smooth,
				     output.pressure_vector_smooth,
				     output.pressure_max,
				     output.phase_max);

		    if(optionsMenu.mainMenu.ionProbe=="true")
		    {
		        cMenu.display_message("Finding ion probe peaks");
			    find_min(output.encoder_vector_smooth,
					     output.ion_vector_smooth,
					     output.ion_min,
					     output.ion_phase);
		    }
		}

		cMenu.display_message("Finished finding peaks");

/*
        output.liquidFuelMassFlowRate = computeLiquidFuelFlowRate( output.combustor_frequency,
                                                                   output.fuel_pressure,
                                                                   output.pressure_vector_smooth);
*/

		if(optionsMenu.dataAnalysisMenu.spectrogram=="true")
		{

		    cMenu.display_message("Building pressure spectrogram");
			build_spectrogram(cInput,
							  output.combustor_frequency,
							  output.spectrogram_f,
						      output.spectrogram_m,
						      optionsMenu.dataAnalysisMenu.spectNcyclesSpace,
						      optionsMenu.dataAnalysisMenu.spectNcyclesWindow,
						      output.spectrogram_n,
							  output.pressure_colmn,
						      output.time);

		    if(optionsMenu.mainMenu.ionProbe=="true")
			{
		        cMenu.display_message("Building ion probe spectrogram");
			    build_spectrogram(cInput,
			    				  output.combustor_frequency,
			    				  output.ion_spectrogram_f,
			    			      output.ion_spectrogram_m,
			    			      optionsMenu.dataAnalysisMenu.spectNcyclesSpace,
			    			      optionsMenu.dataAnalysisMenu.spectNcyclesWindow,
			    			      output.spectrogram_n,
			    				  output.encoder_colmn,
			    			      output.ion_time);
			}
		}
		if(optionsMenu.dataAnalysisMenu.trackPeaks=="true")
		{
		    cMenu.display_message("Tracking peaks");
			track_peaks_active(cInput,
						       output.pressure_vector,
						       output.ion_vector,
						       output.pressure_colmn,
						       output.ion_probe_colmn,
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


		if(optionsMenu.dataAnalysisMenu.computeMassFlowRate=="true")
		{
		    cMenu.display_message("Computing mass flow rates");
	    	computeMassFlowRate(optionsMenu.dataAnalysisMenu.inletTemperature,
	    			            optionsMenu.dataAnalysisMenu.staticPressure,
	    						optionsMenu.dataAnalysisMenu.gasConstant,
	    						output.areaFunction,
	    						output.combustor_frequency,
	    						optionsMenu.dataAnalysisMenu.inletLength,
	    						optionsMenu.dataAnalysisMenu.C,
	    						optionsMenu.dataAnalysisMenu.eta,
	    						output.massFlowRateAir,
	    						output.airMassFlow,
	    						output.airMassFlowTime,
	    						output.pressure_vector_smooth,
	    						output.encoder_vector_smooth);

            if (output.mass_flow_rate_fuel_1 > 10.0)
                output.mass_flow_rate_fuel_1 = output.mass_flow_rate_fuel_1/1000.0;

            if (output.mass_flow_rate_fuel_1 > 0.01)
                output.mass_flow_rate_fuel_1 = output.mass_flow_rate_fuel_1/1000.0;

	
		    output.equivalenceRatio = (output.mass_flow_rate_fuel_1
				                   /output.massFlowRateAir)/(1.0/3.0984);

		    std::cout << std::endl;
		    std::cout << "*** Equivalence Ratio ***" << std::endl;
		    std::cout << "    Mdot_fuel = " << output.mass_flow_rate_fuel_1 << std::endl;
		    std::cout << "    Mdot_air  = " << output.massFlowRateAir << std::endl;
		    std::cout << "    Phi       = " << output.equivalenceRatio << std::endl;
		    std::cout << std::endl;

		    std::cout << "encoder size: " << output.encoder_vector_smooth.size() << '\t'
			          << "air size: "     << output.airMassFlow.size() << std::endl;

		}

		if(optionsMenu.dataAnalysisMenu.bodePlots=="true")
		{
		    cMenu.display_message("Building bode plot");
	    	int check = buildBodePlot (cInput,
	    			                   0,
	    							   output.pressure_colmn,
	    							   output.bodeP,
	    							   output.bodePprime);
		}
				                    


	}
	/* process active valve data */
	else if(optionsMenu.mainMenu.combustorType=="average")
	{
		cMenu.display_message("Processing active valve data");
		cMenu.display_message("Setting valve geometry params...");
		valveProperties valveGeom; 
		valveGeom.spanCount = optionsMenu.dataAnalysisMenu.spanCount;
		valveGeom.inletRadius = optionsMenu.dataAnalysisMenu.inletRadius;
		valveGeom.ballInternalRadius = optionsMenu.dataAnalysisMenu.ballInternalRadius;
		valveGeom.h = optionsMenu.dataAnalysisMenu.lengthBall;

		cMenu.display_message("Computing valve area function");
		computeBallValveAreaFunction(valveGeom,
				                     1000,
									 output.areaFunction);

		if(optionsMenu.dataAnalysisMenu.findPeaks=="true")
		{
		    cMenu.display_message("Finding pressure peaks");
		    find_min(cInput.table_column(3),
				     cInput.table_column(2),
				     output.pressure_min,
				     output.phase_min);

		    find_max(cInput.table_column(3),
				     cInput.table_column(2),
				     output.pressure_max,
				     output.phase_max);

		    if(optionsMenu.mainMenu.ionProbe=="true")
		    {
		        cMenu.display_message("Finding ion probe peaks");
			    find_min(output.time_vector_smooth,
					     output.ion_vector_smooth,
					     output.ion_min,
					     output.ion_phase);
		    }
		}

		if(optionsMenu.dataAnalysisMenu.computeMassFlowRate=="true")
		{
		    cMenu.display_message("Computing mass flow rates");
	    	computeMassFlowRate(optionsMenu.dataAnalysisMenu.inletTemperature,
	    			            optionsMenu.dataAnalysisMenu.staticPressure,
	    						optionsMenu.dataAnalysisMenu.gasConstant,
	    						output.areaFunction,
	    						output.combustor_frequency,
	    						optionsMenu.dataAnalysisMenu.inletLength,
	    						optionsMenu.dataAnalysisMenu.C,
	    						optionsMenu.dataAnalysisMenu.eta,
	    						output.massFlowRateAir,
	    						output.airMassFlow,
	    						output.airMassFlowTime,
	    						output.pressure_vector_smooth,
	    						output.encoder_vector_smooth);

            if (output.mass_flow_rate_fuel_1 > 10.0)
                output.mass_flow_rate_fuel_1 = output.mass_flow_rate_fuel_1/1000.0;

            if (output.mass_flow_rate_fuel_1 > 0.01)
                output.mass_flow_rate_fuel_1 = output.mass_flow_rate_fuel_1/1000.0;

	
		    output.equivalenceRatio = (output.mass_flow_rate_fuel_1
				                   /output.massFlowRateAir)/(1.0/3.0984);

		    std::cout << std::endl;
		    std::cout << "*** Equivalence Ratio ***" << std::endl;
		    std::cout << "    Mdot_fuel = " << output.mass_flow_rate_fuel_1 << std::endl;
		    std::cout << "    Mdot_air  = " << output.massFlowRateAir << std::endl;
		    std::cout << "    Phi       = " << output.equivalenceRatio << std::endl;
		    std::cout << std::endl;

		    std::cout << "encoder size: " << output.encoder_vector_smooth.size() << '\t'
			          << "air size: "     << output.airMassFlow.size() << std::endl;

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
		else if (strncmp(line.c_str(),"Total pressure:",15)==0)
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
		else if (strncmp(line.c_str(),"Gas probe position:",18)==0)
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+3);
			output.gasProbePosition = std::stod(tmp);
		}
		else if (strncmp(line.c_str(),"Air On:",6)==0)
		{
			size_t position = line.find(":");
			std::string tmp = line.substr(position+3);
			output.airOn = std::stod(tmp);
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
				 "Time [s]        Pressure [bar]  Ion Probe [volts]",45)==0 ||
				 (strncmp(line.c_str(),
				 "Time [s]        Pressure [bar]  Position [volts]",45)==0)))
		{

			if ((strncmp(line.c_str(),
					"Time [s]        Pressure [bar]  Ion Probe [volts] encoder [volts]",56)==0))
			{
				encoder_colmn = 4;
				output.pressure_colmn = 1;
			}
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
			if ((strncmp(line.c_str(),
			"Time [s]        Pressure [bar]  Ion Probe [volts] encoder [volts]",64)==0))
			{
				encoder_colmn = 4;
                output.pressure_colmn = 1;
			}
			if ((strncmp(line.c_str(),
			"Time [s]        Pressure [bar]  Position [volts]",40)==0))
			{
				encoder_colmn = 3;
                output.pressure_colmn = 1;
			}
			if ((strncmp(line.c_str(), "Time [s]                Pressure [bar]          Position [volts]",
							60)==0))
			{
				encoder_colmn = 3;
				output.pressure_colmn = 1;
			}
/*
			else
			{
				encoder_colmn = 2;
				output.pressure_colmn = 3;
			}
*/
			output.encoder_colmn = encoder_colmn-1;
			exit_loop = true;
		}

		skip_lines++;

		if (exit_loop)
			break;
	}
}
