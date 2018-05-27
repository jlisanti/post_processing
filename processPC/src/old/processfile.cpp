#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include  <cmath>
#include <cstring>

#include <gsl/gsl_statistics.h>

#include "datafileinput.h"
#include "processfile.h"
#include "fft.h"
#include "averagecycle.h"
#include "curvefit.h"
#include "findroots.h"
#include "scatterplot.h"
#include "computerms.h"
#include "computemass.h"
#include "spectrogram.h"
#include "ssa.h"
#include "trackpeaks.h"

void process_file(std::string file, output_data &output)
{
	int skip_lines = 0;
	int encoder_colmn = 0;

	std::ifstream input_file_stream(file);
	
	//output.file_name = file;
	//std::cout << file << std::endl;

	std::cout << " Reading input data" << std::endl;

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
		         (strncmp(line.c_str(),"Valve frequency:",15)==0))
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
			/*
			size_t position = line.find(":");
			std::string tmp = line.substr(position+2);
			output.fuel_pressure = std::stod(tmp);
			*/
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
        else if (strncmp(line.c_str(),"Data collected on",17)==0)
        {
			/*
            size_t position = line.find(":");
            std::string time = line.substr(position-2);

            std::istringstream ss(time);
            std::string token;

            std::vector<double> t;

            while(std::getline(ss, token, ':')) {
                std::cout << token << std::endl;
                t.push_back(std::stod(token));
            }

            std::cout << "time " << t[0] << ":" << t[1] << ":" << t[2] << std::endl;
            double tmp = t[0]*60*60 + t[1]*60.0 + t[2];

            output.acquisition_time = tmp;
			*/
        }
		else if ((strncmp(line.c_str(),
				"Time [s]            Pressure [bar]",34)==0) ||
				 (strncmp(line.c_str(),
				 "Time [s]                Pressure [bar]",36)==0) ||
				 (strncmp(line.c_str(),
				 "Time [s]        Pressure [bar]  Position [volts]",45)==0))
		{

			if ((strncmp(line.c_str(),
					"Time [s]            Pressure [bar]      Injector [volts]",56)==0))
			{
				encoder_colmn = 3;
			}
			else
				encoder_colmn = 2;
			exit_loop = true;
		}

		skip_lines++;

		if (exit_loop)
			break;
	}

	std::cout << "Reading data file..." << std::endl;
    std::cout << file << std::endl;
	DataFileInput cInput(file,skip_lines);
    std::cout << "file length " << cInput.file_length() << std::endl;
	double dumb;

	//compute_prms(cInput,output.p_rms,1.01325);

	int nn = cInput.file_length();

	std::cout << "Building array from data" << std::endl;
	double * data = new double[nn];

	for (int i = 0; i < nn; i++)
		data[i] = cInput.table_value(i,1);


	std::cout << "Compute RMS" << std::endl;
	output.p_rms = sqrt(gsl_stats_variance(data, 1, nn));
	std::cout << output.p_rms << std::endl;

	delete[] data;

	std::cout << "Compute frequency spectrum" << std::endl;
	compute_frequency_spectrum(cInput,
			                   output.spectrum_magnitude,
							   output.spectrum_frequency,
							   dumb);
							   //output.combustor_frequency);
	std::cout << output.combustor_frequency << std::endl;

	ssa(cInput,
		output.time_ssa,
		output.pressure_ssa,
		output.M,
		output.kmax,
		1);

	ssa(cInput,
		output.time_ssa,
		output.ion_ssa,
		output.M,
		output.kmax,
		3);

	std::cout << "Averaging cycle..." << std::endl;
	std::cout << "endcoder: " << encoder_colmn << std::endl;
	output.encoder_colmn = encoder_colmn;
	average_cycle(cInput,
			      output.time_ssa,
				  output.pressure_ssa,
				  output.ion_ssa,
			      encoder_colmn,
				  1,
				  output.window,
				  output.encoder_vector_smooth,
				  output.pressure_vector_smooth,
				  output.ion_vector_smooth,
				  output.pressure_scatter,
				  output.encoder_scatter,
				  output.ion_scatter);
	/*
		*/
	std::cout << encoder_colmn << std::endl;
	std::cout << "Fitting poly..." << std::endl;

	
	find_min(output.encoder_vector_smooth,
			 output.pressure_vector_smooth,
			 output.pressure_min,
			 output.phase_min);

	find_max(output.encoder_vector_smooth,
			 output.pressure_vector_smooth,
			 output.pressure_max,
			 output.phase_max);

	find_min(output.encoder_vector_smooth,
			 output.ion_vector_smooth,
			 output.ion_min,
			 output.ion_phase);

	build_spectrogram(cInput,
			          output.combustor_frequency,
					  output.spectrogram_f,
					  output.spectrogram_m,
					  output.spectrogram_n,
					  output.dt,
					  output.time);
	/*
	track_peaks(cInput,
			    output.pressure_vector_smooth,
				output.encoder_vector_smooth,
				output.ion_vector_smooth,
			    encoder_colmn,
				1,
				3,
				0,
				output.e_max_p,
				output.e_min_p,
				output.e_max_i,
				output.e_min_i,
				output.t_max_p,
				output.t_min_p,
				output.t_max_i,
				output.t_min_i,
				output.p_max,
				output.p_min,
				output.i_max,
				output.i_min);
				*/


	/*
	for (int i = 0; i < output.encoder_vector_smooth.size(); i++)
		std::cout << output.encoder_vector_smooth[i] << std::endl;
		*/

	std::cout << "Exiting proccessfile..." << std::endl;


	/*
	compute_mass_flow_rate(output.pressure_vector_smooth,
			               output.encoder_vector_smooth,
						   output.mass_flow_rate_fuel_1,
						   output.mass_flow_rate_fuel_2,
						   19.6/2000.0,
						   23.0/2000.0,
						   40.0/1000.0,
						   output.combustor_frequency,
						   output.mdot_air,
						   output.phi,
						   output.AF_ratio,
						   output.mdot_inlet,
						   output.area_inlet,
						   output.area_inlet_2);
						   */
	/*
	curve_fit(output.encoder_vector_smooth, 
			 output.pressure_vector_smooth,
    		 output.coefficients);
			 */
}
