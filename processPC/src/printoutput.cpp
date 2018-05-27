#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream> 

#include "processfile.h"
#include "datafileinput.h"
#include "readinputfiles.h"
#include "printoutput.h"
#include "averagecycle.h"

inline bool file_exists (const std::string& name) {
	std::ifstream f(name.c_str());
	return f.good();
}

int round_value(int number, int multiple)
{
	int remainder = number % multiple;
	if (remainder == 0)
		return number;

	int upper_value = std::round(number + multiple - remainder);
	int lower_value = std::round(number - remainder);

	if (std::abs(upper_value-number) < std::abs(lower_value-number))
		return upper_value;
	else
		return lower_value;
}

void print_output_active (output_data &output)
{

	output.combustor_frequency = round_value(output.combustor_frequency,5);
	std::string file_name;
	std::ostringstream convert;

	std::ostringstream convert1;
	std::ostringstream convert2;
	std::ostringstream convert3;
	std::ostringstream convert4;

	if (output.mass_flow_rate_fuel_1 < 10.0)
		output.mass_flow_rate_fuel_1 = output.mass_flow_rate_fuel_1*1000.0;

	convert << int(output.combustor_frequency);

	convert1 << int(output.fuel_injection_phase);
	convert2 << int(output.fuel_injection_pulse_width);
	convert3 << int(output.mass_flow_rate_fuel_1*1.0);
	convert4 << int(output.mass_flow_rate_fuel_2*1.0);

	std::cout << "printing output..." << std::endl;

	/*
	file_name = convert.str() + "Hz_" + convert3.str() + "mgpsCH4_" + 
		        convert4.str() + "mgpsH2_p" +  
		        convert1.str() + "_pw" + convert2.str() + "_out.dat";
	std::cout << file_exists(file_name) << std::endl;

    std::string averaged_file = convert.str() + "Hz_" + convert3.str() + "mgpsCH4_" + 
                                convert4.str() + "mgpsH2_p" +  
                                convert1.str() + "_pw" + convert2.str() + "_scatter.dat";

    std::string averaged_curve = convert.str() + "Hz_" + convert3.str() + "mgpsCH4_" + 
                                convert4.str() + "mgpsH2_p" +  
                                convert1.str() + "_pw" + convert2.str() + "_averaged.dat";
								*/

	file_name = convert.str() + "Hz_" + convert3.str() + "mgpsC2H4_out.dat";
	std::cout << file_exists(file_name) << std::endl;

    std::string averaged_file = convert.str() + "Hz_" + convert3.str() + "mgpsC2H4_scatter.dat";
    std::string averaged_curve = convert.str() + "Hz_" + convert3.str() + "mgpsC2H4_averaged.dat";

    std::string spectrum_file = convert.str() + "Hz_" + convert3.str() + "mgpsC2H4_spectrum.dat";
    std::string spectrogram_file = convert.str() + "Hz_" + convert3.str() + "mgpsC2H4_spectrogram.dat";

    std::string scatter_file = convert.str() + "Hz_" + convert3.str() + "mgpsC2H4_pressure_scatter.dat";
    std::string p_max_file = convert.str() + "Hz_" + convert3.str() + "mgpsC2H4_p_max_peaks.dat";
    std::string p_min_file = convert.str() + "Hz_" + convert3.str() + "mgpsC2H4_p_min_peaks.dat";

	for (int i = 0; i < 100; i++)
	{

		std::ostringstream convert2;
		convert2 << (i+1);
		if(file_exists(file_name))
		{
			// Add to averaged file
			// Check if averaged file exist
			// Combine
		    if(!file_exists(averaged_file ))
			{
				// read existing file
				int skip_lines = 0;
				DataFileInput cInput(file_name,skip_lines);

				//  create file
				std::ofstream fout_avg(averaged_file);
				
				int count = 0;

				for (int j = 0; j < cInput.file_length(); j++)
				{
					for (int k = count; k < output.encoder_vector_smooth.size(); k++)
					{
						if (cInput.table_value(j,0) > output.encoder_vector_smooth[k])
						{
							fout_avg << output.encoder_vector_smooth[k] << '\t'
								     << output.pressure_vector_smooth[k] << '\t'
								     << output.ion_vector_smooth[k] << '\t'
								     << std::endl;
							count++;
						}
						else
						{
							fout_avg << cInput.table_value(j,0) << '\t' 
								     << cInput.table_value(j,1) << '\t'
								     << cInput.table_value(j,2) << '\t'
								     << cInput.table_value(j,3) << std::endl; 
							break;
						}
					}
				}
			}
			else
			{
				// read existing file
				int skip_lines = 0;
				DataFileInput cInput(averaged_file,skip_lines);

				//  create file
				std::ofstream fout_avg(averaged_file);

				int count = 0;

				for (int j = 0; j < cInput.file_length(); j++)
				{
					for (int k = count; k < output.encoder_vector_smooth.size(); k++)
					{
						if (cInput.table_value(j,0) > output.encoder_vector_smooth[k])
						{
							fout_avg << output.encoder_vector_smooth[k] << '\t'
								     << output.pressure_vector_smooth[k] << '\t'
								     << output.ion_vector_smooth[k] << '\t'
								     << std::endl;
							count++;
						}
						else
						{
							fout_avg << cInput.table_value(j,0) << '\t' 
								     << cInput.table_value(j,1) << '\t'
								     << cInput.table_value(j,2) << '\t'
								     << cInput.table_value(j,3) << std::endl; 
							break;
						}
					}
				}
				DataFileInput cInput_avg(averaged_file,skip_lines);
				std::vector<double> p_vec;
				std::vector<double> e_vec;
				std::vector<double> i_vec;
				std::vector<double> p_scat;
				std::vector<double> e_scat;
				std::vector<double> i_scat;

				average_cycle_active(cInput_avg,
						      output.time_ssa,
							  output.pressure_ssa,
							  output.ion_ssa,
						      output.encoder_colmn,
							  1,
							  output.window,
							  p_vec,
							  e_vec,
							  i_vec,
							  p_scat,
							  e_scat,
							  i_scat);

				std::ofstream fout_curv(averaged_curve);

				for (int i = 0; i < output.encoder_vector_smooth.size(); i++)
				{
					fout_curv << output.encoder_vector_smooth[i] << '\t'
						      << output.pressure_vector_smooth[i] << '\t'
							  << output.ion_vector_smooth[i] << '\t'
					          << std::endl;
				}
			}
			file_name = convert.str() + "Hz_" + convert3.str() + "mgpsC2H4_out_" + convert2.str() + ".dat";
			/*
			file_name = convert.str() + "Hz_" + convert3.str() + "mgpsCH4_" + 
						convert4.str() + "mgpsH2_p" +  
						convert1.str() + "_pw" + convert2.str() + "_out_"
						+ convert2.str() + ".dat";
						*/
		}
		else
			break;

	}


	if (!file_exists("combustor_data.dat"))
	{
		std::ofstream fout3("combustor_data.dat");
		fout3 //<< std::setw(10) << "mdot_air" << '\t'
			  << std::setw(12) << "mdot_C2H4" << '\t'
			  //<< std::setw(10) << "mdot_fuel_2" << '\t'
			  << std::setw(12) << "frequency" << '\t'
			  << std::setw(12) << "max_P" << '\t'
			  << std::setw(12) << "min_P" << '\t'
			  << std::setw(12) << "max_phase" << '\t'
			  << std::setw(12) << "min_phase" << '\t'
			  //<< std::setw(10) << "phase" << '\t'
			  //<< std::setw(10) << "pulse_width" << '\t'
			  << std::setw(12) << "Fuel_pres" << '\t'
			  << std::setw(12) << "Po" << '\t'
			  << std::setw(12) << "To" << '\t'
			  << std::setw(12) << "P_rms" << '\t'
			  << std::setw(12) << "ion_phase" << '\t'
			  << std::setw(12) << "ion_min" << '\t'
			  << std::setw(12) << "UHC" << '\t'
			  << std::setw(12) << "NOx" << '\t'
			  << std::setw(12) << "CO" << '\t'
			  << std::setw(12) << "CO2" << '\t'
			  << std::setw(12) << "O2" << std::endl;  

	}

	
	std::ofstream fout1(file_name);
	for (int i = 0; i < output.encoder_vector_smooth.size(); i++)
	{
		fout1 << output.encoder_vector_smooth[i] << '\t'
			  << output.pressure_vector_smooth[i] << '\t'
			  << output.ion_vector_smooth[i] << '\t'
			  << std::endl;
	}

	std::ofstream fout2("combustor_data.dat", std::ios_base::app);
	fout2 //<< std::setw(12) << output.mass_flow_rate_air << '\t'
		  << std::setw(12) << output.mass_flow_rate_fuel_1 << '\t'
		  //<< std::setw(12) << output.mass_flow_rate_fuel_2 << '\t'
		  << std::setw(12) << output.combustor_frequency << '\t'
		  << std::setw(12) << output.pressure_max << '\t'
		  << std::setw(12) << output.pressure_min << '\t'
		  << std::setw(12) << output.phase_max << '\t'
		  << std::setw(12) << output.phase_min << '\t'
		  //<< std::setw(12) << output.fuel_injection_phase << '\t'
		  //<< std::setw(12) << output.fuel_injection_pulse_width << '\t' 
		  << std::setw(12) << output.fuel_pressure << '\t' 
		  << std::setw(12) << output.total_pressure << '\t' 
		  << std::setw(12) << output.total_temperature << '\t'
		  << std::setw(12) << output.p_rms << '\t' 
		  << std::setw(12) << output.ion_phase << '\t'
		  << std::setw(12) << output.ion_min << '\t'
		  << std::setw(12) << output.UHC << '\t'
		  << std::setw(12) << output.NOx << '\t'
		  << std::setw(12) << output.CO << '\t'
		  << std::setw(12) << output.CO2 << '\t'
		  << std::setw(12) << output.O2 << std::endl;

	std::ofstream fout3(scatter_file, std::ios_base::app);
	for (int i = 0; i < output.encoder_scatter.size(); i++)
		fout3 << output.encoder_scatter[i] << '\t'
			  << output.pressure_scatter[i] << '\t'
			  << output.ion_scatter[i] << std::endl;


	std::ofstream fout4(spectrum_file);
	for (int i = 0; i < output.spectrum_magnitude.size(); i++)
		fout4 << output.spectrum_frequency[i] << '\t'
			  << output.spectrum_magnitude[i] << std::endl;
	std::cout << "about to print spetrogram" << std::endl;

	std::ofstream fout5(spectrogram_file);
	for (int i = 0; i < output.time.size(); i++)
	{
		for (int j = 0; j < output.spectrogram_n; j++)
		{
			fout5 << output.time[i] << '\t'
				  << output.spectrogram_f[i][j] << '\t'
				  << output.spectrogram_m[i][j] << std::endl;
		}
		fout5 << std::endl;
	}

	std::ofstream fout6(p_max_file);
	fout6 << std::setw(10) << "index"   << '\t'
		  << std::setw(10) << "t_max_p" << '\t' 
		  << std::setw(10) << "e_max_p" << '\t'
		  << std::setw(10) << "p_max"   << std::endl;
	for (int i = 0; i < output.e_max_p.size(); i++)
	{
		fout6 << std::setw(10) << i                 << '\t'
              << std::setw(10) << output.t_max_p[i] << '\t'
			  << std::setw(10) << output.e_max_p[i] << '\t'
			  << std::setw(10) << output.p_max[i]   << std::endl;
	}

	/*
	std::ofstream fout6(peak_file);
	fout6 << std::setw(10) << "index" << '\t'
		  << std::setw(10) << "t_max_p" << '\t'
		  << std::setw(10) << "t_min_p" << '\t'
		  //<< std::setw(10) << "t_max_i" << '\t'
		  //<< std::setw(10) << "t_min_i" << '\t'
		  << std::setw(10) << "e_max_p" << '\t'
		  << std::setw(10) << "e_min_p" << '\t'
		  //<< std::setw(10) << "e_max_i" << '\t'
		  //<< std::setw(10) << "e_min_i" << '\t'
		  << std::setw(10) << "p_max"   << '\t'
		  << std::setw(10) << "p_min"   << std::endl;//'\t'
		  //<< std::setw(10) << "i_max"   << '\t'
		  //<< std::setw(10) << "i_min"   << std::endl;
	for (int i = 0; i < output.e_max_p.size(); i++)
	{
		fout6 << std::setw(10) << i << '\t'
              << std::setw(10) << output.t_max_p[i] << '\t'
			  << std::setw(10) << output.t_min_p[i] << '\t'
			  //<< std::setw(10) << output.t_max_i[i] << '\t'
			  //<< std::setw(10) << output.t_min_i[i] << '\t'

              << std::setw(10) << output.e_max_p[i] << '\t'
			  << std::setw(10) << output.e_min_p[i] << '\t'
			  //<< std::setw(10) << output.e_max_i[i] << '\t'
			  //<< std::setw(10) << output.e_min_i[i] << '\t'

              << std::setw(10) << output.p_max[i] << '\t'
			  << std::setw(10) << output.p_min[i] << std::endl; //'\t'
			  //<< std::setw(10) << output.i_max[i] << '\t'
			  //<< std::setw(10) << output.i_min[i] << std::endl; 
	}

			  */
}

void print_output_passive (output_data &output, options &optionsMenu)
{
	double round_comb_freq = round_value(output.combustor_frequency,5);
	std::string file_name;
	std::ostringstream convert;

	std::ostringstream convert1;
	std::ostringstream convert2;
	std::ostringstream convert3;
	std::ostringstream convert4;

	convert << int(round_comb_freq);
	convert1 << int(output.total_temperature);
	file_name = convert.str() + "Hz_" + convert1.str() + "C.dat";
    std::string averaged_file = convert.str() + "Hz_" + convert1.str() + "C_scatter.dat";
    std::string averaged_curve = convert.str() + "Hz_" + convert1.str() + "C_averaged.dat";

    std::string spectrum_file = convert.str() + "Hz_" + convert1.str() + "C_spectrum.dat";
    std::string spectrogram_file = convert.str() + "Hz_" + convert1.str() + "C_spectrogram.dat";


	for (int i = 0; i < 100; i++)
	{

		std::ostringstream convert2;
		convert2 << (i+1);
		if(file_exists(file_name))
		{

			// Add to averaged file
			// Check if averaged file exist
			// Combine
		    if(!file_exists(averaged_file ))
			{
				// read existing file
				int skip_lines = 0;
				DataFileInput cInput(file_name,skip_lines);

				//  create file
				std::ofstream fout_avg(averaged_file);
				
				int count = 0;

				for (int j = 0; j < cInput.file_length(); j++)
				{
					for (int k = count; k < output.position_vector_smooth.size(); k++)
					{
						if (cInput.table_value(j,0) > output.position_vector_smooth[k])
						{
							fout_avg << output.position_vector_smooth[k] << '\t'
								     << output.pressure_vector_smooth[k] << std::endl;
							count++;
						}
						else
						{
							fout_avg << cInput.table_value(j,0) << '\t' 
								     << cInput.table_value(j,1) << std::endl; 
							break;
						}
					}
				}
			}
			else
			{

				// read existing file
				int skip_lines = 0;
				DataFileInput cInput(averaged_file,skip_lines);

				//  create file
				std::ofstream fout_avg(averaged_file);

				int count = 0;

				for (int j = 0; j < cInput.file_length(); j++)
				{
					for (int k = count; k < output.position_vector_smooth.size(); k++)
					{
						if (cInput.table_value(j,0) > output.position_vector_smooth[k])
						{
							fout_avg << output.position_vector_smooth[k] << '\t'
								     << output.pressure_vector_smooth[k] << std::endl;
							count++;
						}
						else
						{
							fout_avg << cInput.table_value(j,0) << '\t' 
								     << cInput.table_value(j,1) << std::endl; 
							break;
						}
					}
				}
				DataFileInput cInput_avg(averaged_file,skip_lines);
				std::vector<double> p_vec;
				std::vector<double> e_vec;
				std::vector<double> p_scat;
				std::vector<double> e_scat;

				std::ofstream fout_curv(averaged_curve);

				for (int i = 0; i < output.position_vector_smooth.size(); i++)
				{
					fout_curv << output.position_vector_smooth[i] << '\t'
						      << output.pressure_vector_smooth[i] << std::endl;
				}
			}

			file_name = convert.str() + "Hz_" + convert1.str() + "C_" + convert2.str() + ".dat";
			/*
			file_name = convert.str() + "Hz_" + convert3.str() + "mgpsCH4_" + 
						convert4.str() + "mgpsH2_p" +  
						convert1.str() + "_pw" + convert2.str() + "_out_"
						+ convert2.str() + ".dat";
						*/
		}
		else
			break;

	}


	if (!file_exists("combustor_data.dat"))
	{
		std::ofstream fout3("combustor_data.dat");
		fout3 //<< std::setw(10) << "mdot_air" << '\t'
			  << std::setw(12) << "time" << '\t'
			  //<< std::setw(10) << "mdot_fuel_2" << '\t'
			  << std::setw(12) << "frequency" << '\t'
			  << std::setw(12) << "max_P" << '\t'
			  << std::setw(12) << "min_P" << '\t'
			  << std::setw(12) << "max_phase" << '\t'
			  << std::setw(12) << "min_phase" << '\t'
			  //<< std::setw(10) << "phase" << '\t'
			  //<< std::setw(10) << "pulse_width" << '\t'
			  << std::setw(12) << "Fuel_pres" << '\t'
			  << std::setw(12) << "Po" << '\t'
			  << std::setw(12) << "To" << '\t'
			  << std::setw(12) << "P_rms" << '\t' 
			  << std::setw(12) << "UHC" << '\t'
			  << std::setw(12) << "NOx" << '\t'
			  << std::setw(12) << "CO" << '\t'
			  << std::setw(12) << "CO2" << '\t'
			  << std::setw(12) << "O2" << std::endl;  

	}
	

	std::ofstream fout1(file_name);
	for (int i = 0; i < output.position_vector_smooth.size(); i++)
	{
		fout1 << output.position_vector_smooth[i] << '\t'
			  << output.pressure_vector_smooth[i] << '\t'
			  << output.ion_vector_smooth[i] << std::endl; 
	}

	std::ofstream fout2("combustor_data.dat", std::ios_base::app);
	fout2 //<< std::setw(12) << output.mass_flow_rate_air << '\t'
		  << std::setw(12) << output.acquisition_time << '\t'
		  //<< std::setw(12) << output.mass_flow_rate_fuel_2 << '\t'
		  << std::setw(12) << output.combustor_frequency << '\t'
		  << std::setw(12) << output.pressure_max << '\t'
		  << std::setw(12) << output.pressure_min << '\t'
		  << std::setw(12) << output.phase_max << '\t'
		  << std::setw(12) << output.phase_min << '\t'
		  //<< std::setw(12) << output.fuel_injection_phase << '\t'
		  //<< std::setw(12) << output.fuel_injection_pulse_width << '\t' 
		  << std::setw(12) << output.fuel_pressure << '\t' 
		  << std::setw(12) << output.total_pressure << '\t' 
		  << std::setw(12) << output.total_temperature << '\t'
		  << std::setw(12) << output.p_rms << '\t' 
		  << std::setw(12) << output.UHC << '\t'
		  << std::setw(12) << output.NOx << '\t'
		  << std::setw(12) << output.CO << '\t'
		  << std::setw(12) << output.CO2 << '\t'
		  << std::setw(12) << output.O2 << std::endl;

	std::ofstream fout4(spectrum_file);
    for (int i = 0; i < output.spectrum_magnitude.size(); i++)
        fout4 << output.spectrum_frequency[i] << '\t'
              << output.spectrum_magnitude[i] << std::endl;

	if(optionsMenu.dataAnalysisMenu.spectrogram=="true")
	{
    std::cout << "about to print spetrogram" << std::endl;

    std::ofstream fout5(spectrogram_file);
    for (int i = 0; i < output.time.size(); i++)
    {
        for (int j = 0; j < output.spectrogram_n; j++)
        {
            fout5 << output.time[i] << '\t'
                  << output.spectrogram_f[i][j] << '\t'
                  << output.spectrogram_m[i][j] << std::endl;
        }
        fout5 << std::endl;
    }
	}



	/*

	std::ofstream fout4("pressure_scatter.dat");
	for (int i = 1; i < output.position_scatter.size()-1; i++)
		fout4 << output.position_scatter[i] << '\t'
			  << output.pressure_scatter[i] << '\t'
			  << output.ion_scatter[i] << std::endl;
			  */
}
