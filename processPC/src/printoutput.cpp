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
#include "computeballvalvearea.h"

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


void print_output_active (output_data &output, options &optionsMenu)
{
	std::cout << " Printing output" << std::endl;
	std::string mainFileOutput;
	std::string fileFront;

	double encoderStep = M_PI/2.5;

	// Round combustor frequency
	output.combustor_frequency = round_value(output.combustor_frequency,5);
	std::ostringstream frequency;
	std::ostringstream massFlowRateA;

	if (output.mass_flow_rate_fuel_1 < 10.0)
		output.mass_flow_rate_fuel_1 = output.mass_flow_rate_fuel_1*1000.0;

	if (output.mass_flow_rate_fuel_1 < 10.0)
		output.mass_flow_rate_fuel_1 = output.mass_flow_rate_fuel_1*1000.0;

	frequency << int(output.combustor_frequency);
	massFlowRateA << int(output.mass_flow_rate_fuel_1); 

	fileFront = frequency.str() + "Hz_" + massFlowRateA.str() + "mgps_" 
		+ optionsMenu.mainMenu.fuelA;

	if(optionsMenu.mainMenu.numberFuels > 1)
	{
		std::ostringstream massFlowRateB;
	    if (output.mass_flow_rate_fuel_2 < 10.0)
	    	output.mass_flow_rate_fuel_2 = output.mass_flow_rate_fuel_2*1000.0;

	    massFlowRateB << int(output.mass_flow_rate_fuel_2); 

	    fileFront = frequency.str() + "Hz_" + massFlowRateA.str() + "mgps_" 
	    	+ optionsMenu.mainMenu.fuelA + massFlowRateB.str() + "mgps_" 
	    	+ optionsMenu.mainMenu.fuelB;
	}

	mainFileOutput = fileFront + "_output.dat";

    std::string averagedFile    = fileFront + "_scatter.dat";
    std::string averagedCurve   = fileFront + "_averaged.dat";



	/* check for duplicate cases and average */
	bool end = false;
	std::string fileNameNext = "none";

	if(file_exists(mainFileOutput))
	{
		/*
		std::cout << "***************" << std::endl;
		std::cout << mainFileOutput << '\t' << fileNameNext << std::endl;
		std::cout << "***************" << std::endl;
		*/

	    for (int i = 0; i < 100; i++)
	    {
		    std::ostringstream fileIndex;
		    std::ostringstream fileIndexi;
		    fileIndex  << (i+1);
		    fileIndexi << (i+2);

		    fileNameNext   = fileFront + "_output_" + fileIndex.str() + ".dat";

			std::cout << std::endl;
			std::cout << "file name " << mainFileOutput << '\t'
				      << "next      " << fileNameNext   << std::endl;
			std::cout << std::endl;

		    if(!file_exists(fileNameNext))
		    {
				end = true;
			    // Add to averaged file
			    // Check if averaged file exist
			    // Combine
		        if(!file_exists(averagedFile))
			    {
				    // read existing file
				    int skip_lines = 0;
				    DataFileInput cInput(mainFileOutput,skip_lines);

				    //  create file
				    std::ofstream fout_avg(averagedFile);
				
				    int count = 0;


	                if(optionsMenu.mainMenu.ionProbe=="true")
				    {
				        for (int j = 0; j < cInput.file_length(); j++)
				        {
					        for (int k = count; k < output.encoder_vector_smooth.size(); k++)
					        {
						        if (cInput.table_value(j,0) > output.encoder_vector_smooth[k])
						        {
							        fout_avg << output.encoder_vector_smooth[k] << '\t'
								             << output.pressure_vector_smooth[k] << '\t'
								             << output.ion_vector_smooth[k] << '\t'
			                                 << ballValveArea(output.encoder_vector_smooth[i]*encoderStep,
					                            output.areaFunction) << '\t'
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
				        for (int j = 0; j < cInput.file_length(); j++)
				        {
					        for (int k = count; k < output.encoder_vector_smooth.size(); k++)
					        {
						        if (cInput.table_value(j,0) > output.encoder_vector_smooth[k])
						        {
							        fout_avg << output.encoder_vector_smooth[k] << '\t'
								             << output.pressure_vector_smooth[k] << '\t'
			                                 << ballValveArea(output.encoder_vector_smooth[i]*encoderStep,
					                            output.areaFunction) << '\t'
								             << std::endl;
							        count++;
						        }
						        else
						        {
							        fout_avg << cInput.table_value(j,0) << '\t' 
								             << cInput.table_value(j,1) << '\t'
								             << cInput.table_value(j,2) << std::endl; 
							        break;
						        }
					        }
				        }
				    }
			    }
		    	else
		    	{
					std::cout << "*********" << std::endl;
					std::cout << "fileNext did not exist" << std::endl;
					std::cout << "*********" << std::endl;
		    		// read existing file
		    		int skip_lines = 0;
		    		DataFileInput cInput(averagedFile,skip_lines);

		    		//  create file
		    		std::ofstream fout_avg(averagedFile);

		    		int count = 0;


	                if(optionsMenu.mainMenu.ionProbe=="true")
			    	{
			    	    for (int j = 0; j < cInput.file_length(); j++)
			    	    {
			    		    for (int k = count; k < output.encoder_vector_smooth.size(); k++)
			    		    {
			    			    if (cInput.table_value(j,0) > output.encoder_vector_smooth[k])
			    			    {
			    				    fout_avg << output.encoder_vector_smooth[k] << '\t'
			    				             << output.pressure_vector_smooth[k] << '\t'
			    					         << output.ion_vector_smooth[k] << '\t'
			                                 << ballValveArea(output.encoder_vector_smooth[i]*encoderStep,
			    		                        output.areaFunction) << '\t'
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
				        for (int j = 0; j < cInput.file_length(); j++)
				        {
					        for (int k = count; k < output.encoder_vector_smooth.size(); k++)
					        {
						        if (cInput.table_value(j,0) > output.encoder_vector_smooth[k])
						        {
							        fout_avg << output.encoder_vector_smooth[k] << '\t'
							                 << output.pressure_vector_smooth[k] << '\t'
			                                 << ballValveArea(output.encoder_vector_smooth[i]*encoderStep,
					                            output.areaFunction) << '\t'
								             << std::endl;
							        count++;
						        }
						        else
						        {
							        fout_avg << cInput.table_value(j,0) << '\t' 
								             << cInput.table_value(j,1) << '\t'
								             << cInput.table_value(j,2) << std::endl; 
							        break;
						        }
					        }
				        }
				    }
			    }

				DataFileInput cInput_avg(averagedFile,0);
				std::vector<double> p_vec;
				std::vector<double> e_vec;
				std::vector<double> i_vec;
				std::vector<double> p_scat;
				std::vector<double> e_scat;
				std::vector<double> i_scat;
		        std::vector<double> tmp1 = cInput_avg.table_column(0);
                std::vector<double> tmp2 = cInput_avg.table_column(1);
                std::vector<double> tmp3 = cInput_avg.table_column(2);

				average_files_active(tmp1,
							         tmp2,
							         tmp3,
							         e_vec,
							         p_vec,
							         i_vec,
							         optionsMenu.mainMenu.ionProbe);

				std::cout << "printing averaged curve" << std::endl;
				std::ofstream fout_curv(averagedCurve);
	            if(optionsMenu.mainMenu.ionProbe=="true")
				{
				    for (int i = 0; i < e_vec.size(); i++)
				    {
					    fout_curv << e_vec[i] << '\t'
						          << p_vec[i] << '\t'
							      << i_vec[i] << '\t'
			                      << ballValveArea(e_vec[i]*encoderStep,
					                 output.areaFunction) << '\t'
					              << std::endl;
				    }
				}
				else
				{
				    for (int i = 0; i < e_vec.size(); i++)
				    {
					    fout_curv << e_vec[i] << '\t'
						          << p_vec[i] << '\t'
			                      << ballValveArea(e_vec[i]*encoderStep,
					                 output.areaFunction) << '\t'
					              << std::endl;
				    }
				}

		        mainFileOutput = fileFront + "_output_" + fileIndex.str()  + ".dat";
		        fileNameNext   = fileFront + "_output_" + fileIndexi.str() + ".dat";
		    }

	        if(end)	
			    break;
		}

	}

    std::string spectrumFile    = fileFront + "_spectrum.dat";
    std::string spectrogramFile = fileFront + "_spectrogram.dat";

    std::string scatterFile     = fileFront + "_pressure_scatter.dat";
    std::string pMaxFile        = fileFront + "_p_max_peaks.dat";
    std::string pMinFile        = fileFront + "_p_min_peaks.dat";
	std::string peaksFile       = fileFront + "_peaks.dat";

	if(optionsMenu.mainMenu.ionProbe=="true")
	{
	    std::ofstream foutPrimary(mainFileOutput);
	    for (int i = 0; i < output.encoder_vector_smooth.size(); i++)
	    {
		    foutPrimary << output.encoder_vector_smooth[i]    << '\t'
		    	        << output.pressure_vector_smooth[i]   << '\t'
		    	        << output.ion_vector_smooth[i]        << '\t'
		    	        << ballValveArea(output.encoder_vector_smooth[i]*encoderStep,
		    			           output.areaFunction) << '\t'
		    	        << std::endl;
	    }
	}
	else
	{
	    std::ofstream foutPrimary(mainFileOutput);
	    for (int i = 0; i < output.encoder_vector_smooth.size(); i++)
	    {
		    foutPrimary << output.encoder_vector_smooth[i]    << '\t'
		    	        << output.pressure_vector_smooth[i]   << '\t'
		    	        << ballValveArea(output.encoder_vector_smooth[i]*encoderStep,
		    			           output.areaFunction) << '\t'
		    	        << std::endl;
	    }
	}

	if (!file_exists("combustor_data.dat"))
	{
		std::ofstream fout3("combustor_data.dat");
		fout3 << std::setw(12) << "mdot_air"  << '\t'
			  << std::setw(12) << "mdot_fuel" << '\t'
			  << std::setw(12) << "equiv_rat" << '\t'
			  << std::setw(12) << "frequency" << '\t'
			  << std::setw(12) << "max_P" << '\t'
			  << std::setw(12) << "min_P" << '\t'
			  << std::setw(12) << "max_phase" << '\t'
			  << std::setw(12) << "min_phase" << '\t'
			  << std::setw(12) << "Fuel_pres" << '\t'
			  << std::setw(12) << "delP" << '\t'
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

	std::ofstream fout2("combustor_data.dat", std::ios_base::app);
	fout2 << std::setw(12) << output.massFlowRateAir << '\t'
		  << std::setw(12) << output.mass_flow_rate_fuel_1 << '\t'
		  << std::setw(12) << output.equivalenceRatio << '\t'
		  << std::setw(12) << output.combustor_frequency << '\t'
		  << std::setw(12) << output.pressure_max << '\t'
		  << std::setw(12) << output.pressure_min << '\t'
		  << std::setw(12) << output.phase_max << '\t'
		  << std::setw(12) << output.phase_min << '\t'
		  //<< std::setw(12) << output.fuel_injection_phase << '\t'
		  //<< std::setw(12) << output.fuel_injection_pulse_width << '\t' 
		  << std::setw(12) << output.fuel_pressure << '\t' 
		  << std::setw(12) << output.total_pressure << '\t' 
		  << std::setw(12) << 1.0+(output.total_pressure/output.p_static) << '\t' 
		  << std::setw(12) << output.total_temperature << '\t'
		  << std::setw(12) << output.p_rms << '\t' 
		  << std::setw(12) << output.ion_phase << '\t'
		  << std::setw(12) << output.ion_min << '\t'
		  << std::setw(12) << output.UHC << '\t'
		  << std::setw(12) << output.NOx << '\t'
		  << std::setw(12) << output.CO << '\t'
		  << std::setw(12) << output.CO2 << '\t'
		  << std::setw(12) << output.O2 << std::endl;


	std::ofstream fout4(spectrumFile);
	if(optionsMenu.mainMenu.ionProbe=="true")
	{
        for (int i = 0; i < output.spectrum_magnitude.size(); i++)
            fout4 << output.spectrum_frequency[i] << '\t'
                  << output.spectrum_magnitude[i] << '\t' 
			      << output.ion_spectrum_frequency[i] << '\t'
			      << output.ion_spectrum_magnitude[i] << std::endl;
	}
	else
	{
        for (int i = 0; i < output.spectrum_magnitude.size(); i++)
            fout4 << output.spectrum_frequency[i] << '\t'
                  << output.spectrum_magnitude[i] << std::endl;
	}

	if(optionsMenu.dataAnalysisMenu.spectrogram=="true")
	{
        std::ofstream fout5(spectrogramFile);
	    if(optionsMenu.mainMenu.ionProbe=="true")
		{
            for (int i = 0; i < output.time.size(); i++)
            {
                for (int j = 0; j < output.spectrogram_n; j++)
                {
                    fout5 << output.time[i] << '\t'
                          << output.spectrogram_f[i][j] << '\t'
                          << output.spectrogram_m[i][j] << '\t' 
                          << output.ion_time[i] << '\t' 
				          << output.ion_spectrogram_f[i][j] << '\t'
				          << output.ion_spectrogram_m[i][j] << std::endl;
                }
                fout5 << std::endl;
            }
		}
		else
		{
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
	}


	if (!file_exists(peaksFile))
	{
	    std::ofstream fout10(peaksFile);
		fout10 << "tMaxP" << '\t'
			   << "pMax" << '\t'
			   << "tMinP" << '\t'
			   << "pMin" << '\t'
			   << "tMinIP" << '\t'
			   << "IPmin" << std::endl;
	    for (int i = 0; i < output.p_max.size(); i++)
		    fout10 << output.t_max_p[i] << '\t'
			       << output.p_max[i]   << '\t'
			       << output.t_min_p[i] << '\t'
			       << output.p_min[i]   << '\t'
			       << output.t_min_i[i] << '\t'
			       << output.i_min[i]   << std::endl;
	}



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
	std::ostringstream convert5;

	convert << int(round_comb_freq);
	convert1 << int(output.total_temperature);
	//file_name = convert.str() + "Hz_" + convert1.str() + "C.dat";
	convert2 << int(optionsMenu.mainMenu.gasoline);
	convert3 << int(optionsMenu.mainMenu.heptane);
	convert4 << int(optionsMenu.mainMenu.ethanol);
	convert5 << int(optionsMenu.mainMenu.diesel);

	file_name = convert2.str() + "G_" + convert3.str() 
		         + "H_" + convert4.str() + "E_" + convert5.str() + "D.dat"; 
    //std::string averaged_file = convert.str() + "Hz_" + convert1.str() + "C_scatter.dat";
    //std::string averaged_curve = convert.str() + "Hz_" + convert1.str() + "C_averaged.dat";

    std::string averaged_file = convert2.str() + "G_" + convert3.str() 
		         + "H_" + convert4.str() + "E_" + convert5.str() + "D_scatter.dat"; 

    std::string averaged_curve = convert2.str() + "G_" + convert3.str() 
		         + "H_" + convert4.str() + "E_" + convert5.str() + "D_average.dat"; 

    std::string spectrum_file = convert.str() + "Hz_" + convert1.str() + "C_spectrum.dat";
    std::string spectrogram_file = convert.str() + "Hz_" + convert1.str() + "C_spectrogram.dat";
	
    std::string p_max_peaks_file = convert.str() + "Hz_" + convert1.str() + "C_p_max_peaks.dat";
    std::string i_max_peaks_file = convert.str() + "Hz_" + convert1.str() + "C_i_max_peaks.dat";

    std::string p_min_peaks_file = convert.str() + "Hz_" + convert1.str() + "C_p_min_peaks.dat";
    std::string i_min_peaks_file = convert.str() + "Hz_" + convert1.str() + "C_i_min_peaks.dat";

    std::string peaks_file = convert.str() + "Hz_" + convert1.str() + "C_peaks.dat";

    bool end = false;
	std::string file_name_next;


	if(file_exists(file_name))
	{

		for (int i = 0; i < 100; i++)
		{
		    std::ostringstream convert6;
		    std::ostringstream convert7;
		    convert6 << (i+1);
		    convert7 << (i+2);

			if(!(file_exists(file_name_next)))
		    {
				end = true;
			    // Add to averaged file
			    // Check if averaged file exist
			    // Combine
		        if(!file_exists(averaged_file ))
			    {
				    // read existing file
				    int skip_lines = 0;
				    DataFileInput cInput(file_name,skip_lines);

				    std::cout << "creating average file" << std::endl;

				    //  create file
				    std::ofstream fout_avg(averaged_file);
				
				    int count = 0;

				    for (int j = 0; j < cInput.file_length(); j++)
			    	{
					    for (int k = count; k < output.time_vector_smooth.size(); k++)
				    	{
						    if (cInput.table_value(j,0) > output.time_vector_smooth[k])
						    {
							    fout_avg << output.time_vector_smooth[k] << '\t'
			                             << output.pressure_vector_smooth[k] << '\t'
			                             << output.ion_vector_smooth[k] << '\t' 
			                             << output.phase_vector_smooth[k] << std::endl; 
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

			        std::cout << "average file already exits. File length: " << 
						cInput.file_length() << std::endl;

				    for (int j = 0; j < cInput.file_length(); j++)
				    {
					    for (int k = count; k < output.time_vector_smooth.size(); k++)
					    {
						    if (cInput.table_value(j,0) > output.time_vector_smooth[k])
						    {
							    fout_avg << output.time_vector_smooth[k] << '\t'
			                             << output.pressure_vector_smooth[k] << '\t'
			                             << output.ion_vector_smooth[k] << '\t' 
			                             << output.phase_vector_smooth[k] << std::endl; 
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
			    	std::vector<double> timeVec;
			    	std::vector<double> pressureVec;
			    	std::vector<double> ionProbeVec;
			    	std::vector<double> phaseVec;

			    	average_files(cInput_avg.table_column(0),
				    		      cInput_avg.table_column(1),
				    			  cInput_avg.table_column(2),
				    			  timeVec,
				    			  pressureVec,
				    			  ionProbeVec,
				    			  phaseVec);

				    std::ofstream fout_curv(averaged_curve);

				    for (int i = 0; i < timeVec.size(); i++)
			    	{
					    fout_curv << timeVec[i]      << '\t'
						          << pressureVec[i]  << '\t'
						    	  << ionProbeVec[i]  << '\t'
						    	  << phaseVec[i]     << std::endl;
				    }
			    }

			    file_name = convert2.str() + "G_" + convert3.str() 
		             + "H_" + convert4.str() + "E_" + convert5.str() + "D_" + convert7.str() + ".dat"; 

				file_name_next = convert2.str() + "G_" + convert3.str() 
		             + "H_" + convert4.str() + "E_" + convert5.str() + "D_" + convert7.str() + ".dat"; 
		    }
			if(end)
				break;
	    }
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
			  << std::setw(12) << "ion_phase" << '\t'
			  << std::setw(12) << "ion_min" << '\t'
			  << std::setw(12) << "UHC" << '\t'
			  << std::setw(12) << "NOx" << '\t'
			  << std::setw(12) << "CO" << '\t'
			  << std::setw(12) << "CO2" << '\t'
			  << std::setw(12) << "O2" << '\t' 
			  << std::setw(12) << "gasoline" << '\t'
			  << std::setw(12) << "heptane" << '\t'
			  << std::setw(12) << "ethanol" << '\t'
			  << std::setw(12) << "diseal" << std::endl;  

	}
	

//	if(!file_exists(file_name))
//	{
	std::ofstream fout1(file_name);
	for (int i = 0; i < output.time_vector_smooth.size(); i++)
	{
		fout1 << output.time_vector_smooth[i] << '\t'
			  << output.pressure_vector_smooth[i] << '\t'
			  << output.ion_vector_smooth[i] << '\t' 
			  << output.phase_vector_smooth[i] << std::endl; 
	}
//	}

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
		  << std::setw(12) << output.ion_phase << '\t'
		  << std::setw(12) << output.ion_min << '\t'
		  << std::setw(12) << output.UHC << '\t'
		  << std::setw(12) << output.NOx << '\t'
		  << std::setw(12) << output.CO << '\t'
		  << std::setw(12) << output.CO2 << '\t'
		  << std::setw(12) << output.O2 << '\t' 
		  << std::setw(12) << optionsMenu.mainMenu.gasoline << '\t'
		  << std::setw(12) << optionsMenu.mainMenu.heptane << '\t'
		  << std::setw(12) << optionsMenu.mainMenu.ethanol << '\t'
		  << std::setw(12) << optionsMenu.mainMenu.diesel << std::endl; 

	std::ofstream fout4(spectrum_file);
    for (int i = 0; i < output.spectrum_magnitude.size(); i++)
        fout4 << output.spectrum_frequency[i] << '\t'
              << output.spectrum_magnitude[i] << '\t' 
			  << output.ion_spectrum_frequency[i] << '\t'
			  << output.ion_spectrum_magnitude[i] << std::endl;

	if(optionsMenu.dataAnalysisMenu.spectrogram=="true")
	{
        std::ofstream fout5(spectrogram_file);
        for (int i = 0; i < output.time.size(); i++)
        {
            for (int j = 0; j < output.spectrogram_n; j++)
            {
                fout5 << output.time[i] << '\t'
                      << output.spectrogram_f[i][j] << '\t'
                      << output.spectrogram_m[i][j] << '\t' 
                      << output.ion_time[i] << '\t' 
				      << output.ion_spectrogram_f[i][j] << '\t'
				      << output.ion_spectrogram_m[i][j] << std::endl;
            }
            fout5 << std::endl;
        }
	}


	if (!file_exists(peaks_file))
	{
	    std::ofstream fout10(peaks_file);
		fout10 << "tMaxP" << '\t'
			   << "pMax" << '\t'
			   << "tMinP" << '\t'
			   << "pMin" << '\t'
			   << "tMinIP" << '\t'
			   << "IPmin" << std::endl;
	    for (int i = 0; i < output.p_max.size(); i++)
		    fout10 << output.t_max_p[i] << '\t'
			       << output.p_max[i]   << '\t'
			       << output.t_min_p[i] << '\t'
			       << output.p_min[i]   << '\t'
			       << output.t_min_i[i] << '\t'
			       << output.i_min[i]   << std::endl;
	}

	/*
	std::ofstream fout6(p_max_peaks_file);
	for (int i = 0; i < output.p_max.size(); i++)
		fout6 << output.t_max_p[i] << '\t'
			  << output.p_max[i] << std::endl;

	std::ofstream fout7(p_min_peaks_file);
	for (int i = 0; i < output.p_min.size(); i++)
		fout7 << output.t_min_p[i] << '\t'
			  << output.p_min[i] << std::endl;

	std::ofstream fout8(i_max_peaks_file);
	for (int i = 0; i < output.i_max.size(); i++)
		fout8 << output.t_max_i[i] << '\t'
			  << output.i_max[i] << std::endl;

	std::ofstream fout9(i_min_peaks_file);
	for (int i = 0; i < output.i_min.size(); i++)
		fout9 << output.t_min_i[i] << '\t'
			  << output.i_min[i] << std::endl;
			  */
}
