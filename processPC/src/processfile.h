#ifndef CODES_POSTPROCESSING_PROCESSPC_SRC_PROCESSFILE_H_
#define CODES_POSTPROCESSING_PROCESSPC_SRC_PROCESSFILE_H_

#include <string>
#include <vector>

#include "readinputfiles.h"

typedef struct {
	double mass_flow_rate_air;
	double mass_flow_rate_fuel_1;
	double mass_flow_rate_fuel_2;
	std::string fuel;
	std::string file_name;
	double combustor_frequency;
	double ion_frequency;
	double mean_combustor_pressure;
	double phase_max;
	double phase_min;
	double pressure_max;
	double pressure_min;
	double ion_min;
	double ion_phase;
	double fuel_injection_phase;
	double fuel_injection_pulse_width;
	double p_rms;
	double p_static;
	double static_pressure;
	std::vector<double> spectrum_magnitude;
	std::vector<double> spectrum_frequency;
	std::vector<double> ion_spectrum_magnitude;
	std::vector<double> ion_spectrum_frequency;
	std::vector<double> encoder_vector_smooth;
	std::vector<double> pressure_vector_smooth;
	std::vector<double> ion_vector_smooth;
	std::vector<double> pressure_scatter;
	std::vector<double> encoder_scatter;
	std::vector<double> position_vector_smooth;
	std::vector<double> position_scatter;
	std::vector<double> ion_scatter;
	std::vector<double> coefficients;
	std::vector<double> ion_vector;
	std::vector<double> position_vector;
	std::vector<double> encoder_vector;
	std::vector<double> pressure_vector;
	std::vector<double> upCrossOver;
	std::vector<double> downCrossOver;
	bool fuel_injection;
	double mdot_air;
	double phi;
	double AF_ratio;
	std::vector<double> mdot_inlet;
	std::vector<double> area_inlet;
	std::vector<double> area_inlet_2;
	double total_pressure;
	double total_temperature;
	double fuel_pressure;
	std::vector<std::vector<double>> ion_spectrogram_m;
	std::vector<std::vector<double>> ion_spectrogram_f;
	std::vector<std::vector<double>> spectrogram_m;
	std::vector<std::vector<double>> spectrogram_f;
	int spectrogram_n;
	double dt;
	std::vector<double> time;
	std::vector<double> ion_time;
	std::vector<double> time_ssa;
	std::vector<double> pressure_ssa;
	std::vector<double> ion_ssa;
	int M;
	int kmax;
	double window;
	double acquisition_time;
	int encoder_colmn;
	std::vector<double> e_max_p;
	std::vector<double> e_min_p;
	std::vector<double> e_max_i;
	std::vector<double> e_min_i;
	std::vector<double> t_max_p;
	std::vector<double> t_min_p;
	std::vector<double> t_max_i;
	std::vector<double> t_min_i;
	std::vector<double> p_max;
	std::vector<double> p_min;
	std::vector<double> i_max;
	std::vector<double> i_min;
	double UHC;
	double NOx;
	double CO;
	double CO2;
	double O2;
} output_data;

void process_file(std::string file, output_data &output, options &optionsMenu);
void read_header_file_data(std::string file, output_data &output, int &skip_lines);

#endif
