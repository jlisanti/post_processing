#include <vector>
#include <iostream>
#include <algorithm>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_statistics.h>

#include "datafileinput.h"
#include "processfile.h"
#include "spectrogram.h"

void build_spectrogram(DataFileInput &cInput,
								double &freq,
								std::vector<std::vector<double>> &spectrogram_f,
								std::vector<std::vector<double>> &spectrogram_m,
								double numberCyclesSpace,
								double numberCyclesWindow,
								int &n,
								std::vector<double> &time)
{
	/* Choose window size
	   Choose step size 
	   Step through file performing fft at each step 
	     accross specified window
	*/

	//double t_cycle = 1.0/freq;
	//int num_cycle = 3;
	//double time_window = num_cycle*t_cycle; 
	//n = ceil(time_window/dt); 

	/* compute acquisition time step - should be better */
	double dt = cInput.table_value(2,0) - cInput.table_value(1,0);
	double spacing = numberCyclesSpace*(1.0/freq);
	double window  = numberCyclesWindow*(1.0/freq);

	std::cout << 1.0/freq << '\t' << spacing << '\t' << window << std::endl;

	n = ceil(window/dt);
	int step = ceil(spacing/dt);

	int current_i = 0;
	//int number_steps = ceil(cInput.file_length()/step);
	int k = 0;

	//int n = cInput.file_length();
	//for (int j = current_i; j < number_steps; j++)
	while ((current_i+n) < cInput.file_length())
	//for (int k = 0; k < 100; k++)
	{
		double * timestep = new double[n];
		double * pressure = new double[n];
		for (int i = 0; i < n; i++)
		{
			//std::cout << i << std::endl;
			timestep[i] = cInput.table_value(current_i+i,0);
			pressure[i] = cInput.table_value(current_i+i,1);
		}
		time.push_back(cInput.table_value(current_i,0));
		//std::cout << "time " << time[k] << std::endl;

	/* NOT SURE WHY THIS DOES NOT WORK - GIVES STRANGE ERROR
	std::copy(cInput.table_column(0).begin(), cInput.table_column(0).end(),
			  timestep);

	std::copy(cInput.table_column(1).begin(), cInput.table_column(1).end(),
			  pressure);

			  */

		gsl_fft_real_wavetable * real;
		gsl_fft_real_workspace * work;

		work = gsl_fft_real_workspace_alloc (n);
		real = gsl_fft_real_wavetable_alloc (n);

		gsl_fft_real_transform (pressure, 1, n, real, work);
		gsl_fft_real_wavetable_free (real);

		
		std::vector<double> frequency;
		std::vector<double> magnitude;

		double Fs = 1.0/(timestep[1] - timestep[0]);

		for (int i = 0; i < n/2; i++)
		{
			frequency.push_back(((Fs*i)/n)/2.0);
			magnitude.push_back(pressure[i]*pressure[i]);
		}

		magnitude[0] = 0.0;
		spectrogram_m.push_back(magnitude);
		spectrogram_f.push_back(frequency);
	
	//auto it = std::max_element(std::begin(magnitude), std::end(magnitude));
	//int index = std::distance(std::begin(magnitude), it);

	//peak_frequency = frequency[index];

		frequency.clear();
		magnitude.clear();

		gsl_fft_real_workspace_free (work);
		delete[] timestep;
		delete[] pressure;

		current_i = current_i+step;
		k++;
	}
	n = floor(n/2.0);
}
