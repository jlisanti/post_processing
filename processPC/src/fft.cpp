#include <vector>
#include <iostream>
#include <algorithm>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_statistics.h>

#include "datafileinput.h"

void compute_frequency_spectrum(DataFileInput &cInput,
		                        std::vector<double> &magnitude,
								std::vector<double> &frequency,
								int col,
								double &peak_frequency)
{

	std::cout << peak_frequency << std::endl;

	int n = cInput.file_length();

	std::cout << n << std::endl;
	double * timestep = new double[n];
	double * pressure = new double[n];
	for (int i = 0; i < n; i++)
	{
		timestep[i] = cInput.table_value(i,0);
		pressure[i] = cInput.table_value(i,col);
	}

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

	double Fs = 1.0/(timestep[1] - timestep[0]);

	for (int i = 0; i < n/2; i++)
	{
		frequency.push_back(((Fs*i)/n)/2.0);
		magnitude.push_back(pressure[i]*pressure[i]);
	}

	magnitude[0] = 0.0;
	magnitude[1] = 0.0;
	magnitude[2] = 0.0;
	
	auto it = std::max_element(std::begin(magnitude), std::end(magnitude));
	int index = std::distance(std::begin(magnitude), it);

	peak_frequency = frequency[index];

	gsl_fft_real_workspace_free (work);
	delete[] timestep;
	delete[] pressure;
}
