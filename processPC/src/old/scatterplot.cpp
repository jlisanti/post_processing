#include <iostream>
#include <fstream>
#include <vector>

#include "datafileinput.h"
#include "smoothencoder.h"
#include "scatterplot.h"

void plot_pressure_scatter(DataFileInput &cInput,
		                   int encoder_colmn,
						   int pressur_colmn)
{
	std::vector<double> encoder_vec = smooth_encoder(cInput,
			                                         encoder_colmn);

	std::ofstream fout1("smoothed_encoder.dat");
	for (int i = 0; i < cInput.file_length(); i++)
	{
		fout1 << cInput.table_value(i,0) << '\t'
			  << cInput.table_value(i,1) << '\t'
			  << cInput.table_value(i,2) << '\t'
			  << encoder_vec[i] << std::endl;
	}
}
