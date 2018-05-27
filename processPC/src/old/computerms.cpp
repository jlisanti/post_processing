#include <cmath>
#include <iostream>

#include "datafileinput.h"

void compute_prms(DataFileInput &cInput, double &p_rms, double mean)
{
	double tmp = 0.0;
	int count = 0;
	for (int i = 0; i < cInput.file_length(); i++)
	{
		tmp = tmp + pow(cInput.table_value(i,1)-mean,2);
		count++;
	}
	p_rms = sqrt(tmp/double(count));
	std::cout << "pRMS " << p_rms << std::endl;
}
