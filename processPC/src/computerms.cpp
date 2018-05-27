#include <cmath>
#include <iostream>

#include <gsl/gsl_statistics.h>

#include "datafileinput.h"

double compute_prms(DataFileInput &cInput, double mean)
{
	int nn = cInput.file_length();
    double * data = new double[nn];

    for (int i = 0; i < nn; i++)
        data[i] = cInput.table_value(i,1);

    delete[] data;

    return sqrt(gsl_stats_variance(data, 1, nn));
}
