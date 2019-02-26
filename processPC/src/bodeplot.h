#ifndef CODES_POSTPROCESSING_PROCESSPC_SRC_BODEPLOT_H_
#define CODES_POSTPROCESSING_PROCESSPC_SRC_BODEPLOT_H_

struct func_param
{
	int N;
	double timestep;
	std::vector<double> p;
};

inline double pinterp(double time,
		              double timestep,
					  std::vector<double> pressure);

inline double func(double tau,
		           void *params);

int buildBodePlot (DataFileInput &cInput,
		           int timeCol,
				   int presCol,
				   std::vector<double> &bodeP,
				   std::vector<double> &bodePprime);

#endif
