#include <iostream>
#include <string>
#include <cstdio>
#include <cmath>
#include <istream>
#include <fstream>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "datafileinput.h"
#include "bodeplot.h"


inline double pinterp(double time, double timestep,
        std::vector<double> pressure)
{
    double p;
    double tol = 1.e-10;
    int iHigh = ceil(time/timestep);
    int iLow = floor(time/timestep);
    if(iHigh == iLow)
       iHigh = iLow + 1;
	/*
    if(std::abs((iLow*timestep)-time) < tol)
        iLow = iLow - 1;
    if(std::abs((iHigh*timestep)-time) < tol)
        iHigh = iHigh + 1;
		*/
	if(iHigh > (pressure.size()-1))
	{
		iHigh = pressure.size()-1;
		iLow  = pressure.size()-2;
	}
	if(((double(iHigh)*timestep)-time) < tol)
		return pressure[iHigh];
	else if(((double(iLow)*timestep)-time) < tol)
		return pressure[iLow];
	else 
	{
        int N = 2;
        double x[2] = {iLow*timestep,iHigh*timestep};
        double y[2] = {pressure[iLow],pressure[iHigh]};

        gsl_interp_accel * acc = gsl_interp_accel_alloc ();
        const gsl_interp_type * t = gsl_interp_linear;
        gsl_spline * spline = gsl_spline_alloc (t,N);

        gsl_spline_init (spline,x,y,N);
        p = gsl_spline_eval (spline,time,acc);
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);

        return p;
	}
}

inline double func(double tau, void *params)
{
    struct func_param input = *(struct func_param *)params;
    double timestep = input.timestep;
    std::vector<double> p = input.p;
    int N = floor(((input.N*timestep)-tau)/timestep);
    double xbar = 0.0;

    for (int i=0; i<N; i++)
        xbar = xbar + p[i];

    xbar = xbar*(1.0/N);

    double c1 = 0.0;
    double c2 = 0.0;
    double c3 = 1.0/N;
    for (int i=0; i<N; i++)
        c1 = c1+((pinterp((i*timestep)+tau,timestep,p)-xbar)*(p[i]-xbar));
    for (int i=0; i<N; i++)
        c2 = c2+pow((p[i]-xbar),2);
    return (c3*c1)/(c3*c2);
}

int buildBodePlot(DataFileInput &cInput,
		          int timeCol,
				  int presCol,
				  std::vector<double> &bodeP,
				  std::vector<double> &bodePprime)
{
    double timeStep = cInput.table_value(1,timeCol) 
		               - cInput.table_value(0,timeCol);

	struct func_param input {
		cInput.file_length(),
		timeStep,
		cInput.table_column(presCol),
	};

	double t = 0.0;
	double dt = timeStep/1.0;
	double dtau = 0.0;
	int status;
	int max_iter = 100;

	const gsl_root_fsolver_type * T;
	gsl_root_fsolver * s;

	double r    = 0.0;
	double x_lo = 0.0;
	double x_hi = 1.0e-6;
	double tau  = 0.0;
	int max_step = 10000;

    for (int i=1; i<=max_step; i++)
    {
        dtau = dt*i;
        double test = func(dtau,&input);
        std::cout << " tau: " << dtau << " C(tau): " << test << std::endl;
        if (i == max_step)
            std::cout << "  No crossover found within max iteration"
                << std::endl;
        if (test < 0.0)
        {
            x_lo = dt*(i-1);
            x_hi = dtau;
            break;
        }
    }

	std::cout << "x_lo = " << x_lo << '\t'
		      << "x_hi = " << x_hi << '\t' 
			  << "dt =   " << dt   << std::endl;

	gsl_function F;

	F.function = &func;
	F.params = &input;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	for (int i = 1; i < max_iter; i++)
	{
		status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi,
                0, 1.0e-10);
        if (status == GSL_SUCCESS)
        {
            std::cout << "Converged: " << std::endl;
            tau = x_lo;
            std::cout << "tau: " << '\t' << tau << std::endl;
            break;
        }

        std::cout << i << '\t' <<  x_lo << '\t' << x_hi << '\t' << r << '\t'
                        << (x_hi - x_lo) << std::endl;
	}

	for (int i = 0; i < floor(((double(cInput.file_length())*timeStep)-tau)/timeStep); i++)
	{
		bodeP.push_back(cInput.table_value(i,presCol));
		bodePprime.push_back((pinterp((double(i)*timeStep)+tau,
					          timeStep,
							  cInput.table_column(presCol))));
	}

	return 0;
}
