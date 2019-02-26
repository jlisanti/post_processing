#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>

#include <gsl/gsl_poly.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>

#include "findroots.h"

/*

void find_max(std::vector<double> coefficients, double &max, double &phase_max,
		                                        double &min, double &phase_min)
{
	int n = coefficients.size()-1;
	double *a = new double[n];
	double *z = new double[(n-1)*2];

	max = 0.0;
	min = 0.0;
	phase_max = 0.0;
	phase_min = 0.0;

	std::cout << "Finding roots..." << std::endl;
	std::cout << "Number coefficients " << n << std::endl;

	for (int i = 1; i < n+1; i++)
	{
		a[i-1] = double(i)*coefficients[i];
	}

	//for (int i = 0; i < n-1; i++)
//		std::cout << i  << '\t' << coefficients[i] << std::endl;

	gsl_poly_complex_workspace * w 
		= gsl_poly_complex_workspace_alloc (n);
	gsl_poly_complex_solve (a, n, w, z);
	gsl_poly_complex_workspace_free (w);

	for (int i = 0; i < n-1; i++)
    {
		if(z[2*i+1] == 0.0)
		{
			if((z[2*i] < 2.5) && (z[2*i] > 0.0))
			{
				double func = 0.0;
				double t = z[2*i];
				for (int j = 0; j < coefficients.size(); j++)
					func = func + coefficients[j]*pow(t,j);

				if((phase_max != 0.0) && (phase_min != 0.0))
				{
					if (func > max)
					{
						max = func;
						phase_max = t;
					}
					if (func < min)
					{
						min = func;
						phase_min = t;
					}

				}

				if(phase_max == 0.0)
				{
					phase_max = z[2*i];
					max = func;
				}
				else
				{
					phase_min = z[2*i];
					min = func;
				}
			}
		}
		//std::cout << i << '\t' << z[2*i] << '\t' << z[2*i+1] << std::endl;
    }

	double t = phase_max;
	double fit = 0.0;
	for (int j = 0; j < coefficients.size(); j++)
		fit = fit + coefficients[j]*pow(t,j);

	max = fit;
	
	t = phase_min;
	fit = 0.0;
	for (int j = 0; j < coefficients.size(); j++)
		fit = fit + coefficients[j]*pow(t,j);

	min = fit;

	for (int i = 0; i < (n-1)*2; i++)
	{
		if ()
		{
		}
	}

	delete[] a;
	delete[] z;
}
*/

void find_max(std::vector<double> x_col,
		      std::vector<double> y_col,
			  double &max,
			  double &phase_max)
{
	auto it = std::max_element(std::begin(y_col), std::end(y_col));
	int index = std::distance(std::begin(y_col), it);

	int window = 10;
	int p = 3;
	double chisq;

	gsl_multifit_linear_workspace * ws;
	gsl_matrix *X, *cov;
	gsl_vector *y, *c;

	int n = 2*window + 1;
	
	X = gsl_matrix_alloc(n,p);
	y = gsl_vector_alloc(n);

	c = gsl_vector_alloc (p);
	cov = gsl_matrix_alloc(p,p);

	for (int i = 0; i < n; i++)
	{
		double xi = x_col[index - window + i];
		double yi = y_col[index - window + i];
		for (int j = 0; j < p; j++)
			gsl_matrix_set (X,i,j,pow(xi,j));
		gsl_vector_set (y,i,yi);
	}
	
	ws = gsl_multifit_linear_alloc(n,p);
	gsl_multifit_linear(X,y,c,cov,&chisq,ws);

	/*
	std::vector<double> coefficients;
	for (int i = 0; i < p; i++)
		coefficients.push_back(gsl_vector_get(c,i));
		*/

	double A = gsl_vector_get (c,0);
	double B = gsl_vector_get (c,1);
	double C = gsl_vector_get (c,2);

	std::ofstream fout("max_fit.dat");
	for (int i = 0; i < n; i++)
	{
		double fit = 0.0;
		double t = x_col[index - window + i];

		fit = gsl_vector_get(c,0) + gsl_vector_get(c,1)*t +
			  gsl_vector_get(c,2)*pow(t,2);
		fout << t << '\t'  << fit << std::endl;

	}

	gsl_multifit_linear_free(ws);
	gsl_matrix_free(X);
	gsl_matrix_free(cov);
	gsl_vector_free(y);
	gsl_vector_free(c);

	phase_max = -B/(2.0*C);
	max = A + phase_max*B + pow(phase_max,2)*C;

	/*
    n = p;
	double *a = new double[n];
	double *z = new double[(n-1)*2];

	for (int i = 1; i < n; i++)
	{
		a[i-1] = double(i)*coefficients[i-1];
		//std::cout << a[i-1] << std::endl;
	}

	std::cout << "shits on fire" << std::endl;
	gsl_poly_complex_workspace * w 
		= gsl_poly_complex_workspace_alloc (n);
	gsl_poly_complex_solve (a, n, w, z);
	gsl_poly_complex_workspace_free (w);

	std::cout << "made it  here" << std::endl;
	for (int i = 0; i < n-1; i++)
    {
		std::cout << z[2*i] << '\t' << z[2*i+1] << std::endl;
		if(z[2*i+1] == 0.0)
		{
			if((z[2*i] < x_col[index+window]) && (z[2*i] > x_col[index-window]))
			{
				double func = 0.0;
				double t = z[2*i];
				for (int j = 0; j < p; j++)
					func = func + gsl_vector_get(c,j)*pow(t,j);

				min = func;
				phase_min = t;

			}
		}
    }


	*/

}

void find_min(std::vector<double> x_col,
		      std::vector<double> y_col,
			  double &min,
			  double &phase_min)
{
	auto it = std::min_element(std::begin(y_col), std::end(y_col));
	int index = std::distance(std::begin(y_col), it);

	int window = 20;
	int p = 3;
	double chisq;

	gsl_multifit_linear_workspace * ws;
	gsl_matrix *X, *cov;
	gsl_vector *y, *c;

	int n = 2*window + 1;
	
	X = gsl_matrix_alloc(n,p);
	y = gsl_vector_alloc(n);

	c = gsl_vector_alloc (p);
	cov = gsl_matrix_alloc(p,p);

	for (int i = 0; i < n; i++)
	{
		double xi = x_col[index - window + i];
		double yi = y_col[index - window + i];
		for (int j = 0; j < p; j++)
			gsl_matrix_set (X,i,j,pow(xi,j));
		gsl_vector_set (y,i,yi);
	}
	
	ws = gsl_multifit_linear_alloc(n,p);
	gsl_multifit_linear(X,y,c,cov,&chisq,ws);

	/*
	std::vector<double> coefficients;
	for (int i = 0; i < p; i++)
		coefficients.push_back(gsl_vector_get(c,i));
		*/

	double A = gsl_vector_get (c,0);
	double B = gsl_vector_get (c,1);
	double C = gsl_vector_get (c,2);

	std::ofstream fout("min_fit.dat");
	for (int i = 0; i < n; i++)
	{
		double fit = 0.0;
		double t = x_col[index - window + i];

		fit = gsl_vector_get(c,0) + gsl_vector_get(c,1)*t +
			  gsl_vector_get(c,2)*pow(t,2);
		fout << t << '\t'  << fit << std::endl;

	}

	gsl_multifit_linear_free(ws);
	gsl_matrix_free(X);
	gsl_matrix_free(cov);
	gsl_vector_free(y);
	gsl_vector_free(c);

	phase_min = -B/(2.0*C);
	min = A + phase_min*B + pow(phase_min,2)*C;

	/*
    n = p;
	double *a = new double[n];
	double *z = new double[(n-1)*2];

	for (int i = 1; i < n; i++)
	{
		a[i-1] = double(i)*coefficients[i-1];
		//std::cout << a[i-1] << std::endl;
	}

	std::cout << "shits on fire" << std::endl;
	gsl_poly_complex_workspace * w 
		= gsl_poly_complex_workspace_alloc (n);
	gsl_poly_complex_solve (a, n, w, z);
	gsl_poly_complex_workspace_free (w);

	std::cout << "made it  here" << std::endl;
	for (int i = 0; i < n-1; i++)
    {
		std::cout << z[2*i] << '\t' << z[2*i+1] << std::endl;
		if(z[2*i+1] == 0.0)
		{
			if((z[2*i] < x_col[index+window]) && (z[2*i] > x_col[index-window]))
			{
				double func = 0.0;
				double t = z[2*i];
				for (int j = 0; j < p; j++)
					func = func + gsl_vector_get(c,j)*pow(t,j);

				min = func;
				phase_min = t;

			}
		}
    }


	*/

}

void find_sub_inlet_pressure()
{
}
