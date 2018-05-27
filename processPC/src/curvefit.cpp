#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
//#include <gsl/gsl_multifit_nlin.h>

#include "datafileinput.h"

/*
struct data {
	size_t n;
	double * t;
	double * y;
	double * sigma;
};

int sin_f ( const gsl_vector * x, void *data, gsl_vector * f)
{
	size_t n = ((struct data *)data)->n;
	double *t = ((struct data *)data)->t;
	double *y = ((struct data *)data)->y;
	double *sigma = ((struct data *)data)->sigma;

	double A1 = gsl_vector_get (x,0);
	double B1 = gsl_vector_get (x,1);
	double C1 = gsl_vector_get (x,2);

	double A2 = gsl_vector_get (x,3);
	double B2 = gsl_vector_get (x,4);
	double C2 = gsl_vector_get (x,5);

	double A3 = gsl_vector_get (x,6);
	double B3 = gsl_vector_get (x,7);
	double C3 = gsl_vector_get (x,8);

	double A4 = gsl_vector_get (x,9);
	double B4 = gsl_vector_get (x,10);
	double C4 = gsl_vector_get (x,11);

	double D = gsl_vector_get (x,12);

	size_t i; 
	
	for (i = 0; i < n; i++)
	{
		double Yi = 0.0;
		//for (int j = 0; j < 3; j++)
	//		Yi = Yi + coef[j];
		gsl_vector_set (f,  i, (Yi - y[i])/sigma[i]);
	}

	return GSL_SUCCESS;
}
*/
/*

int sin_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
	size_t n = ((struct data *)data)->n;
	double *t = ((struct data *)data)->t;
	double *sigma = ((struct data *)data)->sigma;

	double A1 = gsl_vector_get (x,0);
	double B1 = gsl_vector_get (x,1);
	double C1 = gsl_vector_get (x,2);

	double A2 = gsl_vector_get (x,3);
	double B2 = gsl_vector_get (x,4);
	double C2 = gsl_vector_get (x,5);

	double A3 = gsl_vector_get (x,6);
	double B3 = gsl_vector_get (x,7);
	double C3 = gsl_vector_get (x,8);

	double A4 = gsl_vector_get (x,9);
	double B4 = gsl_vector_get (x,10);
	double C4 = gsl_vector_get (x,11);

	double D = gsl_vector_get (x,12);

	size_t i;

	for (i = 0; i < n; i++)
	{
		double s = sigma[i];
		double e1 = sin(B1*t[i] + C1);
		double e2 = A1*cos(B1*t[i] + C1);
		double e3 = cos(B2*t[i] + C2);
		double e4 = -A2*sin(B2*t[i] + C2);
		double e5 = sin(B3*t[i] + C3);
		double e6 = A3*cos(B3*t[i] + C3);
		double e7 = cos(B4*t[i] + C4);
		double e8 = -A4*sin(B4*t[i] + C4);
		gsl_matrix_set (J, i, 0, e1);
		gsl_matrix_set (J, i, 1, e2*t[i]);
		gsl_matrix_set (J, i, 2, e2);
		gsl_matrix_set (J, i, 3, e3); 
		gsl_matrix_set (J, i, 4, e3*t[i]); 
		gsl_matrix_set (J, i, 5, e4); 
		gsl_matrix_set (J, i, 6, e5);
		gsl_matrix_set (J, i, 7, e6*t[i]);
		gsl_matrix_set (J, i, 8, e6);
		gsl_matrix_set (J, i, 9, e7);
		gsl_matrix_set (J, i, 10, e8*t[i]);
		gsl_matrix_set (J, i, 11, e8);
		gsl_matrix_set (J, i, 12, 1.0); 
	}
	return GSL_SUCCESS;
}

int sin_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
	sin_f (x, data, f);
	sin_df (x, data, J);

	return GSL_SUCCESS;
}
*/

//void print_state (size_t iter, gsl_multifit_fdfsolver * s);


void curve_fit(std::vector<double> x_col,
		       std::vector<double> y_col, 
			   std::vector<double> &coefficients)
{
	
	/*
	int len = y_col.size();

	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	unsigned int i, iter = 0;
	const size_t n = len;
	const size_t p = 13;

	gsl_matrix *covar = gsl_matrix_alloc (p,p);
	double t[len], 
		   y[len], 
		   sigma[len];
	struct data d = { n, t, y, sigma };
	gsl_multifit_function_fdf f;
	double x_init[13] = { 0.1, 2.5, 0.0, 0.5, 0.1, 2.5, 0.05, 
	                      0.1, 2.5, 0.0, 2.0, 0.1, 1.5};
	gsl_vector_view x = gsl_vector_view_array (x_init, p);

	f.f = &sin_f;
	f.df = &sin_df;
	f.fdf = &sin_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;

	for (int i = 0; i < n; i++)
	{
		t[i] = x_col[i];
		y[i] = y_col[i];
		sigma[i] = 1.0;
	}

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);

	status = GSL_CONTINUE;

	while ((status == GSL_CONTINUE) && (iter < 500))
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);

		print_state (iter, s);

		if (status)
			break;

		status = gsl_multifit_test_delta (s->dx, s->x, 1.0e-6, 1.0e-6);
	}

	//gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	{
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - p;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));
		std::cout << "chisq/dof = " << pow(chi,2.0)/dof << std::endl;
		std::cout << "A1         = " << FIT(0) << " " << c*ERR(0) << std::endl;
		std::cout << "B1         = " << FIT(1) << " " << c*ERR(1) << std::endl;
		std::cout << "C1         = " << FIT(2) << " " << c*ERR(2) << std::endl;
		std::cout << "A2         = " << FIT(3) << " " << c*ERR(0) << std::endl;
		std::cout << "B2         = " << FIT(4) << " " << c*ERR(1) << std::endl;
		std::cout << "C2         = " << FIT(5) << " " << c*ERR(2) << std::endl;
		std::cout << "A3         = " << FIT(6) << " " << c*ERR(0) << std::endl;
		std::cout << "B3         = " << FIT(7) << " " << c*ERR(1) << std::endl;
		std::cout << "C3         = " << FIT(8) << " " << c*ERR(2) << std::endl;
		std::cout << "A4         = " << FIT(9) << " " << c*ERR(0) << std::endl;
		std::cout << "B4         = " << FIT(10) << " " << c*ERR(1) << std::endl;
		std::cout << "C4         = " << FIT(11) << " " << c*ERR(2) << std::endl;
		std::cout << "D         = " << FIT(12) << " " << c*ERR(3) << std::endl;
	}

	coefficients.push_back(FIT(0));
	coefficients.push_back(FIT(1));
	coefficients.push_back(FIT(2));
	coefficients.push_back(FIT(3));
	coefficients.push_back(FIT(4));
	coefficients.push_back(FIT(5));
	coefficients.push_back(FIT(6));
	coefficients.push_back(FIT(7));
	coefficients.push_back(FIT(8));
	coefficients.push_back(FIT(9));
	coefficients.push_back(FIT(10));
	coefficients.push_back(FIT(11));
	coefficients.push_back(FIT(12));

	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	*/

	//int skip_lines = 0;
	int time_colmn = 0;
	int pres_colmn = 0;

	int p = 30;
	double chisq;

	gsl_multifit_linear_workspace * ws;
	gsl_matrix *X, *cov;
	gsl_vector *y, *c;

	int n = x_col.size();

	X = gsl_matrix_alloc(n,p);
	y = gsl_vector_alloc(n);

	c = gsl_vector_alloc (p);
	cov = gsl_matrix_alloc (p,p);

	for (int i = 0; i < n-1; i++)
	{
		double xi = x_col[i];
		double yi = y_col[i];
		for (int j = 0; j < p; j++)
			gsl_matrix_set (X, i, j, pow(xi,j));
		gsl_vector_set (y, i, yi);
	}

	ws = gsl_multifit_linear_alloc(n,p);
	gsl_multifit_linear(X,y,c,cov, &chisq, ws);

	for (int i = 0; i < p; i++)
	{
		coefficients.push_back(gsl_vector_get(c,i));
	}

	gsl_multifit_linear_free(ws);
	gsl_matrix_free(X);
	gsl_matrix_free(cov);
	gsl_vector_free(y);
	gsl_vector_free(c);
}

/*
void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
	std::cout << "iter: " << iter 
		      << " x = " << gsl_vector_get (s->x, 0)
			  << " " << gsl_vector_get (s->x, 1)
			  << " " << gsl_vector_get (s->x, 2)
			  << " |f(x)| = " << gsl_blas_dnrm2 (s->f)
			  << std::endl;
}
*/
