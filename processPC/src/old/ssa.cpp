#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_math.h>

#include "boost/multi_array.hpp"

#include "datafileinput.h"

#include "ssa.h"


int ssa(DataFileInput &cInput,
		std::vector<double> &time,
		std::vector<double> &pressure,
		int M,
		int kmax,
		int colmn)
{
	//std::vector<double> time;
	//std::vector<double> pressure;

	for (int i = 0; i < cInput.file_length(); i++)
	{
		time.push_back(cInput.table_value(i,0));
		pressure.push_back(cInput.table_value(i,colmn));
	}

	double timestep = time[1] - time[0];
	int N = pressure.size();
	/* Not sure what this value is for */
	//int M = 1000;

	/* Compute eigenvalues */
	std::vector<double> eigenvalues;
	std::vector<double> eigenvalues_unsorted;
	boost::multi_array<double, 2> eigenvectors;
	eigenvectors.resize(boost::extents[M][M]);

	std::cout << " Computing eigenvalues..." << std::endl;
	std::cout << std::endl;

	/* Define MxM covariace matrix */
	gsl_matrix * cij = gsl_matrix_alloc (M,M);

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j  < M; j++)
		{
			gsl_matrix_set (cij,i,j,covariance(i,j,pressure,N));
		}
	}

	gsl_vector_complex *eval = gsl_vector_complex_alloc (M);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (M,M);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (M);
    gsl_eigen_nonsymmv (cij,eval,evec,w);
    gsl_eigen_nonsymmv_free (w);
    gsl_eigen_nonsymmv_sort (eval,evec,GSL_EIGEN_SORT_ABS_ASC);

	/* Store eigenvalues in vector */
	for (int i = 0; i < M; i++)
	{
        gsl_complex eval_i = gsl_vector_complex_get (eval,i);
        eigenvalues.push_back(GSL_REAL(eval_i));
        eigenvalues_unsorted.push_back(GSL_REAL(eval_i));

        // Store eigenvectors
        for (int j = 0; j < M; j++)
        {
            gsl_complex evec_ij = gsl_matrix_complex_get(evec,j,i);
            eigenvectors[i][j] = GSL_REAL(evec_ij);
        }
	}

	gsl_vector_complex_free (eval);
	gsl_matrix_complex_free (evec);

	/* Sort eigenvalues in order of decreasing variance */

	sort(eigenvalues.begin(),eigenvalues.end(),SortType);

	//int kmax;

	/* number of eigenvalues to use */
	//kmax = 10;

	boost::multi_array<double, 2> eigenvectors_sorted;
	eigenvectors_sorted.resize(boost::extents[kmax][M]);

	int step = 0;
	for (int kstep = 0; kstep < kmax; kstep++)
    {
        for (int i = 0; i < eigenvalues.size(); i++)
        {
            if (step >= kmax)
                break;
            if (gsl_fcmp (eigenvalues_unsorted[i],eigenvalues[step],1.e-6)
                    == 0)
            {
                for (int j = 0; j < M; j++)
                    eigenvectors_sorted[step][j] = eigenvectors[i][j];
                step = step + 1;
            }
        }
    }

    // Compute principal components
    boost::multi_array<double, 2> A;
    A.resize(boost::extents[kmax][N]);

    for (int t = 0; t < N; t++)
    {
        for (int k = 0; k < kmax; k++)
        {
            double sum = 0.0;

            for (int j = 0; j < M; j++)
            {
                sum = sum + pressure[t+j]*eigenvectors_sorted[k][j];
            }

            A[k][t] = sum;
        }
    }

    int Lt, Ut;
    double Mt;
    int Nprime = N - M + 1;
    boost::multi_array<double, 1> R;
    R.resize(boost::extents[N]);

    std::cout << "  Reconstructing time series with " << kmax
        << " eigenvectors" << std::endl;
    std::cout << std::endl;

    for (int t = 1; t <= N; t++)
    {
        if ((1 <= t) && (t <= (M-1)))
        {
            Mt = 1.0/t;
            Lt = 1;
            Ut = t;
        }
        else if ((M <= t) && (t <= Nprime))
        {
            Mt = 1.0/M;
            Lt = 1;
            Ut = M;
        }
        else if (((Nprime + 1) <= t) && (t <= N))
        {
            Mt = 1.0/(N-t+1);
            Lt = t - N + M;
            Ut = M;
        }
        else
        {
			std::cout << "   Error!" << std::endl;
        }

        double sumk = 0.0;
        sumk = 0.0;

        for (int k = 0; k < kmax; k++)
        {
            double sumj = 0.0;
            sumj = 0.0;

            for (int j = (Lt-1); j <= (Ut-1); j++)
            {
                sumj = sumj + A[k][t-j-1]*eigenvectors_sorted[k][j];
            }

            sumk = sumk + sumj;
        }

        R[t-1] = (Mt)*(sumk);
    }

	std::ostringstream file_colmn;
	file_colmn << int(colmn);

	std::string eigen_file = "eigen_" + file_colmn.str() + ".dat"; 
	std::string R_file = "R_" + file_colmn.str() + ".dat"; 

	std::ofstream fout7(eigen_file);
    if (!fout7)
    {
		std::cerr << "   Output file coult not be openned for writing!" << std::endl;
    }
    for (int i = 0; i < eigenvalues.size(); i++)
    {
        fout7 << i << '\t';
        fout7 << eigenvalues[i] << std::endl;
    }

	// Write reconstructed time series to file

	//std::string outfileR;
    //outfileR = std::string("R-") + "test.dat";

	//std::ofstream foutR(outfileR);
	std::ofstream foutR(R_file);
    if (!foutR)
    {
		std::cerr << "   Output file count not be opened for writing!" << std::endl;
    }
    for (int t = 0; t < N; t++)
    {
        foutR << t*timestep << '\t';
        foutR << R[t] << std::endl;
    }

	return 0;
}

bool SortType(double i, double j) { return i > j; }

inline double covariance(int i, int j, std::vector<double> X, int N)
{
    double sum = 0.0;
    for (int t = 0; t < (N-abs(i-j)); t++)
    {
        sum = sum + (X[t]*X[t + abs(i-j)]);
    }
    return ((1.0/(N - abs(i-j))) * sum);
}

inline double mean(std::vector<double> X)
{
    double data[X.size()];
    copy(X.begin(),X.end(),data);
    return gsl_stats_mean(data,1,X.size());
}

inline double variance(std::vector<double> X)
{
    double data[X.size()];
    copy(X.begin(),X.end(),data);
    return gsl_stats_variance(data,1,X.size());
}

inline double sd(std::vector<double> X)
{
    double data[X.size()];
    copy(X.begin(),X.end(),data);
    return gsl_stats_sd(data,1,X.size());
}
