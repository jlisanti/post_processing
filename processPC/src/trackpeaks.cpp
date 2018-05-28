#include <iostream>
#include <vector>
#include <cmath>

#include "datafileinput.h"
#include "trackpeaks.h"
#include "smoothencoder.h"


void track_peaks(DataFileInput &cInput,
		         std::vector<double> &p_vector,
				 std::vector<double> &e_vector,
				 std::vector<double> &i_vector,
		         int e_colmn,
				 int p_colmn,
				 int i_colmn,
		         int t_colmn,
				 std::vector<double> &ve_max_p,
				 std::vector<double> &ve_min_p,
				 std::vector<double> &ve_max_i,
				 std::vector<double> &ve_min_i,
				 std::vector<double> &vt_max_p,
				 std::vector<double> &vt_min_p,
				 std::vector<double> &vt_max_i,
				 std::vector<double> &vt_min_i,
				 std::vector<double> &vp_max,
				 std::vector<double> &vp_min,
				 std::vector<double> &vi_max,
				 std::vector<double> &vi_min)
{

	// loop through file
	// define start and end regions
	// look for mins and max
	//std::vector<double> encoder_vec = smooth_encoder(cInput,
//			                                         e_colmn);
	bool found_p_max=false; 
	bool found_p_min=false;

	int search_distance = 50;
	int offset = 1;
	int check = 0;
	int Ntol = 1;


	for (int i = search_distance; i < p_vector.size()-search_distance; i++)
	{
		check=0;

		for (int j = offset; j < search_distance; j++)
		{
			std::cout << p_vector[i+j] << '\t'
				      << p_vector[i]   << '\t'
					  << p_vector[i-j] << std::endl;
			if (p_vector[i+j] > p_vector[i])
				check++;

			if (p_vector[i-j] > p_vector[i])
				check++;
		}
		if(check<Ntol)
		{
			found_p_max=true;
			check=0;
		}

		check=0;

		for (int j = offset; j < search_distance; j++)
		{
			if (p_vector[i+j] < p_vector[i])
				check++;

			if (p_vector[i-j] < p_vector[i])
				check++;
		}
		if(check<Ntol)
		{
			found_p_min=true;
			check=0;
		}
		
		if(found_p_max==true)
		{
			ve_max_p.push_back(e_vector[i]);
			vt_max_p.push_back(cInput.table_value(i,t_colmn));
			vp_max.push_back(p_vector[i]);
			found_p_max=false;
		}
		if(found_p_min==true)
		{
			ve_min_p.push_back(e_vector[i]);
			vt_min_p.push_back(cInput.table_value(i,t_colmn));
			vp_min.push_back(p_vector[i]);
			found_p_min=false;
		}
	}
	// Check for duplicates
	// Determine mean spacing between points
	double sum = 0.0;
	for (int i = 0; i < vp_max.size(); i++)
		sum = sum + vt_max_p[i];

	double tol = 0.1;
	
	std::vector<int> replace_index;

	double mean_max_p_spacing = sum/double(vp_max.size());
	for (int i = 0; i < vp_max.size()-1; i++)
		if (std::abs(vp_max[i]-vp_max[i+1]) < tol*mean_max_p_spacing)
		{
			replace_index.push_back(i);
			replace_index.push_back(i+1);
		}

	for (int i = 0; i < replace_index.size(); i++)
	{
		double sum = vp_max[replace_index[i]];
		for (int j = i; j < replace_index.size()-i; j++)
		{
			if(j==i)
				sum = sum + vp_max[replace_index[i+1]];
			else if(replace_index[i]==replace_index[j])
				sum = sum + vp_max[replace_index[i+1]];
		}
	}
}

		/*
		if((cInput.table_value(i,p_colmn) > cInput.table_value(i-4,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) > cInput.table_value(i-5,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) > cInput.table_value(i-6,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) > cInput.table_value(i-7,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) > cInput.table_value(i-8,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) > cInput.table_value(i+4,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) > cInput.table_value(i+5,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) > cInput.table_value(i+6,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) > cInput.table_value(i+7,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) > cInput.table_value(i+8,p_colmn))
		  )
		{
			found_p_max=true;
		}

		if((cInput.table_value(i,p_colmn) < cInput.table_value(i-4,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) < cInput.table_value(i-5,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) < cInput.table_value(i-6,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) < cInput.table_value(i-7,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) < cInput.table_value(i-8,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) < cInput.table_value(i+4,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) < cInput.table_value(i+5,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) < cInput.table_value(i+6,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) < cInput.table_value(i+7,p_colmn)) &&
		   (cInput.table_value(i,p_colmn) < cInput.table_value(i+8,p_colmn))
		  )
		{
			found_p_min=true;
		}
		*/

		/*
		std::cout << cInput.table_value(i-10,p_colmn) << '\t'
			      << cInput.table_value(i-5,p_colmn)  << '\t'
				  << cInput.table_value(i-4,p_colmn)  << '\t'
				  << cInput.table_value(i-2,p_colmn)  << '\t'
				  */

	/*
	double e_max_p = 0.0;
	double e_min_p = 0.0;
	double p_max = 0.0; 
	double p_min = 1.0e10;
	double e_max_i = 0.0;
	double e_min_i = 0.0;
	double i_max = 0.0;
	double i_min = 1.0e10;

	double t_max_p = 0.0;
	double t_min_p = 0.0;
	double t_max_i = 0.0;
	double t_min_i = 0.0;

	bool new_cycle=false;
	bool rising=false;

	int peak_count = 0;
	int check_count = 0;

	int peak_check = 10;
	*/

	/*
	int max_peak_count = 0;
	int min_peak_count = 0;

	int max_check_count = 0;
	int min_check_count = 0;
	*/

	/* dangerous way of determining encoder direction */

	/*
	if(encoder_vec[10] < encoder_vec[11])
		rising=true;
		*/

		/*
		if(cInput.table_value(i,p_colmn) > p_max)
		{
			p_max = cInput.table_value(i,p_colmn);
			e_max_p = encoder_vec[i];
			t_max_p = cInput.table_value(i,t_colmn);
		}
		if(cInput.table_value(i,p_colmn) < p_min)
		{
			p_min = cInput.table_value(i,p_colmn);
			e_min_p = encoder_vec[i];
			t_min_p = cInput.table_value(i,t_colmn);
		}

		if(cInput.table_value(i,i_colmn) > i_max)
		{
			i_max = cInput.table_value(i,i_colmn);
			e_max_i = encoder_vec[i];
			t_max_i = cInput.table_value(i,t_colmn);
		}
		if(cInput.table_value(i,i_colmn) < i_min)
		{
			i_min = cInput.table_value(i,i_colmn);
			e_min_i = encoder_vec[i];
			t_min_i = cInput.table_value(i,t_colmn);
		}
		*/

		/*
		if(peak_count==2)
			new_cycle=true;
			*/

		/*
		if(new_cycle==true)
		{
			ve_max_i.push_back(e_max_i);
			ve_min_i.push_back(e_min_i);

			ve_max_p.push_back(e_max_p);
			ve_min_p.push_back(e_min_p);

			vt_max_i.push_back(t_max_i);
			vt_min_i.push_back(t_min_i);

			vt_max_p.push_back(t_max_p);
			vt_min_p.push_back(t_min_p);

			vp_max.push_back(p_max);
			vp_min.push_back(p_min);

			vi_max.push_back(i_max);
			vi_min.push_back(i_min);
			
			p_max = 0.0;
			p_min = 1.0e10;
			i_max = 0.0;
			i_min = 1.0e10;
			peak_count = 0;
			check_count = 0;
			max_peak_count=0;
			min_peak_count=0;
			max_check_count = 0;
			min_check_count = 0;

			new_cycle=false;
			std::cout << "***new cycle***" << std::endl;
		}
			*/
	
/*

		int dir_1 = 0;
		int dir_2 = 0;
		int dir_3 = 0;
		int dir_4 = 0;
		int dir_5 = 0;

		if(cInput.table_value(i,p_colmn) > cInput.table_value(i+1,p_colmn))
			dir_1 = 1;
		if(cInput.table_value(i+1,p_colmn) > cInput.table_value(i+2,p_colmn))
			dir_2 = 1;
		if(cInput.table_value(i+2,p_colmn) > cInput.table_value(i+3,p_colmn))
			dir_3 = 1;
		if(cInput.table_value(i+3,p_colmn) > cInput.table_value(i+4,p_colmn))
			dir_4 = 1;
		if(cInput.table_value(i+4,p_colmn) > cInput.table_value(i+5,p_colmn))
			dir_5 = 1;
		if((dir_1==dir_2) &&
		   (dir_2==dir_3) &&
		   (dir_3==dir_4) &&
		   (dir_4==dir_5) &&
		   (dir_5==dir_1))
			check_count++;

        if(check_count>=peak_check)
		{
			peak_count++;
			check_count = 0;
		}
		
		std::cout << i << '\t' << check_count << '\t' << peak_count << std::endl;
		*/

		/*
		if(cInput.table_value(i,p_colmn) < p_max)
			max_check_count++;
		else
			if(max_check_count > 0)
				max_check_count--;

		if(cInput.table_value(i,p_colmn) > p_min)
			min_check_count++;
		else
			if(min_check_count > 0)
				min_check_count--;
				*/

		/*
        if(max_check_count==peak_check)
			max_peak_count++;

        if(min_check_count==peak_check)
			min_peak_count++;
			*/


		/*
		if(rising==true)
			//if((encoder_vec[i] > 5.0*encoder_vec[i+1]) || ((encoder_vec[i] < 2.5) && (encoder_vec[i+1] > 2.5)))
			if((en
			{
				if((max_peak_count==2) && (min_peak_count==2))
					new_cycle=true;
				else
				{
					p_max = 0.0;
					p_min = 1.0e10;
					i_max = 0.0;
					i_min = 1.0e10;
				}
			}
			else
				new_cycle=false;
		else
			if((encoder_vec[i] < 5.0*encoder_vec[i+1]) || ((encoder_vec[i] > 2.5) && (encoder_vec[i+1] < 2.5)))
			{
				new_cycle=true;
			}
			else
				new_cycle=false; 
				*/
