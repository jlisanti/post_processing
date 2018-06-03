#include <iostream>
#include <vector>
#include <cmath>

#include "datafileinput.h"
#include "trackpeaks.h"
#include "smoothencoder.h"


void track_peaks(DataFileInput &cInput,
		         std::vector<double> &p_vector,
				 std::vector<double> &i_vector,
				 int p_colmn,
				 int i_colmn,
		         int t_colmn,
				 std::vector<double> &vt_max_p,
				 std::vector<double> &vt_min_p,
				 std::vector<double> &vt_max_i,
				 std::vector<double> &vt_min_i,
				 std::vector<double> &vp_max,
				 std::vector<double> &vp_min,
				 std::vector<double> &vi_max,
				 std::vector<double> &vi_min,
				 std::vector<double> &downCrossOver,
				 std::vector<double> &upCrossOver)
{
	bool upFirst = false;
	int crossIndxSize = 0;

	/* determine which is first, up or down */
	if(cInput.table_value(0,downCrossOver[0]) > 
			cInput.table_value(0,upCrossOver[0]))
		upFirst=true;
	else
		upFirst=false;

	if(downCrossOver.size() > upCrossOver.size())
		crossIndxSize = downCrossOver.size();
	else
		crossIndxSize = upCrossOver.size();

	for(int i = 0; i < crossIndxSize; i++)
	{
		int strtIndx = downCrossOver[i];
		int endIndx  = downCrossOver[i+1];

		//if(upFirst==true)
		//if(strtIndx > endIndx)
	//		endIndx = upCrossOver[i+1];

		double pMax = 0.0;
		double pMin = 100.0;
		double iMin = 100.0;
		int pMaxIndx = 0;
		int pMinIndx = 0;
		int iMinIndx = 0;
		/* loop through cycle and find min/max */
		for(int j = strtIndx; j < endIndx; j++)
		{
			if(cInput.table_value(j,p_colmn) > pMax)
			{
				pMax = cInput.table_value(j,p_colmn);
				pMaxIndx = j;
			}
			if(cInput.table_value(j,p_colmn) < pMin)
			{
				pMin = cInput.table_value(j,p_colmn);
				pMinIndx = j;
			}
			if(cInput.table_value(j,i_colmn) < iMin)
			{
				iMin = cInput.table_value(j,i_colmn);
				iMinIndx = j;
			}
		}

		vt_max_p.push_back(cInput.table_value(pMaxIndx,t_colmn));
		vt_min_p.push_back(cInput.table_value(pMinIndx,t_colmn));
		vt_min_i.push_back(cInput.table_value(iMinIndx,t_colmn));
		vp_max.push_back(pMax);
		vp_min.push_back(pMin);
		vi_min.push_back(iMin);
	}

	std::cout << "finished tracking peaks" << std::endl;
}


	// loop through file
	// define start and end regions
	// look for mins and max
	//std::vector<double> encoder_vec = smooth_encoder(cInput,
//			                                         e_colmn);
	/*
	bool found_p_max=false; 
	bool found_p_min=false;

	int search_distance = 50;
	int offset = 1;
	int check = 0;
	int Ntol = 1;

	std::cout << "Tracking peaks..." << std::endl;

	for (int i = 0; i < cInput.file_length(); i++)
	{

		p_vector.push_back(cInput.table_value(i,p_colmn));
		i_vector.push_back(cInput.table_value(i,i_colmn));
	}


	for (int i = search_distance; i < p_vector.size()-search_distance; i++)
	{
		check=0;

		for (int j = offset; j < search_distance; j++)
		{
//			std::cout << p_vector[i+j] << '\t'
//				      << p_vector[i]   << '\t'
//					  << p_vector[i-j] << std::endl;
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
			vt_max_p.push_back(cInput.table_value(i,t_colmn));
			vp_max.push_back(p_vector[i]);
			found_p_max=false;
		}
		if(found_p_min==true)
		{
			vt_min_p.push_back(cInput.table_value(i,t_colmn));
			vp_min.push_back(p_vector[i]);
			found_p_min=false;
		}
	}
	std::cout << "checking for duplicates" << std::endl;
	// Check for duplicates
	// Determine mean spacing between points
	double sum = 0.0;
	for (int i = 0; i < vp_max.size(); i++)
		sum = sum + vt_max_p[i];

	double tol = 0.1;
	
	std::vector<int> replace_index;

	std::cout << "compute mean spacing" << '\t' << vp_max.size() << std::endl;
	double mean_max_p_spacing = sum/double(vp_max.size());
	for (int i = 0; i < vp_max.size()-1; i++)
	{
		if (std::abs(vp_max[i]-vp_max[i+1]) < tol*mean_max_p_spacing)
		{
			replace_index.push_back(i);
			replace_index.push_back(i+1);
		}
	}

	std::cout << "itendified points to remove" << std::endl; 
	*/
	/*
	for (int i = 0; i < replace_index.size()-1; i++)
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
	*/
		/*
	std::cout << "finished tracking peaks" << std::endl;
}
*/
