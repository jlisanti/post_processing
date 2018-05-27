#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "datafileinput.h"
#include "smoothencoder.h"

std::vector<double> smooth_encoder(DataFileInput &cInput,
		                           int enco_colmn)
{
	//DataFileInput cInput(input_file,skip_lines);

	double gain = 5.0;
	std::vector<int> end_index;
	std::vector<int> start_index;
	std::vector<double> encoder_value;

	start_index.push_back(0);
	
	int count1 = 0;
	int count2 = 0;

	for (int i = 0; i < cInput.file_length()-11; i++)	
	{

		if(cInput.table_value(i,enco_colmn) 
				> cInput.table_value(i+10,enco_colmn))
			count1++;
		else if(cInput.table_value(i,enco_colmn) 
				< cInput.table_value(i+10,enco_colmn))
			count2++;

	}

	if (count1 < count2)
	{
		for (int i = 0; i < cInput.file_length()-1; i++)
		{
			if(cInput.table_value(i,enco_colmn) 
					> gain*cInput.table_value(i+1,enco_colmn))
			{
				end_index.push_back(i);
				start_index.push_back(i+1);
			}
		}

		end_index.push_back(cInput.file_length()-1);

		double start_encoder_value = cInput.table_value(0,enco_colmn);
		double end_encoder_value 
			= cInput.table_value(cInput.file_length()-1,enco_colmn);

		for (int i = 0; i < end_index.size(); i++)
		{
			for (int j = start_index[i]-1; j < end_index[i]; j++)
			{
				if (i == 0)
				{
					encoder_value.push_back(((5.0-start_encoder_value)
							/(end_index[i]-start_index[i])
							*(j-start_index[i]))+start_encoder_value);
				}
				else if (i == end_index.size()-1)
				{
					encoder_value.push_back(((end_encoder_value)
							/(end_index[i]-start_index[i])
							*(j-start_index[i])));
				}
				else
				{
					encoder_value.push_back((5.0/(end_index[i]-start_index[i])*
						(j-start_index[i])));
				}
			}
		}

	}
	else 
	{
		for (int i = 0; i < cInput.file_length()-1; i++)
		{
			if(cInput.table_value(i+1,enco_colmn) 
					> gain*cInput.table_value(i,enco_colmn))
			{
				end_index.push_back(i);
				start_index.push_back(i+1);
			}
		}

		end_index.push_back(cInput.file_length()-1);

		double start_encoder_value = cInput.table_value(0,enco_colmn);
		double end_encoder_value 
			= cInput.table_value(cInput.file_length()-1,enco_colmn);

		for (int i = 0; i < end_index.size(); i++)
		{
			for (int j = start_index[i]-1; j < end_index[i]; j++)
			{
				if (i == 0)
				{
					encoder_value.push_back(((start_encoder_value)
							/(end_index[i]+start_index[i])
							*(-j+start_index[i]))+start_encoder_value);
				}
				else if (i == end_index.size()-1)
				{
					encoder_value.push_back(((5.0-end_encoder_value)
							/(end_index[i]-start_index[i])
							*(-j+start_index[i]))+5.0);
				}
				else
				{
					encoder_value.push_back((5.0/(end_index[i]-start_index[i])*
						(-j+start_index[i]))+5.0);
				}

				encoder_value[j] = encoder_value[j]*(-1.0) + 5.0;
			}
		}
		/*
		for (int i = 0; i  < encoder_value.size(); i++)
		{
			double old_value = encoder_value[i];
			std::cout << old_value << '\t' << encoder_value[i] << std::endl;
		}
		*/
	}
	return encoder_value;
}
