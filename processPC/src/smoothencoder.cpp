#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "datafileinput.h"
#include "smoothencoder.h"
#include "ssa.h"

std::vector<double> smooth_encoder(DataFileInput &cInput,
		                           int enco_colmn)
{
	double gain = 4.0;
	std::vector<int> end_index;
	std::vector<int> start_index;
	std::vector<double> encoder_value;
    std::vector<double> encoder_smooth;

    std::vector<double> timeEncoder;

	start_index.push_back(0);
	
	int count1 = 0;
	int count2 = 0;

	ssa(cInput,
        timeEncoder,
        encoder_smooth,
        1,
        1,
        enco_colmn);

    /* determine encoder direction */

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
		std::cout << "Encoder moving forwards" << std::endl;
		for (int i = 0; i < encoder_smooth.size()-2; i++)
		{
			if(encoder_smooth[i] > gain*encoder_smooth[i+2])
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
		std::cout << "Encoder moving backwards" << std::endl;
		for (int i = 0; i < encoder_smooth.size()-1; i++)
		{
			if(encoder_smooth[i+1] > gain*encoder_smooth[i])
			{
				end_index.push_back(i);
				start_index.push_back(i+1);
			}
		}
       
        std::cout << cInput.file_length() << '\t' << encoder_smooth.size() << std::endl;
		
		end_index.push_back(encoder_smooth.size()-1);

		double start_encoder_value = cInput.table_value(0,enco_colmn);
		double end_encoder_value 
			= cInput.table_value(cInput.file_length()-1,enco_colmn);


		for (int i = 0; i < end_index.size(); i++)
		{

			std::cout << "start index: " << start_index[i] 
                      << " end index:  " << end_index[i] << std::endl;


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
			    if(j<0)
					j=0;	
				encoder_value[j] = encoder_value[j]*(-1.0) + 5.0;
			}
		}
    }

/*
    std::ofstream foutEnc("Encoderaverage.dat");
    for (int i = 0; i < encoder_value.size(); i++)
    {
        foutEnc << timeEncoder[i]      << '\t'
              << encoder_value[i]      << std::endl;
    }

    exit(1);
*/


/*
	bool negs = true;

	while(negs==true)
	{
		negs=false;

		for(int i = 0; i < encoder_value.size(); i++)
		{
			if(encoder_value[i] < 0.0)
			{
				if(i==0)
					encoder_value[i] = encoder_value[i+1];
                else if(i==encoder_value.size()-1)
					encoder_value[i] = encoder_value[i-1];
				else
                    encoder_value[i] = (encoder_value[i-1] + encoder_value[i+1])/2.0; 
					
				negs=true;
			}	
		}	
	} 
*/
	return encoder_value;
}
