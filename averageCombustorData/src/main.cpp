#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "datafileinput.h"

int main(int argc, char*argv[])
{
	/* read processed input file */
    if (argc <= 1)
    {
		std::cout << " Error : no input file" << std::endl;;
        exit(1);
    }

	int skip_lines = 1;
    DataFileInput cInput(argv[1],skip_lines);
	
	int freqColmn = -1;
	int fuelColmn = -1;

	std::cout << " Enter freq column index: " << std::endl;
	std::cout << " -> ";
	std::cin >> freqColmn;

	std::cout << " Enter fuel column index: " << std::endl;
	std::cout << " -> ";
	std::cin >> freqColmn;

	double freq = -1;
	double fuel = -1;
	int count = 0;

	std::vector<std::vector<double>> averaged_data;

	for (int i = 0; i < cInput.size(); i++)
	{
		bool found = false;
		std::vector<double> tmp;

		for (int j = 0; j < cInput.number_cols(); j++)
		{
			double tableValue = cInput.table_value(i,j);
			std::cout << i << '\t' << j << '\t' <<  tableValue << std::endl;
			tmp.push_back(cInput.table_value(i,j));
		}

		for (int j = 0; j < count; j++) 
		{
		    if((tmp[freqColmn]==averaged_data[j][freqColmn])
					&& (std::abs(tmp[fuelColmn]-averaged_data[j][fuelColmn]) < 5.0))
			    found=true;
		}

		if(found==false)
		{
		    std::vector<double> sums(0.0,cInput.number_cols());
		    for(int j = 0; j < cInput.size(); j++)
		    {
		        if((tmp[freqColmn]==cInput.table_value(j,freqColmn))
					&& (std::abs(tmp[fuelColmn]-cInput.table_value(j,fuelColmn))<5.0))
			    {
			        for(int k = 0; k < cInput.number_cols(); k++)
			        {
					   sums[k] = sums[k] + tmp[k]; 
			        }
			    }
		    }
		    for(int k = 0; k < cInput.number_cols(); k++)
			{
				sums[k] = sums[k]/double(sums.size()); 
			}

		    count++;
		}
	}

	std::ofstream fout("combustor_average_data.dat");
	for (int i = 0; i < count; i++)
	{
		for (int j = 0; j < cInput.number_cols(); j++)
		{
			fout << averaged_data[i][j] << '\t';
		}
		fout << std::endl;
	}


	return 0;
}
