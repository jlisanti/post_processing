#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "datafileinput.h"
#include "smoothencoder.h"
#include "averagecycle.h"

void average_cycle_active(DataFileInput &cInput, 
						  std::vector<double> &time_ssa,
						  std::vector<double> &pressure_ssa,
						  std::vector<double> &ion_ssa,
						  int encoder_colmn,
				          int pressur_colmn,
				          double window,
				          std::vector<double> &e_vector,
				          std::vector<double> &p_vector,
						  std::vector<double> &i_vector,
						  std::vector<double> &pressure_scatter,
						  std::vector<double> &encoder_scatter,
						  std::vector<double> &ion_scatter)
{
	std::vector<double> encoder_vec = smooth_encoder(cInput,
			                                         encoder_colmn);
	//double window = 2.5/50.0;
	//int steps = cInput.column_max(encoder_colmn);
	int steps = ceil(2.5/window);
	int lines = cInput.file_length(); 

	double step_min = 0.0;
	double step_max = window;

	//std::vector<double> p_vector;
	//std::vector<double> e_vector;

	std::string prompt;
	std::cout << " Print from valve openning? (Y/n)" << std::endl;
	std::cin >> prompt;
	if((prompt=="Y") ||
	   (prompt=="y") ||
	   (prompt=="yes") ||
	   (prompt=="Yes") ||
	   (prompt=="YES"))
	{
		double r_valve = 0.0;
		double h_valve = 0.0;
		double r_inlet  = 0.0;
		std::cout << " Enter valve radius" << std::endl;
		std::cin >> r_valve;
		std::cout << " Enter valve length" << std::endl;
		std::cin >> h_valve;
		std::cout << " Enter inlet radius" << std::endl;
		std::cin >> r_inlet;

		double h = sqrt(pow(h_valve/2.0,2) + pow(r_inlet,2));
		double alpha = 2.0*asin(r_inlet/h);
		double encoder_close = alpha*(2.5/(2.0*M_PI));
		double encoder_open = (2.0*M_PI-alpha)*(2.5/(2.0*M_PI));
		int open_index = 0;
		std::vector<double> e_tmp;
		std::vector<double> index;
		std::cout << encoder_close << '\t' << encoder_open << std::endl;

		/*
		for (int i = 0; i < e_vector.size(); i++)
			if(e_vector[i] >= encoder_open)
			{
				open_index = i;
				break;
			}
		for (int i = 0; i < e_vector.size(); i++)
		{
			if(e_vector[i] <= encoder_open)
			{
				e_vector[i] = e_vector[i] + e_vector[e_vector.size()-1];
			}
		}
		double offset = e_vector[open_index];
		for (int i = 0; i < e_vector.size(); i++)
		{
			e_vector[i] = e_vector[i] - offset;
		}
		*/

	for (int i = 0; i < steps; i++)
	{
		double p_sum = 0.0;
		double i_sum = 0.0;
		double e_sum = 0.0;
		double count = 0.0;
		double encoder = 0.0;

		for (int j = 0; j < lines; j++)
		{
			if(encoder_vec[j] > 2.5)
			{
				encoder = encoder_vec[j] - 2.5;	
				if((encoder >= encoder_open) || (encoder <= encoder_close))
				{
					if(encoder <= encoder_close)
					{
						encoder = encoder + (2.5-encoder_open);
					}
					else
					{
						encoder = encoder - encoder_open;
					}
				}
				else
				{
					encoder = encoder + (2.5-encoder_open);
				}
			}
			else
			{
				encoder = encoder_vec[j];
				if((encoder >= encoder_open) || (encoder <= encoder_close))
				{
					if(encoder <= encoder_close)
					{
						encoder = encoder + (2.5-encoder_open);
					}
					else
					{
						encoder = encoder - encoder_open;
					}
				}
				else
				{
					encoder = encoder + (2.5-encoder_open);
				}
			}

			if((encoder > step_min) && (encoder < step_max))
			{
				//pressure_scatter.
			    //		push_back(cInput.table_value(j,pressur_colmn));
				pressure_scatter.push_back(pressure_ssa[j]);
				ion_scatter.push_back(ion_ssa[j]);
				encoder_scatter.push_back(encoder);
				//p_sum = p_sum + cInput.table_value(j,pressur_colmn);
				p_sum = p_sum + pressure_ssa[j];
				i_sum = i_sum + ion_ssa[j];
				e_sum = e_sum + encoder;
				count++;
			}

		}
		p_vector.push_back(p_sum/double(count));
		i_vector.push_back(i_sum/double(count));
		e_vector.push_back(e_sum/double(count));
		step_min = step_min + window; 
		step_max = step_max + window;
	}
	}
	/*
	bool check = false;
	while (check==false)
	{
		check = true;
		for (int i = 0; i < e_vector.size()-1; i++)
		{
			if(e_vector[i] > e_vector[i+1])
			{
				double tmp_a = e_vector[i];
				double tmp_b = e_vector[i+1];
				e_vector[i] = tmp_b;
				e_vector[i+1] = tmp_a;

				tmp_a = p_vector[i];
				tmp_b = p_vector[i+1];
				p_vector[i] = tmp_b;
				p_vector[i+1] = tmp_a;

				tmp_a = i_vector[i];
				tmp_b = i_vector[i+1];
				i_vector[i] = tmp_b;
				i_vector[i+1] = tmp_a;

				check = false;
			}
		}
	}
	*/

	
	//std::ofstream fout("test_out.dat");
	//for (int i = 0; i < p_vector.size(); i++)
	//{
    //	fout << e_vector[i] << '\t' << p_vector[i] << std::endl;
	//}

	/*
	std::string prompt;
	std::cout << " Print from valve openning? (Y/n)" << std::endl;
	std::cin >> prompt;
	if((prompt=="Y") ||
	   (prompt=="y") ||
	   (prompt=="yes") ||
	   (prompt=="Yes") ||
	   (prompt=="YES"))
	{
		double r_valve = 0.0;
		double h_valve = 0.0;
		double r_inlet  = 0.0;
		std::cout << " Enter valve radius" << std::endl;
		std::cin >> r_valve;
		std::cout << " Enter valve length" << std::endl;
		std::cin >> h_valve;
		std::cout << " Enter inlet radius" << std::endl;
		std::cin >> r_inlet;

		double h = sqrt(pow(h_valve/2.0,2) + pow(r_inlet,2));
		double alpha = 2.0*asin(r_inlet/h);
		double encoder_open = alpha*(2.5/(2.0*M_PI));
		int open_index = 0;
		std::vector<double> e_tmp;
		std::vector<double> index;

		for (int i = 0; i < e_vector.size(); i++)
			if(e_vector[i] >= encoder_open)
			{
				open_index = i;
				break;
			}
		for (int i = 0; i < e_vector.size(); i++)
		{
			if(e_vector[i] <= encoder_open)
			{
				e_vector[i] = e_vector[i] + e_vector[e_vector.size()-1];
			}
		}
		double offset = e_vector[open_index];
		for (int i = 0; i < e_vector.size(); i++)
		{
			e_vector[i] = e_vector[i] - offset;
		}
	}
	bool check = false;
	while (check==false)
	{
		check = true;
		for (int i = 0; i < e_vector.size()-1; i++)
		{
			if(e_vector[i] > e_vector[i+1])
			{
				double tmp_a = e_vector[i];
				double tmp_b = e_vector[i+1];
				e_vector[i] = tmp_b;
				e_vector[i+1] = tmp_a;

				tmp_a = p_vector[i];
				tmp_b = p_vector[i+1];
				p_vector[i] = tmp_b;
				p_vector[i+1] = tmp_a;

				tmp_a = i_vector[i];
				tmp_b = i_vector[i+1];
				i_vector[i] = tmp_b;
				i_vector[i+1] = tmp_a;

				check = false;
			}
		}
	}
*/
}


void average_cycle_passive(DataFileInput &cInput, 
                   int pressur_colmn,
                   int ion_colmn,
                   double window,
                   std::vector<double> &position_average,
                   std::vector<double> &pressure_average,
                   std::vector<double> &ion_average,
                   std::vector<double> &position_scatter,
                   std::vector<double> &pressure_scatter,
                   std::vector<double> &ion_scatter,
                   double static_pressure)
{

    // Compute mean of pressure signal
    // Perhaps better to record static pressure - 
    //    then the crossover value can be this. 

    // Proxy
    //double static_pressure = 1.01325;
    /*
    int n = cInput.file_length();
    double * data = new double[n];

    for (int i = 0; i < n; i++)
        data[i] = cInput.table_value(i,pressur_colmn);

    double mean = gsl_stats_mean(data, 1, n);

    std::cout << "mean: " <<  mean << std::endl;

    exit(1);
    */

    // Isolate cycle stop and start positions from noise

    if ((static_pressure < 0.9) || (static_pressure > 1.1))
        static_pressure = 1.01325;

    int buffer = 4;

    std::vector<int> down_cross_over;
    std::vector<int> up_cross_over;

	for (int i = buffer; i < cInput.file_length()-buffer; i ++)
    {
        if(cInput.table_value(i,pressur_colmn)   > static_pressure &&
           cInput.table_value(i-1,pressur_colmn) > static_pressure &&
           cInput.table_value(i-2,pressur_colmn) > static_pressure &&
           cInput.table_value(i+1,pressur_colmn) < static_pressure &&
           cInput.table_value(i+2,pressur_colmn) < static_pressure)
        {
            down_cross_over.push_back(i);
        }
        if(cInput.table_value(i,pressur_colmn)   < static_pressure &&
           cInput.table_value(i-1,pressur_colmn) < static_pressure &&
           cInput.table_value(i-2,pressur_colmn) < static_pressure &&
           cInput.table_value(i+1,pressur_colmn) > static_pressure &&
           cInput.table_value(i+2,pressur_colmn) > static_pressure)
        {
            up_cross_over.push_back(i);
        }
    }

    /*
    std::ofstream fout1("down_cross_over.dat");
    for (int i = 0; i < down_cross_over.size(); i++)
    {
        fout1 << cInput.table_value(down_cross_over[i],0) << '\t' 
              << cInput.table_value(down_cross_over[i],pressur_colmn) << std::endl; 
    }

    std::ofstream fout2("up_cross_over.dat");
    for (int i = 0; i < up_cross_over.size(); i++)
    {
        fout2 << cInput.table_value(up_cross_over[i],0) << '\t' 
              << cInput.table_value(up_cross_over[i],pressur_colmn) << std::endl; 
    }
    */

    // Map pressure signal from each cycle into single space 
   
    // start from first up crossing 
    for (int i = 0; i < up_cross_over.size()-1; i++)
    {
        for (int j = up_cross_over[i]; j < up_cross_over[i+1]; j++)
        {
            // time mapping
            double dt = cInput.table_value(up_cross_over[i+1],0)
                       -cInput.table_value(up_cross_over[i],0);
            double t_initial = cInput.table_value(up_cross_over[i],0);
            double pos = (cInput.table_value(j,0)
                * (1.0/dt)) - (t_initial/dt);
            position_scatter.push_back(pos);
            pressure_scatter.push_back(cInput.table_value(j,pressur_colmn));
            ion_scatter.push_back(cInput.table_value(j,ion_colmn));
        }
    }
 std::ofstream fout1("scatter.dat");
    for (int i = 0; i < position_scatter.size(); i++)
    {
        fout1 << position_scatter[i] << '\t'
              << pressure_scatter[i] << '\t'
              << ion_scatter[i]      << std::endl;
    }

    // Average cycles into single curve

    int steps = ceil(1.0/window);

    double step_min = 0.0;
    double step_max = window;

    for (int i = 0; i < steps; i++)
    {
        double count = 0;
        double press_sum = 0;
        double posit_sum = 0;
        double ion_sum = 0;

        for (int j = 0; j < position_scatter.size(); j++)
        {
            if((position_scatter[j] > step_min) && (position_scatter[j] < step_max))
            {
                press_sum = press_sum + pressure_scatter[j];
                posit_sum = posit_sum + position_scatter[j];
                ion_sum   = ion_sum   + ion_scatter[j];
                count++;
            }
        }
        position_average.push_back(posit_sum/double(count));
        pressure_average.push_back(press_sum/double(count));
        ion_average.push_back(ion_sum/double(count));
        step_min = step_min + window;
        step_max = step_max + window;
    }

    std::ofstream fout2("average.dat");
    for (int i = 0; i < position_average.size(); i++)
    {
        fout2 << position_average[i] << '\t'
              << pressure_average[i] << '\t'
              << ion_average[i]      << std::endl;
    }
}


