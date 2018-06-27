#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "datafileinput.h"
#include "smoothencoder.h"
#include "averagecycle.h"

void average_cycle_active(DataFileInput &cInput, 
		                  double window,
						  double step,
						  std::vector<double> &time_ssa,
						  std::vector<double> &pressure_ssa,
						  std::vector<double> &ion_ssa,
						  int encoder_colmn,
				          int pressur_colmn,
				          std::vector<double> &e_vector,
				          std::vector<double> &p_vector,
						  std::vector<double> &i_vector,
						  std::vector<double> &pressure_scatter,
						  std::vector<double> &encoder_scatter,
						  std::vector<double> &ion_scatter,
						  std::vector<double> &up_cross_over,
						  std::vector<double> &down_cross_over,
						  std::string printFromStart)
{
	std::vector<double> encoder_vec = smooth_encoder(cInput,
			                                         encoder_colmn);
	//double window = 2.5/50.0;
	//int steps = cInput.column_max(encoder_colmn);
	int steps = ceil(2.5/window);
	int lines = cInput.file_length(); 

	double step_min = 0.0;
	double step_max = window;

	if(printFromStart=="true")
	{
		double r_valve = 8.65/1000.0;
		double h_valve = 40.0/1000.0;
		double r_inlet = 8.65/1000.0;

		double h = sqrt(pow(h_valve/2.0,2) + pow(r_inlet,2));
		double alpha = 2.0*asin(r_inlet/h);
		double encoder_close = alpha*(2.5/(2.0*M_PI));
		double encoder_open = (2.0*M_PI-alpha)*(2.5/(2.0*M_PI));
		int open_index = 0;
		std::vector<double> e_tmp;
		std::vector<double> index;
		std::cout << encoder_close << '\t' << encoder_open << std::endl;

		for (int i = 0; i < encoder_vec.size(); i++)
			if(encoder_vec[i] > 2.5)
				encoder_vec[i] = encoder_vec[i] - 2.5;

		/* loop through encoder file */
		for (int i = 0; i < encoder_vec.size(); i++)
		{
			/* check for cross opening value */
			if((encoder_vec[i]   <= encoder_open) &&
			   (encoder_vec[i-1] < encoder_open)  &&
			   (encoder_vec[i-2] < encoder_open)  &&
			   (encoder_vec[i-3] < encoder_open)  && 
			   (encoder_vec[i-4] < encoder_open)  &&
			   (encoder_vec[i+1] > encoder_open)  &&
			   (encoder_vec[i+2] > encoder_open)  &&
			   (encoder_vec[i+3] > encoder_open)  &&
			   (encoder_vec[i+4] > encoder_open))
				up_cross_over.push_back(i);
			if((encoder_vec[i]   <= encoder_close) &&
			   (encoder_vec[i-1] < encoder_close)  &&
			   (encoder_vec[i-2] < encoder_close)  &&
			   (encoder_vec[i-3] < encoder_close)  &&
			   (encoder_vec[i-4] < encoder_close)  &&
			   (encoder_vec[i+1] > encoder_close)  &&
			   (encoder_vec[i+2] > encoder_close)  &&
			   (encoder_vec[i+3] > encoder_close)  &&
			   (encoder_vec[i+4] > encoder_close))
				down_cross_over.push_back(i);
			
		}

    std::ofstream fout_cross("crossover.dat");
	for (int i = 0; i < up_cross_over.size(); i++)
	{
		fout_cross << cInput.table_value(down_cross_over[i],0) << '\t'
			       << cInput.table_value(down_cross_over[i],3) << std::endl;
	}
	exit(1);
	}


		/*

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

		//	if((encoder > step_min) && (encoder < step_max))
		//	{
				//pressure_scatter.
			    //		push_back(cInput.table_value(j,pressur_colmn));
				pressure_scatter.push_back(pressure_ssa[j]);
				ion_scatter.push_back(ion_ssa[j]);
				encoder_scatter.push_back(encoder);
				//p_sum = p_sum + cInput.table_value(j,pressur_colmn);
				//p_sum = p_sum + pressure_ssa[j];
				//i_sum = i_sum + ion_ssa[j];
				//e_sum = e_sum + encoder;
		//		count++;
		//	}

		}
		//p_vector.push_back(p_sum/double(count));
		//i_vector.push_back(i_sum/double(count));
		//e_vector.push_back(e_sum/double(count));
		//step_min = step_min + window; 
		//step_max = step_max + window;
	}

	    std::ofstream fout1("scatter.dat",std::ios_base::app);
        for (int i = 0; i < encoder_scatter.size(); i++)
        {
            fout1 << encoder_scatter[i]    << '\t'
                  << pressure_scatter[i] << '\t'
                  << ion_scatter[i]      << std::endl;
        }
		exit(1);
	}
	*/

	/*
        int buffer = 4;

	    for (int i = buffer; i < cInput.file_length()-buffer; i ++)
        {
            if(((encoder_vec[i+2] > encoder_open) && (encoder_vec[i-1] > encoder_open))
					|| ((encoder_vec[i+2] > (encoder_open+2.5)) 
						&& (encoder_vec[i-1] > (encoder_open+2.5))))
                up_cross_over.push_back(i);
            if(((encoder_vec[i+2] > encoder_close) && (encoder_vec[i-1] > encoder_close))
					|| ((encoder_vec[i+2] > (encoder_close+2.5)) 
						&& (encoder_vec[i-1] > (encoder_close+2.5))))
                down_cross_over.push_back(i);
        }

        // Map pressure signal from each cycle into single space 

	    std::vector<double> cycle_time;
	    int size = up_cross_over.size()-1;

	    if(printFromStart=="true")
	    	size = down_cross_over.size()-1;
   
        // start from first up crossing 
        for (int i = 0; i < size; i++)
        {

	        int strtIndx = up_cross_over[i];
	        int endIndx  = up_cross_over[i+1];

		    if(printFromStart=="true")
		    {
			    strtIndx = down_cross_over[i];
			    endIndx   = down_cross_over[i+1];
			    size = down_cross_over.size()-1;
		    }

            for (int j = strtIndx; j < endIndx; j++)
            {
			    encoder_scatter.push_back(encoder_vec[j]);
                pressure_scatter.push_back(cInput.table_value(j,pressur_colmn));
                ion_scatter.push_back(cInput.table_value(j,2));
            }
        }
	    std::ofstream fout1("scatter.dat",std::ios_base::app);
        for (int i = 0; i < encoder_scatter.size(); i++)
        {
            fout1 << encoder_scatter[i]    << '\t'
                  << pressure_scatter[i] << '\t'
                  << ion_scatter[i]      << std::endl;
        }
		exit(1);
	}

    // Average cycles into single curve
	*/

		/*
	int steps = ceil((1.0/frequency)/step);

    double time_step_min = 0.0;
    double time_step_max = step + window;
	double time_step = 0.0;

    for (int i = 0; i < steps; i++)
    {
        double countTime = 0;
        double press_time_sum = 0;
        double time_sum = 0;
        double ion_time_sum = 0;

        for (int j = 0; j < time_scatter.size(); j++)
        {
            if((time_scatter[j] > time_step_min) && (time_scatter[j] < time_step_max))
            {
                press_time_sum = press_time_sum + pressure_scatter[j];
                time_sum = time_sum + time_scatter[j];
                ion_time_sum   = ion_time_sum   + ion_scatter[j];
                countTime++;
            }
        }

		if(countTime != 0)
		{
			time_average.push_back(time_sum/double(countTime));
			pressure_average.push_back(press_time_sum/double(countTime));
			ion_average.push_back(ion_time_sum/double(countTime));
		}

		time_step = time_step + step;
        time_step_min = time_step - window;
        time_step_max = time_step + window;
    }

	for (int i = 0; i < time_average.size(); i++)
	{
		phase_average.push_back(time_average[i]/time_average[time_average.size()-1]);
	}

    std::ofstream fout2("average.dat");
    for (int i = 0; i < time_average.size(); i++)
    {
        fout2 << time_average[i]     << '\t'
              << pressure_average[i] << '\t'
              << ion_average[i]      << std::endl;
    }


	*/
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

}


void average_cycle_passive(DataFileInput &cInput, 
                   int pressur_colmn,
                   int ion_colmn,
                   double window,
				   double step,
                   std::vector<double> &phase_average,
				   std::vector<double> &time_average,
                   std::vector<double> &pressure_average,
                   std::vector<double> &ion_average,
                   std::vector<double> &phase_scatter,
				   std::vector<double> &time_scatter,
                   std::vector<double> &pressure_scatter,
                   std::vector<double> &ion_scatter,
				   std::vector<double> &down_cross_over,
				   std::vector<double> &up_cross_over,
				   double frequency,
                   double static_pressure,
				   std::string printFromStart)
{


    // Isolate cycle stop and start positions from noise

    if ((static_pressure < 0.9) || (static_pressure > 1.1))
        static_pressure = 1.01325;
	int time_colmn = 0;

    int buffer = 4;

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

    std::ofstream fout_cross("crossover.dat");
	for (int i = 0; i < up_cross_over.size(); i++)
	{
		fout_cross << cInput.table_value(down_cross_over[i],time_colmn) << '\t'
			       << cInput.table_value(down_cross_over[i],pressur_colmn) << std::endl;
	}

    // Map pressure signal from each cycle into single space 

	std::vector<double> cycle_time;


	int size = up_cross_over.size()-1;

	if(printFromStart=="true")
		size = down_cross_over.size()-1;
   
    // start from first up crossing 
    for (int i = 0; i < size; i++)
    {

	    int strtIndx = up_cross_over[i];
	    int endIndx   = up_cross_over[i+1];

		if(printFromStart=="true")
		{
			strtIndx = down_cross_over[i];
			endIndx   = down_cross_over[i+1];
			size = down_cross_over.size()-1;
		}

        for (int j = strtIndx; j < endIndx; j++)
        {
            // time mapping
            double dt = cInput.table_value(endIndx,0)
                       -cInput.table_value(strtIndx,0);
            double t_initial = cInput.table_value(strtIndx,0);
            double posPhase = (cInput.table_value(j,0)
                               * (1.0/dt)) - (t_initial/dt);
			double posTime = cInput.table_value(j,0)
				             - cInput.table_value(strtIndx,0);
            phase_scatter.push_back(posPhase);
			time_scatter.push_back(posTime);
            pressure_scatter.push_back(cInput.table_value(j,pressur_colmn));
            ion_scatter.push_back(cInput.table_value(j,ion_colmn));
        }
    }
	std::ofstream fout1("scatter.dat",std::ios_base::app);
    for (int i = 0; i < phase_scatter.size(); i++)
    {
        fout1 << phase_scatter[i]    << '\t'
			  << time_scatter[i]     << '\t'
              << pressure_scatter[i] << '\t'
              << ion_scatter[i]      << std::endl;
    }

    // Average cycles into single curve

    //int steps = ceil(1.0/window);
//	std::cout << "step " << step << " window " << window <<  std::endl;
//	exit(1);
	int steps = ceil((1.0/frequency)/step);

    //double phase_step_min = 0.0;
    //double phase_step_max = window;
    double time_step_min = 0.0;
    double time_step_max = step + window;
	double time_step = 0.0;

    for (int i = 0; i < steps; i++)
    {
        //double countPhase = 0;
        double countTime = 0;
        //double press_phase_sum = 0;
        double press_time_sum = 0;
        //double phase_sum = 0;
        double time_sum = 0;
        //double ion_phase_sum = 0;
        double ion_time_sum = 0;
		/*

        for (int j = 0; j < phase_scatter.size(); j++)
        {
            if((phase_scatter[j] > phase_step_min) && (phase_scatter[j] < phase_step_max))
            {
                press_phase_sum = press_phase_sum + pressure_scatter[j];
                phase_sum = phase_sum + phase_scatter[j];
                ion_phase_sum   = ion_phase_sum   + ion_scatter[j];
                countPhase++;
            }
        }
		*/
        for (int j = 0; j < time_scatter.size(); j++)
        {
			/*
			std::cout << time_step_min << '\t'
				      << time_scatter[j] << '\t'
				      << time_step_max << std::endl; 
					  */
            if((time_scatter[j] > time_step_min) && (time_scatter[j] < time_step_max))
            {
                press_time_sum = press_time_sum + pressure_scatter[j];
                time_sum = time_sum + time_scatter[j];
                ion_time_sum   = ion_time_sum   + ion_scatter[j];
                countTime++;
            }
        }
        //phase_average.push_back(phase_sum/double(countPhase));
		if(countTime != 0)
		{
			time_average.push_back(time_sum/double(countTime));
			pressure_average.push_back(press_time_sum/double(countTime));
			ion_average.push_back(ion_time_sum/double(countTime));
		}
        //phase_step_min = phase_step_min + window;
        //phase_step_max = phase_step_max + window;
		time_step = time_step + step;
        time_step_min = time_step - window;
        time_step_max = time_step + window;
    }

	for (int i = 0; i < time_average.size(); i++)
	{
		phase_average.push_back(time_average[i]/time_average[time_average.size()-1]);
	}

    std::ofstream fout2("average.dat");
    for (int i = 0; i < time_average.size(); i++)
    {
        fout2 << time_average[i]     << '\t'
              << pressure_average[i] << '\t'
              << ion_average[i]      << std::endl;
    }
}

void average_files(std::vector<double> time_scatter,
		           std::vector<double> pressure_scatter,
				   std::vector<double> ion_scatter,
				   std::vector<double> &time_average,
				   std::vector<double> &pressure_average,
				   std::vector<double> &ion_average,
				   std::vector<double> &phase_average)
{
	double max_time = *max_element(std::begin(time_scatter), std::end(time_scatter));
	int steps = 100;

	double step = max_time/double(steps);
	double window = 1.5*step;
    //int steps = window/step;//ceil((1.0/frequency)/step);

    double time_step_min = 0.0;
    double time_step_max = step + window;
	double time_step = 0.0;

    for (int i = 0; i < steps; i++)
    {
        double countTime = 0;
        double press_time_sum = 0;
        double time_sum = 0;
        double ion_time_sum = 0;

        for (int j = 0; j < time_scatter.size(); j++)
        {
            if((time_scatter[j] > time_step_min) && (time_scatter[j] < time_step_max))
            {
                press_time_sum = press_time_sum + pressure_scatter[j];
                time_sum = time_sum + time_scatter[j];
                ion_time_sum   = ion_time_sum   + ion_scatter[j];
                countTime++;
            }
        }
		if(countTime != 0)
		{
			time_average.push_back(time_sum/double(countTime));
			pressure_average.push_back(press_time_sum/double(countTime));
			ion_average.push_back(ion_time_sum/double(countTime));
		}

		time_step = time_step + step;
        time_step_min = time_step - window;
        time_step_max = time_step + window;
    }

	for (int i = 0; i < time_average.size(); i++)
	{
		phase_average.push_back(time_average[i]/time_average[time_average.size()-1]);
	}
}
