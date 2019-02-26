#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "datafileinput.h"
#include "smoothencoder.h"
#include "averagecycle.h"
#include "findroots.h"

//		std::vector<double> encoder_average;
//		std::vector<double> pressure_average;
//		std::vector<double> ion_probe_average;

void average_cycle_active(DataFileInput &cInput, 
		                  double window,
						  double step,
						  double pressureOffset,
						  std::vector<double> &time_ssa,
						  std::vector<double> &pressure_ssa,
						  std::vector<double> &ion_ssa,
						  int encoder_colmn,
				          int pressur_colmn,
						  int ion_colmn,
				          std::vector<double> &encoder_average,
				          std::vector<double> &pressure_average,
						  std::vector<double> &ion_average,
						  std::vector<double> &pressure_scatter,
						  std::vector<double> &encoder_scatter,
						  std::vector<double> &ion_scatter,
						  std::vector<double> &up_cross_over,
						  std::vector<double> &down_cross_over,
						  std::vector<double> &areaFunction,
						  std::string shiftCurve,
						  std::string printFromStart)
{
	std::cout << "Average active valve data" << std::endl;
	std::cout << "encoder colmn: " << encoder_colmn << std::endl;
	std::vector<double> encoder_vec = smooth_encoder(cInput,
			                                         encoder_colmn);

	std::cout << "ion_comn " << ion_colmn << std::endl;
	std::cout << "averaged encoder" << std::endl;
	int steps = ceil(2.5/window);
	int lines = cInput.file_length(); 

	double step_min = 0.0;
	double step_max = window;

	if(printFromStart=="true")
	{
		std::cout << "printing from start" << std::endl;

		int openIndex = 0;
		int closeIndex = 0;
		for (int i = 1; i < areaFunction.size()-1; i++)
		{
			if((areaFunction[i]==0.0) && (areaFunction[i-1]!=0.0))
				closeIndex = i;
			if((areaFunction[i]==0.0) && (areaFunction[i+1]!=0.0))
				openIndex = i;
		}

		double encoderStep = 2.5/double(areaFunction.size());
		double encoder_open = double(openIndex)*encoderStep;
		double encoder_close = double(closeIndex)*encoderStep;

		/* re-order area function */
		std::vector<double> tmp;
		for (int i = openIndex; i < areaFunction.size(); i++)
			tmp.push_back(areaFunction[i]);
		for (int i = 0; i < openIndex; i++)
			tmp.push_back(areaFunction[i]);

		areaFunction.clear();
		areaFunction = tmp;
		tmp.clear();


		std::vector<double> encoder_clean;

		for (int i = 0; i < encoder_vec.size(); i++)
		{
		    if(encoder_vec[i] > 2.5)
			    encoder_vec[i] = encoder_vec[i] - 2.5;
			if(encoder_vec[i] < 0.0)
				encoder_vec[i] = 0.0;
		}

		std::cout << "encoder_open: " << encoder_open << std::endl;
		std::cout << "encoder_close: " << encoder_close << std::endl;

		/* loop through encoder file */
		for (int i = 4; i < encoder_vec.size()-4; i++)
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
				down_cross_over.push_back(i);
			if((encoder_vec[i]   <= encoder_close) &&
			   (encoder_vec[i-1] < encoder_close)  &&
			   (encoder_vec[i-2] < encoder_close)  &&
			   (encoder_vec[i-3] < encoder_close)  &&
			   (encoder_vec[i-4] < encoder_close)  &&
			   (encoder_vec[i+1] > encoder_close)  &&
			   (encoder_vec[i+2] > encoder_close)  &&
			   (encoder_vec[i+3] > encoder_close)  &&
			   (encoder_vec[i+4] > encoder_close))
				up_cross_over.push_back(i);
			
		}

		std::vector<double> time_scatter;
		std::vector<double> phase_scatter;

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
                            - cInput.table_value(strtIndx,0);
                double t_initial = cInput.table_value(strtIndx,0);
                double posPhase = (cInput.table_value(j,0)
                                   * (1.0/dt)) - (t_initial/dt);
			    double posTime = cInput.table_value(j,0)
				                 - cInput.table_value(strtIndx,0);
                encoder_scatter.push_back(encoder_vec[j]);
				phase_scatter.push_back(posPhase);
			    time_scatter.push_back(posTime);
                pressure_scatter.push_back(cInput.table_value(j,1));
                ion_scatter.push_back(cInput.table_value(j,ion_colmn));
            }
        }

		for (int i = 0; i < encoder_scatter.size(); i++)
		{
			if(encoder_scatter[i] < encoder_close)
			{
				encoder_scatter[i] = 2.5 + encoder_scatter[i];
			}
			encoder_scatter[i] = encoder_scatter[i] - encoder_close;
		}
		
        // Average cycles into single curve
	    int steps = ceil(2.5/step);


        double time_step_min = 0.0;
        double time_step_max = step + window;
	    double time_step = 0.0;

        double encoder_step_min = 0.0;
        double encoder_step_max = step + window;
	    double encoder_step = 0.0;


		std::cout << "Number of steps: " << steps << std::endl;

        for (int i = 0; i < steps; i++)
        {
            double countStep = 0;
            double press_time_sum = 0;
            double encoder_sum = 0;
            double ion_time_sum = 0;

            for (int j = 0; j < time_scatter.size(); j++)
            {
                if((encoder_scatter[j] > encoder_step_min) && (encoder_scatter[j] < encoder_step_max))
                {
                    press_time_sum = press_time_sum + pressure_scatter[j];
                    encoder_sum = encoder_sum + encoder_scatter[j];
                    ion_time_sum   = ion_time_sum   + ion_scatter[j];
                    countStep++;
                }
            }

		    if(countStep != 0)
		    {
			    encoder_average.push_back(encoder_sum/double(countStep));
			    pressure_average.push_back(press_time_sum/double(countStep));
			    ion_average.push_back(ion_time_sum/double(countStep));
		    }
		    encoder_step = encoder_step + step;
            encoder_step_min = encoder_step - window;
            encoder_step_max = encoder_step + window;
        }
		std::cout << "encoder average: " <<
			encoder_average.size() << std::endl;
		encoder_average[0] = 0.0;
		/* compute shift to aline pressure min */
		/* find pressure min */
		/*
		int minPresIndx = -1;
		for (int i = 0; i < pressure_average.size(); i++)
		{
			double minTmp = 1.0e+10;
			if(pressure_average[i] < minTmp)
			{
				minTmp = pressure_average[i];
				minPresIndx = i;
			}
		}
		double shift = encoder_average[minPresIndx] - pressureOffset;
		*/

		double Pmin = 0.0;
		double PminPhase = 0.0;
		std::cout << "about to search for min" << std::endl;
		find_min(encoder_average,
				 pressure_average,
				 Pmin,
				 PminPhase);

		std::cout << "found min" << std::endl;

		double shift = pressureOffset - PminPhase;
		if(shiftCurve=="false")
			shift = 0.0;
		std::cout << std::endl;
		std::cout << "****finding min offset****" << std::endl;
		std::cout << "    pMinPhase = " << PminPhase << std::endl;
		std::cout << "    pOffset   = " << pressureOffset << std::endl;
		std::cout << "    shift     = " << shift << std::endl;
		/* shift for pressure offset */
		for (int i = 0; i < encoder_average.size(); i++)
		{
			//double encoder_tmp = encoder_average[i] - pressureOffset;
			double encoder_tmp = encoder_average[i] + shift;
			if(encoder_tmp < 0.0)
				encoder_tmp = 2.5 + encoder_tmp;
			else if(encoder_tmp > 2.5)
				encoder_tmp = encoder_tmp - 2.5;
			encoder_average[i] = encoder_tmp;
		}
		std::vector<double> encoderTmp;
		std::vector<double> pressureTmp;
		std::vector<double> ionTmp;

		while (encoderTmp.size() < encoder_average.size())
		{
			double max = -1.0;
			double min = 10.0e+10;
			int arrangeIndex = 0;
			/*

			for (int i = 0; i < encoderTmp.size(); i++)
			{
				if(encoderTmp[i] > max) 
					max = encoderTmp[i];
			}
			*/
			if(encoderTmp.size() != 0)
				max = encoderTmp[encoderTmp.size()-1];
		    for (int i = 0; i < encoder_average.size(); i++)
		    {
			    if((encoder_average[i] < min) && (encoder_average[i] > max)) 
				{
					min = encoder_average[i]; 
					arrangeIndex = i;
				}
		    }
			encoderTmp.push_back(encoder_average[arrangeIndex]);
			pressureTmp.push_back(pressure_average[arrangeIndex]);
			ionTmp.push_back(ion_average[arrangeIndex]);
		}

		encoder_average.clear();
		encoder_average = encoderTmp;

		pressure_average.clear();
		pressure_average = pressureTmp;

		ion_average.clear();
		ion_average = ionTmp;
	}
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

void average_files_active(std::vector<double> encoder_scatter,
		                  std::vector<double> pressure_scatter,
				          std::vector<double> ion_scatter,
						  std::vector<double> mass_scatter,
				          std::vector<double> &encoder_average,
				          std::vector<double> &pressure_average,
				          std::vector<double> &ion_average,
						  std::vector<double> &mass_average,
				          std::string ion_probe)
{
	double max_encoder = *max_element(std::begin(encoder_scatter), std::end(encoder_scatter));
	int steps = 60;

	double step = max_encoder/double(steps);
	double window = 5.0*step;
    //int steps = window/step;//ceil((1.0/frequency)/step);

    double encoder_step_min = 0.0;
    double encoder_step_max = step + window;
	double encoder_step = 0.0;

    for (int i = 0; i < steps; i++)
    {
        double countEncoder = 0.0;
        double press_encoder_sum = 0.0;
        double encoder_sum = 0.0;
        double ion_encoder_sum = 0.0;
        double mass_encoder_sum = 0.0;

        for (int j = 0; j < encoder_scatter.size(); j++)
        {
            if((encoder_scatter[j] > encoder_step_min) && (encoder_scatter[j] < encoder_step_max))
            {
                press_encoder_sum = press_encoder_sum + pressure_scatter[j];
                encoder_sum = encoder_sum + encoder_scatter[j];
				mass_encoder_sum = mass_encoder_sum + mass_scatter[j];
				if(ion_probe=="true")
                    ion_encoder_sum   = ion_encoder_sum   + ion_scatter[j];
                countEncoder++;
            }
        }
		if(countEncoder != 0)
		{
			encoder_average.push_back(encoder_sum/double(countEncoder));
			pressure_average.push_back(press_encoder_sum/double(countEncoder));
			mass_average.push_back(mass_encoder_sum/double(countEncoder));
			if(ion_probe=="true")
			    ion_average.push_back(ion_encoder_sum/double(countEncoder));
		}

		encoder_step = encoder_step + step;
        encoder_step_min = encoder_step - window;
        encoder_step_max = encoder_step + window;
    }
}
