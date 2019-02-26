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


void track_peaks_active(DataFileInput &cInput,
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
			/*
			if(cInput.table_value(j,i_colmn) < iMin)
			{
				iMin = cInput.table_value(j,i_colmn);
				iMinIndx = j;
			}
			*/
		}

		vt_max_p.push_back(cInput.table_value(pMaxIndx,t_colmn));
		vt_min_p.push_back(cInput.table_value(pMinIndx,t_colmn));
		//vt_min_i.push_back(cInput.table_value(iMinIndx,t_colmn));
		vt_min_i.push_back(0);
		vp_max.push_back(pMax);
		vp_min.push_back(pMin);
		//vi_min.push_back(iMin);
		vi_min.push_back(0);
	}

	std::cout << "finished tracking peaks" << std::endl;
}
