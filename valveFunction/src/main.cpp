#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include "functions.h"

int main()
{
	double theta = 0;
	double thetaMax = M_PI;
	int thetaCount = 100;

	double thetaStep = thetaMax/double(thetaCount);

	double ballInternalRadius = 21.0/2000.0;
	int xCount = 50;
	double xStep = ballInternalRadius/double(xCount);
	double x = 0;

	/* inlet properties */
	double inletRadius = 21.0/2000.0;


	std::vector<double> inletX;
	std::vector<double> inletY;

	circleCoordinates(inletRadius,
			          50,
					  inletX,
					  inletY);

    std::ofstream foutCircle("circle.dat");
	for (int i = 0; i < inletX.size(); i++)
		foutCircle << inletX[i] << '\t'
			       << inletY[i] << std::endl;
	for (int i = 0; i < inletX.size(); i++)
		foutCircle << inletX[i] << '\t'
			       << -inletY[i] << std::endl;

	int yCount = 50;
	double yStep = ballInternalRadius/double(yCount);

    
	double h = 42.0/1000.0;
	
    std::ofstream foutLine("line.dat");
	for (int i = 0; i < yCount; i++)
	{
		double y = yStep*double(i);
		double l = sqrt(pow(y,2) + pow((h/2),2));
		double alpha = asin(y/l);
		foutLine << xp(l,alpha) << '\t' 
			     << yp(l,alpha) << std::endl;
	}

	double alphaOpen = M_PI/1.5;

    std::ofstream foutLine2("line2.dat");
	for (int i = 0; i < yCount; i++)
	{
		double y = yStep*double(i);
		double l = sqrt(pow(y,2) + pow((h/2),2));
		double alpha = asin(y/l) + alphaOpen;
		foutLine2 << xp(l,alpha) << '\t' 
			      << yp(l,alpha) << std::endl;
	}

	double yOffSet = h*sin(alphaOpen);
	std::vector<double> valveX;
	std::vector<double> valveY;

	std::vector<double> valveOpenX;
	std::vector<double> valveOpenY;

	int count = 0;

    std::ofstream foutEllipse("ellipse.dat");
	for (int i = 0; i < inletX.size(); i++)
	{
		double l = sqrt(pow(inletY[i],2) + pow((h/2),2));
		double alpha = asin(inletY[i]/l) + alphaOpen;

		valveX.push_back(inletX[i]);
		valveY.push_back(yp(l,alpha));
		count++;

		//if((std::abs(valveX[i]) < inletRadius) && (std::abs(valveY[i]) < inletRadius)

		if(sqrt(pow(valveX[i],2) + pow(valveY[i],2)) < inletRadius)
		{
			valveOpenX.push_back(inletX[i]);
			valveOpenY.push_back(yp(l,alpha));
            foutEllipse << inletX[i] << '\t'
		    	        << yp(l,alpha) << std::endl;
		}
	}
	for (int i = 0; i < inletX.size(); i++)
	{
		double l = sqrt(pow(inletY[i],2) + pow((h/2),2));
		double alpha = asin(inletY[i]/l) + alphaOpen;

		valveX.push_back(inletX[i]);
		valveY.push_back(-yp(l,alpha)+yOffSet);

		if(sqrt(pow(valveX[i+count-1],2) + pow(valveY[i+count-1],2)) < inletRadius)
		{
			valveOpenX.push_back(inletX[i]);
			valveOpenY.push_back(-yp(l,alpha)+yOffSet);
            foutEllipse << inletX[i] << '\t'
		    	        << -yp(l,alpha)+yOffSet << std::endl;
		}
	}

    std::ofstream foutInt("integrate.dat");
    std::ofstream foutIntNeg("integrateNeg.dat");
	double areaEllipseNeg = 0.0;
	double areaEllipse    = 0.0; 
	for (int i = 0; i < valveOpenY.size()-1; i++)
	{
		if(valveOpenY[i] > 0.0)
		{
			areaEllipseNeg = areaEllipseNeg 
				+ (((valveOpenY[i+1]+valveOpenY[i])/2.0)*(valveOpenX[i]-valveOpenX[i+1]));
			foutInt << valveOpenX[i]   << '\t'
				    << valveOpenY[i]   << '\n'
					<< valveOpenX[i+1] << '\t'
					<< valveOpenY[i]   << '\n'
					<< valveOpenX[i]   << '\t'
					<< 0.0             << '\n'
					<< valveOpenX[i+1] << '\t'
					<< 0.0             << std::endl;
		}
		if(valveOpenY[i] < 0.0)
		{
			areaEllipse = areaEllipse 
				+ (((valveOpenY[i+1]+valveOpenY[i])/2.0)*(valveOpenX[i]-valveOpenX[i+1]));
			foutIntNeg << valveOpenX[i]   << '\t'
				       << valveOpenY[i]   << '\n'
					   << valveOpenX[i+1] << '\t'
					   << valveOpenY[i]   << '\n'
					   << valveOpenX[i]   << '\t'
					   << 0.0             << '\n'
					   << valveOpenX[i+1] << '\t'
					   << 0.0             << std::endl;
		}
	}
	double areaCircle = 0.0;
	for (int i = 0; i < inletX.size(); i++)
	{
		if((inletX[i] > valveOpenX[0])) 
        {
			if(inletY[i] > 0.0)
			{
				areaCircle = areaCircle 
					+ (((inletY[i+1]+inletY[i])/2.0)*(inletX[i]-inletX[i+1]));

			    foutInt << inletX[i]   << '\t'
				        << inletY[i]   << '\n'
				    	<< inletX[i+1] << '\t'
				    	<< inletY[i]   << '\n'
				    	<< inletX[i]   << '\t'
				    	<< 0.0             << '\n'
				    	<< inletX[i+1] << '\t'
				    	<< 0.0             << std::endl;
			}
		}
	}

	double areaInlet = pow(inletRadius,2)*M_PI;
	double areaOpen = 2.0*(0.25*areaInlet - areaEllipseNeg - areaCircle) 
		                   - 2.0*areaEllipse;

	std::cout << "Fully Open: " << areaInlet << '\n' 
		      << "Area  Open: " << areaOpen  << std::endl; 

	std::cout << "areaNeg: " << areaEllipseNeg << '\n'
		      << "areaEllipse: " << areaEllipse << '\n'
		      << "area_circle: " << areaCircle  << '\n'
		      << "circle: " << pow(inletRadius,2)*M_PI << std::endl; 
	
	return 0;
}

double xp(double h,
		  double alpha)
{
	return h*sin((M_PI/2.0)-alpha);
}

double yp(double h,
		  double alpha)
{
	return h*sin(alpha);
}

void circleCoordinates(double radius,
		               int angleSteps,
		               std::vector<double> &x,
			           std::vector<double> &y)
{
	double angleStep = (M_PI/2.0)/double(angleSteps);
	for (int i = 0; i < angleSteps+1; i++)
	{
		double angle = angleStep*double(i);
		x.push_back(cos(angle)*radius);
		y.push_back(sin(angle)*radius);
	}
}

