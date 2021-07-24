#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include "computeballvalvearea.h"

double ballValveArea(double theta,
		             std::vector<double> &areaFunction)
{
	int thetaCount = areaFunction.size();
	double thetaStep = M_PI/double(thetaCount);

	if (theta==0.0)
	{
		return areaFunction[0];
    }
	else if((theta > 0.0) && (theta < M_PI))
	{
	    /* if theta between 0 and Pi/2 
	         valve closing */

		int thetaIndxHigh = floor(theta/thetaStep);
		int thetaIndxLow  = ceil(theta/thetaStep);
		if(thetaIndxHigh > (areaFunction.size()-1))
			thetaIndxHigh = areaFunction.size()-1;
		if(thetaIndxLow > (areaFunction.size()-1))
			thetaIndxLow = areaFunction.size()-2;
		if(thetaIndxHigh==thetaIndxLow)
		{
			if(thetaIndxHigh+1 < (areaFunction.size()-1))
				thetaIndxHigh = thetaIndxHigh + 1;
		    else
		    	thetaIndxLow = thetaIndxLow - 1;
		}
		return interpolate(areaFunction[thetaIndxLow],
				           areaFunction[thetaIndxHigh],
						   double(thetaIndxLow)*thetaStep,
						   double(thetaIndxHigh)*thetaStep,
						   theta);
	}
	else if((theta >= M_PI) && (theta < 2.0*M_PI))
	{
	    /* if theta between Pi/2 and Pi
	         valve openning */

		int thetaIndxHigh = floor(theta/thetaStep);
		int thetaIndxLow  = ceil(theta/thetaStep);
		if(thetaIndxHigh > (areaFunction.size()-1))
			thetaIndxHigh = areaFunction.size()-1;
		if(thetaIndxLow > (areaFunction.size()-1))
			thetaIndxLow = areaFunction.size()-2;
		if(thetaIndxHigh==thetaIndxLow)
		{
			if(thetaIndxHigh+1 < (areaFunction.size()-1))
				thetaIndxHigh = thetaIndxHigh + 1;
		    else
		    	thetaIndxLow = thetaIndxLow - 1;
		}
		return interpolate(areaFunction[thetaIndxLow],
				           areaFunction[thetaIndxHigh],
						   double(thetaIndxLow)*thetaStep,
						   double(thetaIndxHigh)*thetaStep,
						   theta);
	}
	else if(theta==2.0*M_PI)
	{
		return areaFunction[areaFunction.size()-1];
	}
	else
	{
		/* theta out of range */
	}
}

double interpolate(double y0,
		           double y1,
				   double x0,
				   double x1,
				   double x)
{
	return y0 + (x - x0)*((y1-y0)/(x1-x0));
}

void computeBallValveAreaFunction(valveProperties &valveGeom,
		                          int thetaCount,
		                          std::vector<double> &areaFunction)
{
	double thetaStep = (M_PI/2.0)/double(thetaCount);

	for (int i = 0; i < thetaCount; i++)
	{
		double theta = double(i)*thetaStep;
		areaFunction.push_back(computeBallValveArea(valveGeom,theta));
	}
	int halfCount = areaFunction.size();
	for (int i = halfCount; i > 0; i--)
		areaFunction.push_back(areaFunction[i-1]);
}

double computeBallValveArea(valveProperties &valveGeom,
		                    double alphaOpen)
{
//	double ballInternalRadius = valveGeom.ballInternalRadius;
	int xCount = valveGeom.spanCount;
	//double xStep = ballInternalRadius/double(xCount);
	//double x = 0;

	/* inlet properties */
	double inletRadius = valveGeom.inletRadius;


	std::vector<double> inletX;
	std::vector<double> inletY;

	circleCoordinates(inletRadius,
			          valveGeom.spanCount,
					  inletX,
					  inletY);

	int yCount = valveGeom.spanCount;
	//double yStep = ballInternalRadius/double(yCount);

    
	double h = valveGeom.h;
	

	double yOffSet = h*sin(alphaOpen);
	std::vector<double> valveX;
	std::vector<double> valveY;

	std::vector<double> valveOpenX;
	std::vector<double> valveOpenY;

	int count = 0;
	
	/* need to not used "inletY" for projection if
	      the diameter of ball internal diameter does
		  not equal the inlet diameter - change in future */

	for (int i = 0; i < inletX.size(); i++)
	{
		double l = sqrt(pow(inletY[i],2) + pow((h/2),2));
		double alpha = asin(inletY[i]/l) + alphaOpen;

		valveX.push_back(inletX[i]);
		valveY.push_back(yp(l,alpha));
		count++;

		if(sqrt(pow(valveX[i],2) + pow(valveY[i],2)) < inletRadius)
		{
			valveOpenX.push_back(inletX[i]);
			valveOpenY.push_back(yp(l,alpha));
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
		}
	}

	bool valveClosed = true;

	if(valveOpenY.size()!=0) 
		valveClosed = false;

	if(valveClosed)
		return 0.0;

	double areaEllipseNeg = 0.0;
	double areaEllipse    = 0.0; 
	for (int i = 0; i < valveOpenY.size()-1; i++)
	{
		if(valveOpenY[i] > 0.0)
			areaEllipseNeg = areaEllipseNeg 
				+ (((valveOpenY[i+1]+valveOpenY[i])/2.0)*(valveOpenX[i]-valveOpenX[i+1]));
		if(valveOpenY[i] < 0.0)
			areaEllipse = areaEllipse 
				+ (((valveOpenY[i+1]+valveOpenY[i])/2.0)*(valveOpenX[i]-valveOpenX[i+1]));
	}
	double areaCircle = 0.0;
	for (int i = 0; i < inletX.size(); i++)
	{
		if((inletX[i] > valveOpenX[0])) 
        {
			if(inletY[i] > 0.0)
				areaCircle = areaCircle 
					+ (((inletY[i+1]+inletY[i])/2.0)*(inletX[i]-inletX[i+1]));
		}
	}

	double areaInlet = pow(inletRadius,2)*M_PI;
	double areaOpen = 2.0*(0.25*areaInlet - areaEllipseNeg - areaCircle) 
		                   - 2.0*areaEllipse;
	if(areaOpen < 0.0)
		areaOpen = 0.0;
	else if(alphaOpen == 0.0)
		areaOpen = areaInlet;
	else if(alphaOpen == 2.0*M_PI)
		areaOpen = areaInlet;

	return areaOpen;
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

