#ifndef CODES_POSTPROCESSING_PROCESSPC_SRC_COMPUTEVALVEAREA_H_
#define CODES_POSTPROCESSING_PROCESSPC_SRC_COMPUTEVALVEAREA_H_

#include <vector>

typedef struct {
	int spanCount;
	double inletRadius;
	double ballInternalRadius;
	double h;
} valveProperties;

double ballValveArea(double theta,
		             std::vector<double> &areaFunction);

void computeBallValveAreaFunction(valveProperties &valveGeom,
		                          int thetaCount,
								  std::vector<double> &areaFunction);

double computeBallValveArea(valveProperties &valveGeom,
		                    double alphaOpen);


double interpolate(double y0,
		           double y1,
				   double x0,
				   double x1,
				   double x);

double xp(double h,
          double alpha);

double yp(double h,
          double alpha);

double projectedEllipseY(double x,
                         double Radius,
                         double theta);

void circleCoordinates(double radius,
		               int angleSteps,
					   std::vector<double> &x,
					   std::vector<double> &y);


#endif
