#ifndef CODES_POSTPROCESSING_VALVEFUNCTION_SRC_FUNCTIONS_H_
#define CODES_POSTPROCESSING_VALVEFUNCTION_SRC_FUNCTIONS_H_

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
