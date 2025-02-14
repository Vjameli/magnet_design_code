//*****************************************************************************
//        C++      					                                          *
//        Iman Dayarian  	                                                  *
//        Toronto - CANADA            		                                  *
//                                                                            *
//        June 2015                   				                          *
//*****************************************************************************

#include "Ellipse.h"

Ellipse::Ellipse(double x, double y, double a, double b, double angle, int nbPoints, vector<double>& X, vector<double>& Y, vector<double>& PolarR, vector<double>& PolarTheta)
{
#define PI (3.141592654)
	
	double beta = -angle * (PI / 180);
    double sinbeta = sin(beta);
    double cosbeta = cos(beta);

	double step = 180.0/(nbPoints-1);
	int cmp = 0;
	for (double ang = 270; ang < 451; cmp++){
		double alpha = ang * (PI / 180);
		double sinalpha = sin(alpha);
		double cosalpha = cos(alpha);
		X[cmp] = x + (a * cosalpha * cosbeta - b * sinalpha * sinbeta);
		Y[cmp] = y + (a * cosalpha * sinbeta + b * sinalpha * cosbeta);
		
		double r = -1;
		double theta = -1;
		
		Polar(X[cmp], Y[cmp], r, theta);
		
		PolarR[cmp] = r;
		PolarTheta[cmp] = theta;
		
		ang = ang + step;
	}
}

void Ellipse::Polar(double x, double y, double& r, double& theta)
{
	r = sqrt((pow(x,2))+(pow(y,2)));

    theta = atan(y/x);
	
	theta = (theta*180)/3.141592654;
}
