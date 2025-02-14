#ifndef ELLIPSE_H
#define ELLIPSE_H


#include <list>
#include<vector>
#include <string>
#include <stdlib.h>
#include <math.h>
using namespace std;


class Ellipse
{
  public :
  	Ellipse(double x, double y, double a, double b, double angle, int steps, vector<double>& X, vector<double>& Y, vector<double>& PolarR, vector<double>& PolarTheta);
	void Polar(double x, double y, double& r, double& theta);
};

#endif

