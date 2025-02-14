#ifndef BFIELD_H
#define BFIELD_H


#include <list>
#include<vector>
#include <string>
#include <stdlib.h>
#include <cmath>


using namespace std;


class BField
{
  public :
  	BField(double i, double a, double z0, double r0, double zmin, double zmax, double rmin, double rmax, double step, vector<vector<double> > &Bz, vector<vector<double> > &Br, vector<vector<double> > &Z, vector<vector<double> > &R);
	BField(double i, double a, double Zc, double Rc, double Zt, double Rt, double &Bz, double &Br);

    double drf(double x, double y, double z, int* piErr);
	double drd(double x, double y, double z, int* piErr);

};

#endif

