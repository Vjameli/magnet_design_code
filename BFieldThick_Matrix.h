#ifndef BFIELDTHICKMATRIX_H
#define BFIELDTHICKMATRIX_H


#include <list>
#include<vector>
#include <string>
#include <stdlib.h>
#include <math.h>
#include "Integral.h"

using namespace std;


class BFieldThick_Matrix
{
  public :
	BFieldThick_Matrix(double J0, double r0, double dr, double dz, double xmin, double zmin, double xmax, double zmax, double xshift, double zshift, double s, vector<vector<double> > &BzMatrix, vector<vector<double> > &BxMatrix);
	void Polar(double x, double y, double z, double& r, double& theta, double& z1);

};

#endif

