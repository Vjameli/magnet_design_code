#ifndef BFIELDTHICK_H
#define BFIELDTHICK_H

#include "Integral.h"
#include <list>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

class BFieldThick {
public:
  BFieldThick(double J0, double r0, double dr, double dz, double Xt, double Zt,
              double xshift, double zshift, double &Bz, double &Br);
  void Polar(double x, double y, double z, double &r, double &theta,
             double &z1);
};

#endif
