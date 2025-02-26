#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <list>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>
// #include "stdafx.h"
// #include "integration.h"
#include <cmath>    //Header file for mathematical operartions
#include <iomanip>  //Header file for precession
#include <iostream> //Header file for cin & cout

using namespace std;

class Integral {
public:
  Integral(double _t1, double _t2, double _t3, double _t4);

  double t1;
  double t2;
  double t3;
  double t4;
  double N1(double c, double d);
  double N2(double c, double d);
  // double fn1(double x);
  // double fn2(double x);

private:
  long double f1(long double x);
  long double f2(long double x);
  long double pn(long double a[], int n, int m, long double x);
  long double dn(long double a[], int n, int m, long double x);
  long double fact(int n);
};

#endif
