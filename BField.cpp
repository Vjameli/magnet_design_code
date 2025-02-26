//*****************************************************************************
//        C++ * Iman Dayarian * Toronto - CANADA *
//                                                                            *
//        July 2015 *
//*****************************************************************************

#include "BField.h"
#include "float.h"
#include "math.h"
#include "stdio.h"
#include <iostream>

#define MAX(x, y) (x > y ? x : y)
#define MAX3(x, y, z) MAX(MAX(x, y), z)
#define MIN(x, y) (x > y ? y : x)
#define MIN3(x, y, z) MIN(MIN(x, y), z)

double BField::drf(double x, double y, double z, int *piErr) {
  int iErr = 0;
  double mu, xn, yn, zn, xndev, yndev, zndev, xnroot, ynroot, znroot, lambda,
      epsilon, e2, e3, result, s;

  const double c1 = 1.0 / 24.0;
  const double c2 = 3.0 / 44.0;
  const double c3 = 1.0 / 14.0;
  const double errtol = pow(DBL_EPSILON * 4.0, 1.0 / 6.0);
  const double lolim = 5.0 * DBL_MIN;
  const double hilim = DBL_MAX / 5.0;

  if (piErr) {
    if (MIN3(x, y, z) < 0.0) {
      iErr = 1;
    } else if (MIN3(x + y, x + z, y + z) < lolim) {
      iErr = 2;
    } else if (MAX3(x, y, z) > hilim) {
      iErr = 3;
    }
  }
  if (iErr) {
    if (piErr) {
      *piErr = iErr;
    }
    result = 0.0;
  } else {
    xn = x;
    yn = y;
    zn = z;

    while (1) {
      mu = (xn + yn + zn) / 3.0;
      xndev = 2.0 - (mu + xn) / mu;
      yndev = 2.0 - (mu + yn) / mu;
      zndev = 2.0 - (mu + zn) / mu;
      epsilon = MAX3(fabs(xndev), fabs(yndev), fabs(zndev));
      if (epsilon < errtol)
        break;
      xnroot = sqrt(xn);
      ynroot = sqrt(yn);
      znroot = sqrt(zn);
      lambda = xnroot * (ynroot + znroot) + ynroot * znroot;
      xn = (xn + lambda) * 0.25;
      yn = (yn + lambda) * 0.25;
      zn = (zn + lambda) * 0.25;
    }
    e2 = xndev * yndev - pow(zndev, 2);
    e3 = xndev * yndev * zndev;
    s = 1.0 + (c1 * e2 - 0.1 - c2 * e3) * e2 + c3 * e3;

    if (piErr) {
      *piErr = 0;
    }
    result = s / sqrt(mu);
  }
  return result;
}
/* END drf() */

/*
 * E(k): complete elliptic integral of the second kind
 */

#define E(k, ierr)                                                             \
  (drf(0.0, 1.0 - pow(k, 2), 1.0, &ierr) -                                     \
   (pow(k, 2) / 3.0) * drd(0.0, 1.0 - pow(k, 2), 1.0, &ierr))

/*
 * FastE(k): fast, complete elliptic integral of the second kind.
 *           Use this macro if the complete elliptic integral of
 *           the first kind was previously computed for the same
 *           value of k.
 */

#define FastE(F, k, ierr)                                                      \
  ((F) - (pow(k, 2) / 3.0) * drd(0.0, 1.0 - pow(k, 2), 1.0, &ierr))

/*
 * drd.c:   Compute the complete or incomplete elliptic integral of the
 *          second kind.
 *
 * Description:
 *
 *  For x and y non-negative, x+y and z positive, drf(x,y,z) = integral
 *  from zero to infinity of
 *
 *
 *            -1/2     -1/2     -3/2
 *  (3/2)(t+x)    (t+y)    (t+z)    dt.
 *
 *  If x or y is zero, the integral is complete.
 *
 *  *piErr returns non-zero if any arguments are invalid.
 *
 * Credits:
 *
 *  This function is adapted by E. Dennison from FORTRAN code written by:
 *
 *  Carlson, B.C.
 *  Notis, E.M
 *  Pexton, R.L.
 *
 */
double BField::drd(double x, double y, double z, int *piErr) {
  int iErr = 0;
  double mu, xn, yn, zn, xndev, yndev, zndev, xnroot, ynroot, znroot, lambda,
      epsilon, ea, eb, ec, ed, ef, sigma, power4, result, s1, s2;

  const double c1 = 3.0 / 14.0;
  const double c2 = 1.0 / 6.0;
  const double c3 = 9.0 / 22.0;
  const double c4 = 3.0 / 26.0;
  const double errtol = pow(DBL_EPSILON / 3.0, 1.0 / 6.0);
  double uplim;
  const double lolim = 2.0 / pow(DBL_MAX, 2.0 / 3.0);
  double tuplim = pow(DBL_MIN, 1.0 / 3.0);
  tuplim = pow(0.1 * errtol, 1.0 / 3.0) / tuplim;
  uplim = pow(tuplim, 2.0);

  if (piErr) {
    if (MIN(x, y) < 0.0) {
      iErr = 1;
    } else if (MAX3(x, y, z) > uplim) {
      iErr = 2;
    } else if (MIN(x + y, z) < lolim) {
      iErr = 3;
    }
  }
  if (iErr) {
    if (piErr) {
      *piErr = iErr;
    }
    result = 0.0;
  } else {
    xn = x;
    yn = y;
    zn = z;
    sigma = 0.0;
    power4 = 1.0;
    while (1) {
      mu = (xn + yn + 3.0 * zn) * 0.2;
      xndev = (mu - xn) / mu;
      yndev = (mu - yn) / mu;
      zndev = (mu - zn) / mu;
      epsilon = MAX3(fabs(xndev), fabs(yndev), fabs(zndev));
      if (epsilon < errtol)
        break;
      xnroot = sqrt(xn);
      ynroot = sqrt(yn);
      znroot = sqrt(zn);
      lambda = xnroot * (ynroot + znroot) + ynroot * znroot;
      sigma = sigma + power4 / (znroot * (zn + lambda));
      power4 = power4 * 0.25;
      xn = (xn + lambda) * 0.25;
      yn = (yn + lambda) * 0.25;
      zn = (zn + lambda) * 0.25;
    }
    ea = xndev * yndev;
    eb = zndev * zndev;
    ec = ea - eb;
    ed = ea - 6.0 * eb;
    ef = ed + ec + ec;
    s1 = ed * (-c1 + 0.25 * c3 * ed - 1.5 * c4 * zndev * ef);
    s2 = zndev * (c2 * ef + zndev * (-c3 * ec + zndev * c4 * ea));
    if (piErr) {
      *piErr = 0;
    }
    result = 3.0 * sigma + power4 * (1.0 + s1 + s2) / (mu * sqrt(mu));
  }
  return result;
}
/* END drd() */

/*
 * MAIN:    Sample program using drd and drf
 *
 * Description:
 *
 *          Demonstrate use of drd and drf in computing the magnetic
 *          field at any point in space due to a circular current
 *          filament (one turn coil).
 *
 *
 */

BField::BField(double i, double a, double z0, double r0, double zmin,
               double zmax, double rmin, double rmax, double step,
               vector<vector<double>> &Bz, vector<vector<double>> &Br,
               vector<vector<double>> &Z, vector<vector<double>> &R) {
#define PI (3.141592654)
#define MU0 (PI * 4.0E-7)

  double k, q, rq, fk, ek, al, be, ga, al2, be2, alt4, Ht;
  int ierr;

  /* begin computation here */

  int in1 = 0;
  for (double z = zmin - z0; z <= zmax - z0; in1++) {
    int in2 = 0;
    for (double r = rmin - r0; r <= rmax - r0; in2++) {
      al = abs(r / a);
      alt4 = al * 4.0;
      be = abs(z / a);
      ga = z / r;
      al2 = al * al;
      be2 = be * be;
      q = pow(1 + al, 2) + be2;
      rq = sqrt(q);
      k = sqrt(alt4 / q);
      fk = F(k, ierr);
      ek = FastE(fk, k, ierr);
      Ht = i / (2.0 * a * PI * rq);
      Bz[in1][in2] = MU0 * Ht * (ek * (1 - al2 - be2) / (q - alt4) + fk);
      Br[in1][in2] =
          (r == 0.0)
              ? (0.0)
              : MU0 * (Ht * ga * (ek * (1 + al2 + be2) / (q - alt4) - fk));

      Z[in1][in2] = z;
      R[in1][in2] = r;
      r = r + step;
    }
    z = z + step;
  }
}

// BField::BField(double i, double a, double Zc, double Rc, double Zt, double
// Rt, double &Bz, double &Br)
// {
// 		//field = new BField(1, CoilsR[n], CoilsZ[n], 0, Zt[m], Rt[m],
// Bz, Br);
//
// #define PI (3.141592654)
// #define MU0 (PI*4.0E-7)
//
//     double k,q,rq,fk,ek,al,be,ga,al2,be2,alt4,Ht;
//     int    ierr;
//
// 	/* begin computation here */
//
//
// 			double r = Rt - Rc;
// 			double z = Zt - Zc;
//
// 			al = abs(r/a);
// 			alt4 = al*4.0;
// 			be = abs(z/a);
// 			ga = z/r;
// 			al2 = al*al;
// 			be2 = be*be;
// 			q = pow(1+al,2)+be2;
// 			rq = sqrt(q);
// 			k = sqrt(alt4/q);
// 			fk = F(k,ierr);
// 			ek = FastE(fk,k,ierr);
// 			Ht = i/(2.0*a*PI*rq);
// 			Bz = MU0 * Ht*(ek*(1-al2-be2)/(q-alt4)+fk);
// 			Br = (r==0.0)?(0.0):MU0
// *(Ht*ga*(ek*(1+al2+be2)/(q-alt4)-fk)); 			std::cout<<"CoilsR[n]: "<<a<<",
// CoilsZ[n] : "<<Zc<<", Zt[m] : "<<Zt<<", Rt[m] : "<<Rt<<", AMN : "<< Bz<<",
// First: "<<fk<<", Second: "<< ek <<std::endl;
//
//
//
// }

BField::BField(double i, double a, double Zc, double Rc, double Zt, double Rt,
               double &Bz, double &Br) {
  // field = new BField(1, CoilsR[n], CoilsZ[n], 0, Zt[m], Rt[m], Bz, Br);

#define pi (3.141592654)

  double miu0 = 4 * pi * 1e-7;
  double B0 = (i * miu0) / (2 * a);

  double r = Rt - Rc;
  double z = Zt - Zc;
  if (r != 0) {
    double al = abs(r / a);                     // variable
    double be = abs(z / a);                     // variable
    double ga = (z / r);                        // variable
    double q = ((1 + al) * (1 + al) + be * be); // variable
    double k = sqrt(4 * al / q); // elliptical integrals coefficient
    double K = alglib::ellipticintegralk(k * k);
    double E = alglib::ellipticintegrale(k * k);
    //[K,E] = ellipke(k^2);    % elliptical integrals
    // display(K);
    // display(E);
    Bz = B0 * (1 / (pi * sqrt(q))) *
         (E * (1 - al * al - be * be) / (q - 4 * al) +
          K); // axial component, e.g. along z-axis
    Br = B0 * (ga / (pi * sqrt(q))) *
         (E * (1 + al * al + be * be) / (q - 4 * al) - K); // radial component
    //        	std::cout<<"CoilsR[n]: "<<a<<", CoilsZ[n] : "<<Zc<<", Zt[m] :
    //        "<<Zt<<", Rt[m] : "<<Rt<<", AMN : "<< Bz<<", First: "<<K<<",
    //        Second: "<< E <<std::endl;

  } else // if r = 0 which means that the magnetic field is computed along the
         // main axis, e.g. z-axis, only
  {
    Bz = B0 * (a / sqrt(a * a + z * z)) * (a / sqrt(a * a + z * z)) *
         (a / sqrt(a * a + z * z));
    Br = 0;
  }
}

// BField::BField(double i, double a, double Zc, double Rc, double Zt, double
// Rt, double &Bz, double &Br)
// {
// #define PI (3.141592654)
// #define MU0 (PI*4.0E-7)
//
//     /* input parameters */
//     // float a; /* loop radius */
// //     float r; /* measurement point radius */
// //     float x; /* measurement point axial position */
// //     float i; /* loop current */
//     /* output parameters */
//     //double Hx,Hr; /* axial, radial field components, A/m */
//     /* working vars */
//     double k,q,rq,fk,ek,al,be,ga,al2,be2,alt4,Ht;
//     int    ierr;
//
//     /* gather input parameters */
// //     printf("Loop radius (meters): ");
// //     scanf("%f",&a);
// //     printf("Radius at measurement point (meters): ");
// //     scanf("%f",&r);
// //     printf("Axial position at measurement point (meters): ");
// //     scanf("%f",&x);
// //     printf("Loop current (amperes): ");
// //     scanf("%f",&i);
//     /* begin computation here */
//
// 	double r = abs(Rt - Rc);
// 	double z = abs(Zt - Zc);
//     al = r/a;
//     alt4 = al*4.0;
//     be = z/a;
//     ga = z/r;
//     al2 = al*al;
//     be2 = be*be;
//     q = pow(1+al,2)+be2;
//     rq = sqrt(q);
//     k = sqrt(alt4/q);
//     fk = alglib::ellipticintegralk(k);
//     ek = FastE(fk,k,ierr);
//     Ht = i/(2.0*a*PI*rq);
//     Bz = MU0*Ht*(ek*(1-al2-be2)/(q-alt4)+fk);
//     Br = (r==0.0)?(0.0):MU0*(Ht*ga*(ek*(1+al2+be2)/(q-alt4)-fk));
//
// 	std::cout<<"CoilsR[n]: "<<a<<", CoilsZ[n] : "<<Zc<<", Zt[m] : "<<Zt<<",
// Rt[m] : "<<Rt<<", AMN : "<< Bz<<", First: "<<fk<<", Second: "<<ek <<", k
// :"<<k<<std::endl;
//
//     /* display the output */
//     // printf("Axial field,  Bx: %e (Tesla)\n",(float)Hx*MU0);
// //     printf("Radial field, Br: %e (Tesla)\n",(float)Hr*MU0);
//
// }
