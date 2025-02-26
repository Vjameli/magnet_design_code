//*****************************************************************************
//        Iman Dayarian * Toronto, Canada *
//                                                                            *
//        Generated June 2015                                                 *
//*****************************************************************************

#include <ilcplex/ilocplex.h>

ILOSTLBEGIN

#include "BField.h"
#include "BFieldThick.h"
#include "BFieldThick_Matrix.h"
#include "Ellipse.h"
#include "Model2.h"
#include "tempsC++.h"
#include "supplementary_functions.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip> // std::setprecision
#include <list>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <values.h>
#include <vector>

using namespace std;

int main() {

  cout << "Main is working" << endl;

  double B_zero = 0.5; // [T] - Magnetic field to achieve
  double tol = 10e-6;  // Magnetic field homogeneity tolerance

  double Gamma = 1.4; // Cross-section ratio
  double CCD = 21e+7; // Critical current density

  // Define candidate coils for North pole: Initial grid
  double zmin = 34; // [cm] z-domain
  double zmax = 42;
  double rmin = 10; // [cm] r-domain
  double rmax = 48;
  double step = 10; // BIGGER STEP -> LESS COILS (it was 2 before)

  // Define and initialize vector of coordinates for the ellipse
  int nbPoints = 5; // (it was 40 before)
  vector<double> Xtarget;
  vector<double> Ytarget;
  vector<double> PolarR;
  vector<double> PolarTheta;

  initializeTargetPoints(nbPoints, Xtarget, Ytarget, PolarR, PolarTheta);

  // Sample optimization points on the ellipse's border
  sampleEllipseBorder(nbPoints, Xtarget, Ytarget, PolarR, PolarTheta);

  vector<double> Zt;
  vector<double> Rt;
  vector<double> Xt;
  Zt.resize(nbPoints);
  Rt.resize(nbPoints);
  Xt.resize(nbPoints);

  // Cilyndrical coordinates of the points
  Zt = Ytarget;
  Rt = PolarR;

  // Initialiize Xt with the absolute value of Xtarget
  for (int i = 0; i < nbPoints; i++) {
    Xt[i] = abs(Xtarget[i]);
  }

  int Nt = Zt.size();

  // Sampling a set of points inside the target volume
  Ellipse *ellipse2;
  vector<double> Zt_in; // Stores the y target coordinates
  vector<double> Xt_in; // Stores the x target coordinates
  Zt_in.push_back(0);
  Xt_in.push_back(0);

  sampleEllipseInside(nbPoints, Xtarget, Ytarget, PolarR, PolarTheta, Zt_in,
                      Xt_in);

  int Nt_in = Zt_in.size();

  // Coordinates of the cente of the coils based on the step (IN METERS!)
  vector<double> NPZ; // North pole Z
  vector<double> NPR;
  vector<double> SPZ; // South pole Z
  vector<double> SPR;

  computeCoilCenters(NPZ, NPR, SPZ, SPR, zmin, zmax, rmin, rmax, step);

  // CoilsZ = concatenation of (NPZ + SPZ), analogous for CoilsR
  vector<double> CoilsZ;
  vector<double> CoilsR;
  CoilsZ = NPZ;
  CoilsR = NPR;
  CoilsZ.insert(CoilsZ.end(), SPZ.begin(), SPZ.end());
  CoilsR.insert(CoilsR.end(), SPR.begin(), SPR.end());
  int Nc = CoilsZ.size();

  // Set of potential cross sections (Potential widths)
  // VALUE WAS ORIGINALLY EQUAL TO 5!!!!!!!!!! <-------
  vector<double> widths;
  vector<double> cross_sections;

  computeCrossSections(cross_sections, widths, Gamma);

  int Ns = cross_sections.size();

  // Compute BField and AMNS matrix
  vector<vector<vector<double>>> AMNS;
  computeAMNS(AMNS, Nt, Nc, Ns, Zt, Xt, CoilsZ, CoilsR, cross_sections, widths,
              Gamma);

  AMNSToExcel(AMNS, Nt, Nc, Ns);

  // Initialize vectors for optimization model and cpu time
  Model2 *design;
  double duration = 0;

  int NbCoils = 7;

  // Initialize 2D SOL vector
  vector<vector<double>> SOL;

  initialize2DSolVector(SOL, Nc, Ns);

  vector<vector<vector<double>>> SOL_VECTOR;

  // Optimization environment starts
  duration = 0;
  ChronoCPU chrono;
  chrono.start();

  design = new Model2(IloEnv(), AMNS, CoilsR, CoilsZ, B_zero, tol, widths,
                      cross_sections, CCD, Gamma, zmin, zmax, rmin, rmax, step,
                      NbCoils);

  design->solve();

  duration = chrono.getTime();

  if (design->getCplexStatus() != IloAlgorithm::Infeasible &&
      design->getCplexStatus() != IloAlgorithm::InfeasibleOrUnbounded) {

    // Just in case we run the optimization many times and need to store more
    // than one result
    design->getSol(SOL);

    double B[Nt];
    evaluateSolutionBorder(SOL, AMNS, Nt, Nc, Ns, B);

    // Compute the resultant magnetic field from the solution inside the FOV.
    double B_in[Nt_in];
    evaluateSolutionInside(SOL, Nt_in, Nc, Ns, Xt_in, Zt_in, CoilsR, CoilsZ,
                           widths, cross_sections, Gamma, B_in);

    // Peak to peak on the border and inside of the FOV
    // RMS on the border and inside of the FOV
    double Peak2Peak_in = computePeak2Peak(B_in, Nt_in, B_zero);
    double RMS_in = computeRMS(B_in, Nt_in);
    double Peak2Peak = computePeak2Peak(B, Nt, B_zero);
    double RMS = computeRMS(B, Nt);

    ofstream outFile;
    string filename;

    filename = "Coils_" + IntToStr(NbCoils) + ".txt";

    writeOutfile(filename, design->getObjValue(), duration, Peak2Peak, RMS,
                 Peak2Peak_in, RMS_in, SOL, CoilsR, CoilsZ, widths, Gamma, Ns);

    printOutfile(SOL, Nc, Ns, CoilsR, CoilsZ, widths, Gamma,
                 design->getObjValue(), duration);

    // Compute the magnetic field matrix
    double B_MATRIX[61][66];
    computeMagneticFieldMatrix(SOL, Nc, Ns, CoilsR, CoilsZ, widths, Gamma,
                               cross_sections, 0.0, -0.3, 0.65, 0.3, 0.01,
                               B_MATRIX);

    // Write the magnetic field matrix to a CSV file
    matrixToExcel(B_MATRIX, "FieldTotal_Matrix.csv");

    cout << "CODE IS DONE!!" << endl;
  } else {
    cout << "Infeasible or unbounded!" << endl;
  }
  cout << "CODE IS DONE!!" << endl;
  return 0;
}
