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

std::string IntToStr(int n) {
  std::stringstream result;
  result << n;
  return result.str();
}

void initializeTargetPoints(int nbPoints, std::vector<double> &Xtarget,
                            std::vector<double> &Ytarget,
                            std::vector<double> &PolarR,
                            std::vector<double> &PolarTheta) {
  // Resize vectors
  Xtarget.resize(nbPoints);
  Ytarget.resize(nbPoints);
  PolarR.resize(nbPoints);
  PolarTheta.resize(nbPoints);

  // Initialize arrays with zeros
  for (int i = 0; i < nbPoints; i++) {
    Xtarget[i] = 0;
    Ytarget[i] = 0;
    PolarR[i] = 0;
    PolarTheta[i] = 0;
  }
}

void sample_ellipse_border(int nbPoints, std::vector<double> &Xtarget,
                           std::vector<double> &Ytarget,
                           std::vector<double> &PolarR,
                           std::vector<double> &PolarTheta) {

  Ellipse *ellipse;
  ellipse = new Ellipse(0, 0, 0.2, 0.2, 0, nbPoints, Xtarget, Ytarget, PolarR,
                        PolarTheta);
}

void sample_ellipse_inside(int nbPoints, std::vector<double> &Xtarget,
                           std::vector<double> &Ytarget,
                           std::vector<double> &PolarR,
                           std::vector<double> &PolarTheta,
                           std::vector<double> &Zt_in,
                           std::vector<double> &Xt_in) {

  Ellipse *ellipse2;
  ellipse2 = new Ellipse(0, 0, 0.2, 0.2, 0, nbPoints, Xtarget, Ytarget, PolarR,
                         PolarTheta);
  Zt_in.insert(Zt_in.end(), Ytarget.begin(), Ytarget.end());
  Xt_in.insert(Xt_in.end(), Xtarget.begin(), Xtarget.end());
}

void computeAMNS(vector<vector<vector<double>>> &AMNS, int Nt, int Nc, int Ns,
                 vector<double> &Zt, vector<double> &Xt, vector<double> &CoilsZ,
                 vector<double> &CoilsR, vector<double> &cross_sections,
                 vector<double> &widths, double Gamma) {

  // Initialize AMNS matrix
  AMNS.resize(Nt);
  for (int i = 0; i < Nt; i++) {
    AMNS[i].resize(Nc);
    for (int j = 0; j < Nc; j++) {
      AMNS[i][j].resize(Ns);
    }
  }

  // Compute BField and fill AMNS matrix
  BFieldThick *field;
  for (int m = 0; m < Nt; m++) {
    for (int n = 0; n < Nc; n++) {
      for (int s = 0; s < Ns; s++) {
        double Bz = 0;
        double Br = 0;
        double J = (1e7) / cross_sections[s];
        field = new BFieldThick(J, CoilsR[n], widths[s], widths[s] / Gamma,
                                Xt[m], Zt[m], 0, CoilsZ[n], Bz, Br);

        AMNS[m][n][s] = Bz;
      }
    }
  }
}

void AMNS_to_excel(const vector<vector<vector<double>>> &AMNS, int Nt, int Nc,
                   int Ns) {
  ofstream ExcelFile3;
  ExcelFile3.open("AMNS.csv");

  for (int m = 0; m < Nt; m++) {
    for (int n = 0; n < Nc; n++) {
      for (int s = 0; s < Ns; s++) {
        ExcelFile3 << AMNS[m][n][s] << ",";
      }
    }
    ExcelFile3 << endl;
  }
  ExcelFile3.close();
}

void evaluate_solution_border(const vector<vector<double>> &SOL,
                              const vector<vector<vector<double>>> &AMNS,
                              int Nt, int Nc, int Ns, double *B) {
  // Initialize B array with zeros
  for (int i = 0; i < Nt; i++) {
    B[i] = 0;
  }

  // Compute the resultant magnetic field from the solution on the border of the
  // FOV
  for (int n = 0; n < Nc; n++) {
    for (int s = 0; s < Ns; s++) {
      if (SOL[n][s] * 1e7 > 1 || -SOL[n][s] * 1e7 > 1) {
        for (int m = 0; m < Nt; m++) {
          double Bz = AMNS[m][n][s] * SOL[n][s];
          B[m] = B[m] + Bz;
        }
      }
    }
  }
}

void evaluate_solution_inside(const vector<vector<double>> &SOL, int Nt_in,
                              int Nc, int Ns, const vector<double> &Xt_in,
                              const vector<double> &Zt_in,
                              const vector<double> &CoilsR,
                              const vector<double> &CoilsZ,
                              const vector<double> &widths,
                              const vector<double> &cross_sections,
                              double Gamma, double *B_in) {
  // Initialize B_in array
  for (int i = 0; i < Nt_in; i++)
    B_in[i] = 0;

  // Compute the resultant magnetic field from the solution inside the FOV
  BFieldThick *field2;
  for (int m = 0; m < Nt_in; m++) {
    for (int n = 0; n < Nc; n++) {
      for (int s = 0; s < Ns; s++) {
        if (SOL[n][s] * 1e7 > 1 || -SOL[n][s] * 1e7 > 1) {
          double Bz = 0;
          double Br = 0;
          double J = SOL[n][s] * 1e7 / cross_sections[s];
          field2 = new BFieldThick(J, CoilsR[n], widths[s], widths[s] / Gamma,
                                   Xt_in[m], Zt_in[m], 0, CoilsZ[n], Bz, Br);
          B_in[m] = B_in[m] + Bz;
        }
      }
    }
  }
}

double compute_Peak2peak(const double *B, int Nt, double B_zero) {
  double min_B = 10;
  double max_B = -10;
  for (int m = 0; m < Nt; m++) {
    if (B[m] < min_B)
      min_B = B[m];
    if (B[m] > max_B)
      max_B = B[m];
  }
  return 1e6 * ((max_B - B_zero) + (B_zero - min_B));
}

double compute_RMS(const double *B, int Nt) {
  double sum_square = 0;
  for (int m = 0; m < Nt; m++) {
    sum_square += B[m] * B[m];
  }
  return sqrt(sum_square / Nt);
}

void write_outFile(const string &filename, double objValue, double duration,
                   double Peak2Peak, double RMS, double Peak2Peak_in,
                   double RMS_in, const vector<vector<double>> &SOL,
                   const vector<double> &CoilsR, const vector<double> &CoilsZ,
                   const vector<double> &widths, double Gamma, int Ns) {
  ofstream outFile;
  outFile.open(filename.c_str());
  outFile << "OBJECTIVE VALUE: " << objValue << endl;
  outFile << "Total Duration: " << duration << endl;
  outFile << "Peak2Peak on the surface of the volums: " << Peak2Peak << endl;
  outFile << "RMS on the surface of the volums: " << RMS << endl;
  outFile
      << "Peak2Peak over the sampling inside and on the surface of the volums: "
      << Peak2Peak_in << endl;
  outFile << "RMS over the sampling inside and on the surface of the volums: "
          << RMS_in << endl;
  for (unsigned int n = 0; n < CoilsR.size(); n++) {
    for (int s = 0; s < Ns; s++) {
      if (SOL[n][s] * 1e7 > 1 || -SOL[n][s] * 1e7 > 1) {
        outFile << n << ", " << CoilsR[n] << ", " << CoilsZ[n] << ", "
                << setprecision(9) << SOL[n][s] * 1e7 << ", " << widths[s]
                << ", " << widths[s] / Gamma << endl;
      }
    }
  }
  outFile.close();
}

void print_outFile(const vector<vector<double>> &SOL, int Nc, int Ns,
                   const vector<double> &CoilsR, const vector<double> &CoilsZ,
                   const vector<double> &widths, double Gamma, double objValue,
                   double duration) {
  cout << "OBJECTIVE VALUE: " << objValue << endl;
  cout << "*** Total Duration: " << duration << endl;
  for (int n = 0; n < Nc; n++) {
    for (int s = 0; s < Ns; s++) {
      if (SOL[n][s] * 1e7 > 1 || -SOL[n][s] * 1e7 > 1) {
        cout << "N: " << n << ", Zn: " << CoilsZ[n] << ", Rn: " << CoilsR[n]
             << ", Current: " << SOL[n][s] * 1e7 << ", W: " << widths[s]
             << ", H: " << widths[s] / Gamma << endl;
      }
    }
  }
}

void compute_magnetic_field_matrix(const vector<vector<double>> &SOL, int Nc,
                                   int Ns, const vector<double> &CoilsR,
                                   const vector<double> &CoilsZ,
                                   const vector<double> &widths, double Gamma,
                                   const vector<double> &cross_sections,
                                   double XMIN, double ZMIN, double XMAX,
                                   double ZMAX, double STEP,
                                   double B_MATRIX[61][66]) {
  vector<vector<double>> BzMatrix(61, vector<double>(66, 0));
  vector<vector<double>> BxMatrix(61, vector<double>(66, 0));

  // Initialize B_MATRIX
  for (int i = 0; i < 61; i++)
    for (int j = 0; j < 66; j++)
      B_MATRIX[i][j] = 0;

  for (int n = 0; n < Nc; n++) {
    for (int s = 0; s < Ns; s++) {
      if (SOL[n][s] * 1e7 > 1 || -SOL[n][s] * 1e7 > 1) {
        // Reset BzMatrix and BxMatrix
        for (int i = 0; i < 61; i++)
          for (int j = 0; j < 66; j++) {
            BzMatrix[i][j] = 0;
            BxMatrix[i][j] = 0;
          }

        double J = SOL[n][s] * 1e7 / cross_sections[s];

        BFieldThick_Matrix *field_Matrix = new BFieldThick_Matrix(
            J, CoilsR[n], widths[s], widths[s] / Gamma, XMIN, ZMIN, XMAX, ZMAX,
            0, CoilsZ[n], STEP, BzMatrix, BxMatrix);

        // Accumulate results in B_MATRIX
        for (int i = 0; i < 61; i++)
          for (int j = 0; j < 66; j++)
            B_MATRIX[i][j] += BzMatrix[i][j];

        delete field_Matrix;
      }
    }
  }
}

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
  sample_ellipse_border(nbPoints, Xtarget, Ytarget, PolarR, PolarTheta);

  vector<double> Zt;
  vector<double> Rt;
  vector<double> Xt;
  Zt.resize(nbPoints);
  Rt.resize(nbPoints);
  Xt.resize(nbPoints);

  // Coordinates of the points in cilyndrical coordinates
  Zt = Ytarget;
  Rt = PolarR;

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

  sample_ellipse_inside(nbPoints, Xtarget, Ytarget, PolarR, PolarTheta, Zt_in,
                        Xt_in);

  int Nt_in = Zt_in.size();

  // Coordinates of the cente of the coils based on the step (IN METERS!)
  vector<double> NPZ; // North pole Z
  vector<double> NPR; // North pole R
  vector<double> SPZ; // South pole Z
  vector<double> SPR; // South pole R

  for (int zc0 = zmin; zc0 <= zmax;) {
    for (int rc0 = rmin; rc0 <= rmax;) {
      NPZ.push_back(zc0 / 100.0);
      NPR.push_back(rc0 / 100.0);
      SPZ.push_back(-zc0 / 100.0);
      SPR.push_back(rc0 / 100.0);
      rc0 = rc0 + step;
    }
    zc0 = zc0 + step;
  }

  // CoilsZ = concatenation of (NPZ + SPZ)
  // CoilsR = concatenation of (NPR + SPR)
  vector<double> CoilsZ;
  vector<double> CoilsR;
  CoilsZ = NPZ;
  CoilsZ.insert(CoilsZ.end(), SPZ.begin(), SPZ.end());
  CoilsR = NPR;
  CoilsR.insert(CoilsR.end(), SPR.begin(), SPR.end());
  int Nc = CoilsZ.size();

  // Set of potential cross sections (Potential widths)
  // VALUE WAS ORIGINALLY EQUAL TO 5!!!!!!!!!! <-------
  vector<double> widths;
  for (int value = 50; value <= 100; value += 5)
    widths.push_back(value / 1000.0);

  vector<double> cross_sections;

  for (unsigned int i = 0; i < widths.size(); i++) {
    cross_sections.push_back(widths[i] * widths[i] / Gamma);
  }

  int Ns = cross_sections.size();

  // Compute BField and AMNS matrix
  vector<vector<vector<double>>> AMNS;
  computeAMNS(AMNS, Nt, Nc, Ns, Zt, Xt, CoilsZ, CoilsR, cross_sections, widths,
              Gamma);
  AMNS_to_excel(AMNS, Nt, Nc, Ns);

  // Initialize vectors for optimization model and cpu time
  Model2 *design;
  double duration = 0;

  int NbCoils = 7;
  vector<double> CPUs;
  vector<double> OBJECTIVEs;
  vector<int> NBCOILS;

  // Initialize 2D SOL vector
  vector<vector<double>> SOL;
  SOL.resize(Nc);
  for (int n = 0; n < Nc; n++) {
    SOL[n].resize(Ns);
    for (int s = 0; s < Ns; s++) {
      SOL[n][s] = 0;
    }
  }

  vector<vector<vector<double>>> SOL_VECTOR;

  // Here, I removed the while loop for now.
  cout << "************** NBCOILS = " << NbCoils << " **************" << endl;
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
    SOL_VECTOR.push_back(SOL);
    CPUs.push_back(duration);
    OBJECTIVEs.push_back(design->getObjValue());
    NBCOILS.push_back(NbCoils);

    ofstream outFile;
    string filename;

    filename = "Coils_" + IntToStr(NbCoils) + ".txt";
    cout << filename << "  \n";

    double B[Nt];
    evaluate_solution_border(SOL, AMNS, Nt, Nc, Ns, B);

    for (int m = 0; m < Nt; m++) {
      cout << m << " " << B[m] << endl;
      if (B[m] > B_zero * (1 + tol) || B[m] < B_zero * (1 - tol))
        cout << "******" << endl;
    }

    // Compute the resultant magnetic field from the solution inside the FOV.
    double B_in[Nt_in];
    evaluate_solution_inside(SOL, Nt_in, Nc, Ns, Xt_in, Zt_in, CoilsR, CoilsZ,
                             widths, cross_sections, Gamma, B_in);

    // Peak to peak on the border and inside of the FOV
    // RMS on the border and inside of the FOV
    double Peak2Peak_in = compute_Peak2peak(B_in, Nt_in, B_zero);
    double RMS_in = compute_RMS(B_in, Nt_in);
    double Peak2Peak = compute_Peak2peak(B, Nt, B_zero);
    double RMS = compute_RMS(B, Nt);

    write_outFile(filename, design->getObjValue(), duration, Peak2Peak, RMS,
                  Peak2Peak_in, RMS_in, SOL, CoilsR, CoilsZ, widths, Gamma, Ns);
  }

  print_outFile(SOL, Nc, Ns, CoilsR, CoilsZ, widths, Gamma,
                design->getObjValue(), duration);

  // Compute the magnetic field matrix
  double B_MATRIX[61][66];
  compute_magnetic_field_matrix(SOL, Nc, Ns, CoilsR, CoilsZ, widths, Gamma,
                                cross_sections, 0.0, -0.3, 0.65, 0.3, 0.01,
                                B_MATRIX);

  // Write the magnetic field matrix to a CSV file
  ofstream ExcelFile2;
  ExcelFile2.open("FieldTotal_Matrix.csv");
  for (int i = 0; i < 61; i++) {
    for (int j = 0; j < 66; j++)
      ExcelFile2 << setprecision(9) << B_MATRIX[i][j] << ", ";
    ExcelFile2 << endl;
  }
  ExcelFile2.close();

  // Do the same thing over again for the final solution of the optimization
  cout << "CODE IS DONE!!" << endl;
  return 0;
}
