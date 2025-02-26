#include "supplementary_functions.h"
#include "BField.h"
#include "BFieldThick.h"
#include "BFieldThick_Matrix.h"
#include "Ellipse.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
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

void sampleEllipseBorder(int nbPoints, std::vector<double> &Xtarget,
                         std::vector<double> &Ytarget,
                         std::vector<double> &PolarR,
                         std::vector<double> &PolarTheta) {

  Ellipse *ellipse;
  ellipse = new Ellipse(0, 0, 0.2, 0.2, 0, nbPoints, Xtarget, Ytarget, PolarR,
                        PolarTheta);
}

void sampleEllipseInside(int nbPoints, std::vector<double> &Xtarget,
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

void computeCoilCenters(vector<double> &NPZ, vector<double> &NPR,
                        vector<double> &SPZ, vector<double> &SPR, double zmin,
                        double zmax, double rmin, double rmax, double step) {
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
}

void computeCrossSections(vector<double> &cross_sections,
                          vector<double> &widths, double Gamma) {
  for (int value = 50; value <= 100; value += 5) {
    widths.push_back(value / 1000.0);
  }

  for (unsigned int i = 0; i < widths.size(); i++) {
    cross_sections.push_back(widths[i] * widths[i] / Gamma);
  }
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

void AMNSToExcel(const vector<vector<vector<double>>> &AMNS, int Nt, int Nc,
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

void initialize2DSolVector(vector<vector<double>> &SOL, int Nc, int Ns) {
  SOL.resize(Nc);
  for (int n = 0; n < Nc; n++) {
    SOL[n].resize(Ns);
    for (int s = 0; s < Ns; s++) {
      SOL[n][s] = 0;
    }
  }
}

void evaluateSolutionBorder(const vector<vector<double>> &SOL,
                            const vector<vector<vector<double>>> &AMNS, int Nt,
                            int Nc, int Ns, double *B) {
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

void evaluateSolutionInside(const vector<vector<double>> &SOL, int Nt_in,
                            int Nc, int Ns, const vector<double> &Xt_in,
                            const vector<double> &Zt_in,
                            const vector<double> &CoilsR,
                            const vector<double> &CoilsZ,
                            const vector<double> &widths,
                            const vector<double> &cross_sections, double Gamma,
                            double *B_in) {
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

double computePeak2Peak(const double *B, int Nt, double B_zero) {
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

double computeRMS(const double *B, int Nt) {
  double sum_square = 0;
  for (int m = 0; m < Nt; m++) {
    sum_square += B[m] * B[m];
  }
  return sqrt(sum_square / Nt);
}

void writeOutfile(const string &filename, double objValue, double duration,
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

void printOutfile(const vector<vector<double>> &SOL, int Nc, int Ns,
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

void computeMagneticFieldMatrix(const vector<vector<double>> &SOL, int Nc,
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

void matrixToExcel(const double B_MATRIX[61][66], const string &filename) {
  ofstream ExcelFile;
  ExcelFile.open(filename);
  for (int i = 0; i < 61; i++) {
    for (int j = 0; j < 66; j++)
      ExcelFile << setprecision(9) << B_MATRIX[i][j] << ", ";
    ExcelFile << endl;
  }
  ExcelFile.close();
}
