#ifndef SUPPLEMENTARY_FUNCTIONS_H
#define SUPPLEMENTARY_FUNCTIONS_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

class Ellipse;
class BFieldThick;
class BFieldThick_Matrix;

std::string IntToStr(int n);

void initializeTargetPoints(int nbPoints, std::vector<double> &Xtarget,
                            std::vector<double> &Ytarget,
                            std::vector<double> &PolarR,
                            std::vector<double> &PolarTheta);

void sampleEllipseBorder(int nbPoints, std::vector<double> &Xtarget,
                         std::vector<double> &Ytarget,
                         std::vector<double> &PolarR,
                         std::vector<double> &PolarTheta);

void sampleEllipseInside(int nbPoints, std::vector<double> &Xtarget,
                         std::vector<double> &Ytarget,
                         std::vector<double> &PolarR,
                         std::vector<double> &PolarTheta,
                         std::vector<double> &Zt_in,
                         std::vector<double> &Xt_in);

void computeCoilCenters(std::vector<double> &NPZ, std::vector<double> &NPR,
                        std::vector<double> &SPZ, std::vector<double> &SPR,
                        double zmin, double zmax, double rmin, double rmax,
                        double step);

void computeCrossSections(std::vector<double> &cross_sections,
                          std::vector<double> &widths, double Gamma);

void computeAMNS(std::vector<std::vector<std::vector<double>>> &AMNS, int Nt,
                 int Nc, int Ns, std::vector<double> &Zt,
                 std::vector<double> &Xt, std::vector<double> &CoilsZ,
                 std::vector<double> &CoilsR,
                 std::vector<double> &cross_sections,
                 std::vector<double> &widths, double Gamma);

void AMNSToExcel(const std::vector<std::vector<std::vector<double>>> &AMNS,
                 int Nt, int Nc, int Ns);

void initialize2DSolVector(std::vector<std::vector<double>> &SOL, int Nc,
                           int Ns);

void evaluateSolutionBorder(
    const std::vector<std::vector<double>> &SOL,
    const std::vector<std::vector<std::vector<double>>> &AMNS, int Nt, int Nc,
    int Ns, double *B);

void evaluateSolutionInside(
    const std::vector<std::vector<double>> &SOL, int Nt_in, int Nc, int Ns,
    const std::vector<double> &Xt_in, const std::vector<double> &Zt_in,
    const std::vector<double> &CoilsR, const std::vector<double> &CoilsZ,
    const std::vector<double> &widths,
    const std::vector<double> &cross_sections, double Gamma, double *B_in);

double computePeak2Peak(const double *B, int Nt, double B_zero);

double computeRMS(const double *B, int Nt);

void writeOutfile(const std::string &filename, double objValue, double duration,
                  double Peak2Peak, double RMS, double Peak2Peak_in,
                  double RMS_in, const std::vector<std::vector<double>> &SOL,
                  const std::vector<double> &CoilsR,
                  const std::vector<double> &CoilsZ,
                  const std::vector<double> &widths, double Gamma, int Ns);

void printOutfile(const std::vector<std::vector<double>> &SOL, int Nc, int Ns,
                  const std::vector<double> &CoilsR,
                  const std::vector<double> &CoilsZ,
                  const std::vector<double> &widths, double Gamma,
                  double objValue, double duration);

void computeMagneticFieldMatrix(
    const std::vector<std::vector<double>> &SOL, int Nc, int Ns,
    const std::vector<double> &CoilsR, const std::vector<double> &CoilsZ,
    const std::vector<double> &widths, double Gamma,
    const std::vector<double> &cross_sections, double XMIN, double ZMIN,
    double XMAX, double ZMAX, double STEP, double B_MATRIX[61][66]);

void matrixToExcel(const double B_MATRIX[61][66], const std::string &filename);

#endif // SUPPLEMENTARY_FUNCTIONS_H
