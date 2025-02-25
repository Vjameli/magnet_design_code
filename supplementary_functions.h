#ifndef SUPPLEMENTARY_FUNCTIONS_H
#define SUPPLEMENTARY_FUNCTIONS_H

#include <sstream>
#include <string>
#include <vector>

// String conversion utility
std::string IntToStr(int n);

// Target point initialization
void initializeTargetPoints(int nbPoints, std::vector<double> &Xtarget,
                          std::vector<double> &Ytarget,
                          std::vector<double> &PolarR,
                          std::vector<double> &PolarTheta);

// Area calculations
void initializeAreaArrays(std::vector<double> &AREA_X,
                         std::vector<double> &AREA_Z, double X_init = 0.01,
                         double Z_init = -0.3);

#endif // SUPPLEMENTARY_FUNCTIONS_H
