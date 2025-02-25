#include "supplementary_functions.h"

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

void initializeAreaArrays(std::vector<double> &AREA_X,
                          std::vector<double> &AREA_Z, double X_init,
                          double Z_init) {
  // Initialize X area points
  AREA_X.clear();
  for (int i = 1; i <= 40; i++) {
    AREA_X.push_back(X_init * i);
  }

  // Initialize Z area points
  AREA_Z.clear();
  for (int i = 0; i < 61; i++) {
    AREA_Z.push_back(Z_init + i * 0.01);
  }
}
