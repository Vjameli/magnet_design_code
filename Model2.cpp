//*****************************************************************************
//        C++ * Class Model *
//                                                                            *
//        Iman Dayarian *
//                                                                            *
//        June 2015 *
//*****************************************************************************

#include "Model2.h"

Model2::Model2(bool surpress_cplex, IloEnv envir,
               vector<vector<vector<double>>> AMNS, vector<double> CoilsR,
               vector<double> CoilsZ, double B_zero, double tol,
               vector<double> widths, vector<double> cross_sections, double CCD,
               double Gamma, double zmin, double zmax, double rmin, double rmax,
               double step, int NbCoils) {
  env = envir;
  Nt = AMNS.size();
  Nc = AMNS[0].size();
  NbCS = widths.size();

  double L = 0.01;

  x = IloNumVarArray2(env, Nc);
  for (int n = 0; n < Nc; n++) {
    x[n] = IloNumVarArray(env, NbCS, -IloInfinity, IloInfinity, ILOFLOAT);
  }

  y = IloNumVarArray2(env, Nc);
  for (int n = 0; n < Nc; n++) {
    y[n] = IloNumVarArray(env, NbCS, 0, 1, ILOBOOL);
  }

  model = IloModel(env);
  cplex = IloCplex(model);
  obj = IloMinimize(env);

  for (int m = 0; m < Nt; m++) {
    IloExpr exp(env);
    for (int n = 0; n < Nc; n++) {
      for (int s = 0; s < NbCS; s++) {
        exp += AMNS[m][n][s] * x[n][s];
      }
    }
    model.add(exp <= B_zero * (1 + tol));
    exp.end();
  }

  for (int m = 0; m < Nt; m++) {
    IloExpr exp(env);
    for (int n = 0; n < Nc; n++) {
      for (int s = 0; s < NbCS; s++) {
        exp += AMNS[m][n][s] * x[n][s];
      }
    }
    model.add(-exp <= -B_zero * (1 - tol));
    exp.end();
  }

  // //Symmetry constraints
  for (int n = 0; n < Nc / 2; n++) {
    for (int s = 0; s < NbCS; s++) {
      IloExpr exp(env);
      exp = y[n][s] - y[n + Nc / 2][s];
      model.add(exp == 0);
      exp.end();
    }
  }

  // Maximum Nb Coils
  IloExpr expT(env);
  for (int n = 0; n < Nc / 2; n++) {
    for (int s = 0; s < NbCS; s++) {
      expT += y[n][s];
    }
  }
  model.add(expT <= NbCoils);
  expT.end();
  // //********************************

  for (int n = 0; n < Nc; n++) {
    IloExpr exp(env);

    for (int s = 0; s < NbCS; s++) {
      exp += y[n][s];
    }
    model.add(exp <= 1);
    exp.end();
  }

  //  // //*******************************

  for (int n = 0; n < Nc; n++) {
    for (int s = 0; s < NbCS; s++) {
      IloExpr exp(env);
      exp += x[n][s];
      exp -= CCD * 1e-7 * y[n][s] * cross_sections[s];
      model.add(exp <= 0);
      exp.end();
    }
  }

  for (int n = 0; n < Nc; n++) {
    for (int s = 0; s < NbCS; s++) {
      IloExpr exp(env);
      exp -= x[n][s];
      exp -= CCD * 1e-7 * y[n][s] * cross_sections[s];
      model.add(exp <= 0);
      exp.end();
    }
  }

  // //********************************

  for (int n1 = 0; n1 < Nc / 2 - 1; n1++) {
    for (int n2 = n1 + 1; n2 < Nc / 2; n2++) {
      for (int s1 = 0; s1 < NbCS; s1++) {
        for (int s2 = 0; s2 < NbCS; s2++) {
          if (abs(CoilsZ[n1] - CoilsZ[n2]) <
                  widths[s1] / (Gamma * 2) + widths[s2] / (2 * Gamma) + L &&
              abs(CoilsR[n1] - CoilsR[n2]) <
                  widths[s1] / 2 + widths[s2] / 2 + L) {
            IloExpr exp1(env);
            IloExpr exp2(env);
            exp1 += y[n1][s1] + y[n2][s2];
            exp2 += y[n1][s2] + y[n2][s1];
            model.add(exp1 <= 1);
            exp1.end();
            model.add(exp2 <= 1);
            exp2.end();
          }
        }
      }
    }
  }

  for (int n1 = Nc / 2; n1 < Nc - 1; n1++) {
    for (int n2 = n1 + 1; n2 < Nc; n2++) {
      for (int s1 = 0; s1 < NbCS; s1++) {
        for (int s2 = 0; s2 < NbCS; s2++) {
          if (abs(CoilsZ[n1] - CoilsZ[n2]) <
                  widths[s1] / (2 * Gamma) + widths[s2] / (2 * Gamma) + L &&
              abs(CoilsR[n1] - CoilsR[n2]) <
                  widths[s1] / 2 + widths[s2] / 2 + L) {
            IloExpr exp1(env);
            IloExpr exp2(env);
            exp1 += y[n1][s1] + y[n2][s2];
            exp2 += y[n1][s2] + y[n2][s1];
            model.add(exp1 <= 1);
            model.add(exp2 <= 1);
            exp1.end();
            exp2.end();
          }
        }
      }
    }
  }

  //********************************

  // Objectif
  IloExpr exp(env);
  for (int n = 0; n < Nc; n++) {
    for (int s = 0; s < NbCS; s++) {
      exp += 1e4 * CoilsR[n] * y[n][s] * cross_sections[s];
    }
  }

  obj.setExpr(exp);
  exp.end();
  model.add(obj);

  //    cplex.setParam(IloCplex::NodeFileInd,1);
  cplex.setParam(IloCplex::EpGap, 10e-9);
  cplex.setParam(IloCplex::EpRHS, 10e-9);
  cplex.setParam(IloCplex::EpInt, 10e-9);
  cplex.setParam(IloCplex::EpMrk, 0.99);
  cplex.setParam(IloCplex::NumericalEmphasis, 1);
  cplex.setParam(IloCplex::SimDisplay, 0);
  cplex.setParam(IloCplex::Threads, 1);

  if (surpress_cplex) {
    cplex.setOut(env.getNullStream());
  }
}

Model2::~Model2() {
  cplex.exportModel("model.sav");
  env.end();
}

void Model2::exportModel(char *file) { cplex.exportModel(file); }

void Model2::AddOverlappingRow(int coil1, int coil2, int corss1, int cross2) {
  IloExpr exp1(env);
  IloExpr exp2(env);
  exp1 += y[coil1][corss1] + y[coil2][cross2] + z[coil1] + z[coil2];
  exp2 += y[coil1][cross2] + y[coil2][corss1] + z[coil1] + z[coil2];
  model.add(exp1 <= 3);
  exp1.end();
  model.add(exp2 <= 3);
  exp2.end();
}

void Model2::getSol(vector<vector<double>> &sol) {
  for (int n = 0; n < Nc; n++) {
    for (int s = 0; s < NbCS; s++) {
      sol[n][s] = cplex.getValue(x[n][s]);
    }
  }
}

int Model2::getCrossIndex(int coil) {

  for (int s = 0; s < NbCS; s++) {
    if (cplex.getValue(y[coil][s])) {
      return s;
    }
  }
  return -1;
}
