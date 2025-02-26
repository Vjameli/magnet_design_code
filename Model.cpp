//*****************************************************************************
//        C++ * Class Model *
//                                                                            *
//        Iman Dayarian *
//                                                                            *
//        June 2015 *
//*****************************************************************************

#include "Model.h"

Model::Model(IloEnv envir, vector<vector<double>> AMN, vector<double> CoilsR,
             double B_zero, double tol) {
  env = envir;
  Nt = AMN.size();
  Nc = AMN[0].size();

  x = IloNumVarArray(env, Nc, -IloInfinity, IloInfinity, ILOFLOAT);
  t = IloNumVarArray(env, Nc, 0.0, IloInfinity, ILOFLOAT);

  model = IloModel(env);
  cplex = IloCplex(model);
  obj = IloMinimize(env);

  for (int m = 0; m < Nt; m++) {
    IloExpr exp(env);
    for (int n = 0; n < Nc; n++) {
      exp += AMN[m][n] * x[n];
    }
    model.add(exp <= B_zero * (1 + tol));
  }

  for (int m = 0; m < Nt; m++) {
    IloExpr exp(env);
    for (int n = 0; n < Nc; n++) {
      exp += AMN[m][n] * x[n];
    }
    model.add(-exp <= -B_zero * (1 - tol));
  }

  for (int n = 0; n < Nc; n++) {
    IloExpr exp(env);
    exp += x[n] - t[n];
    model.add(exp <= 0);
  }

  for (int n = 0; n < Nc; n++) {
    IloExpr exp(env);
    exp += -x[n] - t[n];
    model.add(exp <= 0);
  }

  // Objectif
  IloExpr exp(env);
  for (int n = 0; n < Nc; n++) {
    exp += CoilsR[n] * t[n];
  }

  obj.setExpr(exp);
  model.add(obj);

  cplex.setParam(IloCplex::SimDisplay, 0);
  cplex.setParam(IloCplex::Threads, 1);
}

Model::~Model() {
  cplex.exportModel("model.lp");
  env.end();
}

void Model::exportModel(char *file) { cplex.exportModel(file); }
