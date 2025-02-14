#ifndef MODEL_H
#define MODEL_H

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <list>
#include <string>
#include<vector>


using namespace std;

typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloArray<IloNumVarArray> > IloNumVarArray3;


class Model
{
  public :
	  Model(IloEnv env, vector<vector<double> > AMN, vector<double> CoilsR, double B_zero, double tol);
	  ~Model();
	  void exportModel(char *file);
	  int getCplexStatus() {return cplex.getCplexStatus();}
	  string nom(string deb,int num);
	  IloBool solve() {return cplex.solve();}
	  double getX(int num) {return cplex.getValue(x[num]);}	
	  double getT(int num) {return cplex.getValue(t[num]);}	  
  
	  double getObjValue() {return cplex.getObjValue();}

	  IloEnv env;
	  IloObjective obj;
	  IloNumVarArray  x;
	  IloNumVarArray  t;
	  IloNumVarArray2  y;
	  IloModel  model;
	  IloCplex cplex;

	  int Nt;
	  int Nc;

};

#endif
