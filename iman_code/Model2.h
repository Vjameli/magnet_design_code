#ifndef MODEL2_H
#define MODEL2_H

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <list>
#include <string>
#include<vector>
#include <cmath>



using namespace std;

typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloArray<IloNumVarArray> > IloNumVarArray3;


class Model2
{
  public :
	  Model2(IloEnv envir, vector<vector<vector<double> > > AMNS, 
vector<double> CoilsR, vector<double> CoilsZ, double B_zero, 
double tol, vector<double> widths, vector<double> cross_sections, 
double CCD, double Gamma, double zmin, double zmax, double rmin, double rmax, double step, int NbCoils);
	  ~Model2();
	  void exportModel(char *file);
	  void AddOverlappingRow(int coil1, int coil2, int corss1, int cross2);
	  int getCplexStatus() {return cplex.getCplexStatus();}
	  string nom(string deb,int num);
	  IloBool solve() {return cplex.solve();}
	  double getX(int num1, int num2) {return cplex.getValue(x[num1][num2]);}	
	  //double getT(int num) {return cplex.getValue(t[num]);}
	  double getH(int num) {return cplex.getValue(h[num]);}
	  double getW(int num) {return cplex.getValue(w[num]);}	 
	  double getZ(int num) {return cplex.getValue(z[num]);}
	  int getCrossIndex(int coil); 
	  void getSol(vector<vector <double> >& sol);
  
	  double getObjValue() {return cplex.getObjValue();}

	  IloEnv env;
	  IloObjective obj;
	  IloNumVarArray2  x;
	  IloNumVarArray  w;
	  IloNumVarArray  h;
	  IloNumVarArray  z;
	  IloNumVarArray2  t;
	  IloNumVarArray2  y;
	  IloNumVarArray2  rnn;
	  IloNumVarArray2  znn;
	  IloModel  model;
	  IloCplex cplex;

	  int Nt;
	  int Nc;
	  int NbCS;	//Nb Cross-sections

};

#endif
