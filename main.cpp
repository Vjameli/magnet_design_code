//*****************************************************************************
//        Iman Dayarian  	                                                  *
//        Toronto, Canada       	                                          *
//                                                                            *
//        Generated June 2015                                                 *
//*****************************************************************************

#include<ilcplex/ilocplex.h>

ILOSTLBEGIN

#include <math.h>
#include "BField.h"
#include "BFieldThick.h"
#include "BFieldThick_Matrix.h"
#include "Ellipse.h"
#include "Model2.h"
#include <cstdlib>
#include <list>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <values.h>
#include <vector>
#include "tempsC++.h"
#include <cmath>
#include <iomanip>      // std::setprecision



using namespace std;

string IntToStr(int n) 
{
    stringstream result;
    result << n;
    return result.str();
}


int main()
{
	

	double B_zero = 0.5;    // [T] - Magnetic field to achieve
	double tol = 10e-6;    // Magnetic field homogeneity tolerance

	double Gamma = 1.4;
	double CCD = 21e+7;

	
	// ***** Define candidate coils: Initial grid *****
	// North pole
	double zmin = 34;	// [cm] z-domain
	double zmax = 42;
	double rmin = 10;	// [cm] r-domain
	double rmax = 48;
	double step = 2; // BIGGER STEP -> LESS COILS (it was 2 before)

	
	
	// ***** Define the Nt target points *****
	int nbPoints = 40; //IM PUTTING A SMALL VALUE HERE JUST FOR TESTING! (it was 40 before)
	vector<double> Xtarget;
	vector<double> Ytarget;
	vector<double> PolarR;
	vector<double> PolarTheta;
	Xtarget.resize(nbPoints);
	Ytarget.resize(nbPoints);
	PolarR.resize(nbPoints);
	PolarTheta.resize(nbPoints);
	
	for (int i = 0; i< nbPoints; i++){
		Xtarget[i] = 0;
		Ytarget[i] = 0;
		PolarR[i] = 0;
		PolarTheta[i] = 0;
	}

	vector<double> AREA_X;
	vector<double> AREA_Z;
	double X_init = 0.01;
	for (int i = 1; i<=40 ; i++)
		AREA_X.push_back(X_init*i);
	
	double Z_init = - 0.3;
	for (int i = 0; i<61 ; i++)
		AREA_Z.push_back(Z_init + i*0.01);

	
	Ellipse *ellipse;
	ellipse = new Ellipse(0, 0, 0.2, 0.2, 0, nbPoints, Xtarget, Ytarget, PolarR, PolarTheta);
	
	for (int i = 0; i< nbPoints; i++){
//		cout<< i <<", X: "<<Xtarget[i]<<", Y: "<<Ytarget[i]<<", PolarR: "<<PolarR[i]<< ", PolarTheta: "<<PolarTheta[i]<<endl;
	}
	vector<double> Zt;
	vector<double> Rt;
	vector<double> Xt;
	Zt.resize(nbPoints);
	Rt.resize(nbPoints);
	Xt.resize(nbPoints);
	
	Zt = Ytarget;
	Rt = PolarR;
	
	for(int i = 0; i< nbPoints; i++)
		Xt[i]= abs(Xtarget[i]);
	
	int Nt = Zt.size();	
    //cout << "Nt: " << Nt << endl;
	

	// Sampling a set of point inside the target volume
	Ellipse *ellipse2;
	vector<double> Zt_in;
	vector<double> Xt_in;
	Zt_in.push_back(0);
	Xt_in.push_back(0);
	for (int k = 1; k < 5; k++){
        ellipse2 = new Ellipse(0, 0, k * 0.25 * 0.2, k * 0.25 * 0.2, 0, nbPoints, Xtarget, Ytarget, PolarR, PolarTheta);	
        Zt_in.insert(Zt_in.end(), Ytarget.begin(), Ytarget.end());
        Xt_in.insert(Xt_in.end(), Xtarget.begin(), Xtarget.end());	
	}
	int Nt_in = Zt_in.size();

	
	vector<double> NPZ;		//North pole Z
	vector<double> NPR;		//North pole R
	
	vector<double> SPZ;		//South pole Z
	vector<double> SPR;		//South pole R
	
	for(int zc0 = zmin; zc0 <= zmax; ){
		for(int rc0 = rmin; rc0 <= rmax; ){
		NPZ.push_back(zc0/100.0);
		NPR.push_back(rc0/100.0);
		SPZ.push_back(-zc0/100.0);
		SPR.push_back(rc0/100.0);
		rc0 = rc0 + step;		
		}
	zc0 = zc0 + step;	
	}


	vector<double> CoilsZ;
	vector<double> CoilsR;
	CoilsZ = NPZ;
	CoilsZ.insert( CoilsZ.end(), SPZ.begin(), SPZ.end() );
	CoilsR = NPR;
	CoilsR.insert( CoilsR.end(), SPR.begin(), SPR.end() );
	int Nc = CoilsZ.size();
	
	
	//Set of potential cross sections (Potential widths)
	vector<double> widths;
	for(int value = 5; value <= 100; value += 5)
		widths.push_back(value/1000.0);
    

	vector<double> cross_sections;
	
	for(unsigned int i = 0; i< widths.size(); i++){
	cross_sections.push_back(widths[i]*widths[i]/Gamma);
	}
	
	int Ns = cross_sections.size();
	

	vector<vector<vector<double> > > AMNS;
	AMNS.resize(Nt);
	for(int i = 0; i < Nt; i++){
		AMNS[i].resize(Nc);
		for(int j = 0; j<Nc; j++){
			AMNS[i][j].resize(Ns);
			for(int k = 0; k<Ns; k++){
				AMNS[i][j][k] = 0;
			}
		}
	}
	
	
	// ifstream file ( "csvlist-0.005-0.1.csv" ); 
	// string value;
	// int Mcmp = 0;
	// int Ncmp = 0;
	// int Scmp = 0;
	// int cmp = 0;
	// while ( file.good() )
	// {
     // getline ( file, value, ',' ); 
	 // double value_new = atof(value.c_str());
		// if(cmp == Nc * Ns ) 
		// {
		// Mcmp++;
		// cmp = 0;
		// Ncmp = 0;
		// Scmp = 0;
		// }
		// if(cmp == (Ncmp+1) * Ns)
		// {
		// Ncmp++;
		// Scmp = 0;
		// }
	// AMNS[Mcmp][Ncmp][Scmp] = value_new;
	// cmp++;
	// Scmp++;
	// }


	
 	BFieldThick *field;
	for(int m = 0; m < Nt; m++){
		for(int n = 0; n < Nc; n++){
			for(int s = 0; s < Ns; s++){
				double Bz = 0;
				double Br = 0;
				double J = (1e7)/cross_sections[s];
				field = new BFieldThick(J, CoilsR[n], widths[s], widths[s]/Gamma, Xt[m], Zt[m], 0, CoilsZ[n], Bz, Br);

				AMNS[m][n][s] = Bz;

                if (m == 1 ){
                }
			}
		}
	} 
	
	delete field; 	
	
    cout << "Excel stuff" << endl;

	ofstream ExcelFile3;
    ExcelFile3.open("AMNS.csv");
	 
	for(int m = 0; m < Nt; m++){
		for(int n = 0; n < Nc; n++){
			for(int s = 0; s < Ns; s++){
				ExcelFile3 <<AMNS[m][n][s] <<", ";
			}
		}
		ExcelFile3 <<endl;
	}
      ExcelFile3.close();



	Model2* design;
	double duration = 0;


//bool isFeasible = 1;
int NbCoils = 7;
vector<double> CPUs;
vector<double> OBJECTIVEs;
vector<int> NBCOILS;


vector<vector<double> > SOL;
SOL.resize(Nc);
for(int n = 0; n<Nc; n++){
	SOL[n].resize(Ns);
		for(int s = 0; s<Ns; s++){
			SOL[n][s] = 0;
		}
	}
		
vector<vector<vector<double> > > SOL_VECTOR;

		
		
do {
    cout <<"************** NBCOILS = "<< NbCoils <<" **************"<<endl;
	duration = 0;
    ChronoCPU  chrono;
    chrono.start();
    
	design = new Model2(IloEnv(), AMNS, CoilsR, CoilsZ, B_zero, tol, widths, cross_sections, CCD, Gamma, zmin, zmax, rmin, rmax, step, NbCoils);


	//design->exportModel("bbb11.lp");
	design->solve();
	
	//chrono.stop();
    duration = chrono.getTime();
    
	if(design->getCplexStatus() != IloAlgorithm::Infeasible  &&
            design->getCplexStatus() != IloAlgorithm::InfeasibleOrUnbounded){
            
            design->getSol(SOL);
            SOL_VECTOR.push_back(SOL);
            CPUs.push_back(duration);
    		OBJECTIVEs.push_back(design->getObjValue());
    		NBCOILS.push_back(NbCoils);
		
		
		ofstream outFile;
    	string filename;

        filename="Coils_" + IntToStr(NbCoils) +".txt";
        cout<< filename << "  \n";


 		double B[Nt];
		for(int i = 0; i < Nt; i++)
				B[i] = 0;
	 
	 	 for(int n = 0; n < Nc; n++){
		 for(int s = 0; s < Ns; s++){
			 if (SOL[n][s] * 1e7 > 1 || -SOL[n][s] * 1e7> 1){
				 for(int m = 0; m < Nt; m++){
					 
				double Bz = AMNS[m][n][s] * SOL[n][s] ;
				// double Br = 0;
				
				// double J = SOL[n][s] * 10e7/cross_sections[s];
				// field = new BFieldThick(J, CoilsR[n], widths[s], widths[s]/Gamma, Xt[m], Zt[m], 0, CoilsZ[n], Bz, Br);
				B[m] = B[m] + Bz;
			 }
			 }
		 }
	 } 

	for(int m = 0; m < Nt; m++){
		cout<<m<<" "<<B[m]<<endl;
		if(B[m] > B_zero * (1 + tol) || B[m] < B_zero * (1 - tol))
			cout<<"******"<<endl;
	}


	BFieldThick *field2;
	double B_in[Nt_in];
	for(int i = 0; i < Nt_in; i++)
				B_in[i] = 0;
	for(int m = 0; m < Nt_in; m++){
		for(int n = 0; n < Nc; n++){
			for(int s = 0; s < Ns; s++){
				if (SOL[n][s]* 1e7> 1 || -SOL[n][s]* 1e7> 1){
				double Bz = 0;
				double Br = 0;
				double J = SOL[n][s] * 1e7/cross_sections[s];
				field2 = new BFieldThick(J, CoilsR[n], widths[s], widths[s]/Gamma, Xt_in[m], Zt_in[m], 0, CoilsZ[n], Bz, Br); 
				B_in[m] = B_in[m] + Bz;
				}
			}
		}
	}


	double min_B = 10;
	double max_B = -10;
	double sum_square = 0;
	for(int m = 0; m < Nt; m++){	
		if(B[m] < min_B) min_B = B[m];
		if(B[m] > max_B) max_B = B[m];	
		sum_square += B[m]*B[m];	
	}

	//Peak to Peak over the sample points inside the volume
	//and also the computation of RMS
	double min_B_in = 10;
	double max_B_in = -10;
	double sum_square_in = 0;
	for(int m = 0; m < Nt_in; m++){	
		if(B_in[m] < min_B_in) min_B_in = B_in[m];
		if(B_in[m] > max_B_in) max_B_in = B_in[m];	
		sum_square_in += B_in[m]*B_in[m];	
	}
	
	double Peak2Peak_in = 1e6*((max_B_in - 0.5) + (0.5 - min_B_in));
	double RMS_in = sqrt(sum_square_in/Nt_in);
	double Peak2Peak = 1e6*((max_B - 0.5) + (0.5 - min_B));
	double RMS = sqrt(sum_square/Nt);
	
	
    outFile.open(filename.c_str());
	outFile <<"OBJECTIVE VALUE: "<<design->getObjValue()<<endl;
	outFile<<"Total Duration: "<<duration<<endl;
	outFile<<"Peak2Peak on the surface of the volums: "<<Peak2Peak<<endl; 
	outFile<<"RMS on the surface of the volums: "<<RMS<<endl;
	outFile<<"Peak2Peak over the sampling inside and on the surface of the volums: "<<Peak2Peak_in<<endl; 
	outFile<<"RMS over the sampling inside and on the surface of the volums: "<<RMS_in<<endl;
     for(unsigned int n = 0; n<CoilsR.size(); n++)
		 for(int s = 0; s<Ns; s++){
			 if (SOL[n][s] * 1e7 > 1 || -SOL[n][s] * 1e7 > 1){
				 outFile <<n<<", "<<CoilsR[n]<<", "<<CoilsZ[n]<<", "<<setprecision(9)<<SOL[n][s] * 1e7 <<", "<<widths[s]<<", " << widths[s]/Gamma<<endl;
			 }
       	 }         
     outFile.close();


//         for(int n = 0; n<CoilsR.size(); n++)
// 		 for(int s = 0; s<Ns; s++){
// 			 if (SOL[n][s] * 1e7 > 1 || -SOL[n][s] * 1e7 > 1){
// 				 outFile <<n<<", "<<CoilsR[n]<<", "<<CoilsZ[n]<<", "<<setprecision(9)<<SOL[n][s] * 1e7 <<", "<<widths[s]<<", " << widths[s]/Gamma<<endl;
// 			 }
//        	 } 
//         outFile.close();
    
             
    }
	
	
    
    
    NbCoils--;

}	while ((design->getCplexStatus() != IloAlgorithm::Infeasible  &&
            design->getCplexStatus() != IloAlgorithm::InfeasibleOrUnbounded) && NbCoils > 6);
            
            

	

	 
	 cout<<"OBJECTIVE VALUE: "<<design->getObjValue()<<endl;
	 cout<<"*** Total Duration: "<<duration<<endl;
	 for(int n = 0; n < Nc; n++){
		 for(int s = 0; s < Ns; s++){
			 if (SOL[n][s]* 1e7> 1 || -SOL[n][s]* 1e7> 1){
				 cout << "N: "<<n<<", Zn: "<<CoilsZ[n]<<", Rn: "<<CoilsR[n]<<", Current: "<<SOL[n][s] * 1e7 <<", W: "<<widths[s]<<", H: " <<   widths[s]/Gamma<<endl;
			 }
		 }
	 }

    //cout << "Defining some more vectors" << endl;

		vector<vector<double> > BzMatrix;
		vector<vector<double> > BxMatrix;
		BzMatrix.resize(61); 
		BxMatrix.resize(61); 
		for(int i = 0; i< 61; i++){
			BzMatrix[i].resize(66); 
			BxMatrix[i].resize(66);
			for(int j = 0; j< 66; j++){
				BzMatrix[i][j] = 0; 
				BxMatrix[i][j] = 0;
			}
		}
		
    //cout << "Some other stuff" << endl;
	 
	  	double XMIN = 0.0;
		double ZMIN = -0.3;
		double XMAX = 0.65;
		double ZMAX = 0.3;
		double STEP = 0.01;
		
 		double B_MATRIX[61][66];
		for(int i = 0; i < 61; i++)
			for(int j = 0; j < 66; j++)
				B_MATRIX[i][j] = 0;
		
    //cout << "More other stuff with the matrix" << endl;
    //cout << "Nc = " << Nc << endl;
	 	 for(int n = 0; n < Nc; n++){
		 for(int s = 0; s < Ns; s++){
			 if (SOL[n][s]* 1e7> 1 || -SOL[n][s]* 1e7> 1){
				
				for(int i = 0; i< 61; i++)
					for(int j = 0; j< 66; j++){
						BzMatrix[i][j] = 0; 
						BxMatrix[i][j] = 0; 
						}
				
					 
				//double Bz = 0;
				//double Br = 0;
				double J = SOL[n][s] * 1e7/cross_sections[s];

				BFieldThick_Matrix *field_Matrix;
                //cout << "maybe the problem is in here???????" << endl;
				field_Matrix = new BFieldThick_Matrix(J, CoilsR[n], widths[s], widths[s]/Gamma, XMIN, ZMIN, XMAX, ZMAX, 0, CoilsZ[n], STEP, BzMatrix, BxMatrix);
                //cout << "Probably here" << endl;
					for(int i = 0; i < 61; i++)
						for(int j = 0; j < 66; j++)
							B_MATRIX[i][j] = B_MATRIX[i][j] + BzMatrix[i][j];
				delete field_Matrix;
			 }
			 }
                cout << "n = " << n << endl;
		 } 
	 
	 
    //cout << "Now the dense excel starts again" << endl;
	 ofstream ExcelFile2;
     ExcelFile2.open("FieldTotal_Matrix.csv");
		for(int i = 0; i < 61; i++){
			for(int j = 0; j < 66; j++)
				ExcelFile2 <<setprecision(9)<< B_MATRIX[i][j]<<", ";
				ExcelFile2 << endl;
		}
     ExcelFile2.close(); 
	 

	 
 		double B[Nt];
		for(int i = 0; i < Nt; i++)
				B[i] = 0;
	 
	 	 for(int n = 0; n < Nc; n++){
		 for(int s = 0; s < Ns; s++){
			 if (SOL[n][s] * 1e7 > 1 || -SOL[n][s] * 1e7> 1){
				 for(int m = 0; m < Nt; m++){
					 
				double Bz = AMNS[m][n][s] * SOL[n][s] ;
				// double Br = 0;
				
				// double J = SOL[n][s] * 10e7/cross_sections[s];
				// field = new BFieldThick(J, CoilsR[n], widths[s], widths[s]/Gamma, Xt[m], Zt[m], 0, CoilsZ[n], Bz, Br);
				B[m] = B[m] + Bz;
			 }
			 }
		 }
	 } 

	for(int m = 0; m < Nt; m++){
		cout<<m<<" "<<B[m]<<endl;
		if(B[m] > B_zero * (1 + tol) || B[m] < B_zero * (1 - tol))
			cout<<"******"<<endl;
	}


	BFieldThick *field2;
	double B_in[Nt_in];
	for(int i = 0; i < Nt_in; i++)
				B_in[i] = 0;
	for(int m = 0; m < Nt_in; m++){
		for(int n = 0; n < Nc; n++){
			for(int s = 0; s < Ns; s++){
				if (SOL[n][s]* 1e7> 1 || -SOL[n][s]* 1e7> 1){
				double Bz = 0;
				double Br = 0;
				double J = SOL[n][s] * 1e7/cross_sections[s];
				field2 = new BFieldThick(J, CoilsR[n], widths[s], widths[s]/Gamma, Xt_in[m], Zt_in[m], 0, CoilsZ[n], Bz, Br); 
				B_in[m] = B_in[m] + Bz;
				}
			}
		}
	}


	double min_B = 10;
	double max_B = -10;
	double sum_square = 0;
	for(int m = 0; m < Nt; m++){	
		if(B[m] < min_B) min_B = B[m];
		if(B[m] > max_B) max_B = B[m];	
		sum_square += B[m]*B[m];	
	}

	//Peak to Peak over the sample points inside the volume
	//and also the computation of RMS
	double min_B_in = 10;
	double max_B_in = -10;
	double sum_square_in = 0;
	for(int m = 0; m < Nt_in; m++){	
		if(B_in[m] < min_B_in) min_B_in = B_in[m];
		if(B_in[m] > max_B_in) max_B_in = B_in[m];	
		sum_square_in += B_in[m]*B_in[m];	
	}
	
	double Peak2Peak_in = 1e6*((max_B_in - 0.5) + (0.5 - min_B_in));
	double RMS_in = sqrt(sum_square_in/Nt_in);
	double Peak2Peak = 1e6*((max_B - 0.5) + (0.5 - min_B));
	double RMS = sqrt(sum_square/Nt);
		
		
 ofstream ExcelFile;
     ExcelFile.open("Coils.csv");
    
     for(unsigned int n = 0; n<CoilsR.size(); n++)
		 for(int s = 0; s<Ns; s++){
			 if (SOL[n][s] * 1e7 > 1 || -SOL[n][s] * 1e7 > 1){
				 ExcelFile <<n<<", "<<CoilsR[n]<<", "<<CoilsZ[n]<<", "<<setprecision(9)<<SOL[n][s] * 1e7 <<", "<<widths[s]<<", " << widths[s]/Gamma<<endl;
			 }
       	 }         
     ExcelFile.close();
    

	ofstream ExcelFile4;
    ExcelFile4.open("Results.csv");
	ExcelFile4 <<"OBJECTIVE VALUE: "<<design->getObjValue()<<endl;
	ExcelFile4<<"Total Duration: "<<duration<<endl;
	ExcelFile4<<"Peak2Peak on the surface of the volums: "<<Peak2Peak<<endl; 
	ExcelFile4<<"RMS on the surface of the volums: "<<RMS<<endl;
	ExcelFile4<<"Peak2Peak over the sampling inside and on the surface of the volums: "<<Peak2Peak_in<<endl; 
	ExcelFile4<<"RMS over the sampling inside and on the surface of the volums: "<<RMS_in<<endl;
     for(unsigned int n = 0; n<CoilsR.size(); n++)
		 for(int s = 0; s<Ns; s++){
			 if (SOL[n][s] * 1e7 > 1 || -SOL[n][s] * 1e7 > 1){
				 ExcelFile4 <<n<<", "<<CoilsR[n]<<", "<<CoilsZ[n]<<", "<<setprecision(9)<<SOL[n][s] * 1e7 <<", "<<widths[s]<<", " << widths[s]/Gamma<<endl;
			 }
       	 }         
     ExcelFile4.close();
	 
 ofstream ExcelFile5;
     ExcelFile5.open("SmallestCoils.csv");
    
     for(unsigned int i = 0; i<CPUs.size(); i++)
		ExcelFile5 <<NBCOILS[i]<<", ";  
		ExcelFile5 <<endl;
     for(unsigned int i = 0; i<CPUs.size(); i++)
		ExcelFile5 <<CPUs[i]<<", " ; 
		ExcelFile5 <<endl;
	for(unsigned int i = 0; i<CPUs.size(); i++)
		ExcelFile5 <<OBJECTIVEs[i]<<", " ; 
		ExcelFile5 <<endl;
	
		
     ExcelFile5.close();	
     
 
 
 	 vector<vector<vector<double> > >::iterator IterSOL;
 	  
    ofstream ExcelFile6;
     ExcelFile6.open("AllCoils.csv");
    	int cmp = 0;
     for( IterSOL = SOL_VECTOR.begin(); IterSOL != SOL_VECTOR.end(); IterSOL++, cmp++){
     	ExcelFile6 << "NBCOIL: "<<NBCOILS[cmp]<<endl;
     	for(unsigned int n = 0; n<CoilsR.size(); n++)
		 for(int s = 0; s<Ns; s++){
			 if ((*IterSOL)[n][s] * 1e7 > 1 || -(*IterSOL)[n][s] * 1e7 > 1){
				 ExcelFile6 <<n<<", "<<CoilsR[n]<<", "<<CoilsZ[n]<<", "<<setprecision(9)<<(*IterSOL)[n][s] * 1e7 <<", "<<widths[s]<<", " << widths[s]/Gamma<<endl;
			 }
       	 } 
       	 
       	ExcelFile6 << "****************************"<<endl;   
     }
	
     ExcelFile6.close();	 
    cout << "CODE IS DONE!!" << endl;
	return 0;
}


