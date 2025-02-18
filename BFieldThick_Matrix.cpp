//*****************************************************************************
//        C++      					                                          *
//        Iman Dayarian  	                                                  *
//        Toronto - CANADA            		                                  *
//                                                                            *
//        August 2015                   				                          *
//*****************************************************************************

#include "BFieldThick_Matrix.h"

BFieldThick_Matrix::BFieldThick_Matrix(double J0, double r0, double dr, double dz, double xmin, double zmin, double xmax, double zmax, double xshift, double zshift, double s, vector<vector<double> > &BzMatrix, vector<vector<double> > &BxMatrix)
										//(J0,r0,dr,dz,xmin,zmin,xmax,zmax,xshift,zshift,s)
{		
#define PI (3.14159265358979311599796346854418516)

// *************************
// Coil located at CAX, (0,0,0)
// J0 = current density - the same for all coils
// r0 = coil radius at mid-position
// dr = coil width
// dz = coil height
// xmin, xmax, zmin, zmax - the domain where the magnetic field is computed
// xshift, zshift - traslation shifts of the coil wrt z-axis 
// s = sampling of the points where magnetic field is computed
// *************************






		double u0 = 4*PI*1e-7;    // Permeability of free space is a global variable !!!!

		double C0 = u0*J0/(2*PI);


double z1 = zmin-zshift;
cout << BzMatrix.size() << endl;
cout << BzMatrix[0].size() << endl;
for (unsigned int in1 = 0; in1 <  BzMatrix.size(); in1++){   // defines the relative distance b/w coil and target point
    double x = xmin-xshift;
	for(unsigned int in2 = 0; in2 < BzMatrix[0].size(); in2++){	// defines the relative lateral shift from z-axis b/w coil and target point
		double theta = 0;
		double z = 0;
		double r = 0;
        Polar(x, 0, z1, r, theta, z);
		
		       if (r >= 0.000001){ 
			Integral* int1;
			Integral* int2;
			Integral* int3;
			Integral* int4;
			
			int1 = new Integral(r0+dr/2,-dz/2,r,z);
			int2 = new Integral(r0-dr/2,-dz/2,r,z);
			int3 = new Integral(r0+dr/2,dz/2,r,z);
			int4 = new Integral(r0-dr/2,dz/2,r,z);
			
            //cout << "Bfieldthick matrix dAdZ" << endl;

			double dAdz = int1->N1(0.0, PI) - int2->N1(0.0, PI) - int3->N1(0.0, PI) + int4->N1(0.0, PI);

            //cout << "Bfieldthick matrix dAdZ after" << endl;
            //dAdz = N1Integral(r0+dr/2,-dz/2,r,z) - N1Integral(r0-dr/2,-dz/2,r,z) - N1Integral(r0+dr/2,dz/2,r,z) + N1Integral(r0-dr/2,dz/2,r,z);
			//double a1 = int3->N2(0.0, PI);
			//double a2 = int4->N2(0.0, PI);
			//double a3 = int1->N2(0.0, PI);
			//double a4 = int2->N2(0.0, PI);
			
			
			double dAdr = int3->N2(0.0, PI) - int4->N2(0.0, PI) - int1->N2(0.0, PI) + int2->N2(0.0, PI);
			//dAdr = N2Integral(r0+dr/2,dz/2,r,z) - N2Integral(r0-dr/2,dz/2,r,z) - N2Integral(r0+dr/2,-dz/2,r,z) + N2Integral(r0-dr/2,-dz/2,r,z);
            
            BxMatrix[in1][in2] = -C0*(x/r)*dAdz;
//             By(in1,in2) = -C0*(y/r)*dAdz;
            BzMatrix[in1][in2] = (C0/r)*dAdr;
		}
        else
			if (r < 0.000001){
            BxMatrix[in1][in2] = 0;
//             By(in1,in2) = 0;
            
            double term1 = (r0 + dr/2 + sqrt((r0+dr/2)* (r0+dr/2) +(dz/2-z) * (dz/2-z)))/(r0 - dr/2 + sqrt((r0-dr/2) * (r0-dr/2) + (dz/2-z) * (dz/2-z)));
            double term2 = (r0 + dr/2 + sqrt((r0+dr/2) * (r0+dr/2) + (dz/2+z) * (dz/2+z)))/(r0 - dr/2 + sqrt((r0-dr/2) * (r0-dr/2) + (dz/2+z) * (dz/2+z)));
            BzMatrix[in1][in2] = C0*PI*((dz/2-z)*log(term1) + (dz/2+z)*log(term2));
        }          
        x = x + s;
    }
    z1 = z1 + s;
}


}

void BFieldThick_Matrix::Polar(double x, double y, double z, double& r, double& theta, double& z1)
{
    //cout << "maybe polar?" << endl;
	r = sqrt((pow(x,2))+(pow(y,2)));

    theta = atan2(y,x);
	
	theta = (theta*180)/3.14159265358979311599796346854418516;
	
	z1 = z;
}
