//*****************************************************************************
//        C++      					                                          *
//        Iman Dayarian  	                                                  *
//        Toronto - CANADA            		                                  *
//                                                                            *
//        July 2015                   				                          *
//*****************************************************************************

#include "BFieldThick.h"

BFieldThick::BFieldThick(double J0, double r0, double dr, double dz, double Xt, double Zt, double xshift, double zshift,  double &Bz, double &Br)
{
#define PI (3.14159265358979311599796346854418516)

        //cout << "BFieldThick constructor parameters:" << endl;
        //cout << "J0=" << J0 << " r0=" << r0 << " dr=" << dr << " dz=" << dz << endl;
        //cout << "Xt=" << Xt << " Zt=" << Zt << " xshift=" << xshift << " zshift=" << zshift << endl;

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

		double z1 = Zt - zshift;
		double x = Xt - xshift;

		double theta = 0;
		double z = 0;
		double r = 0;
		Polar(x, 0, z1, r, theta, z);

        if (r >= 0.000001){ 
            //cout << "I am bigger than r: defining objects" << endl;
			Integral* int1;
			Integral* int2;
			Integral* int3;
			Integral* int4;
			
//            cout << "Defining more objects" << endl;
			int1 = new Integral(r0+dr/2,-dz/2,r,z);
			int2 = new Integral(r0-dr/2,-dz/2,r,z);
			int3 = new Integral(r0+dr/2,dz/2,r,z);
			int4 = new Integral(r0-dr/2,dz/2,r,z);
			
  //          cout << "First double" << endl;
			double dAdz = int1->N1(0.0, PI) - int2->N1(0.0, PI) - int3->N1(0.0, PI) + int4->N1(0.0, PI);
            //dAdz = N1Integral(r0+dr/2,-dz/2,r,z) - N1Integral(r0-dr/2,-dz/2,r,z) - N1Integral(r0+dr/2,dz/2,r,z) + N1Integral(r0-dr/2,dz/2,r,z);
			//double a1 = int3->N2(0.0, PI);
			//double a2 = int4->N2(0.0, PI);
			//double a3 = int1->N2(0.0, PI);
			//double a4 = int2->N2(0.0, PI);
			
			
            //cout << "Second double" << endl;
			double dAdr = int3->N2(0.0, PI) - int4->N2(0.0, PI) - int1->N2(0.0, PI) + int2->N2(0.0, PI);
			//dAdr = N2Integral(r0+dr/2,dz/2,r,z) - N2Integral(r0-dr/2,dz/2,r,z) - N2Integral(r0+dr/2,-dz/2,r,z) + N2Integral(r0-dr/2,-dz/2,r,z);
            
            //cout << "Br" << endl;
            Br = -C0*(x/r)*dAdz;
//             By(in1,in2) = -C0*(y/r)*dAdz;

            //cout << "Bz" << endl;
            Bz = (C0/r)*dAdr;
		}
        else
			if (r < 0.000001){
            //cout << "I am smaller than r" << endl;
            Br = 0;
//             By(in1,in2) = 0;
            
            double term1 = (r0 + dr/2 + sqrt((r0+dr/2)* (r0+dr/2) +(dz/2-z) * (dz/2-z)))/(r0 - dr/2 + sqrt((r0-dr/2) * (r0-dr/2) + (dz/2-z) * (dz/2-z)));
            double term2 = (r0 + dr/2 + sqrt((r0+dr/2) * (r0+dr/2) + (dz/2+z) * (dz/2+z)))/(r0 - dr/2 + sqrt((r0-dr/2) * (r0-dr/2) + (dz/2+z) * (dz/2+z)));
            Bz = C0*PI*((dz/2-z)*log(term1) + (dz/2+z)*log(term2));
        }  
 














// global u0    // Permeability of free space is a global variable !!!!

// double C0 = u0*J0/(2*pi);

	// sz1 = (int)(((zmax-zshift) - (zmin-zshift))/s);
	// sz2 = (int)(((xmax-xshift) - (xmin-xshift))/s);
	
	// vector<vector<double> > Bx;
	// Bx.resize(sz1);
	// for(int j = 0; j < sz1; j++){
		// Bx[j].resize(sz2);
		// for(int i = 0; i < sz1; i++)
			// Bx[j][i] = 0;	
	// }

	// int in1 = 1;
	// for (int z1 = zmin-zshift; z1 <= zmax-zshift; z1 = z1 + s)    // defines the relative distance b/w coil and target point
	// {
		// int in2 = 1;
		// for (int x = xmin-xshift; x <= xmax-xshift; x = x + s)			// defines the relative lateral shift from z-axis b/w coil and target point
		// {
		// double theta = 0;
		// double z = 0;
		// double r = 0;
		// Polar(x, 0, z1, r, theta, z);

        // if (r ~= 0){ 
			// Integral* int1, int2, int3, int4;
			
			// int1 = new Integral(r0+dr/2,-dz/2,r,z);
			// int2 = new Integral(r0-dr/2,-dz/2,r,z);
			// int3 = new Integral(r0+dr/2,dz/2,r,z);
			// int4 = new Integral(r0-dr/2,dz/2,r,z);
			
			// double dAdz = int1->N1(0, PI) - int2->N1(0, PI) - int3->N1(0, PI) + int4->N1(0, PI);
            // //dAdz = N1Integral(r0+dr/2,-dz/2,r,z) - N1Integral(r0-dr/2,-dz/2,r,z) - N1Integral(r0+dr/2,dz/2,r,z) + N1Integral(r0-dr/2,dz/2,r,z);
			// double dAdr = int3->N2(0, PI) - int4->N2(0, PI) - int1->N2(0, PI) + int2->N2(0, PI);
			// //dAdr = N2Integral(r0+dr/2,dz/2,r,z) - N2Integral(r0-dr/2,dz/2,r,z) - N2Integral(r0+dr/2,-dz/2,r,z) + N2Integral(r0-dr/2,-dz/2,r,z);
            
            // Bx[in1][in2] = -C0*(x/r)*dAdz;
// //             By(in1,in2) = -C0*(y/r)*dAdz;
            // Bz[in1][in2] = (C0/r)*dAdr;
		// }
        // else
			// if (r == 0){
            // Bx(in1,in2) = 0;
// //             By(in1,in2) = 0;
            
            // double term1 = (r0 + dr/2 + sqrt((r0+dr/2)* (r0+dr/2) +(dz/2-z) * (dz/2-z)))/(r0 - dr/2 + sqrt((r0-dr/2) * (r0-dr/2) + (dz/2-z) * (dz/2-z)));
            // double term2 = (r0 + dr/2 + sqrt((r0+dr/2) * (r0+dr/2) + (dz/2+z) * (dz/2+z)))/(r0 - dr/2 + sqrt((r0-dr/2) * (r0-dr/2) + (dz/2+z) * (dz/2+z)));
            // Bz[in1][in2] = C0*PI*((dz/2-z)*log(term1) + (dz/2+z)*log(term2));
        // }  
        // in2 = in2 + 1;
		// }
    // in1 = in1 + 1;
	// }

}

void BFieldThick::Polar(double x, double y, double z, double& r, double& theta, double& z1)
{
	r = sqrt((pow(x,2))+(pow(y,2)));

    theta = atan2(y,x);
	
	theta = (theta*180)/3.14159265358979311599796346854418516;
	
	z1 = z;
}
