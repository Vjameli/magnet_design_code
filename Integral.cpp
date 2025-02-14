//*****************************************************************************
//        C++      					                                          *
//        Iman Dayarian  	                                                  *
//        Toronto - CANADA            		                                  *
//                                                                            *
//        July 2015                   				                          *
//*****************************************************************************

#include "Integral.h"

Integral::Integral(double _t1, double _t2, double _t3, double _t4)
{
#define PI (3.141592654)
	
t1 = _t1;
t2 = _t2;
t3 = _t3;
t4 = _t4;
}


//Given Function of Integration
long double Integral::f1(long double x)
{
	long double d;
	double g1 = t1-t3*cos(x);
	double g2 = t3*t3*(sin(x)*sin(x)) + (t2-t4)*(t2-t4);
	double g3 = t3*cos(x);
	//double g4 = t3*sin(x);
	//double g5 = t2 - t4;
	double gf = sqrt(g1*g1 + g2);
	//double M1 = 2*g5*g3*log(g1+gf) + (1/4.0)*(3*(g3*g3-g4*g4) - (g1+g3)*(g1+g3))*log((gf-g5)/(gf+g5));
	//double M2 = -2*atan((g1*(g1+gf)+g4*g4)/(g4*g5)) + atan((g5+g1+gf)/g4) - atan((-g5+g1+gf)/g4);

	d = cos(x)*(gf + g3*log(g1+gf));
	

	return d;
}


long double Integral::f2(long double x)
{
	long double d;
	double g1 = t1-t3*cos(x);
	double g2 = t3*t3*(sin(x)*sin(x)) + (t2-t4)*(t2-t4);
	double g3 = t3*cos(x);
	double g4 = t3*sin(x);
	double g5 = t2 - t4;
	double gf = sqrt(g1*g1 + g2);
	double M1 = 2*g5*g3*log(g1+gf) + (0.25)*(3*(g3*g3-g4*g4) - (g1+g3)*(g1+g3))*log((gf-g5)/(gf+g5));
	double M2 = -2*atan((g1*(g1+gf)+g4*g4)/(g4*g5)) + atan((g5+g1+gf)/g4) - atan((-g5+g1+gf)/g4);
	
	d = cos(x)*(M1 + g3*g4*M2 + (0.5)*g5*gf);
	

	return d;
}


//For Legendre's Polynomial Pn(x)
long double Integral::pn(long double a[],int n,int m,long double x)
{
	int i;
	long double p=0;
	if(m==0){
		for(i=0;i<=n;i=i+2){
			
			if(x==0)
			break;
			p+=a[i]*pow(x,i);
		}
	}
	else
	{
		for(i=1;i<=n;i=i+2){ 
			p+=a[i]*pow(x,i);
		} 
	}
	return p;
}

//Derivative of Pn(x)
long double Integral::dn(long double a[],int n,int m,long double x)
{
	int i;
	long double p=0;
	if(m==0)
	{
		for(i=0;i<=n;i=i+2)
		{
			if(x==0)
			break;
			p+=i*a[i]*pow(x,i-1);
		}
	}
	else
	{
		for(i=1;i<=n;i=i+2)
		{
			p+=i*a[i]*pow(x,i-1);
		}
	}
	return p;
}

//Factorial Function
long double Integral::fact(int n)
{
	int i;
	long double f=1;
	for(i=2;i<=n;i++)
	{
		f*=i;
	}
	return f;
}

//Main Function
double Integral::N1(double c, double d)
{
	int n,m,i,N;
	
	n = 18;
	// double c,d;
	// cout<<"Enter the value of n for Pn(x) : \n";
	// cin>>n;
	// cout<<"Enter the lower limit a of integration : \n";
	// cin>>c;
	// cout<<"Enter the upper limit b of integration : \n";
	// cin>>d;

	if(n<=0)
	return 0;


	long double a[n],y[n],z[n],w[n],l,v,s,g=0,u[n];
	m = n%2;
	if(m == 0)
	{
		N=n/2;
	}
	else
	{
		N=(n-1)/2;
	}


	for(i=0;i<=N;i++)
	{
		a[n-2*i]=(pow(-1,i)*fact(2*n-2*i))/(pow(2,n)*fact(i)*fact(n-i)*fact(n-2*i));
	}

	// if(m==0)
	// {
		// cout<<"\nThe Legendre's Polynomial is : ";
		// cout<<a[0];
		// for(i=2;i<=n;i=i+2)
		// cout<<" + ("<<setprecision(10)<< a[i]<<") X^"<<i;
	// }
	// else
	// {
		// cout<<"\nThe Legendre's Polynomial is : ";
		// cout<<"("<<a[1]<<") X";
		// for(i=3;i<=n;i=i+2)
		// cout<<" + ("<<a[i]<<") X^"<<i;
	// }
	// cout<<endl;

	//Roots of Pn(x)
	for(i=0;i<n;i++)
	{
		z[i]=cos(3.14*(i+0.75)/(n+0.5));
		l=z[i];
		do
		{
			s=l-(pn(a,n,m,l)/dn(a,n,m,l));
			v=l;
			l=s;
		}
		while(fabs(l-v)>0.0000000000000001);
		y[i]=l;
		w[i]=2/((1-pow(l,2))*(dn(a,n,m,l)*dn(a,n,m,l)));
	}

	for(i=0;i<n;i++)
	{
		u[i]=((d-c)*y[i]/2)+(c+d)/2;
	}
	// cout<<"Roots\t\t\t\t"<<"Weights\n";
	// for(i=0;i<n;i++)
	// {
		// cout<<setprecision(15)<<y[i]<<"\t\t"<<setprecision(15)<<w[i]<<endl;
	// }
	for(i=0;i<n;i++)
	g+=w[i]*f1(u[i]);
	g=g*(d-c)/2;
//	cout<<"The Value of Integration is = "<<setprecision(10)<<g<<endl;
	return g;
}


double Integral::N2(double c, double d)
{
	int n,m,i,N;
	
	n = 10;
	// double c,d;
	// cout<<"Enter the value of n for Pn(x) : \n";
	// cin>>n;
	// cout<<"Enter the lower limit a of integration : \n";
	// cin>>c;
	// cout<<"Enter the upper limit b of integration : \n";
	// cin>>d;

	if(n<=0)
	return 0;


	long double a[n],y[n],z[n],w[n],l,v,s,g=0,u[n];
	m = n%2;
	if(m == 0)
	{
		N=n/2;
	}
	else
	{
		N=(n-1)/2;
	}


	for(i=0;i<=N;i++)
	{
		a[n-2*i]=(pow(-1,i)*fact(2*n-2*i))/(pow(2,n)*fact(i)*fact(n-i)*fact(n-2*i));
	}

	// if(m==0)
	// {
		// cout<<"\nThe Legendre's Polynomial is : ";
		// cout<<a[0];
		// for(i=2;i<=n;i=i+2)
		// cout<<" + ("<<setprecision(10)<< a[i]<<") X^"<<i;
	// }
	// else
	// {
		// cout<<"\nThe Legendre's Polynomial is : ";
		// cout<<"("<<a[1]<<") X";
		// for(i=3;i<=n;i=i+2)
		// cout<<" + ("<<a[i]<<") X^"<<i;
	// }
	// cout<<endl;

	//Roots of Pn(x)
	for(i=0;i<n;i++)
	{
		z[i]=cos(3.14*(i+0.75)/(n+0.5));
		l=z[i];
		do
		{
			s=l-(pn(a,n,m,l)/dn(a,n,m,l));
			v=l;
			l=s;
		}
		while(fabs(l-v)>0.0000000000000001);
		y[i]=l;
		w[i]=2/((1-pow(l,2))*(dn(a,n,m,l)*dn(a,n,m,l)));
	}

	for(i=0;i<n;i++)
	{
		u[i]=((d-c)*y[i]/2)+(c+d)/2;
	}
	// cout<<"Roots\t\t\t\t"<<"Weights\n";
	// for(i=0;i<n;i++)
	// {
		// cout<<setprecision(15)<<y[i]<<"\t\t"<<setprecision(15)<<w[i]<<endl;
	// }
	for(i=0;i<n;i++)
	g+=w[i]*f2(u[i]);
	g=g*(d-c)/2;
//	cout<<"The Value of Integration is = "<<setprecision(10)<<g<<endl;
	return g;
}



// void int_function_1_func(double x, double xminusa, double bminusx, double &y, void *ptr) 
// {

// double g1 = t1-t3*cos(x);
// double g2 = t3*t3*(sin(x)*sin(x)) + (t2-t4)*(t2-t4);
// double g3 = t3*cos(x);
// double g4 = t3*sin(x);
// double g5 = t2 - t4;
// double gf = sqrt(g1()*g1() + g2());
// double M1 = 2*g5()*g3()*log(g1()+gf()) + (1/4)*(3*(g3()*g3()-g4()*g4()) - (g1()+g3())*(g1()+g3()))*log((gf()-g5())/(gf()+g5()));
// double M2 = -2*atan((g1()*(g1()+gf())+g4()*g4())/(g4()*g5())) + atan((g5()+g1()+gf())/g4()) - atan((-g5()+g1()+gf())/g4());

// y = cos(x)*(gf() + g3()*log(g1()+gf()));

// double Integral::fn2(double x)
// {return cos(x)*(M1() + g3()*g4()*M2() + (1/2)*g5()*gf());}










 // // this callback calculates f(x)=fn1(x)
   // y = integ->fn1(x);
// }

// void int_function_2_func(double x, double xminusa, double bminusx, double &y, void *ptr) 
// {
    // // this callback calculates f(x)=fn2(x)
 // //   y = fn2(x);
// }

// double Integral::g1()
// {return t1-t3*cos(b);}

// double Integral::g2()
// {return t3*t3*(sin(b)*sin(b)) + (t2-t4)*(t2-t4);}

// double Integral::g3()
// {return t3*cos(b);}

// double Integral::g4()
// {return t3*sin(b);}

// double Integral::g5()
// {return t2 - t4;}

// double Integral::gf()
// {return sqrt(g1()*g1() + g2());}


// double Integral::M1()
// {return 2*g5()*g3()*log(g1()+gf()) + (1/4)*(3*(g3()*g3()-g4()*g4()) - (g1()+g3())*(g1()+g3()))*log((gf()-g5())/(gf()+g5()));}

// double Integral::M2()
// {return -2*atan((g1()*(g1()+gf())+g4()*g4())/(g4()*g5())) + atan((g5()+g1()+gf())/g4()) - atan((-g5()+g1()+gf())/g4());}

// double Integral::fn1(double x)
// {return cos(x)*(gf() + g3()*log(g1()+gf()));}

// double Integral::fn2(double x)
// {return cos(x)*(M1() + g3()*g4()*M2() + (1/2)*g5()*gf());}

// double Integral::N1()
// {
	// double a = 0;
    // double b = PI;
    // autogkstate s;
    // double v;
    // autogkreport rep;

    // autogksmooth(a, b, s);
    // alglib::autogkintegrate(s, int_function_1_func);
    // autogkresults(s, v, rep);

    // return v;

// }

// // double Integral::N2()
// // {
	// // double a = 0;
    // // double b = PI;
    // // autogkstate s;
    // // double v;
    // // autogkreport rep;

    // // autogksmooth(a, b, s);
    // // alglib::autogkintegrate(s, int_function_2_func);
    // // autogkresults(s, v, rep);

    // // return v;

// // }

