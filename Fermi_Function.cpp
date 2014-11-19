//This macro can be loaded into another macro
// When the function Fermi_Function is called, it will return the fermi function

//gROOT->Reset();

#include "fstream.h"
#include "iostream.h"
#include "string.h"
#include "TMath.h"
#include <complex>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif



/*void main()
{
		
}*/

const Double_t electron_mass_c2 = 510.998; //in keV
const Double_t fine_structure_const = 0.007297352533;


Double_t Fermi_Function(Double_t MassNumber, Double_t DaughterProtonNumber, Double_t Energy) //Energy in keV 
{
	
	
	Double_t W = Energy / electron_mass_c2 + 1.;
	Double_t p_e = sqrt(pow(W, 2) - 1.);
	Double_t beta = sqrt(1. - pow(1. / W, 2));
    Double_t nu = fine_structure_const * DaughterProtonNumber / beta;
    Double_t r_0 = (1.15 + 1.80*pow(MassNumber, -2./3.) - 1.20 * pow(MassNumber, -4./3.)) / 386.1592642;
    Double_t R = r_0 * pow(MassNumber, 1./3.);
    Double_t gamma1 = sqrt(1. - pow(fine_structure_const*DaughterProtonNumber, 2) );
	
	Double_t F_0 = FermiFunction_0(gamma1, p_e, R, nu);

       // cout << "F_O " << F_0 << "   Energy" << Energy <<  endl;
	return F_0;

}



Double_t FermiFunction_0(double g1, double pe, double R, double nu)
{
    complex<double> carg;
    
    carg = complex<double>(g1,nu);
    //cout << "complex_abs(cgamma) = " <<  complex_abs(cgamma(carg,0)) << endl;
    return 2.*(g1+1)*pow(2*pe*R, -2.*(1.-g1))*exp(M_PI*nu)*pow(complex_abs(cgamma(carg,0)),2);
    //return pow(complex_abs(cgamma(carg,0)),2);
}


Double_t complex_abs(complex<double> number)
{ 
  Double_t y,re,im;

  re = real(number);
  im = imag(number);

  y = sqrt(re*re + im*im);
  
  return y;
}


//  cgamma.cpp -- Complex gamma function.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns gamma function or log(gamma) for complex argument 'z'.
//
//  OPT value       function
//  ---------       --------
//      0           complex gamma
//      1           complex log(gamma)
//
//  Returns (1e308,0) if the real part of the argument is a negative integer
//  or 0 or exceeds 171.
//
complex<double> cgamma(complex<double> z,int OPT)
{
    complex<double> g,z0,z1;
    double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
    double na,t,x1,y1,sr,si;
    int i,j,k;
    
    x1 = 0;
    na = 0;

    static double a[] = 
	{
        8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590
       };

    x = real(z);
    y = imag(z);
    if (x > 171) return complex<double>(1e308,0);
    if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
        return complex<double>(1e308,1e-308);
    else if (x < 0.0) 
    {
        x1 = x;
        y1 = y;
        x = -x;
        y = -y;
    }
    x0 = x;
    if (x <= 7.0) 
    {
        na = (int)(7.0-x);
        x0 = x+na;
    }
    q1 = sqrt(x0*x0+y*y);
    th = atan(y/x0);
    gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(2.0*M_PI);
    gi = th*(x0-0.5)+y*log(q1)-y;
    for (k=0;k<10;k++)
    {
        t = pow(q1,-1.0-2.0*k);
        gr += (a[k]*t*cos((2.0*k+1.0)*th));
        gi -= (a[k]*t*sin((2.0*k+1.0)*th));
    }
    if (x <= 7.0) 
    {
        gr1 = 0.0;
        gi1 = 0.0;
        for (j=0;j<na;j++) 
	{
            gr1 += (0.5*log((x+j)*(x+j)+y*y));
            gi1 += atan(y/(x+j));
        }
        gr -= gr1;
        gi -= gi1;
    }
    if (x1 < 0.0) 
    {
        q1 = sqrt(x*x+y*y);
        th1 = atan(y/x);
        sr = -sin(M_PI*x)*cosh(M_PI*y);
        si = -cos(M_PI*x)*sinh(M_PI*y);
        q2 = sqrt(sr*sr+si*si);
	th2 = atan(si/sr);
        if (sr < 0.0) th2 += M_PI;
        gr = log(M_PI/(q1*q2))-gr;
        gi = -th1-th2-gi;
        x = x1;
        y = y1;
    }
    if (OPT == 0) 
    {
        g0 = exp(gr);
        gr = g0*cos(gi);
        gi = g0*sin(gi);
    }
    g = complex<double>(gr,gi);
    return g;
}