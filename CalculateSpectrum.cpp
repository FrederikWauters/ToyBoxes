//

gROOT->Reset();

#include "fstream.h"
#include "iostream.h"
#include "string.h"
#include "TMath.h"
#include "TGraph.h"
#include "TFile.h"
#include <TH1.h>

#include <math.h>
#include <stdlib.h>


//#include <complex>
//#include <iostream>
#include <cctype>
#include <cstdio>


  //initialize
 

 double MassNumber = 8;
  double DaughterProtonNumber = 4;
  double EndPointEnergy = 14500;
  int nsteps = 2000;

TH1D* hSpectrum = new TH1D("hSpectrum","Beta- Spectrum",nsteps,0.,EndPointEnergy);


int main()
{


  double fine_structure_const = 1/137.03599907;
  double electron_mass_c2 = 510.999;

  double EnergyOfMaximum;
  double gamma1 = sqrt(1. - pow(fine_structure_const* DaughterProtonNumber, 2) );
  double r_0 = (1.15 + 1.80*pow(MassNumber, -2./3.)- 1.20 * pow(MassNumber, -4./3.)) / 386.1592642;
  double R = r_0 * pow(MassNumber, 1./3.);


  
  double W_0 = EndPointEnergy / electron_mass_c2 + 1.;
 double a[7];

  float b[7][6];
  b[0][0] = 0.115;	b[0][1] = -1.8123;	b[0][2] = 8.2498;	b[0][3] = -11.223;	b[0][4] = -14.854;	b[0][5] = 32.086; 
  b[1][0] = -0.00062;	b[1][1] = 0.007165;	b[1][2] = 0.01841;	b[1][3] = -0.53736;	b[1][4] = 1.2691;	b[1][5] = -1.5467;
  b[2][0] = 0.02482;	b[2][1] = -0.5975;	b[2][2] = 4.84199; 	b[2][3] = -15.3374;	b[2][4] = 23.9774;	b[2][5] = -12.6534;
  b[3][0] = -0.14038;	b[3][1] = 3.64953;	b[3][2] = -38.8143; 	b[3][3] = 172.1368;	b[3][4] = -346.708;	b[3][5] = 288.7873;
  b[4][0] = 0.008152;	b[4][1] = -1.15664;	b[4][2] = 49.9663; 	b[4][3] = -273.711;	b[4][4] = 657.6292;	b[4][5] = -603.7033;
  b[5][0] = 1.2145;	b[5][1] = -23.9931;	b[5][2] =  149.9718;	b[5][3] = -471.2985;	b[5][4] = 662.1909;	b[5][5] = -305.6804;
  b[6][0] = -1.5632;	b[6][1] = 33.4192;	b[6][2] = -255.1333;	b[6][3] = 938.5297;	b[6][4] = -1641.2845;	b[6][5] = 1095.358;
    
 
  for(int i=0; i<7; i++)
  {
    a[i] = 0;
      for(int j=0; j<6; j++){  a[i] += b[i][j]*pow(fine_structure_const*DaughterProtonNumber,j);} 
  }


 TFile f1("Spectrum.dat","RECREATE");

 double Energy, Result;

 gROOT->ProcessLine(".L Fermi_Function.cpp");

  Double_t End = EndPointEnergy;
  

 for(int i = 1; i<nsteps; i++)
 {
   Energy = i*(EndPointEnergy-0.01)/(nsteps*1.);
   
    double Result;
    double W = Energy / electron_mass_c2 + 1.;
    double p_e = sqrt(pow(W, 2) - 1.);
    double beta = sqrt(1. - pow(1. / W, 2));
    double nu = fine_structure_const * DaughterProtonNumber / beta;
    
    double sum = 0;
    
    for(int j=1; j < 7; j++)
	sum += a[j]*pow(W*R, j);
    

    double F_0 = Get_Fermi_Function(MassNumber, DaughterProtonNumber, Energy);
    
    double L_0 = 1 + 13./60.*pow(fine_structure_const*DaughterProtonNumber,2) - W*R*fine_structure_const*DaughterProtonNumber*(41.-26.*gamma1)/15./(2.*gamma1-1) - 
		    fine_structure_const*DaughterProtonNumber*R*gamma1*(17.-2.*gamma1)/30./W/(2.*gamma1-1)+a[0]*R/W + 
		    0.41*(R-0.0164)*pow(fine_structure_const*DaughterProtonNumber,4.5)+ sum;
		    
    double C, C0, C1, C2;
    
    C0 = -233.*pow(fine_structure_const*DaughterProtonNumber,2)/630. - pow(W_0*R,2)/5.+2*W_0*R*fine_structure_const*DaughterProtonNumber/35;
    C1 = -21.*R*fine_structure_const*DaughterProtonNumber/35.+4./9.*W_0*R*R;
    C2 = -4./9.*R*R;
    
    C = 1 + C0 + C1*W + C2*W*W;
   
   Result =  C * L_0 * F_0 * p_e * W * pow(W_0 - W, 2)* ( pow(W_0 - W, 2) + (W*W-1) );
   
   int bin = hSpectrum->FindBin(Energy);
   hSpectrum->SetBinContent(bin,Result);
  //cout << Energy << " " << Result <<  " " << bin << endl;

 }

 

 hSpectrum->Draw();
 //cout << hSpectrum->Integral() << endl;
 f1.Close();
 return 1;

}


double InvertedSpectrum(double Energy)
{
    return -Spectrum((double) Energy);
}

/*double Spectrum(double Energy, int n)
{
    

    double Result;
    double W = Energy / electron_mass_c2 + 1.;
    double p_e = sqrt(pow(W, 2) - 1.);
    double beta = sqrt(1. - pow(1. / W, 2));
    double nu = fine_structure_const * DaughterProtonNumber / beta;
    
    double sum = 0;
    
    for(int j=1; j < 7; j++)
	sum += a[j]*pow(W*R, j);
    

    double F_0 = Get_Fermi_Function(MassNumber, DaughterProtonNumber, Energy);
    
    double L_0 = 1 + 13./60.*pow(fine_structure_const*DaughterProtonNumber,2) - W*R*fine_structure_const*DaughterProtonNumber*(41.-26.*gamma1)/15./(2.*gamma1-1) - 
		    fine_structure_const*DaughterProtonNumber*R*gamma1*(17.-2.*gamma1)/30./W/(2.*gamma1-1)+a[0]*R/W + 
		    0.41*(R-0.0164)*pow(fine_structure_const*DaughterProtonNumber,4.5)+ sum;
		    
    double C, C0, C1, C2;
    
    C0 = -233.*pow(fine_structure_const*DaughterProtonNumber,2)/630. - pow(W_0*R,2)/5.+2*W_0*R*fine_structure_const*DaughterProtonNumber/35;
    C1 = -21.*R*fine_structure_const*DaughterProtonNumber/35.+4./9.*W_0*R*R;
    C2 = -4./9.*R*R;
    
    C = 1 + C0 + C1*W + C2*W*W;
    retrun 1.;
}*/

double Get_Fermi_Function(double M, double D, double E)
{
  

  return Fermi_Function(M, D, E);
}

double Average_mE(double Emin)
{
  //cout << "Give min energy" << endl;

  //double Emin;
  //cin >> Emin;

  if(Emin > EndPointEnergy){cout << "Energy to large !!!" << endl; return;}

  int startBin = hSpectrum->FindBin(Emin);
  int stopBin = hSpectrum->FindBin(EndPointEnergy);

  double weightedSum = 0;
  double normSum = 0;

  for(int i = startBin;i <= stopBin;i++)
  {
    weightedSum += hSpectrum->GetBinContent(i)*(1/(hSpectrum->GetBinCenter(i)+509.99));
    normSum += hSpectrum->GetBinContent(i);
  }

  cout << "Average mE = " << (weightedSum / normSum)*509.999 << endl;
  
  return (weightedSum / normSum)*509.999;

}

double Average_mE(double Emin,double Emax)
{
  //cout << "Give min energy" << endl;

  //double Emin;
  //cin >> Emin;

  if(Emin > EndPointEnergy && Emax > EndPointEnergy){cout << "Energy to large !!!" << endl; return;}

  int startBin = hSpectrum->FindBin(Emin);
  int stopBin = hSpectrum->FindBin(Emax);

  double weightedSum = 0;
  double normSum = 0;

  for(int i = startBin;i <= stopBin;i++)
  {
    weightedSum += hSpectrum->GetBinContent(i)*(1/(hSpectrum->GetBinCenter(i)+509.99));
    normSum += hSpectrum->GetBinContent(i);
  }

  cout << "Average mE = " << (weightedSum / normSum)*509.999 << endl;
  
  return (weightedSum / normSum)*509.999;

}

void Plot_mE_trend()
{
  const Int_t nSteps = 50;
  Double_t energy[nSteps];
  Double_t mE[nSteps];

  for(Int_t i = 1; i <= nSteps; i++)
  {
    energy[i-1] = (EndPointEnergy / nSteps*1.0) * i * 0.99 ;
   cout << "i = " << i << "  energy : " << energy[i] << endl;
    mE[i-1] = Average_mE(energy[i-1]);
   
  }

  TCanvas *c2 = new TCanvas("c2","mE Values for 6He",800,600);
  TGraph* gmE = new TGraph(nSteps,energy,mE);
  gmE->Draw("alp"); 

}


