/*
  FW - Sep 2014
  MC exercise to generate LT spectra with a e- like oscillating background
  Try out different algorithms to get rid of the oscillation and check how the fit results and chi sqr is affected

  Justins Island Filter
  https://muon.npl.washington.edu/elog/musun/analysis-run4/292
  https://muon.npl.washington.edu/elog/musun/analysis-run4/247

*/


#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TROOT.h"
#include <TStyle.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TRandom3.h>


//gRoot->Reset();

//electron spectrum
TF1* fElSpectrum;
TF1* fTest;
TH1D* hInitElSpectrumFine;
TH1D* hInitElSpectrumWide;
TF1* fitFunction; 
TF1* fitFunctionFine;
TH1D* hInitElSpectrumFine_Res;
TH1D* hInitElSpectrumWide_Res;
TH1D* hInitElSpectrumFine_Sigma;
TH1D* hInitElSpectrumWide_Sigma;

//After Gaussian filtering
TH1D* hUniformElSpectrumFine;
TH1D* hUniformElSpectrumWide;
TH1D* hUniformElSpectrumFine_Res;
TH1D* hUniformElSpectrumWide_Res;
TH1D* hUniformElSpectrumFine_Sigma;
TH1D* hUniformElSpectrumWide_Sigma;
TF1* fitFunctionUniform; 
TF1* fitFunctionFineUniform;

//rates in 1
const double decayRate = 455000.;
const double eRate = 0.003; //relative to mu rate
const double eModRate = 50500000.;
const double constBGRate = 0.001;

const double pi = 3.14159265359;


//range
double tStart = 0.;
double tStop = 0.000025;
double tFitStart = 0.00000016;
double tFitStop = 0.00002;

//bin width in ns
int narrowBinWidth = 2.;
int wideBinWidth = 40.;

//MC choice
const int MC = 1;

//Correction parameters
double uniformJitter = 20; //uniform jitter in ns

double value(double t)
{
  return exp(455000*t) + 1000.*t + 5000*sin(50500000.*2*pi*t);
  
}

void BuildFunction(int N)
{
  //fElSpectrum = new TF1("fElSpectrum","exp(-1*[0]*x)+[1]+[2]*sin([3]*x)",tStart,tStop);
  fElSpectrum = new TF1("fElSpectrum","[4]*(exp(-1*[0]*x)+[1]+[2]*sin([3]*x)+[2])",tStart,tStop);
  fElSpectrum->SetParameter(0,decayRate);
  fElSpectrum->SetParameter(1,constBGRate);
  fElSpectrum->SetParameter(2,eRate);
  fElSpectrum->SetParameter(3,2*pi*eModRate);
  fElSpectrum->SetParameter(4,1);
  fElSpectrum->SetNpx(100000);
  cout << "fElSpectrum max : " << fElSpectrum->GetMaximum() << " fElSpectrum min : " << fElSpectrum->GetMinimum() << endl;
  
  fTest = new TF1("fTest","exp(-1*[0]*x)",tStart,tStop);
  fTest->SetParameter(0,decayRate);

  fitFunction = new TF1("fitFunction","[2]*exp([0]*x*-1)+[1]",tFitStart,tFitStop);
  fitFunction->SetParName(0,"rate");
  fitFunction->SetParName(1,"BG");
  fitFunction->SetParName(2,"N");
  fitFunction->SetParLimits(0,350000.,550000.);
  //fitFunction->SetParLimits(1,1.,N);
  fitFunction->SetParLimits(2,1.,N);


  fitFunctionFine = new TF1("fitFunctionFine","[2]*exp([0]*x*-1)+[1]",tFitStart,tFitStop);
  fitFunctionFine->SetParName(0,"rate");
  fitFunctionFine->SetParName(1,"BG");
  fitFunctionFine->SetParName(2,"N");
  fitFunctionFine->SetParLimits(0,350000.,550000.);
  //fitFunction->SetParLimits(1,1.,N);
  fitFunctionFine->SetParLimits(2,1.,N);

  fitFunctionUniform = (TF1*)fitFunction->Clone(); 
  fitFunctionFineUniform = (TF1*)fitFunctionFine->Clone();

}

void InitHistograms()
{
  int nBins = (int)1e09*(tStop - tStart)/(narrowBinWidth); cout << "nBins " << nBins << endl;
  hInitElSpectrumFine = new TH1D("hInitElSpectrumFine","Uncorrected Electron spectrum; decay time (s)",nBins,tStart,tStop);  
  hInitElSpectrumFine_Res = new TH1D("hInitElSpectrumFine_Res","Uncorrected Electron spectrum Residuals; decay time (s); residuals (sigma)",nBins,tStart,tStop);
  hUniformElSpectrumFine = new TH1D("hUniformElSpectrumFine","Corrected Electron spectrum; decay time (s)",nBins,tStart,tStop);  
  hUniformElSpectrumFine_Res = new TH1D("hUniformElSpectrumFine_Res","Corrected Electron spectrum Residuals; decay time (s); residuals (sigma)",nBins,tStart,tStop);  

  nBins = (int)1e09*(tStop - tStart)/(wideBinWidth);
  hInitElSpectrumWide = new TH1D("hInitElSpectrumWide","Uncorrected Electron spectrum; decay time (s)",nBins,tStart,tStop);
  hInitElSpectrumWide_Res = new TH1D("hInitElSpectrumWide_Res","Uncorrected Electron spectrum Residuals; decay time (s); residuals (sigma)",nBins,tStart,tStop);
  hUniformElSpectrumWide = new TH1D("hUniformElSpectrumWide","Corrected Electron spectrum; decay time (s)",nBins,tStart,tStop);
  hUniformElSpectrumWide_Res = new TH1D("hUniformElSpectrumWide_Res","Corrected Electron spectrum Residuals; decay time (s); residuals (sigma)",nBins,tStart,tStop);

  hInitElSpectrumWide_Sigma = new TH1D("hInitElSpectrumWide_Sigma","Residual distribution; residual (sigma)",200,-10,10); 
  hInitElSpectrumFine_Sigma = new TH1D("hInitElSpectrumFine_Sigma","Residual distribution; residual (sigma)",200,-10,10); 

  hUniformElSpectrumWide_Sigma = new TH1D("hUniformElSpectrumWide_Sigma","Residual distribution; residual (sigma)",200,-10,10); 
  hUniformElSpectrumFine_Sigma = new TH1D("hUniformElSpectrumFine_Sigma","Residual distribution; residual (sigma)",200,-10,10); 

  cout << "Histograms initialzed ... " << endl;
}

void DrawStuff(int N)
{
  
  TCanvas* c1 = new TCanvas("c1","Before corrections",1400,800);
  c1->Divide(3,2);
  c1->cd(1);
  gPad->SetLogy();
  fElSpectrum->Draw();
  c1->cd(2);
  gPad->SetLogy();
  hInitElSpectrumFine->Draw();
  fitFunctionFine->Draw("SAME");
  c1->cd(3);
  gPad->SetLogy();
  hInitElSpectrumWide->Draw();
  fitFunction->Draw("SAME");
  c1->cd(4);
  hInitElSpectrumFine_Res->Draw();
  c1->cd(5);
  hInitElSpectrumWide_Res->Draw();
  c1->cd(6);
  gPad->SetLogy();
  hInitElSpectrumFine_Sigma->Draw();

  TCanvas* c2 = new TCanvas("c2","Uniform correction",1400,800);
  c2->Divide(3,2);
  c2->cd(1);
  gPad->SetLogy();
  hUniformElSpectrumFine->Draw();
  fitFunctionFineUniform->Draw("SAME");
  c2->cd(2);
  gPad->SetLogy();
  hUniformElSpectrumWide->Draw();
  fitFunctionUniform->Draw("SAME");
  c2->cd(3);
  hUniformElSpectrumFine_Res->Draw();
  c2->cd(4);
  hUniformElSpectrumWide_Res->Draw();
  c2->cd(5);
  gPad->SetLogy();
  hUniformElSpectrumFine_Sigma->Draw();
  c2->cd(6);
  gPad->SetLogy();
  hUniformElSpectrumWide_Sigma->Draw();
}


int FillHistograms(int N)
{
  cout << "Start filling histograms ... " << endl;
  //TRandom3 should be the best one according to http://root.cern.ch/root/html/TRandom.html
  //initilazing with "0" should give a different seed everytime I run it http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=13253
  TRandom3 random(0); 
  for(int i = 0; i<N ; i++)
  {
    double value = 0;
      
    if(MC==1) value = fElSpectrum->GetRandom(tStart,tStop);
       
    //alternative method, faster?
    if(MC==2)
    {
      bool accept = false;
      double max = fElSpectrum->GetMaximum();
      while(!accept)
      {
        double x = random.Uniform(tStart,tStop);
        double functionValue = fElSpectrum->Eval(x);
        double y = random.Uniform(0,max);
        //cout << " x " << x << " function value " << functionValue << " y " << y << endl;
        if(y < functionValue)
        {
          accept = true;
          value = x;

          //cout << " accept ! "<< endl;
        }      
      }
    }
    hInitElSpectrumFine->Fill(value);
    hInitElSpectrumWide->Fill(value);
   
    double jitter = 1e-09*gRandom->Uniform(uniformJitter);
    

    hUniformElSpectrumFine->Fill(value+jitter);
    hUniformElSpectrumWide->Fill(value+jitter);

  }
  
  return N;
}

int FitHistograms(int N)
{
  cout << "**********************************" << endl;
  cout << "********* fit results ************" << endl;
  cout << "**********************************" << endl << endl;

  cout << "hInitElSpectrumWide fit results : " << endl;

  fitFunction->SetParLimits(2,0.,2.);
  hInitElSpectrumWide->Fit(fitFunction,"BQR");
  fitFunction->SetParLimits(2,1.,N);
  hInitElSpectrumWide->Fit(fitFunction,"BRQ");
  int ndf = fitFunction->GetNDF();
  cout << "Rate: " << fitFunction->GetParameter(0) << "(" << fitFunction->GetParError(0) << ")" << "  BG: " << fitFunction->GetParameter(1)  <<  "(" << fitFunction->GetParError(1) << ")" << "  N: " << fitFunction->GetParameter(2)  <<  "(" << fitFunction->GetParError(2) << ")" << " ChiSqr: " << fitFunction->GetChisquare()/ndf << endl;
  

  cout << "hInitElSpectrumFine fit results : " << endl;

  fitFunctionFine->SetParLimits(2,0.,1.);
  hInitElSpectrumFine->Fit(fitFunctionFine,"BRQ");
  fitFunctionFine->SetParLimits(2,1.,N);
  hInitElSpectrumFine->Fit(fitFunctionFine,"BRQ");
  ndf = fitFunctionFine->GetNDF();
  cout << "Rate: " << fitFunctionFine->GetParameter(0) << "(" << fitFunctionFine->GetParError(0) << ")" << "  BG: " << fitFunctionFine->GetParameter(1)  <<  "(" << fitFunctionFine->GetParError(1) << ")" << "  N: " << fitFunctionFine->GetParameter(2)  <<  "(" << fitFunctionFine->GetParError(2) << ")" << " ChiSqr: " << fitFunctionFine->GetChisquare()/ndf << endl;

  cout << "hUniformElSpectrumWide fit results : " << endl;

  fitFunctionUniform->SetParLimits(2,0.,2.);
  hUniformElSpectrumWide->Fit(fitFunctionUniform,"BQR");
  fitFunctionUniform->SetParLimits(2,1.,N);
  hUniformElSpectrumWide->Fit(fitFunctionUniform,"BRQ");
  ndf = fitFunction->GetNDF();
  cout << "Rate: " << fitFunctionUniform->GetParameter(0) << "(" << fitFunctionUniform->GetParError(0) << ")" << "  BG: " << fitFunctionUniform->GetParameter(1)  <<  "(" << fitFunctionUniform->GetParError(1) << ")" << "  N: " << fitFunctionUniform->GetParameter(2)  <<  "(" << fitFunctionUniform->GetParError(2) << ")" << " ChiSqr: " << fitFunctionUniform->GetChisquare()/ndf << endl;

  cout << "hUniformElSpectrumFine fit results : " << endl;

  fitFunctionFineUniform->SetParLimits(2,0.,1.);
  hUniformElSpectrumFine->Fit(fitFunctionFineUniform,"BRQ");
  fitFunctionFineUniform->SetParLimits(2,1.,N);
  hUniformElSpectrumFine->Fit(fitFunctionFineUniform,"BRQ");
  ndf = fitFunctionFineUniform->GetNDF();
  cout << "Rate: " << fitFunctionFineUniform->GetParameter(0) << "(" << fitFunctionFineUniform->GetParError(0) << ")" << "  BG: " << fitFunctionFineUniform->GetParameter(1)  <<  "(" << fitFunctionFineUniform->GetParError(1) << ")" << "  N: " << fitFunctionFineUniform->GetParameter(2)  <<  "(" << fitFunctionFineUniform->GetParError(2) << ")" << " ChiSqr: " << fitFunctionFineUniform->GetChisquare()/ndf << endl;

  
  for(int i = 1; i <= hInitElSpectrumWide->GetXaxis()->GetNbins(); i++)
  {
    double x = hInitElSpectrumWide->GetBinCenter(i);
    double value = ( hInitElSpectrumWide->GetBinContent(i) - fitFunction->Eval(x) ) / hInitElSpectrumWide->GetBinError(i);
    hInitElSpectrumWide_Res->SetBinContent(i,value); 
  }

  for(int i = 1; i <= hInitElSpectrumFine->GetXaxis()->GetNbins(); i++)
  {
    double x = hInitElSpectrumFine->GetBinCenter(i);
    double value = ( hInitElSpectrumFine->GetBinContent(i) - fitFunctionFine->Eval(x) ) / hInitElSpectrumFine->GetBinError(i);
    hInitElSpectrumFine_Res->SetBinContent(i,value); 
  }

  for(int i = 1; i <= hUniformElSpectrumWide->GetXaxis()->GetNbins(); i++)
  {
    double x = hUniformElSpectrumWide->GetBinCenter(i);
    double value = ( hUniformElSpectrumWide->GetBinContent(i) - fitFunctionUniform->Eval(x) ) / hUniformElSpectrumWide->GetBinError(i);
    hUniformElSpectrumWide_Res->SetBinContent(i,value); 
  }

  for(int i = 1; i <= hUniformElSpectrumFine->GetXaxis()->GetNbins(); i++)
  {
    double x = hUniformElSpectrumFine->GetBinCenter(i);
    double value = ( hUniformElSpectrumFine->GetBinContent(i) - fitFunctionFineUniform->Eval(x) ) / hUniformElSpectrumFine->GetBinError(i);
    hUniformElSpectrumFine_Res->SetBinContent(i,value); 
  }
    
  int start = hInitElSpectrumWide->FindBin(tFitStart);
  int stop = hInitElSpectrumWide->FindBin(tFitStop);
  for(int i = start; i < stop; i ++)
  {
    hInitElSpectrumWide_Sigma->Fill(hInitElSpectrumWide_Res->GetBinContent(i));
    hUniformElSpectrumWide_Sigma->Fill(hUniformElSpectrumWide_Res->GetBinContent(i));

  }

  start = hInitElSpectrumFine->FindBin(tFitStart);
  stop = hInitElSpectrumFine->FindBin(tFitStop);
  for(int i = start; i < stop; i ++)
  {
    hInitElSpectrumFine_Sigma->Fill(hInitElSpectrumFine_Res->GetBinContent(i));
    hUniformElSpectrumFine_Sigma->Fill(hUniformElSpectrumFine_Res->GetBinContent(i));
  }

  TF1 *f = new TF1("f", "gaus", -4, 4);
  hInitElSpectrumFine_Sigma->Fit("f","RQ");
  hUniformElSpectrumFine_Sigma->Fit("f","RQ");
  
  return 1;
}

int GenerateLtSpectrum(int nEvents)
{
  gStyle->SetOptFit();
  BuildFunction(nEvents);
  InitHistograms();
  if(FillHistograms(nEvents)!=nEvents) cout << "Histograms not properly filled" << endl;;
  if(FitHistograms(nEvents)!=1) cout << "Fitting failed" << endl;
  DrawStuff(nEvents);
  return 1;
}



