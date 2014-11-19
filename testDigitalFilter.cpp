//Try some digital filter/pulse shaping algoritm
//Firs, we start with a simple step function with exponential decay

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>

#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <TH1.h>
#include "TF1.h"
#include "TStyle.h"

using namespace std;

//Input pulse, start with a very long island

TH1F* hStepInput; //A unit step function at time t = 0
TH1F* hExpInput; //A unit step with exponential decay
TH1F* hSpiceInput; //A spice simulated preamp output

//First delay and subtract filter
TH1F* hStepDS; 
TH1F* hExpDS; 
TH1F* hSpiceDS; 

//High pass deconvolutions
TH1F* hStepHPD; 
TH1F* hExpHPD;
TH1F* hSpiceHPD; 

//High pass deconvolutions
TH1F* hStepTrap; 
TH1F* hExpTrap;
TH1F* hSpiceTrap;


const Double_t tStart = -20.; //start time in us
const Double_t tStop = 100.; //stop time in us
const Double_t sRate = 50.; //sample rate in MHz

Double_t decayTime; //decay time of the hExpInput in us;
Double_t decayTimeSpice; //decay time of the spice simulated pulse

//Filter param in us
const Double_t k_us = 1;
const Double_t l_us = 1.5; 

//filter functions
Double_t d_kl(TH1* h, int n)
{
  Int_t k = (int)(sRate*k_us);
  Int_t l = (int)(sRate*l_us);
 
  return h->GetBinContent(n) - h->GetBinContent(n-k) - h->GetBinContent(n-l) + h->GetBinContent(n-k-l);
}

Double_t HPD(TH1* h, int n, double tau)
{
  Double_t p_n = 0.;
  for(Int_t i = 1; i <= n; i ++) {  p_n = h->GetBinContent(i) + p_n; }

  Double_t M = 1./(exp(1/(sRate*tau))-1);// cout << "M = " << M << endl;

  return (1/M)*(p_n + M*h->GetBinContent(n));
}

Double_t Acc(TH1* h, int n)
{
  Double_t s_n = 0;
  for(Int_t i = 1; i <= n; i ++) {  s_n = h->GetBinContent(i) + s_n; }

  return s_n/(sRate*k_us);
}

void SetInputPulse()
{
  

  //set step Function**********************************************
  Int_t nBins = (int) (sRate*(tStop-tStart));
  hStepInput = new TH1F("hStepInput","Digitized step input signal",nBins,tStart,tStop);

  for(Int_t i = 1; i <= hStepInput->GetNbinsX(); i++)
  {
    if(hStepInput->GetBinCenter(i)>0){ hStepInput->SetBinContent(i,1); }
    else{ hStepInput->SetBinContent(i,0);} 
    if(hStepInput->GetBinCenter(i)>1.5){ hStepInput->SetBinContent(i,3); }
    if(hStepInput->GetBinCenter(i)>20){ hStepInput->SetBinContent(i,5); }
  }  



  //set exp decay function*****************************************
  hExpInput = new TH1F("hExpInput","Digitized exponential input signal",nBins,tStart,tStop);

  TF1* fExp = new TF1("fExp","exp(-x/[0])",0,100);
  decayTime;
  cout << " Give decay time (us) : " << endl;
  cin >> decayTime;
  fExp->SetParameter(0,decayTime);

  for(Int_t i = 1; i <= hExpInput->GetNbinsX(); i++)
  {
    if(hExpInput->GetBinCenter(i)>0){ hExpInput->SetBinContent(i,fExp->Eval(hExpInput->GetBinCenter(i))); }
    else{ hExpInput->SetBinContent(i,0);} 
  }  


  //Spice pulse**************************************************
  hSpiceInput = new TH1F("hSpiceInput","Digitized spice simulated input signal",nBins,tStart,tStop);

  cout << "Give filename inputpulse" << endl;
  char fname[1024];
  cin >> fname;
  string fname_full = "../NoiseAnalysis/";
  fname_full.append(fname);

  ifstream *is = new ifstream((char*)fname_full.c_str());
  if(is!=NULL) { cout << "file " << fname << " open ... " << endl; }
  else { cout << "Could not open file " << fname << " !!!" << endl; }

  Double_t par1, par2;
  vector<Double_t> time;
  vector<Double_t> value;

  while(is->good())
  {
    *is >> par1;
    *is >> par2;
    //cout << "par1 : " << par1 << "   par2 " << par2 << endl;
    time.push_back(par1*1000000);
    value.push_back(-1*par2);
  }

  Double_t binCenter, valueSpice;

  for(Int_t i = 1; i <= hSpiceInput->GetNbinsX(); i++)
  {
    binCenter = hSpiceInput->GetBinCenter(i)+5.;

    if(binCenter>0)
    {
      Int_t j = 0;
      //cout << " time " << time.at(j) << endl;
      while( time.at(j) < binCenter && j+1 < time.size() ) { j++; }
      //cout << "j = " << j << "  size " << time.size() << endl;
      if( j < time.size() ) valueSpice = value.at(j-1) + (value.at(j)-value.at(j-1))*(binCenter - time.at(j-1))/(time.at(j) - time.at(j-1)); //interpolate
      else valueSpice = 0;
      //valueSpice = 0;
      hSpiceInput->SetBinContent(i,valueSpice/0.975);
 
    }
    else{ hSpiceInput->SetBinContent(i,0);} 
  }

  delete is;
  
  cout << " ... finished reading file ... give time constant :" << endl;
  cin >> decayTimeSpice;

}


void DelaySubtract()
{
  Int_t nBins = (int) (sRate*(tStop-tStart));


  hStepDS = new TH1F("hStepDS","Step Fuction after 'delay-subtract filter'",nBins,tStart,tStop);
  hExpDS = new TH1F("hExpDS","Exponential function after 'delay-subtract filter'",nBins,tStart,tStop);
  hSpiceDS = new TH1F("hSpiceDS","Spice simulated pulse after 'delay-subtract filter'",nBins,tStart,tStop);

  for(Int_t i = 1; i <= hStepDS->GetNbinsX(); i++) { hStepDS->SetBinContent(i,d_kl(hStepInput,i)); }
  for(Int_t i = 1; i <= hExpDS->GetNbinsX(); i++) { hExpDS->SetBinContent(i,d_kl(hExpInput,i)); }
  for(Int_t i = 1; i <= hSpiceDS->GetNbinsX(); i++) { hSpiceDS->SetBinContent(i,d_kl(hSpiceInput,i)); }
}

void HighPassFilterDeconvolver()
{
  Int_t nBins = (int) (sRate*(tStop-tStart));

  hStepHPD = new TH1F("hStepHPD","Step Fuction after 'high-pass deconvolver'",nBins,tStart,tStop);
  hExpHPD = new TH1F("hExpHPD","Exponential function after 'high-pass deconvolver'",nBins,tStart,tStop);
  hSpiceHPD = new TH1F("hSpiceHPD","Spice simulated pulse after 'high-pass deconvolver'",nBins,tStart,tStop);

  for(Int_t i = 1; i <= hStepHPD->GetNbinsX(); i++) { hStepHPD->SetBinContent(i,HPD(hStepDS,i,1000000.)); }
  for(Int_t i = 1; i <= hExpHPD->GetNbinsX(); i++) { hExpHPD->SetBinContent(i,HPD(hExpDS,i,decayTime)); }
  for(Int_t i = 1; i <= hSpiceHPD->GetNbinsX(); i++) { hSpiceHPD->SetBinContent(i,HPD(hSpiceDS,i,decayTimeSpice)); }
}

void ConstructTrap()
{
  Int_t nBins = (int) (sRate*(tStop-tStart));

  hStepTrap = new TH1F("hStepTrap","Step Fuction converted into a trapezoid'",nBins,tStart,tStop);
  hExpTrap = new TH1F("hExpTrap","Exponential function converted into a trapezoid'",nBins,tStart,tStop);
  hSpiceTrap = new TH1F("hSpiceTrap","Spice simulated pulse converted into a trapezoid'",nBins,tStart,tStop);

  for(Int_t i = 1; i <= hStepTrap->GetNbinsX(); i++) { hStepTrap->SetBinContent(i,Acc(hStepHPD,i)); }
  for(Int_t i = 1; i <= hExpTrap->GetNbinsX(); i++) { hExpTrap->SetBinContent(i,Acc(hExpHPD,i)); }
  for(Int_t i = 1; i <= hSpiceTrap->GetNbinsX(); i++) { hSpiceTrap->SetBinContent(i,Acc(hSpiceHPD,i)); }
}

void DrawHistos()
{

  TCanvas *c1 = new TCanvas("c1","Step input response",1000,1000);
  //c1->Divide(2,2);
  c1->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hStepInput->GetXaxis()->SetTitle("time (#mus)");
  hStepInput->SetLineWidth(4);
  hStepInput->SetLineColor(12);
  hStepInput->Draw();
  hStepInput->GetYaxis()->SetRangeUser(-3.5,3.5);
  //c1->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hStepDS->SetLineColor(kRed);
  hStepDS->GetXaxis()->SetTitle("time (#mus)");
  hStepDS->SetLineWidth(3);
  hStepDS->SetLineStyle(2);
  hStepDS->Draw("SAME");
  hStepHPD->SetLineColor(kBlue);
  hStepHPD->GetXaxis()->SetTitle("time (#mus)");
  hStepHPD->SetLineWidth(1);
  hStepHPD->Draw("SAME");
  hStepTrap->SetLineColor(46);
  hStepTrap->GetXaxis()->SetTitle("time (#mus)");
  hStepTrap->SetLineWidth(3);
 // hStepTrap->SetLineStyle(3);
  hStepTrap->Draw("SAME");



  TCanvas *c2 = new TCanvas("c2","Exponential input response",1000,1000);
  //c2->Divide(2,2);
  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hExpInput->GetXaxis()->SetTitle("time (#mus)");
  hExpInput->SetLineWidth(4);
  hExpInput->SetLineColor(12);
  hExpInput->Draw();
  hExpInput->GetYaxis()->SetRangeUser(-1.5,1.5);
  //c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hExpDS->GetXaxis()->SetTitle("time (#mus)");
  hExpDS->SetLineColor(kRed);
  hExpDS->SetLineWidth(3);
  hExpDS->SetLineStyle(2);
  hExpDS->Draw("SAME");
  hExpHPD->SetLineColor(kBlue);
  hExpHPD->GetXaxis()->SetTitle("time (#mus)");
  hExpHPD->SetLineWidth(1);
  hExpHPD->Draw("SAME");
  hExpTrap->SetLineColor(46);
  hExpTrap->GetXaxis()->SetTitle("time (#mus)");
  hExpTrap->SetLineWidth(3);
  //hExpTrap->SetLineStyle(3);
  hExpTrap->Draw("SAME");
  
  TCanvas *c3 = new TCanvas("c3","Spice input response",1000,1000);
 // c3->Divide(2,2);
  c3->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hSpiceInput->GetXaxis()->SetTitle("time (#mus)");
  hSpiceInput->SetLineWidth(4);
  hExpInput->SetLineColor(12);
  hSpiceInput->Draw();
  hSpiceInput->GetYaxis()->SetRangeUser(-1.5,1.5);
  //c3->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hSpiceDS->GetXaxis()->SetTitle("time (#mus)");
  hSpiceDS->SetLineWidth(3);
  hSpiceDS->SetLineColor(kRed);
  hSpiceDS->SetLineStyle(2);
  hSpiceDS->Draw("SAME");
  hSpiceHPD->SetLineColor(kBlue);
  hSpiceHPD->GetXaxis()->SetTitle("time (#mus)");
  hSpiceHPD->SetLineWidth(1);
  hSpiceHPD->Draw("SAME");
  hSpiceTrap->SetLineColor(46);
  hSpiceTrap->GetXaxis()->SetTitle("time (#mus)");
  hSpiceTrap->SetLineWidth(3);
//  hSpiceTrap->SetLineStyle(3);
  hSpiceTrap->Draw("SAME");
}



int Filter()
{
  gStyle->SetPalette(1);

  SetInputPulse();

  DelaySubtract();
  HighPassFilterDeconvolver();
  ConstructTrap();

  DrawHistos();

  return 1;
}
