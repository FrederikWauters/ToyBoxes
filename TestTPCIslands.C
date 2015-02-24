#include "TRandom2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TDatime.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraphSmooth.h"



using namespace std;

//global var
TFile* file;
TTree* tree;
Int_t nentries;
vector<Short_t> samples;
vector<double> filteredsamples;
Int_t fBlockTime;
TTPCIsland *island;
static TBranch *branch;
TH1F* hIsland;
TH1F* hFilteredIsland;
TCanvas* cIsland;
TTPCMiniPulseFinder* pulsefinder;
TF1* fPulse1 ;
TF1* fPulse2;
TLine *pedestalLine;

int loopcounter = 1;


class TPCMINIPF_PARAM  {
  public:
  Int_t kMiniPedestalChoice;
  Int_t kMiniTPedestalChoice;
  Int_t kLibraryPedestalChoice;
  Int_t kMiniTemplateChoice;
  Int_t kMiniTTemplateChoice;
  Int_t kLibraryTemplateChoice;
  Bool_t kMiniFloatingPed;
  Bool_t kMiniTFloatingPed;
  Bool_t kLibraryFloatingPed;
  Int_t kMiniFilterChoice;
  Int_t kMiniTFilterChoice;
  Int_t kLibraryFilterChoice;
  Int_t kMiniRunChoice;
  double kBadFitChi2Threshold;

};

TPCMINIPF_PARAM tpcminipf_parameters;

class TPC_PARAM {
 public:
   Double_t   kClockPeriod;
   Double_t  fDefaultFixedPedestalValues[48];// = { 18.2961,21.600,17.6631,18.9222,18.4078,20.0988,18.357,26.0079,  19.8125,22.2841,17.8005,19.6552,22.9703,23.128,20.8571,22.0123, 19.0476,18.9536,18.086,19.2679	,23.4196,21.9723,20.7179,19.2699,    12.2999,19.8472,16.4372,18.3173,17.9134,20.5568,20.3518,20.1024,    18.0075,20.0779,20.4841,17.6332,20.5541,17.0828,17.4266,17.5412,   21.4699,20.6677,19.5624,19.8547,17.0872,17.7549	,20.1049,17.8361};
};

TPC_PARAM tpc_parameters;

void Init_Param()
{
  tpcminipf_parameters.kMiniPedestalChoice = 1;
  tpcminipf_parameters.kMiniTemplateChoice = 2;
  tpcminipf_parameters.kMiniFloatingPed = false;
  tpcminipf_parameters.kMiniFilterChoice=2;
  tpcminipf_parameters.kMiniRunChoice = 1;
  tpcminipf_parameters.kBadFitChi2Threshold = 9;

  tpc_parameters.kClockPeriod = 41.;

  tpc_parameters.fDefaultFixedPedestalValues = new double[48];
  double defaultRun7Values_H[48] = { 18.2961,21.600,17.6631,18.9222,18.4078,20.0988,18.357,26.0079,  19.8125,22.2841,17.8005,19.6552,22.9703,23.128,20.8571,22.0123, 19.0476,18.9536,18.086,19.2679	,23.4196,21.9723,20.7179,19.2699,    12.2999,19.8472,16.4372,18.3173,17.9134,20.5568,20.3518,20.1024,    18.0075,20.0779,20.4841,17.6332,20.5541,17.0828,17.4266,17.5412,   21.4699,20.6677,19.5624,19.8547,17.0872,17.7549	,20.1049,17.8361};
  for(int i=0; i<48; i++) { tpc_parameters.fDefaultFixedPedestalValues[i]=defaultRun7Values_H[i];}

  pulsefinder->SetPSIRunNumber(7);
  pulsefinder->SetTPCParameters(tpc_parameters);
  pulsefinder->SetTPCMiniParameters(tpcminipf_parameters);
  pulsefinder->MakeAllHistos();

}

void DrawIsland(int entry)
{
  if(hIsland!=NULL) delete hIsland;
  if(hFilteredIsland!=NULL) delete hFilteredIsland;
  if(cIsland!=NULL) delete cIsland;

  cIsland = new TCanvas("cIsland","Islands",600,600);
  if(entry>=nentries) { cout << " entry not in tree " << endl; return; }
  tree->GetEntry(entry);

  samples = island->GetSampleVector();
  filteredsamples = island->GetFilteredSampleVector();

  if(samples.size() > 0)
  {
    hIsland = new TH1F("hIsland","TPC Island",samples.size(),0.5,samples.size()+1);
    hFilteredIsland = new TH1F("hFilteredIsland","TPC Island",filteredsamples.size(),0.5,filteredsamples.size()+1);
   
    hIsland->SetLineColor(26);
    hFilteredIsland->SetLineColor(1);    

    for(int i = 0; i < samples.size(); i++)
    {
      hIsland->SetBinContent(i,samples.at(i));
      hFilteredIsland->SetBinContent(i,filteredsamples.at(i));
    }

    vector<TTPCMiniPulse*> pulses = ProcessIsland(entry,true);
    int nPulses = 0;    

    hIsland->Draw();
    hFilteredIsland->Draw("SAME");
    
    
    int padNr = island->GetPadZeroIndexed();
    double pedestal = tpc_parameters.fDefaultFixedPedestalValues[padNr];
    pedestalLine = new TLine(0,pedestal,samples.size()+1,pedestal);
    pedestalLine->SetLineColor(kRed);
    pedestalLine->Draw("SAME");
    vector<double> fitParams;


    if(pulses.size()>0)
    {
      for(int i = 0; i < pulses.size(); i++)
      {
        if(i==0) 
        {  
           
           fitParams.push_back( pulses.at(0)->GetAmplitude() ); 
           fitParams.push_back( pulses.at(0)->GetCenterTimeInIsland()/tpc_parameters.kClockPeriod  -1); 
           fitParams.push_back( pulses.at(0)->GetSigma()/tpc_parameters.kClockPeriod );
           fitParams.push_back(pulses.at(0)->GetPedestal()); for(int j = 0; j < 4; j++) cout << fitParams.at(j) << endl;
           fPulse1=pulsefinder->GetFunction(fitParams);
           //fPulse1->SetRange(1,100);
           fPulse1->SetLineColor(3);
           fPulse1->Draw("SAME"); 
         }
        if(i==1) 
        {  
           
           fitParams.push_back( pulses.at(0)->GetAmplitude() ); 
           fitParams.push_back( pulses.at(0)->GetCenterTimeInIsland()/tpc_parameters.kClockPeriod  -1); 
           fitParams.push_back( pulses.at(0)->GetSigma()/tpc_parameters.kClockPeriod );
           fitParams.push_back(pulses.at(0)->GetPedestal()); for(int j = 0; j < 4; j++) cout << fitParams.at(j) << endl;
           fPulse2=pulsefinder->GetFunction(fitParams);
           fPulse2->SetLineColor(4);
           fPulse2->Draw("SAME"); 
         }
         nPulses++;
      }
    }
  }
  else { cout << "no island to draw" << endl; return; }  

}

void Loop()
{
  DrawIsland(loopcounter);
  loopcounter++;
  return;
}

vector<TTPCMiniPulse*>  ProcessIsland(int entry, bool print = false, bool addIsland = false)
{
  cout << "Process island " << endl;
  tree->GetEntry(entry);
  vector<TTPCGenericPulse*> pulses =  pulsefinder->Process(island);
  //AddIsland(
  vector<TTPCMiniPulse*> vOut;
  cout << "size pulses " << pulses.size() << endl;
  if(pulses.size()>0)
  {
    for(int i = 0; i < pulses.size(); i++)
    {
      if(print) (TTPCMiniPulse*)pulses.at(i)->Print();
      vOut.push_back((TTPCMiniPulse*)pulses.at(i));
      
    }
  }
  return vOut;
}

/*void AddPulse(double time, double amplitude, entry)
{
  int padNr = ile->GetPadZeroIndexed();
  double pedestal = tpc_parameters.fDefaultFixedPedestalValues(padNr);
  double pedestal = ile
}*/

int Setup()
{

  

  file = TFile::Open("dump88040.root", "READ"); if ( file->IsZombie() ) { cout << " could not open file " << endl; return 1; }
  tree = (TTree*)file->FindObjectAny("TTPCIslandDumpTree");
  nentries = tree->GetEntries();
  cout << " File open, tree open, number of entries: " << nentries << endl;

  island = new TTPCIsland();

  Int_t entry = 1;
  br = tree->GetBranch("dumped_objects");
  br->SetAddress(&island);

  //test tree input
  tree->GetEntry(entry); cout << "First Island has " << island->GetNSamples() << " samples " << endl;


  pulsefinder = new TTPCMiniPulseFinder(); //before init_param!
  Init_Param();

  

  //test pulses finder
  vector<TTPCMiniPulse*> pulses = ProcessIsland(3);
  //cout << "size pulses " << pulses.size() << endl;
  //pulses.at(0)->Print();

  return 0;
}

