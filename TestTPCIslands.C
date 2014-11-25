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

//gROOT->Reset();

using namespace std;

int TestIslands()
{
  TFile* file = TFile::Open("dump88040.root", "READ"); if ( file->IsZombie() ) { cout << " could not open file " << endl; return 1; }
  TTree* tree = (TTree*)file->FindObjectAny("TTPCIslandDumpTree");
  Int_t nentries = tree->GetEntries();
  cout << " File open, tree open, number of entries: " << nentries << endl;

  return 0;
}


