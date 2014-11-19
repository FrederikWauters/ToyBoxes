gROOT->Reset();
gROOT->SetStyle(1);

#include<iostream.h>
#include<TCanvas.h>
#include<TFile.h>
#include<TH1.h>


void main()
{


const Int_t nP4 = 10;
Double_t x4[nP4]={1, 2, 3,4,5,6,7,8,9,10};
Double_t y4[nP4]={10,11.5,12,14,15.3,15,15.5,14.4,13,11};
TGraph *g4 = new TGraph(nP4,x4,y4);

const Int_t nP5 = 5;
Double_t x5[nP5]={1,2,3,4,5};
Double_t y5[nP5]={10,11,16,13,8};
TGraph *g5 = new TGraph(nP5,x5,y5);

Plot(g5);

}

void Plot(TGraph * g)
{
  g->SetMarkerStyle(8);
  g->Draw("ap");

  TSpline3 *s = new TSpline3("grs",g);
  s->SetLineColor(kRed);
  s->Draw("same");
}
