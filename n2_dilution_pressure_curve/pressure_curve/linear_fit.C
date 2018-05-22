// N. Ton
// 05/2018
// Input: data file: pressure (atm), Yield (nb)
// Assume that at zero pressure, yield is zero. Because y_tg reconstruction is not great, so we can't subtract the glass to extract vacuum
// Output: Linear fitting: pressure (vertical) vs yield (nb)
//  
#include "TMinuit.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
const int dim =5;
Float_t t[dim],v[dim],errorz[dim];

//______________________________________________________________________________
Double_t func(float x,Double_t *par)
{
  Double_t value=par[0]*x+par[1];
  return value;
}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  const Int_t nbins = dim;
  Int_t i;

  //calculate chisquare
  Double_t chisq = 0;
  Double_t delta;
  for (i=0;i<nbins; i++) {
    delta  = (v[i]-func(t[i],par))/errorz[i];
    chisq += delta*delta;
  }
  f = chisq;
}

//______________________________________________________________________________
void linear_fit()
{
  gStyle->SetOptFit(0);//1111111);
  TGraphErrors *gr_f1 = new TGraphErrors();
  gr_f1->SetMarkerSize(1.5);
  gr_f1->SetMarkerColor(kBlue);
  gr_f1->SetMarkerStyle(20);
  TGraphErrors *gr_f2 = new TGraphErrors();
  gr_f2->SetMarkerSize(1.5);
  gr_f2->SetMarkerColor(kBlack);
  gr_f2->SetMarkerStyle(20);
  TMultiGraph *mgr= new TMultiGraph();
  mgr->SetTitle(";Yield(nb);Pressure(atm);");
  mgr->Add(gr_f1);
  mgr->Add(gr_f2);
 ifstream data("yield_atm.txt",ios_base::in);
  double atm, yield;
  int count=0;
  while(!data.eof()){
    data>>yield>>atm;
    v[count]=atm;
    t[count]=yield;
    errorz[count]=0.1;
    gr_f1->SetPoint(count,yield,atm);
    gr_f1->SetPointError(count,0,0.1);
    count++;
      }
  gr_f2->SetPoint(0,37.8, 0.122110);

    TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 5 params
    gMinuit->SetFCN(fcn);

    Double_t arglist[10];
    Int_t ierflg = 0;

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);


 
    // Set starting values and step sizes for parameters
    static Double_t vstart[2] ={0.004, 0.0};//{233.0,9.0}; //
    static Double_t step[2] = {0.001,0.001};//{1.0, 0.1};//
    gMinuit->mnparm(0, "slope", vstart[0], step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "bkgd", vstart[1], step[1], 0,0,ierflg);
 

    // Now ready for minimization step
    arglist[0] = 500000;
    arglist[1] =0.1;

    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    double ParErr[2];
    for (int ii=0;ii<2;ii++){
    gMinuit->GetParameter(ii,vstart[ii],ParErr[ii]);
    cout<<vstart[ii]<<" \t"<<ParErr[ii]<<endl;
    }
    TF1 * line0 = new TF1("line0","[0]*x+[1]",0,50);
    line0->SetParName(0,"slope");
    line0->SetParName(1,"bkgd");
    line0->SetParameters(vstart[0],vstart[1]);
    line0->SetNpx(10000);
    gr_f1->Fit("line0");

    // Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    cout<<" minimum chi2 ====="<<amin<<endl;
    cout<<" err definition ="<<errdef<<endl;
    cout<<" number of variable parameter ="<<nvpar<<endl;
    cout<<" number of parameters ="<<nparx<<endl;
    cout<<icstat<<endl;
    delete gMinuit;
  
    TCanvas *cc1 = new TCanvas("cc1","",1200,900);
    mgr->Draw("AP");
  }
