// N. Ton
// 05/2018
// Input: Nitrogen and empty runs
// Output: N2 yield in arbitrary unit (nb)
#include <iomanip>
#include <vector>
#include "TMath.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include <cstdlib>
const int NPOINTS=6;
void get_yield()
{
  gStyle->SetOptStat(0);
  int nbin=200;

  // Histogram limit
  double hy_min(-0.03),hy_max(0.03);
  double hphitg_min(-0.03),hphitg_max(0.03),hthetatg_min(0.01),hthetatg_max(0.06);
  double hw_min(-0.01),hw_max(0.02);

  int runnum[NPOINTS];
  double p0[NPOINTS];
  double eff_vdc[NPOINTS],eff_pid[NPOINTS],eff_sci[NPOINTS];
  double charge[NPOINTS],ps1[NPOINTS],lt[NPOINTS],cer_lo[NPOINTS],ps_lo[NPOINTS],ep_lo[NPOINTS];

  ifstream scaler("n2_2200MeV_0dp.txt",ios_base::in);
  int count=0;
  if(!scaler.is_open())
    {
      cerr<<"Check input file. File not found!"<<endl;
      exit(EXIT_FAILURE);
    }
  int n2run=2087;
  int n2_indx;
  while(!scaler.eof()){
    scaler>>runnum[count]>>p0[count]>>ps1[count]>>lt[count]>>charge[count]>>ep_lo[count]>>ps_lo[count]>>cer_lo[count]>>eff_sci[count]>>eff_pid[count]>>eff_vdc[count];
    if(runnum[count]==n2run)
      n2_indx=count;
    count++;
  }
  // 0 glass
  // 1 2081
  // 2 2082
  // 3 2083
  // 4 2085

  double scl = ps1[0]*eff_pid[n2_indx]*eff_vdc[n2_indx]*eff_sci[n2_indx]*charge[n2_indx]*lt[n2_indx]
    /(ps1[n2_indx]*eff_pid[0]*eff_vdc[0]*eff_sci[0]*charge[0]*lt[0]);

  /* To shift w peak */
  double w_off_2081= 0.00041;
  double w_off_2082= 0.00047;//GeV n2 shift
  double w_off_2083= 0.00075;//GeV n2 shift
  double w_off_2085= 0.00084;//GeV n2 shift
  double w_off_2087= 10.0;//0.00011;//GeV n2 shift
  double w_off = 0.0;


  /************  Glass *********************/
  TFile* fMT = new TFile("2080_n2_yield.root");
  TH1F * MTw2     = new TH1F("MTw2","",nbin,hw_min,hw_max);
  TTree *tMT=(TTree*)fMT->Get("t1");
  tMT->Draw(Form("rec_w_he3+%f>>MTw2",w_off),Form("rec_w_he3+%f<=%f",w_off,w_off_2087));
  double N_glass = MTw2->GetEntries();

  TFile* fN2 = new TFile("2087_n2_yield.root");
  TH1F * Nw2     = new TH1F("Nw2","",nbin,hw_min,hw_max);
  TTree *tN2=(TTree*)fN2->Get("t1");
  tN2->Draw(Form("rec_w_he3+%f>>Nw2",w_off),Form("rec_w_he3+%f<=%f",w_off,w_off_2087));
  double N_n2 = Nw2->GetEntries();
  cout<<N_glass<<"\t N2 count= "<<N_n2<<endl;
  /* Calculate cross section */
  double SA = 227.0*1e-6;
  double Y_n2 = N_n2*ps1[n2_indx]*1.6*1e-19/(eff_sci[n2_indx]*eff_vdc[n2_indx]*eff_pid[n2_indx]*lt[n2_indx]*charge[n2_indx]*1e-6);
  double Y_glass = N_glass*ps1[0]*1.6*1e-19/(eff_sci[0]*eff_vdc[0]*eff_pid[0]*lt[0]*charge[0]*1e-6);

  double t_target = 34.3;//in cm. I get +/-11 cm from simulation

  /*  1e15 arbitrary unit */
  double xs_n2  = (Y_n2-Y_glass)*1e15/(t_target);
  cout<<"Yield: "<<Y_glass<<"\t N2 yield= "<<Y_n2<<endl;
  cout<<"Total cross section= "<<xs_n2<<" [arbitrary nb]"<<endl;



  TFile* f3 = new TFile("/home/ton/Desktop/MC_He_2200MeV_0dp/cross_section/2087_elastic_he3_ana_cut.root");
  TTree *tsim=(TTree*)f3->Get("h1");
  TH1F * Sw2     = new TH1F("Sw2","",nbin,hw_min,hw_max);
  tsim->Draw("wmm+0.0009>>Sw2","xs");//*(wmm<=0.007&&wmm>=-0.001)");
  TCanvas * c_w2 = new TCanvas("c_w2","",800,600);
  c_w2->Clear();
  Nw2->SetLineColor(kRed);
  Sw2->SetLineColor(kBlack);
  MTw2->Scale(scl);
  Nw2->Add(MTw2,-1);
  Nw2->Draw();
  double scl_sim =Nw2->GetMaximum()/Sw2->GetMaximum();
  Sw2->Scale(scl_sim);
  Sw2->Draw("same");
  MTw2->Draw("same");
}
