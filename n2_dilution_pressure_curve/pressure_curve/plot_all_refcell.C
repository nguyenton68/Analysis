// Nguyen
// 05/2018
// Input: root file for all Nitrogen reference cells, empty run
// Output: plot of all reference cell with glass subtracted

#include <iomanip>
#include "TMath.h"
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
void plot_all_refcell()
{
  gStyle->SetOptStat(0);
  int nbin=150;

  // Histogram limit
  double hy_min(-0.03),hy_max(0.03), hdp_min(-0.02),hdp_max(0.005);
  double hphitg_min(-0.03),hphitg_max(0.03),hthetatg_min(0.01),hthetatg_max(0.06);
  double hw_min(-0.01),hw_max(0.02);
  // Normalization
  // index 0: MT; index 1: N2, index 2: He3
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
  while(!scaler.eof()){
    scaler>>runnum[count]>>p0[count]>>ps1[count]>>lt[count]>>charge[count]>>ep_lo[count]>>ps_lo[count]>>cer_lo[count]>>eff_sci[count]>>eff_pid[count]>>eff_vdc[count];
    count++;
  }
  // 0 glass
  // 1 2081
  // 2 2082
  // 3 2083
  // 4 2085
  double scl_2081 = ps1[0]*eff_pid[1]*eff_vdc[1]*eff_sci[1]*charge[1]*lt[1]
    /(ps1[1]*eff_pid[0]*eff_vdc[0]*eff_sci[0]*charge[0]*lt[0]);
  double scl_2082 = ps1[0]*eff_pid[2]*eff_vdc[2]*eff_sci[2]*charge[2]*lt[2]
    /(ps1[2]*eff_pid[0]*eff_vdc[0]*eff_sci[0]*charge[0]*lt[0]);
  double scl_2083 = ps1[0]*eff_pid[3]*eff_vdc[3]*eff_sci[3]*charge[3]*lt[3]
    /(ps1[3]*eff_pid[0]*eff_vdc[0]*eff_sci[0]*charge[0]*lt[0]);
  double scl_2085 = ps1[0]*eff_pid[4]*eff_vdc[4]*eff_sci[4]*charge[4]*lt[4]
    /(ps1[4]*eff_pid[0]*eff_vdc[0]*eff_sci[0]*charge[0]*lt[0]);
  double scl_2087 = ps1[0]*eff_pid[5]*eff_vdc[5]*eff_sci[5]*charge[5]*lt[5]
    /(ps1[5]*eff_pid[0]*eff_vdc[0]*eff_sci[0]*charge[0]*lt[0]);

  /* To shift w peak */
  double w_off_2081= -0.0003;//GeV n2 shift
  double w_off_2082= -0.00036;//GeV n2 shift
  double w_off_2083= -0.00064;//GeV n2 shift
  double w_off_2085= -0.00073;//GeV n2 shift
  double w_off_2087= 0.0;
    /*********************************/
  TFile* fMT = new TFile("2080_n2_fp_cut.root");
  TH1F * MTw2     = new TH1F("MTw2","",nbin,hw_min,hw_max);
  TTree *tMT=(TTree*)fMT->Get("t1");
  tMT->Draw("rec_w_he3>>MTw2");
  
  TFile* fN2_2081 = new TFile("2081_n2_fp_cut.root");
  TH1F * Nw2_2081     = new TH1F("Nw2_2081","",nbin,hw_min,hw_max);
  TTree *tN2_2081=(TTree*)fN2_2081->Get("t1");
  tN2_2081->Draw(Form("rec_w_he3+%f>>Nw2_2081",w_off_2081));

  TFile* fN2_2082 = new TFile("2082_n2_fp_cut.root");
  TH1F * Nw2_2082     = new TH1F("Nw2_2082","",nbin,hw_min,hw_max);
  TTree *tN2_2082=(TTree*)fN2_2082->Get("t1");
  tN2_2082->Draw(Form("rec_w_he3+%f>>Nw2_2082",w_off_2082)); 

  TFile* fN2_2083 = new TFile("2083_n2_fp_cut.root");
  TH1F * Nw2_2083     = new TH1F("Nw2_2083","",nbin,hw_min,hw_max);
  TTree *tN2_2083=(TTree*)fN2_2083->Get("t1");
  tN2_2083->Draw(Form("rec_w_he3+%f>>Nw2_2083",w_off_2083)); 
  
  TFile* fN2_2085 = new TFile("2085_n2_fp_cut.root");
  TH1F * Nw2_2085     = new TH1F("Nw2_2085","",nbin,hw_min,hw_max);
  TTree *tN2_2085=(TTree*)fN2_2085->Get("t1");
  tN2_2085->Draw(Form("rec_w_he3+%f>>Nw2_2085",w_off_2085)); 

  TFile* fN2_2087 = new TFile("2087_n2_fp_cut.root");
  TH1F * Nw2_2087     = new TH1F("Nw2_2087","",nbin,hw_min,hw_max);
  TTree *tN2_2087=(TTree*)fN2_2087->Get("t1");
  tN2_2087->Draw(Form("rec_w_he3+%f>>Nw2_2087",w_off_2087)); 

  TH1F * glass_2081 =(TH1F*)MTw2->Clone("glass_2081");
  TH1F * glass_2082 =(TH1F*)MTw2->Clone("glass_2082");
  TH1F * glass_2083 =(TH1F*)MTw2->Clone("glass_2083");
  TH1F * glass_2085 =(TH1F*)MTw2->Clone("glass_2085");
  TH1F * glass_2087 =(TH1F*)MTw2->Clone("glass_2087");
  /*************************************/
  TCanvas * c_w2 = new TCanvas("c_w2","",800,600);
  c_w2->Clear();
  Nw2_2081->SetLineColor(kRed);
  Nw2_2082->SetLineColor(kBlue);
  Nw2_2083->SetLineColor(kMagenta);
  Nw2_2085->SetLineColor(kGreen);
  Nw2_2087->SetLineColor(kBlack);

  glass_2081->Scale(scl_2081);
  glass_2082->Scale(scl_2082);
  glass_2083->Scale(scl_2083);
  glass_2085->Scale(scl_2085);
  glass_2087->Scale(scl_2087);

  Nw2_2081->Add(glass_2081,-1);
  Nw2_2082->Add(glass_2082,-1);
  Nw2_2083->Add(glass_2083,-1);
  Nw2_2085->Add(glass_2085,-1);
  Nw2_2087->Add(glass_2087,-1);

  Nw2_2081->Draw();
  Nw2_2082->Draw("same");
  Nw2_2083->Draw("same");
  Nw2_2085->Draw("same");
  Nw2_2087->Draw("same");
  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->AddEntry(Nw2_2081,"N2 29 psi","l");
  legend->AddEntry(Nw2_2082,"N2 73 psi","l");
  legend->AddEntry(Nw2_2083,"N2 110 psi","l");
  legend->AddEntry(Nw2_2085,"N2 147 psi","l");
  legend->Draw();
}
