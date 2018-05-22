// Nguyen
// 05/2018
// Input: Nitrogen, empty runs
// Simulate gaussian shape to take into account the N2 excited states contribution
// Output: N2 cross section
// No beam trip cut
// No software collimator cuts
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
void get_xs()
{
  gStyle->SetOptStat(0);
  int nbin=200;

  // Histogram limit
  double hy_min(-0.03),hy_max(0.03);
  double hphitg_min(-0.03),hphitg_max(0.03),hthetatg_min(0.01),hthetatg_max(0.06);
  double hw_min(-0.01),hw_max(0.02);
  // Normalization
  // index 0: MT; index 1: N2, index 2: He3
  int runnum[3];
  double p0[3];
  double eff_vdc[3],eff_pid[3],eff_sci[3];
  double charge[3],ps1[3],lt[3],cer_lo[3],ps_lo[3],ep_lo[3];
  ifstream scaler("../he3_2200MeV_0dp.txt",ios_base::in);
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


  double scl_N_MT = ps1[0]*eff_pid[1]*eff_vdc[1]*eff_sci[1]*charge[1]*lt[1]
    /(ps1[1]*eff_pid[0]*eff_vdc[0]*eff_sci[0]*charge[0]*lt[0]);

  double w_off = -0.0028;// GeV


  /************  Glass *********************/
  TFile* fMT = new TFile("2080_n2_xs.root");
  TH1F * MTyt     = new TH1F("MTyt","",nbin,hy_min,hy_max);
  TH1F * MTw2     = new TH1F("MTw2","",nbin,hw_min,hw_max);
  TTree *tMT=(TTree*)fMT->Get("t1");
  tMT->Draw("rec_y>>MTyt");
  tMT->Draw(Form("rec_w2+%f>>MTw2",w_off));
  double N_glass = MTyt->GetEntries();

  /************  N2 *********************/
  TFile* fN2 = new TFile("2082_n2_xs.root");
  TH1F * Nyt     = new TH1F("Nyt","",nbin,hy_min,hy_max);
  TH1F * Nw2     = new TH1F("Nw2","",nbin,hw_min,hw_max);
  TTree *tN2=(TTree*)fN2->Get("t1");
  tN2->Draw("rec_y>>Nyt");
  tN2->Draw(Form("rec_w2+%f>>Nw2",w_off));
  double N_n2 = Nyt->GetEntries();

 

  /* Calculate cross section */
  double SA = 227.0*1e-6;
  double Y_n2 =       N_n2*ps1[1]*1.6*1e-19/(eff_sci[1]*eff_vdc[1]*eff_pid[1]*lt[1]*charge[1]*1e-6);
  double Y_glass = N_glass*ps1[0]*1.6*1e-19/(eff_sci[0]*eff_vdc[0]*eff_pid[0]*lt[0]*charge[0]*1e-6);

  double den_n2 =       2.0*5.356*2.6894*1e19;//1amg =2.687*1e19  (#/cm3)
  double den_n2_in_he3 =2.0*0.105*2.6894*1e19;
  double t_target = 34.3;//in cm. I get +/-11 cm from simulation
  double den_ratio = 0.026;
  /*  1cm2 = 1e33 nb */
  double xs_n2  = (Y_n2-Y_glass)*1e33/(SA*den_n2*t_target);
  cout<<"n2 cross section= "<<xs_n2<<" [nb]"<<endl;

  
  TH1F * OrMTyt =(TH1F*)MTyt->Clone("OrMTyt");
  TH1F * OrMTw2 =(TH1F*)MTw2->Clone("OrMTw2");

  TFile* f3 = new TFile("/home/ton/Desktop/MC_N2_2200MeV_0dp/cross_section/2082_n2_elastic_w_ana_cut.root");
  TTree *tsim=(TTree*)f3->Get("h1");
  TH1F * Sw2     = new TH1F("Sw2","",nbin,hw_min,hw_max);
  tsim->Draw("wmm-0.002>>Sw2","xs");
  /***************   Generate gaussian distribution ************/
  TF1 *mygaus = new TF1("mygaus","[0]*exp(-0.5*((x-[1])/[2])**2)",hw_min, hw_max);
  mygaus->SetParameter(0, 1);
  mygaus->SetParameter(1, 0.031);
  mygaus->SetParameter(2, 0.014); 
  TH1F *hGaus = new TH1F("hGaus","",nbin, hw_min, hw_max); 
  hGaus->FillRandom("mygaus",1000000); 
  hGaus->Scale(1.0);
  hGaus->Draw();

  /*************************************/
  TCanvas * c_w2 = new TCanvas("c_w2","",800,600);
  c_w2->Clear();
  Nw2->SetLineColor(kRed);
  Sw2->SetLineColor(kBlack);
  MTw2->Scale(scl_N_MT);
  Nw2->Add(MTw2,-1);
  Nw2->Draw();
  MTw2->Draw("same");
  double scl_sim =Nw2->GetMaximum()/Sw2->GetMaximum();
  Sw2->Scale(scl_sim);
  Sw2->Draw("same");
  TH1F * residual =(TH1F*)Nw2->Clone("residual");
  residual->Add(Sw2,-1);
  residual->SetLineColor(kBlue);
  residual->Draw("same");
  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->AddEntry(Nw2,"Data","l");
  legend->AddEntry(Sw2,"Simulation","l");
  legend->Draw();

  /*************************************/
  TAxis *axis = Nw2->GetXaxis();
  int bmin = axis->FindBin(-0.0093); //in your case xmin=-1.5
  int bmax = axis->FindBin(0.0087); //in your case xmax=0.8
  double integral = Nw2->Integral(bmin,bmax);
  cout<<"integral of n2= "<<integral<<endl;

  TAxis *axisG = residual->GetXaxis();
  int bminG = axisG->FindBin(-0.0093); //in your case xmin=-1.5
  int bmaxG = axisG->FindBin(0.0087); //in your case xmax=0.8
  double integralG = residual->Integral(bminG,bmaxG);
  cout<<"Gaussian contribution= "<<integralG<<endl;
}
