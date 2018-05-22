// Nguyen
// 05/05/2017
// Input: root file for He3, Nitrogen, empty runs
// Output: histogram for (+) and (-) helicity to calculate asymmetry
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
void plot_all()
{
  gStyle->SetOptStat(0);
  int nbin=150;

  // Histogram limit
  double hy_min(-0.03),hy_max(0.03), hdp_min(-0.02),hdp_max(0.005);
  double hphitg_min(-0.03),hphitg_max(0.03),hthetatg_min(0.01),hthetatg_max(0.06);
  double hw_min(-0.005),hw_max(0.02);
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
  cout<<"scaling factor= "<<scl_N_MT<<endl;
  /* To shift w peak */
  double w_off = -0.0029;// GeV
  double w_glass_off=0.00;//GeV glass shift
  double w_n2_off= -0.000;//GeV n2 shift

    /*********************************/
  TFile* fMT = new TFile("2080_n2_fp_cut.root");
  TH1F * MTw2     = new TH1F("MTw2","",nbin,hw_min,hw_max);
  TTree *tMT=(TTree*)fMT->Get("t1");
  tMT->Draw(Form("rec_w2+%f+%f>>MTw2",w_off,w_glass_off));
  
  TFile* fN2 = new TFile("2082_n2_fp_cut.root");
  TH1F * Nw2     = new TH1F("Nw2","",nbin,hw_min,hw_max);
  TTree *tN2=(TTree*)fN2->Get("t1");
  tN2->Draw(Form("rec_w2+%f+%f>>Nw2",w_off,w_n2_off));//,Form("rec_y>=%f&&rec_y<=%f&&(rec_w2+%f)<=%f",y_min,y_max,w_off,w_max));

  

  TFile* f3 = new TFile("/home/ton/Desktop/MC_N2_2200MeV_0dp/cross_section/2082_n2_elastic.root");
  TTree *tsim=(TTree*)f3->Get("h1");
  TH1F * Sw2     = new TH1F("Sw2","",nbin,hw_min,hw_max);
  tsim->Draw("wmm>>Sw2","xs");//*(wmm<=0.007&&wmm>=-0.001)");
  /***************   Generate gaussian distribution ************/
  TF1 *mygaus = new TF1("mygaus","[0]*exp(-0.5*((x-[1])/[2])**2)",hw_min, hw_max);
  mygaus->SetParameter(0, 1);
  mygaus->SetParameter(1, 0.009);
  mygaus->SetParameter(2, 0.0018); 
  TH1F *hGaus = new TH1F("hGaus","",nbin, hw_min, hw_max); 
  hGaus->FillRandom("mygaus",10000); 
  hGaus->Scale(2.6);
  hGaus->Draw("same");
  /*************************************/
  TAxis *axis = Nw2->GetXaxis();
  int bmin = axis->FindBin(-0.001); //in your case xmin=-1.5
  int bmax = axis->FindBin(0.007); //in your case xmax=0.8
  double integral = Nw2->Integral(bmin,bmax);
  cout<<"integral of n2= "<<integral<<endl;

  TAxis *axisG = hGaus->GetXaxis();
  int bminG = axisG->FindBin(-0.001); //in your case xmin=-1.5
  int bmaxG = axisG->FindBin(0.007); //in your case xmax=0.8
  double integralG = hGaus->Integral(bminG,bmaxG);
  cout<<"Gaussian contribution= "<<integralG<<endl;
  /*************************************/
  TCanvas * c_w2 = new TCanvas("c_w2","",800,600);
  c_w2->Clear();
  Nw2->SetLineColor(kRed);
  Sw2->SetLineColor(kBlack);
  hGaus->SetLineColor(kMagenta);
  MTw2->Scale(scl_N_MT);
  Nw2->Add(MTw2,-1);
  Nw2->Draw();
  MTw2->Draw("same");
  hGaus->Draw("same");
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
  legend->AddEntry(Nw2,"Clean N2","l");
  legend->AddEntry(Sw2,"Simulation","l");
  legend->Draw();


}
