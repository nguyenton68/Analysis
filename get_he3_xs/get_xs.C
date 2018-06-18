// Nguyen
// 10/04/2017
// Input: root file for He3, Nitrogen, empty runs
// Output: histogram for He3 cross section
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
  double hw_min(-0.005),hw_max(0.02);
  // Normalization
  // index 0: MT; index 1: N2, index 2: He3
  int runnum[3];
  double p0[3];
  double eff_vdc[3],eff_pid[3],eff_sci[3];
  double charge[3],ps1[3],lt[3],cer_lo[3],ps_lo[3],ep_lo[3];
  ifstream scaler("he3_0dp_runs.txt",ios_base::in);
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

  double scl_He3_MT = ps1[0]*eff_pid[2]*eff_vdc[2]*eff_sci[2]*charge[2]*lt[2]
    /(ps1[2]*eff_pid[0]*eff_vdc[0]*eff_sci[0]*charge[0]*lt[0]);

  /*Check this number NiDens */
  double NiDens = 0.026;//1.751;//0.6667; ratio=2*dentar_n2/(dentar+2*dentar_n2)
  double scl_He3_N = NiDens*ps1[1]*eff_pid[2]*eff_vdc[2]*eff_sci[2]*charge[2]*lt[2]
    /(ps1[2]*eff_pid[1]*eff_vdc[1]*eff_sci[1]*charge[1]*lt[1]);

  cout<<"Two scaling factor= "<<scl_N_MT<<"\t"<<scl_He3_MT<<"\t"<<scl_He3_N<<endl;
  double w_off = 0.0016;// GeV
  double w_n2_off= -0.0003;//GeV n2 shift

    /************  Glass *********************/
  TFile* fMT = new TFile("1901_xs_155data_51pars.root");
  TH1F * MTyt     = new TH1F("MTyt","",nbin,hy_min,hy_max);
  TH1F * MTw2     = new TH1F("MTw2","",nbin,hw_min,hw_max);
  TTree *tMT=(TTree*)fMT->Get("t1");
  tMT->Draw("rec_y>>MTyt");
  tMT->Draw(Form("rec_w2+%f>>MTw2",w_off));
  double N_glass = MTyt->GetEntries();

  /************  N2 *********************/
  TFile* fN2 = new TFile("1903_xs_155data_51pars.root");
  TH1F * Nyt     = new TH1F("Nyt","",nbin,hy_min,hy_max);
  TH1F * Nw2     = new TH1F("Nw2","",nbin,hw_min,hw_max);
  TTree *tN2=(TTree*)fN2->Get("t1");
  tN2->Draw("rec_y>>Nyt");
  tN2->Draw(Form("rec_w2+%f+%f>>Nw2",w_off,w_n2_off));
  double N_n2 = Nyt->GetEntries();

  /************  He3 *********************/
  TFile* fHe3 = new TFile("1904_xs_155data_51pars.root");
  TH1F * HPyt     = new TH1F("HPyt","",nbin,hy_min,hy_max);
  TH1F * HMyt     = new TH1F("HMyt","",nbin,hy_min,hy_max);
  TH1F * Hyt     = new TH1F("Hyt","",nbin,hy_min,hy_max);
  TH1F * Hw2     = new TH1F("Hw2","",nbin,hw_min,hw_max);
  TTree *tHe3=(TTree*)fHe3->Get("t1");
  tHe3->Draw("rec_y>>Hyt");
  tHe3->Draw(Form("rec_w2+%f>>Hw2",w_off));
  double N_he3 = Hyt->GetEntries();
  cout<<N_glass<<" n2= "<<N_n2<<" he3= "<<N_he3<<endl;
  /* Calculate cross section */
  double SA = 234.0*1e-6;
  double Y_he3 =     N_he3*ps1[2]*1.6*1e-19/(eff_sci[2]*eff_vdc[2]*eff_pid[2]*lt[2]*charge[2]*1e-6);
  double Y_n2 =       N_n2*ps1[1]*1.6*1e-19/(eff_sci[1]*eff_vdc[1]*eff_pid[1]*lt[1]*charge[1]*1e-6);
  double Y_glass = N_glass*ps1[0]*1.6*1e-19/(eff_sci[0]*eff_vdc[0]*eff_pid[0]*lt[0]*charge[0]*1e-6);

  double den_he3 =    7.885*2.6894*1e19;//1amg =2.687*1e19  (#/cm3) 
  double den_n2 = 2.0*7.087*2.6894*1e19;//1amg =2.687*1e19  (#/cm3)
  double den_n2_in_he3 =2.0*0.105*2.6894*1e19;
  double t_target = 34.3;//in cm. I get +/-11 cm from simulation
  double den_ratio = 0.026;
  double n2_correction =1.8;
  /*  1cm2 = 1e33 nb */
  double xs_he3 = (Y_he3-Y_glass)*1e33/(SA*den_he3*t_target);
  double xs_n2  = (Y_n2-Y_glass)*1e33/(SA*den_n2*t_target);
  double xs_tot = xs_he3 - n2_correction*den_ratio*xs_n2;
  cout<<"Total cross section= "<<xs_tot<<" [nb]"<<endl;
  cout<<"************************************"<<endl;
  cout<<"(Y_he-Y_n2)/den_he"<<endl;
  double ni_den_ratio = den_n2_in_he3/den_n2;
  double Y_clean_he3 = Y_he3 - Y_glass - n2_correction*(Y_n2-Y_glass)*ni_den_ratio;
  double xs_2 = Y_clean_he3*1e33/(SA*den_he3*t_target);
  cout<<"2nd method= "<<xs_2<<"\t [nb]"<<endl;

  TH1F * OrMTyt =(TH1F*)MTyt->Clone("OrMTyt");
  TH1F * OrMTw2 =(TH1F*)MTw2->Clone("OrMTw2");

  TFile* f3 = new TFile("/home/ton/Desktop/MC_He_1100MeV_0dp/cross_section/1904_elastic.root");
  TTree *tsim=(TTree*)f3->Get("h1");
  TH1F * Sphi = new TH1F("Sphi","",nbin,hphitg_min,hphitg_max);
  TH1F * Stheta = new TH1F("Stheta","",nbin,hthetatg_min,hthetatg_max);
  TH1F * Syt     = new TH1F("Syt","",nbin,hy_min,hy_max);
  TH1F * Sw2     = new TH1F("Sw2","",nbin,hw_min,hw_max);
  tsim->Draw("wmm>>Sw2","xs");

  /************** Plot yt ************/
  /*
  TCanvas * c_yt = new TCanvas("c_yt","",800,600);
  c_yt->Clear();

  TH1F * allHe3 =(TH1F*)Hyt->Clone("allHe3");
  OrMTyt->SetLineColor(kBlue);
  Hyt->SetLineColor(kRed);
  Nyt->SetLineColor(kCyan);
  allHe3->SetLineColor(kGreen);
  // clean N2
  MTyt->Scale(scl_N_MT);
  Nyt->Add(MTyt,-1);

  OrMTyt->Scale(scl_He3_MT);
  Nyt->Scale(scl_He3_N);
  Hyt->Add(OrMTyt,-1);
  TH1F * cleanHe3 =(TH1F*)Hyt->Clone("cleanHe3");
  cleanHe3->SetLineColor(kMagenta);
  cleanHe3->Add(Nyt,-1);
  allHe3->Draw();
  Hyt->Draw("same");
  Nyt->Draw("same");
  OrMTyt->Draw("same");
  cleanHe3->Draw("same");
  // double scl_sim =Syt->GetMaximum()/cleanHe3->GetMaximum();
  //Syt->Scale(scl_sim);
  //Syt->Draw("same");
  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->AddEntry(cleanHe3,"clean He3","l");
  legend->AddEntry(OrMTyt,"Empty","l");
  legend->AddEntry(Nyt,"Clean N2","l");
  legend->Draw();

  */
  /*************Plot phi target  ********/
  /*  TCanvas * c_phi = new TCanvas("c_phi","",800,600);
  c_phi->Clear();

  TH1F * allHe3_phi =(TH1F*)Hphi->Clone("allHe3_phi");
  OrMTphi->SetLineColor(kBlue);
  Hphi->SetLineColor(kRed);
  Nphi->SetLineColor(kCyan);
  allHe3_phi->SetLineColor(kGreen);
  // clean N2
  MTphi->Scale(scl_N_MT);
  Nphi->Add(MTphi,-1);

  OrMTphi->Scale(scl_He3_MT);
  Nphi->Scale(scl_He3_N);
  Hphi->Add(OrMTphi,-1);
  TH1F * cleanHe3_phi =(TH1F*)Hphi->Clone("cleanHe3_phi");
  cleanHe3_phi->SetLineColor(kMagenta);
  cleanHe3_phi->Add(Nphi,-1);
  allHe3_phi->Draw();
  Hphi->Draw("same");
  Nphi->Draw("same");
  OrMTphi->Draw("same");
  cleanHe3_phi->Draw("same");
  // double scl_sim =Sphi->GetMaximum()/cleanHe3->GetMaximum();
  //Sphi->Scale(scl_sim);
  //Sphi->Draw("same");
  */
  /********************************/
  /*  TCanvas * c_theta = new TCanvas("c_theta","",800,600);
  c_theta->Clear();

  TH1F * allHe3_theta =(TH1F*)Htheta->Clone("allHe3_theta");
  OrMTtheta->SetLineColor(kBlue);
  Htheta->SetLineColor(kRed);
  Ntheta->SetLineColor(kCyan);
  allHe3_theta->SetLineColor(kGreen);
  // clean N2
  MTtheta->Scale(scl_N_MT);
  Ntheta->Add(MTtheta,-1);

  OrMTtheta->Scale(scl_He3_MT);
  Ntheta->Scale(scl_He3_N);
  Htheta->Add(OrMTtheta,-1);
  TH1F * cleanHe3_theta =(TH1F*)Htheta->Clone("cleanHe3_theta");
  cleanHe3_theta->SetLineColor(kMagenta);
  cleanHe3_theta->Add(Ntheta,-1);
  allHe3_theta->Draw();
  Htheta->Draw("same");
  Ntheta->Draw("same");
  OrMTtheta->Draw("same");
  cleanHe3_theta->Draw("same");
  // double scl_sim =Stheta->GetMaximum()/cleanHe3->GetMaximum();
  //Stheta->Scale(scl_sim);
  //Stheta->Draw("same");
  */
  /*************************************/
  TCanvas * c_w2 = new TCanvas("c_w2","",800,600);
  c_w2->Clear();

  TH1F * allHe3_w2 =(TH1F*)Hw2->Clone("allHe3_w2");
  OrMTw2->SetLineColor(kBlue);
  Hw2->SetLineColor(kRed);
  Nw2->SetLineColor(kCyan);
  allHe3_w2->SetLineColor(kGreen);
  Sw2->SetLineColor(kBlack);
  // clean N2
  MTw2->Scale(scl_N_MT);
  Nw2->Add(MTw2,-1);

  OrMTw2->Scale(scl_He3_MT);
  Nw2->Scale(scl_He3_N);
  Hw2->Add(OrMTw2,-1);
  TH1F * cleanHe3_w2 =(TH1F*)Hw2->Clone("cleanHe3_w2");
  cleanHe3_w2->SetLineColor(kMagenta);
  cleanHe3_w2->Add(Nw2,-1);
  allHe3_w2->Draw();
  Hw2->Draw("same");
  Nw2->Draw("same");
  OrMTw2->Draw("same");
  cleanHe3_w2->Draw("same");
  double scl_sim =cleanHe3_w2->GetMaximum()/Sw2->GetMaximum();
  Sw2->Scale(scl_sim);
  Sw2->Draw("same");
  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->AddEntry(cleanHe3_w2,"clean He3","l");
  legend->AddEntry(OrMTw2,"Glass","l");
  legend->AddEntry(Nw2,"Clean N2","l");
  legend->AddEntry(Sw2,"Simulation","l");
  legend->Draw();


}
