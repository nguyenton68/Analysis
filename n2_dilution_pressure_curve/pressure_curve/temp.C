#define temp_cxx
#include "temp.h"
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


void temp::Loop()
{
  if (fChain == 0) return;
  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("R.tr.r_ph",1);
  fChain->SetBranchStatus("R.tr.r_th",1);
  fChain->SetBranchStatus("R.tr.r_x",1);
  fChain->SetBranchStatus("R.tr.r_y",1);
  fChain->SetBranchStatus("R.tr.n",1);
  fChain->SetBranchStatus("R.ps.e",1);
  fChain->SetBranchStatus("R.sh.e",1);
  fChain->SetBranchStatus("R.gold.dp",1);
  fChain->SetBranchStatus("R.gold.y",1);
  fChain->SetBranchStatus("R.gold.th",1);
  fChain->SetBranchStatus("R.gold.ph",1);
  fChain->SetBranchStatus("R.cer.asum_c",1);
  fChain->SetBranchStatus("DR.evtypebits",1);
  fChain->SetBranchStatus("beam.H.helicity",1);
  int runnum=2087;

  TFile *fnew = new TFile(Form("%d_n2_fp_cut.root",runnum),"RECREATE");
  TTree * tree = new TTree("t1","new tree"); 
  Float_t xfocal,yfocal,phfocal,thfocal,gold_dp,rec_phi,rec_theta,rec_w_he3,rec_w_n2,rec_y,helicity;
  tree->Branch("xfocal",&xfocal,"xfocal/F");
  tree->Branch("yfocal",&yfocal,"yfocal/F");
  tree->Branch("thfocal",&thfocal,"thfocal/F");
  tree->Branch("phfocal",&phfocal,"phfocal/F");
  tree->Branch("gold_dp",&gold_dp,"gold_dp/F");
  tree->Branch("rec_phi",&rec_phi,"rec_phi/F");
  tree->Branch("rec_theta",&rec_theta,"rec_theta/F");
  tree->Branch("rec_y",&rec_y,"rec_y/F");
  tree->Branch("rec_w_he3",&rec_w_he3,"rec_w_he3/F");
  tree->Branch("rec_w_n2",&rec_w_n2,"rec_w_n2/F");
  tree->Branch("helicity",&helicity,"helicity/F");
  /*  Spectrometer central angle */
  double p0 = 2226.76;// MeV
  double ang_0 = 5.91*TMath::Pi()/180.0;
  double cos_HRS = cos(ang_0);
  double sin_HRS = sin(ang_0);
  double Mhe3 = 2.808;//mass GeV
  double Mn2 = 13.044;//mass GeV
  double Ebeam = 2.23714 ;//beam energy GeV
  double w_min,w_max;
  if(runnum==2080)
    {
      w_min=-0.0096; w_max=0.0084;
    }
  else 
    {
      w_min =-0.0093; w_max = 0.0087;
    }





  double th_max(0.048),ph_min(-0.015),ph_max(0.016),y_min(-0.006),y_max(0.006);
  /* Acceptance cut */
  const int ACC = 6;
  double y_corner[ACC],ph_corner[ACC],acc_edge[ACC];
  /* Analysis cut */
  const int N=4;
  double ana_ph[N],ana_th[N],ana_edge[N];



  int count = 0;//count number of event inside 4D cut
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  double totsh_e, deltat, y_rec, th_rec, ph_rec, hel, momentum, E_pr, ang, w_he3_rec,w_n2_rec;
  double x_foc, y_foc, th_foc, ph_foc;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);// GetEntry(i)
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    /* root variables */
    totsh_e = R_ps_e+R_sh_e;
    deltat = R_gold_dp;
    // double y_rec  = R_gold_y;
    // double th_rec = R_gold_th;
    // double ph_rec = R_gold_ph;
    hel = beam_H_helicity;
    /* Kinematic variables */
    momentum = p0*(1+deltat);
    E_pr = momentum/1000.;// convert E' to GeV
    // double ang = TMath::ACos((cos_HRS+abs(sin_HRS)*ph_rec)/(1+ph_rec*ph_rec+th_rec*th_rec));
    // double w_rec = sqrt(Mhe3**2+2*Mhe3*(Ebeam-E_pr)-4*Ebeam*E_pr*sin(ang/2.)*sin(ang/2.))-Mhe3;
    if((R_tr_n==1)&&(R_cer_asum_c>350)&&(R_ps_e/momentum>0.085)&&(totsh_e/momentum>0.77)&&((evtype&(1<<1))>0)
       ){

      x_foc =   R_tr_r_x[0];
      y_foc =   R_tr_r_y[0];
      ph_foc = R_tr_r_ph[0];
      th_foc = R_tr_r_th[0];

      y_corner[0] = -0.227484*deltat + 0.00965004;
      y_corner[1] = -0.162115*deltat - 0.017178;
      y_corner[2] =  0.203004*deltat - 0.0157547;
      y_corner[3] =  0.524234*deltat + 0.0129056;
      y_corner[4] = -0.260981*deltat + 0.024819;
      y_corner[5] = -0.394433*deltat + 0.0205985;

      ph_corner[0] = -0.120625*deltat + 0.0271181;
      ph_corner[1] = -0.0885365*deltat + 0.012184;
      ph_corner[2] = 0.209205*deltat - 0.0203365;
      ph_corner[3] = 0.677652*deltat - 0.0290798;
      ph_corner[4] = 0.105192*deltat + 0.0445469;
      ph_corner[5] = -0.0579263*deltat + 0.0470789;
      double slp, xpoint;
      for(int ii=0; ii<ACC; ii++)
        {
          if(ii==1||ii==3)
            {
              slp =(y_corner[ii+1]-y_corner[ii])/(ph_corner[ii+1]-ph_corner[ii]);
              xpoint = y_corner[ii]-slp*ph_corner[ii];
              acc_edge[ii]=slp*ph_foc+xpoint;
            }
          else
            {
              if(ii==5)
                {
                  slp=(ph_corner[0]-ph_corner[ii])/(y_corner[0]-y_corner[ii]);
                }
              else
                {
                  slp=(ph_corner[ii+1]-ph_corner[ii])/(y_corner[ii+1]-y_corner[ii]);
                }
              xpoint = ph_corner[ii] - slp*y_corner[ii];
              acc_edge[ii] = slp*y_foc + xpoint;
            }
        }
      /* Acceptance cut */
      if((ph_foc<=acc_edge[0]||ph_foc<=acc_edge[5])&&ph_foc>=acc_edge[2]&&ph_foc<=acc_edge[4]&&y_foc<=acc_edge[3]&&y_foc>=acc_edge[1])
	{
            y_rec = 0.00150805
              +0.38421*y_foc
              -0.39521*ph_foc
              +0.00557369*x_foc
              +0.513705*th_foc
              +3.64968*y_foc*ph_foc
              +7.22965*y_foc*y_foc
              -16.369*ph_foc*th_foc;

            th_rec = 0.0378725
              -0.00376131*x_foc
              +0.693092*y_foc
              -0.269014*th_foc
              +0.0412648*ph_foc;

            ph_rec = 0.00237242
              -0.00698389*x_foc
              +0.234642*y_foc
              -0.0392215*th_foc
              +0.40894*ph_foc
              -16.3558*y_foc*ph_foc
              +16.1314*y_foc*th_foc
              -1.7407*y_foc*y_foc
              +5.57496*ph_foc*ph_foc;
            ang = TMath::ACos((cos_HRS+abs(sin_HRS)*ph_rec)/(1+ph_rec*ph_rec+th_rec*th_rec));
            w_he3_rec = sqrt(Mhe3**2+2*Mhe3*(Ebeam-E_pr)-4*Ebeam*E_pr*sin(ang/2.)*sin(ang/2.))-Mhe3;
            w_n2_rec = sqrt(Mn2**2+2*Mn2*(Ebeam-E_pr)-4*Ebeam*E_pr*sin(ang/2.)*sin(ang/2.))-Mn2;
	  /* Analysis cut: angles, ytg, w*/
	    ana_ph[0] = 0.013583;//0.013065;
	    ana_ph[1] = 0.010217;//0.008647;
	    ana_ph[2] = -0.003288;//-0.005693;
	    ana_ph[3] = 0.000470;//0.0;

	    ana_th[0] = 0.044397;//0.044567;
	    ana_th[1] = 0.031207;//0.032660;
	    ana_th[2] = 0.035808;//0.037486;
	    ana_th[3] = 0.046288;//0.044567;
	  double ana_slope, ana_intercept;
	  for(int ii=0; ii<N; ii++)
	    {
	      if(ii==3)
		ana_slope=(ana_th[0]-ana_th[ii])/(ana_ph[0]-ana_ph[ii]);
	      else
		ana_slope=(ana_th[ii+1]-ana_th[ii])/(ana_ph[ii+1]-ana_ph[ii]);

	      ana_intercept = ana_th[ii]-ana_slope*ana_ph[ii];
	      ana_edge[ii]=ana_slope*ph_rec + ana_intercept;
	    }

	  //if(th_rec>=ana_edge[0]&&th_rec>=ana_edge[1]&&th_rec<=ana_edge[2]&&th_rec<=ana_th[3])
	  //&&w_rec<=w_max && w_rec>=w_min)
	  //{
	      xfocal   = x_foc;
	      yfocal   = y_foc;
	      phfocal  = ph_foc;
	      thfocal  = th_foc;

	      gold_dp  = deltat;

	      rec_phi  = ph_rec;
	      rec_theta= th_rec;
	      rec_y    = y_rec;
	      rec_w_he3= w_he3_rec;
	      rec_w_n2 = w_n2_rec;
	      helicity = hel;
	      count++;
	      tree->Fill();
	      //}//analysis cuts
	}// acceptance cut
    }//cut PID
  }//each entry

  cout<<"Number of event inside 4D cut at focal plane= "<<count<<endl;
  tree->Write();
  fnew->Write();
}
