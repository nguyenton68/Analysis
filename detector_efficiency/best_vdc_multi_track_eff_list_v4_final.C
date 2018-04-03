//  firstly created by Jie Liu, to study multi-track for g2p,   June 2013
//  Applied to Small GDH (E97-110)  Jie Liu, Oct 2013
//  Applied to SCR (E08-015) Jie Liu, Apr 2014
//  any question, contact Jie Liu <jie@jlab.org>
//  version 1.0  Jan 2013 //the same as origin
//  version 2.0  Jan 2014 //updated
//  version 3.0  Aug 2014 //combined left arm and right arm in one script.
//  version 4.0  Jan 2015 // select tg/fp acceptance
// final version: for 2track case more
//  *****note****  just change the part in "best_vdc_multi_track_eff" subroutine
// Nguyen
// - 04/2018: change fabs --> abs
// - Number of argument in vdc_multi_track_eff

#include <iostream>
#include <fstream>
#include <TFile.h> 
#include <TTree.h>
#include <stdlib.h>
#include <TROOT.h>
#include "TMath.h"
#include <cstring>
#include <TChain.h>
#include <TDatime.h>
#include <TFormula.h>
#include </work/halla/gdh/disk1/ton/analyzer-1.4.0-cs65/src/THaAnalysisObject.h>
#include </work/halla/gdh/disk1/ton/analyzer-1.4.0-cs65/include/lib_mysql_support_for_sagdh.h>

using namespace std;
void getdb(Char_t* detname, Float_t* detpara, TDatime date);
void vdc_multi_track_eff(Int_t runnum, Int_t eng, TString path, TDatime date, TString filename, TFormula *acceptance, TString coords); 
void best_vdc_multi_track_eff_list_v4_final(){
//gSystem->Load("/work/halla/gdh/disk1/analyzer-1.4.0-cs65/include/lib_mysql_support_for_sagdh.h");
  //fllows need to be modified to match your experiment
  //**********************************

  Int_t energy=1149; //beam energy setting for run list 
  TString path = "./"; //rootfile directory to analyze  TDatime *date = new TDatime(2011,4,1,0,0,0);
  TDatime date("2003-05-01 00:00:00"); //date for database to use, you need check the lead glass geometry respect to VDC in database
  TString filename = "./output";// where to save the result
  TString runlist = "list.txt";  // a list of runnumber to analyze

  TString coords[2]={"tg", "fp"};
  TFormula *acceptance = new TFormula ("acceptance", "abs([0])<0.6 && abs([1])<0.06 && abs([2])<0.12 && abs([3])<0.05"); //acceptance cut in focal plane [0],[1], [2],[3] are for x, y, th, ph


 //TFormula *acceptance = new TFormula ("acceptance", "fabs([0])<0.05 && fabs([1])<0.05 && fabs([2])<0.05 && fabs([3])<0.05");  //acceptance cut in target plane [0],[1], [2],[3] are for dp, y, th, ph


  // just need to change the part above
  //**********************************



  Int_t runnum[3000],i,count;
  ifstream file_list(runlist.Data(),ios::in);
  if(!file_list){
    cout << "Whoops, can't find the file list!" << endl;
    cout << "Exiting the program." << endl;
    exit(1);
   }

  else {
    i=0;
    file_list >> runnum[i];

    while(!file_list.eof()){
      cout <<" i= "<<i<< "   runnum = " << runnum[i]<< endl;
    i++;
    file_list >> runnum[i];
    }
    file_list.close();
    count=i;
    for (Int_t j=0;j<count;j++)
      {
	cout<<" count "<<count<<"  now run :"<<runnum[j]<<endl;
	vdc_multi_track_eff(runnum[j], energy, path, date, filename, acceptance, coords[1]);
      }
  }
}

void vdc_multi_track_eff(Int_t runnum, Int_t eng, TString path, TDatime date, TString filename, TFormula *acceptance, TString coords) 
{

  Int_t i, j, k, l, m, n, s, t, u, kk, ll, nsplit, column[2];
  Double_t momentum,cer_cut,prl1_cut,prl_total_cut;
  Double_t L_tr_x_array[100],L_tr_y_array[100],L_tr_th_array[100],L_tr_ph_array[100],L_prl1_a_c[100],L_prl2_a_c[100],L_tg_dp_array[100],L_tg_y_array[100],L_tg_th_array[100],L_tg_ph_array[100];
  Double_t *tr_x=L_tr_x_array,*tr_y=L_tr_y_array,*tr_th=L_tr_th_array,*tr_ph=L_tr_ph_array,*prl1_a_c=L_prl1_a_c,*prl2_a_c=L_prl2_a_c,*tg_dp=L_tg_dp_array,*tg_y=L_tg_y_array,*tg_th=L_tg_th_array,*tg_ph=L_tg_ph_array;

  Double_t L_tr_n, prl_cut, prl_cut_total, prl1_e_main, prl2_e_main, L_cer_asum_c,track_tg_dp[7],track_tg_y[7],track_tg_th[7],track_tg_ph[7],track_vdc_x[7],track_vdc_y[7],track_vdc_th[7],track_vdc_ph[7], track_prl1_x[7],track_prl1_y[7],track_prl2_x[7],track_prl2_y[7], prl1_x[7],prl1_y[7],prl2_x[7],prl2_y[7],prl1_e[7], prl2_e[7],prl1_e_bak[7], prl2_e_bak[7];

  Bool_t pidgood[7], pidaccgood[7], pidaccgood2[7], pidaccgood3[7];
  Int_t dxblock[7][7], dx2block[7][7], dxnum[3];
  Bool_t dblock[7][7][5];
  Bool_t pidaccgoodtmp[7][7], pidaccgoodbak[7][7],  pidaccgoodtmp3[7][7][7],  pidaccgoodbak3[7][7][7];
  Double_t prl1_x_main,prl2_x_main,prl1_y_main,prl2_y_main, prl1_trx_main,prl2_trx_main,prl1_try_main,prl2_try_main;
  // mark1[7],mark2[7],mark3[7],mark4[7],prl1_x_main,prl2_x_main,prl1_y_main,prl2_y_main,prl1_e_main,prl2_e_main,sxe_prl1,sxe_prl2,sye_prl1,sye_prl2,sxe_prl1_bak,sxe_prl2_bak,sye_prl1_bak,sye_prl2_bak,E_track[7],E_track_bak[7],r_prl1,r_prl2,d_x[7],d_y[7],d_r[7],d_x_d[7],d_y_d[7],d_r_d[7],aa,bb,cc,prl1_trx_main,prl2_trx_main,prl1_try_main,prl2_try_main;

  Int_t x_center_prl1,y_center_prl1; 
  Int_t x_center_prl2,y_center_prl2; 
  Double_t E_track[7],E_track_bak[7], track[9][4], track_tmp[9][20];
  Int_t event_total, sample_total;
  Double_t sample_pid_total;
  Double_t flag_acc_total, count_tmp;
  bool flag_acc[7], mark[9][4];
  bool flag_tmp[12];
  Char_t detname[2][100]={"L.prl1","L.prl2"};
  Float_t prl1Para[10], prl2Para[10];
  Double_t vdc_eff, vdc_eff_low, vdc_eff_up;

  momentum = get_my_rspec_momentum(runnum); //momentum
  cer_cut = get_my_pid_cut(runnum,"cs","e","cer","lo"); //location of L.cer.asum_c cut
  Bool_t flag = my_result_is_null();
  prl_cut= get_my_pid_cut(runnum,"cs","e","ps","lo") * momentum;   // cut on E/p
  prl_total_cut = get_my_pid_cut(runnum,"cs","e","ep","lo") * momentum;   // cut on preshower
  cout << "momentum = " << momentum<< "runnum = " << runnum << "  cer cut = " << cer_cut<< "  prl1 cut = " << prl_cut<< "  prl1 total cut = " << prl_total_cut<< endl;
  nsplit=2;
  TChain *T = new TChain("T");
  for(i=0;i<=nsplit;i++)
   {
      cout << "i=" << i << endl;
      if(i==0){
	TString input_filename = path + "gdh_";
	input_filename += runnum;
	input_filename += ".root";
	cout << "Input File = " << input_filename << endl;          
	T->Add(input_filename);
	if(T->GetEntries()==0) { 
	  cout << "Whoops, can't find the file list!" << endl;
	  cout << "Exiting the program." << endl;
	  exit(1);
	}
      }
      if(i!=0){
	TString input_filename = path + "gdh_";
	input_filename += runnum;
	input_filename += "_";
	input_filename += i;
	input_filename += ".root";
	cout << "Input File = " << input_filename << endl;
	T->Add(input_filename);
      }
      
   } // end of for loop  
  if(runnum > 20000)
    {
      T->SetBranchAddress("L.tr.n",&L_tr_n);
      T->SetBranchAddress("L.tr.x",tr_x);
      T->SetBranchAddress("L.tr.y",tr_y);   
      T->SetBranchAddress("L.tr.th",tr_th);
      T->SetBranchAddress("L.tr.ph",tr_ph);
      T->SetBranchAddress("L.prl1.e",&prl1_e_main);
      T->SetBranchAddress("L.prl2.e",&prl2_e_main);
      T->SetBranchAddress("L.cer.asum_c",&L_cer_asum_c);
      T->SetBranchAddress("L.prl1.x",&prl1_x_main);
      T->SetBranchAddress("L.prl2.x",&prl2_x_main);
      T->SetBranchAddress("L.prl1.y",&prl1_y_main);
      T->SetBranchAddress("L.prl2.y",&prl2_y_main);
      T->SetBranchAddress("L.prl1.trx",&prl1_trx_main);
      T->SetBranchAddress("L.prl1.try",&prl1_try_main);
      T->SetBranchAddress("L.prl2.trx",&prl2_trx_main);
      T->SetBranchAddress("L.prl2.try",&prl2_try_main);
      T->SetBranchAddress("L.prl1.a_c",prl1_a_c);
      T->SetBranchAddress("L.prl2.a_c",prl2_a_c);
      T->SetBranchAddress("L.tr.tg_dp",tg_dp);
      T->SetBranchAddress("L.tr.tg_y",tg_y);   
      T->SetBranchAddress("L.tr.tg_th",tg_th);
      T->SetBranchAddress("L.tr.tg_ph",tg_ph);
    }
  else
    {
      T->SetBranchAddress("R.tr.n",&L_tr_n);
      T->SetBranchAddress("R.tr.x",tr_x);
      T->SetBranchAddress("R.tr.y",tr_y);   
      T->SetBranchAddress("R.tr.th",tr_th);
      T->SetBranchAddress("R.tr.ph",tr_ph);
      T->SetBranchAddress("R.ps.e",&prl1_e_main);
      T->SetBranchAddress("R.sh.e",&prl2_e_main);
      T->SetBranchAddress("R.cer.asum_c",&L_cer_asum_c);
      T->SetBranchAddress("R.ps.x",&prl1_x_main);
      T->SetBranchAddress("R.sh.x",&prl2_x_main);
      T->SetBranchAddress("R.ps.y",&prl1_y_main);
      T->SetBranchAddress("R.sh.y",&prl2_y_main);
      T->SetBranchAddress("R.ps.trx",&prl1_trx_main);
      T->SetBranchAddress("R.sh.trx",&prl2_trx_main);
      T->SetBranchAddress("R.ps.try",&prl1_try_main);
      T->SetBranchAddress("R.sh.try",&prl2_try_main);
      T->SetBranchAddress("R.ps.a_c",prl1_a_c);
      T->SetBranchAddress("R.sh.a_c",prl2_a_c);
      T->SetBranchAddress("R.tr.tg_dp",tg_dp);
      T->SetBranchAddress("R.tr.tg_y",tg_y);   
      T->SetBranchAddress("R.tr.tg_th",tg_th);
      T->SetBranchAddress("R.tr.tg_ph",tg_ph);
      strcpy(detname[0], "R.ps");
      strcpy(detname[1], "R.sh");
    }
  //read database
  getdb(detname[0], prl1Para, date);
  getdb(detname[1], prl2Para, date);

  //initialize 
  sample_total = 0;
  event_total = 0;
  sample_pid_total = 0;
  for(i=0; i<9;i++)
    {
      for(j=0;j<4;j++)
	{
	  track[i][j]=0;
	}
      for(j=0;j<20;j++)
	{
	  track_tmp[i][j]=0;
	}
    }

  //loop for events
  Int_t nevent = T->GetEntries();
  cout << "number of entries = " << nevent << endl;
  Int_t event_want = 20000;
  if(nevent < event_want)  event_want =  nevent;
  //for (k = 0;k < event_want; k++) 
  for (k = 0;k < nevent; k++) 
    {
    //cout << i << "\r";
      T->GetEntry(k);
    if(k%10000 == 0) cout<<"event  "<<k<<endl;
    if((L_cer_asum_c > cer_cut)&&(prl1_e_main > prl_cut)&&((prl1_e_main + prl2_e_main) > prl_total_cut))
      {
	sample_pid_total++;
	if(L_tr_n==0) 
	  {
	  track[0][0]++;
	  track_tmp[0][0]++;
	  }
	else if (L_tr_n >7)
	  {
	  track[8][0]++;
	  sample_total++;
	  }
	else
	  {
	    flag_acc_total = 0;
	    for(m = 0; m < L_tr_n; m++)
	      {
		flag_acc[m] = false;
		for(i = 0;i < 4;i++)
		  {
		    mark[m][i] = false;//mark for track not locate on shower
		  }
		track_vdc_x[m] = *(tr_x+m);
		track_vdc_y[m] = *(tr_y+m);
		track_vdc_th[m] = *(tr_th+m);
		track_vdc_ph[m] = *(tr_ph+m);
		track_tg_dp[m] = *(tg_dp+m);
		track_tg_y[m] = *(tg_y+m);
		track_tg_th[m] = *(tg_th+m);
		track_tg_ph[m] = *(tg_ph+m);

		acceptance->SetParameters(track_vdc_x[m], track_vdc_y[m],  track_vdc_th[m], track_vdc_ph[m]);
		if(coords == "tg")	acceptance->SetParameters(track_tg_dp[m], track_tg_y[m],  track_tg_th[m], track_tg_ph[m]);
		if(acceptance->Eval(1))
		  {
		    flag_acc_total++;
		    flag_acc[m] = true;
		  }
	      }
	    if(flag_acc_total > 0)
	      {
		sample_total++;
		track[int(L_tr_n)][0]++;
		for(i = 0; i < prl1Para[0]*prl1Para[1]; i++)
		  {
		    L_prl1_a_c[i] = *(prl1_a_c+i);
		  }
		for(i = 0; i < prl2Para[0]*prl2Para[1]; i++)
		  {
		    L_prl2_a_c[i] = *(prl2_a_c+i);
		  }
		for(m = 0; m < L_tr_n; m++)
		  {
		    track_prl1_x[m] = track_vdc_x[m] + track_vdc_th[m]*prl1Para[2];
		    track_prl1_y[m] = track_vdc_y[m] + track_vdc_ph[m]*prl1Para[2];
		    track_prl2_x[m] = track_vdc_x[m] + track_vdc_th[m]*prl2Para[2];
		    track_prl2_y[m] = track_vdc_y[m] + track_vdc_ph[m]*prl2Para[2];
		    x_center_prl1 = int((track_prl1_x[m] - prl1Para[3] + prl1Para[5]/2.0)/ prl1Para[5]);
		    if((track_prl1_x[m] < prl1Para[3] - prl1Para[5]/2.0 + prl1Para[1]*prl1Para[5])^(prl1Para[5]>0))
		      {
			mark[m][0] = true;
		      }
		    if((track_prl1_x[m] > prl1Para[3] - prl1Para[5]/2.0)^(prl1Para[5]>0)) 
		      {
			mark[m][0] = true;
		      }

		    y_center_prl1 = int((track_prl1_y[m] - prl1Para[4] + prl1Para[6]/2.0)/ prl1Para[6]);
		    if((track_prl1_y[m] < prl1Para[4] - prl1Para[6]/2.0 + prl1Para[0]*prl1Para[6])^(prl1Para[6]>0))
		      {
			mark[m][1] = true;
		      }
		    if((track_prl1_y[m] > prl1Para[4] - prl1Para[6]/2.0)^(prl1Para[6]>0)) 
		      {
			mark[m][1] = true;
		      }


		    x_center_prl2 = int((track_prl2_x[m] - prl2Para[3] + prl2Para[5]/2.0)/ prl2Para[5]);
		    if((track_prl2_x[m] < prl2Para[3] - prl2Para[5]/2.0 + prl2Para[1]*prl2Para[5])^(prl2Para[5]>0))
		      {
			mark[m][2] = true;
		      }
		    if((track_prl2_x[m] > prl2Para[3] - prl2Para[5]/2.0)^(prl2Para[5]>0)) 
		      {
			mark[m][2] = true;
		      }

		    y_center_prl2 = int((track_prl2_y[m] - prl2Para[4] + prl2Para[6]/2.0)/ prl2Para[6]);
		    if((track_prl2_y[m] < prl2Para[4] - prl2Para[6]/2.0 + prl2Para[0]*prl2Para[6])^(prl2Para[6]>0))
		      {
			mark[m][3] = true;
		      }
		    if((track_prl2_y[m] > prl2Para[4] - prl2Para[6]/2.0)^(prl2Para[6]>0)) 
		      {
			mark[m][3] = true;
		      }

		    prl1_e[m] = 0;
		    prl2_e[m] = 0;
		    column[0]=prl1Para[1];
		    column[1]=prl2Para[1];
		    for(l = 0;l < prl1Para[0]*prl1Para[1];l++)
		      {
			ll = l / column[0];
			kk = l % column[0];
			if(kk < (x_center_prl1 + 2)&&kk > (x_center_prl1 - 2)&& ll < (y_center_prl1 + 2) && ll > (y_center_prl1 - 2)&&L_prl1_a_c[l]>0)
			  {
			    prl1_e[m] = prl1_e[m] + L_prl1_a_c[l];
			  }
		      }

		    for(l = 0;l < prl2Para[0]*prl2Para[1];l++)
		      {
			ll = l / column[1];
			kk = l % column[1];
			if(kk < (x_center_prl2 + 2) && kk >(x_center_prl2 - 2)&& ll <(y_center_prl2 + 2) && ll > (y_center_prl2 - 2)&&L_prl2_a_c[l] > 0)
			  {
			  prl2_e[m] = prl2_e[m] + L_prl2_a_c[l]; 
			}

		      }

		    if(x_center_prl1 % column[0] == prl1Para[1] - 1)
		      {
			prl1_e_bak[m] = 0;
			for(l = 0;l < prl1Para[0]*prl1Para[1];l++){
			  ll = l/ column[0];
			  kk = l% column[0];
			  if(kk < (x_center_prl1 + 1) && kk > (x_center_prl1 - 2) && ll < (y_center_prl1 + 2)&& ll>(y_center_prl1 - 2)&&L_prl1_a_c[l] > 0)
			    {

			      prl1_e_bak[m] = prl1_e_bak[m] + L_prl1_a_c[l];
			    }
			}
			prl1_e[m] = prl1_e_bak[m];
		      }
		    if(x_center_prl1 % column[0] == prl1Para[1])
		      {
			prl1_e_bak[m] = 0;
			for(l = 0;l < prl1Para[0]*prl1Para[1];l++){
			  ll = l/ column[0];
			  kk = l% column[0];
			  if(kk < (x_center_prl2 + 2)&& kk > (x_center_prl2 - 1)&& ll < (y_center_prl2 + 2)&&ll > (y_center_prl2 - 2)&&L_prl2_a_c[l] > 0)
		 
			    {
			      prl1_e_bak[m] = prl1_e_bak[m] + L_prl1_a_c[l];
			    }
			}
			prl1_e[m] = prl1_e_bak[m];
		      }

		    if(x_center_prl2 % column[1] == prl2Para[1]-1)
		      {
			prl2_e_bak[m] = 0;
			for(l = 0;l < prl2Para[0]*prl2Para[1];l++){
			  ll =l/ column[1];
			  kk =l% column[1];
			  if(kk < (x_center_prl2 + 1)&&kk > (x_center_prl2 - 2)&& ll < (y_center_prl2+2) && ll > (y_center_prl2 - 2)&&L_prl2_a_c[l] > 0)
			    {
			      prl2_e_bak[m] = prl2_e_bak[m] + L_prl2_a_c[l];
			    }
			}
			prl2_e[m] = prl2_e_bak[m];
		      }
		    if(x_center_prl2 % column[1] == prl2Para[1])
		      {
			prl2_e_bak[m] = 0;
			for(l = 0;l < prl2Para[0]*prl2Para[1];l++){
			  ll = l/ column[1];
			  kk = l% column[1];
			  if(kk < (x_center_prl2 + 2) && kk >(x_center_prl2 - 1)&&ll<(y_center_prl2 + 2)&& ll > (y_center_prl2 - 2) && L_prl2_a_c[l] > 0)
		 
			    {
			      prl2_e_bak[m] = prl2_e_bak[m] + L_prl2_a_c[l];
			    }
			}
			prl2_e[m] = prl2_e_bak[m];
		      }

		    if(mark[m][0]||mark[m][1])
		      {
			prl1_e[m]=0;
		      }

		    if(mark[m][2]||mark[m][3])
		      {
			prl2_e[m]=0;
		      }

		    track_prl1_x[m] = int((track_prl1_x[m] - prl1Para[3] + prl1Para[5]/2.0)/ prl1Para[5]);
		    track_prl2_x[m] = int((track_prl2_x[m] - prl2Para[3] + prl2Para[5]/2.0)/ prl2Para[5]);

		    track_prl1_y[m] = int((track_prl1_y[m] - prl1Para[4] + prl1Para[6]/2.0)/ prl1Para[6]);
		    track_prl2_y[m] = int((track_prl2_y[m] - prl2Para[4] + prl2Para[6]/2.0)/ prl2Para[6]);
		    prl1_e_bak[m] = 0;
		    prl2_e_bak[m] = 0;
	
                    for( i =0; i< 3; i++)
		      {
			j=track_prl1_x[m]+ (i-1) * prl1Para[1] + track_prl1_y[m]* prl1Para[1] ;
			n=track_prl2_x[m]+ (i-1) * prl2Para[1] + track_prl2_y[m]* prl2Para[1] ;
			if(j>-1 && j < prl1Para[0]*prl1Para[1] && !mark[m][0] && !mark[m][1])
			  {
			    prl1_e_bak[m]=prl1_e_bak[m] + L_prl1_a_c[j];
			  }
			if(n>-1 && n < prl2Para[0]*prl2Para[1] && !mark[m][2] && !mark[m][3])
			  {
			    prl2_e_bak[m]=prl2_e_bak[m] + L_prl2_a_c[n];
			  }
		      }
		    // cout<<"  *******m  "<<m<<"   prl1 e  "<<prl1_e[m]<<"   prl2  e "<<prl2_e[m]<<endl;
		    E_track[m] = prl1_e[m] + prl2_e[m];
		    E_track_bak[m] = prl1_e_bak[m] + prl2_e_bak[m];
		    pidgood[m] = E_track[m] > prl_total_cut && prl1_e[m] > prl_cut;
		    pidaccgood[m] = E_track[m] > prl_total_cut && prl1_e[m] > prl_cut && flag_acc[m];
		    pidaccgood2[m] = E_track[m] > 2.0 * prl_total_cut && prl1_e[m] > 2.0 * prl_cut && flag_acc[m];
		    pidaccgood3[m] = E_track[m] > 3.0 * prl_total_cut && prl1_e[m] > 3.0 * prl_cut && flag_acc[m];


		  }

		for(i = 0;i < 3;i++)
		  {
		    dxnum[i] = 0; 
		  }
		for(i = 0;i < L_tr_n; i++)
		  {
		    for(j = 0;j < L_tr_n; j++)
		      {
			if(i != j)
			  {
			    dxblock[i][j] = TMath::Abs(track_prl1_x[i]-track_prl1_x[j]);
			    dx2block[i][j] = TMath::Abs(track_prl2_x[i]-track_prl2_x[j]);
			    dblock[i][j][0] = dxblock[i][j] == 0 && (flag_acc[i]||flag_acc[j]);
			    dblock[i][j][1] = dxblock[i][j] == 1 && (flag_acc[i]||flag_acc[j]);
			    dblock[i][j][2] = dxblock[i][j] >0 && (flag_acc[i]||flag_acc[j]);
			    dblock[i][j][3] = dxblock[i][j] >1 && (flag_acc[i]||flag_acc[j]);
			    dblock[i][j][4] = dxblock[i][j] <2 && (flag_acc[i]||flag_acc[j]);

			    if( dxblock[i][j] >0 )  dxnum[0]++;
			    if( dxblock[i][j] >1 )  dxnum[1]++;
			    if( dxblock[i][j] <2 )  dxnum[2]++;

			    pidaccgoodtmp[i][j] = (E_track[i] - E_track_bak[j]>prl_total_cut + E_track_bak[j])&&(prl1_e[i]-prl1_e_bak[j]>prl_cut+prl1_e_bak[j])&&flag_acc[i];
			    if(dx2block[i][j] > 1)
			      {
				pidaccgoodtmp[i][j] = (E_track[i]- prl1_e_bak[j]>prl_total_cut + prl1_e_bak[j])&&(prl1_e[i]-prl1_e_bak[j]>prl_cut+prl1_e_bak[j])&&flag_acc[i];
			      }
			    pidaccgoodbak[i][j] = (E_track_bak[i] > prl_total_cut + E_track_bak[j])&&(prl1_e_bak[i] > prl_cut + prl1_e_bak[j])&&flag_acc[i];
			  }

		      }
		  }

		for(i = 0;i < L_tr_n; i++)
		  {
		    for(j = 0;j < L_tr_n; j++)
		      {
			for(l = 0;l < L_tr_n; l++)
			  {
			    if(i != j && i != l && l != j)
			      {
				//only be used for 3 track events, pidaccgoodtmp3[i][j][l] losse the cuts, compare with pidaccgoodtmp[i][j]
				flag_tmp[0] =  (E_track[i] - E_track_bak[j] > prl_total_cut + E_track_bak[j])&&(prl1_e[i] - prl1_e_bak[j] > prl_cut + prl1_e_bak[j]);
				flag_tmp[1] =  (E_track[i] - E_track_bak[l] > prl_total_cut + E_track_bak[l])&&(prl1_e[i] - prl1_e_bak[l] > prl_cut + prl1_e_bak[l]);
				pidaccgoodtmp3[i][j][l] = flag_tmp[0] && flag_tmp[1] && flag_acc[i];  
				flag_tmp[0] = (E_track_bak[i] > prl_total_cut + E_track_bak[j])&&(prl1_e_bak[i] > prl_cut + prl1_e_bak[j]);
				flag_tmp[1] = (E_track_bak[i] > prl_total_cut + E_track_bak[l])&&(prl1_e_bak[i] > prl_cut + prl1_e_bak[l]);
				pidaccgoodbak3[i][j][l] = flag_tmp[0] && flag_tmp[1] && flag_acc[i];  
			      }
			  }

		      }
		  }
		if(L_tr_n == 1 && pidaccgood[0])
		  {
		    track[1][1]++;
		    track_tmp[1][0]++;
		  }
		if(L_tr_n==2)
		  {

		    if(!pidgood[0] && !pidgood[1])
		      {
			track_tmp[2][0]++;
		      }
		    if(pidaccgood[0] && !pidgood[1])
		      {
			track_tmp[2][1]++;
			track[2][1]++;
		      }
		    if(pidaccgood[1] && !pidgood[0])
		      {
			track_tmp[2][2]++;
			track[2][1]++;
		      }
		    if(pidgood[0] && pidgood[1] && dblock[1][0][4])
		      {
		
			track_tmp[2][3]++;
			track[2][1]++;
		      }
		    if((pidaccgoodtmp[0][1] || pidaccgoodtmp[1][0] || pidaccgoodbak[0][1] || pidaccgoodbak[1][0] || pidaccgood2[0] || pidaccgood2[1]) && pidgood[0] && pidgood[1] && dblock[1][0][0])
		      {
			track_tmp[2][4]++;
		      } 
		    if((pidaccgoodtmp[0][1] || pidaccgoodtmp[1][0] || pidaccgoodbak[0][1] || pidaccgoodbak[1][0] || pidaccgood2[0] || pidaccgood2[1]) && pidgood[0] && pidgood[1] && dblock[1][0][1])
		      {
			track_tmp[2][5]++;
		      } 
		    if(pidgood[0] && pidgood[1] && dblock[1][0][3])
		      {
			track_tmp[2][6]++;
			track[2][1]++;
		      } 
		    if((pidaccgood2[0] || pidaccgood2[1]) && pidgood[0] && pidgood[1] && dblock[1][0][4])
		      {
			track_tmp[2][7]++;
		      } 
		    if((pidaccgood2[0] || pidaccgood2[1]) && pidgood[0] && pidgood[1])
		      {
			track_tmp[2][8]++;
		      } 
		    if((pidaccgoodtmp[0][1] || pidaccgoodtmp[1][0] || pidaccgoodbak[0][1] || pidaccgoodbak[1][0]) && pidgood[0] && pidgood[1] && dblock[1][0][4])
		      {
			track_tmp[2][9]++;
		      } 

		  }

		if(L_tr_n==3)
		  {
		    if(!pidgood[0] && !pidgood[1] && !pidgood[2] )
		      {
			track_tmp[3][0]++;
		      }

		    for(i = 0;i < L_tr_n; i++)
		      {
			for(j = 0;j < L_tr_n; j++)
			  {
			    for(l = 0;l < j; l++) 
			      {
				if(i != j && i != l && l != j)
				  {
				    if(pidaccgood[i] && !pidgood[j] && !pidgood[l] )
				      {
					track_tmp[3][i+1]++;
					track[3][1]++;
				      }	
				  }		      
			      }
			  }
		      }
		    for(i = 0; i<12 ;i++)
		      {
			flag_tmp[i] = false;
		      }
 
		    for(i = 0;i <L_tr_n; i++)
		      {
			for(j = 0;j <L_tr_n; j++)
			  {
			    for(l = 0;l <L_tr_n; l++)
			      {
				if(i != j && i != l && l != j)
				  {
				    if(pidgood[i] && pidgood[j] && ( flag_acc[i] || flag_acc[j]) && !pidgood[l])
				      {
					//2E0
					flag_tmp[0] = true;
				      }
				    if(pidgood[i] && pidgood[j] && dblock[i][j][3] && !pidgood[l])
				      {
					flag_tmp[1] = true;
				      }
				    if((pidaccgoodtmp[i][j] || pidaccgoodbak[i][j] || pidaccgood2[i] || pidaccgood2[j]) && pidgood[i] && pidgood[j] && !pidgood[l] && dblock[i][j][0])
				      {
					flag_tmp[2] = true;
				      }

				    if((pidaccgoodtmp[i][j] || pidaccgoodbak[i][j] || pidaccgood2[i] || pidaccgood2[j]) && pidgood[i] && pidgood[j] && !pidgood[l] && dblock[i][j][1])
				      {
					flag_tmp[3] = true;
				      }
				    if(pidgood[i] && pidgood[j] && pidgood[l])
				      {
					//3E0
					flag_tmp[4] = true;
				      }
				    if(pidaccgood[i] && pidgood[j] && pidgood[l] && dblock[i][j][3] && dblock[i][l][3])
				      {
					flag_tmp[5] = true;
				      }
				    if((pidaccgoodtmp3[i][j][l] || pidaccgoodbak3[i][j][l] || pidaccgood3[i] || pidaccgood3[j] || pidaccgood3[l]) && pidgood[i] && pidgood[j] && pidgood[l] && dxnum[0] == 0)
				      {
					flag_tmp[6] = true;
				      }
				    if((pidaccgoodtmp3[i][j][l] || pidaccgoodbak3[i][j][l] || pidaccgood3[i] || pidaccgood3[j] || pidaccgood3[l]) && pidgood[i] && pidgood[j] && pidgood[l] && dxnum[0]/2.0 == 2  && dxnum[2]/2.0 == 3)
				      {
					flag_tmp[7] = true;
				      }
				  }

			      }

			  }
		      }
		    if(flag_tmp[0])
		      {
			track_tmp[3][4]++;
			track[3][1]++;
		      }
		    if(flag_tmp[1])
		      {
			track_tmp[3][5]++;
		      }
		    if(flag_tmp[2])
		      {
			track_tmp[3][6]++;
		      }
		    if(flag_tmp[3])
		      {
			track_tmp[3][7]++;
		      }
		    if(flag_tmp[4])
		      {
			track_tmp[3][8]++;
			track[3][1]++;
		      } 
		    if(flag_tmp[5])
		      {
			track_tmp[3][9]++;
		      }
		    if(flag_tmp[6])
		      {
			track_tmp[3][10]++;
		      }

		    if(flag_tmp[7])
		      {
			track_tmp[3][11]++;
		      }
		  }

		if(L_tr_n==4)
		  {
		    if(!pidgood[0] && !pidgood[1] && !pidgood[2] && !pidgood[3] )
		      {
			track_tmp[4][0]++;
		      }

		    for(i = 0;i < L_tr_n; i++)
		      {
			for(j = 0;j < L_tr_n;j++)
			  {
			    for(l = 0;l < j;l++)
			      {
				for(n = 0;n < l;n++)
				  {
				    if(i != j && i != l && i != n)
				      {
					if(pidaccgood[i] && !pidgood[j] && !pidgood[l] && !pidgood[n] )
					  {
					    track_tmp[4][i+1]++;
					    track[4][1]++;
					  }
				      }
				  }
			      }
			  }
		      }
		    for(i = 0; i<12 ;i++)
		      {
			flag_tmp[i] = false;
		      }
		    for(i = 0;i <L_tr_n; i++)
		      {
			for(j = 0;j <L_tr_n; j++)
			  {
			    for(l = 0;l <L_tr_n; l++)
			      {
				for(n = 0;n <L_tr_n; n++)
				  {
				    if(i != j && i != l && i != n && j != l && j != n && l != n)
				      {
					if(pidgood[i] && pidgood[j] && ( flag_acc[i] || flag_acc[j]) && !pidgood[l] && !pidgood[n])
					  {
					    //2E
					    flag_tmp[0] = true;
					  }
					if(pidgood[i] && pidgood[j] && dblock[i][j][3] && !pidgood[l] && !pidgood[n])
					  {
					    flag_tmp[1] = true;
					  }
				
					if((pidaccgoodtmp[i][j]||pidaccgoodbak[i][j]||pidaccgood2[i]||pidaccgood2[j])&&pidgood[i]&&pidgood[j]&&!pidgood[l]&&!pidgood[n]&&dblock[i][j][0])
					  {
					    flag_tmp[2] = true;
					  }
					if((pidaccgoodtmp[i][j]||pidaccgoodbak[i][j]||pidaccgood2[i]||pidaccgood2[j])&& pidgood[i] && pidgood[j] && !pidgood[l] && !pidgood[n] && dblock[i][j][1])
					  { //cut stronger
					    flag_tmp[3] = true;
					  }
					if(pidgood[i] && pidgood[j] && pidgood[l] && !pidgood[n] && (flag_acc[i] || flag_acc[j] || flag_acc[l]))
					  {
					    //3E
					    flag_tmp[4] = true;
					  }
			
					if(pidaccgood[i] && pidgood[j] && pidgood[l] && !pidgood[n] && (dblock[i][j][3] && dblock[i][l][3]))
					  {
					    flag_tmp[5] = true;
					  }
					if(pidgood[i] && pidgood[j] && pidgood[l] && pidgood[n])
					  {
					    //4E
					    flag_tmp[6] = true;
					  }
					if(pidaccgood[i] && pidgood[j] && pidgood[l] && pidgood[n] && dblock[i][j][3] && dblock[i][l][3] && dblock[i][n][3])
					  {  
					    flag_tmp[7] = true;
					  }
					
				      }
				  }
			      }
			  }
		      }
		    if(flag_tmp[0])
		      {
			track_tmp[4][5]++;
			track[4][1]++;
		      }
		    if(flag_tmp[1])
		      {
			track_tmp[4][6]++;
		      }
		    if(flag_tmp[2])
		      {
			track_tmp[4][7]++;
		      }
		    if(flag_tmp[3])
		      {
			track_tmp[4][8]++;
		      }
		    if(flag_tmp[4])
		      {
			track_tmp[4][9]++;
			track[4][1]++;
		      } 
		    if(flag_tmp[5])
		      {
			track_tmp[4][10]++;
		      }
		    if(flag_tmp[6])
		      {
			track[4][1]++;
			track_tmp[4][11]++;
		      }
		    if(flag_tmp[7])
		      {
			track_tmp[4][12]++;
		      }
		  }

		if(L_tr_n==5)
		  {
		    if(!pidgood[0] && !pidgood[1] && !pidgood[2] && !pidgood[3] && !pidgood[4] )
		      {
			track_tmp[5][0]++;
		      }
		    for(i = 0;i < L_tr_n; i++)
		      {
			for(j = 0;j < L_tr_n;j++)
			  {
			    for(l = 0;l < j;l++)
			      {
				for(n = 0;n < l;n++)
				  {
				    for(s = 0;s < n;s++)
				      {
					if(i != j && i != l && i != n  && i != s)
					  {
					    if(pidaccgood[i] && !pidgood[j] && !pidgood[l] && !pidgood[n] && !pidgood[s] )
					      {
						track_tmp[5][i+1]++;
						track[5][1]++;
					      }
					  }
				      }
				  }
			      }
			  }
		      }
		    for(i = 0; i<12 ;i++)
		      {
			flag_tmp[i] = false;
		      }
		    for(i = 0;i <L_tr_n; i++)
		      {
			for(j = 0;j <L_tr_n; j++)
			  {
			    for(l = 0;l <L_tr_n; l++)
			      {
				for(n = 0;n <L_tr_n; n++)
				  {
				    for(s = 0;s <L_tr_n; s++)
				      {
					if(i != j && i != l && i != n  && i != s && j != l && j != n && j != s && l != n  && l != s && n != s)
					  {
					    if(pidgood[i] && pidgood[j] && (flag_acc[i] || flag_acc[j]) && !pidgood[l] && !pidgood[n] && !pidgood[s])
					      {
						//2E
						flag_tmp[0] = true;
					      }
					    if(pidgood[i] && pidgood[j] && dblock[i][j][3] && !pidgood[l] && !pidgood[n] && !pidgood[s])
					      {
						flag_tmp[1] = true;
					      }
					    if(pidgood[i] && pidgood[j] && pidgood[l] && (flag_acc[i] || flag_acc[j]|| flag_acc[l])  && !pidgood[n] && !pidgood[s])
					      {
						//3E
						flag_tmp[2] = true;
					      }
					    if(pidaccgood[i] && pidgood[j] && pidgood[l] && !pidgood[n] && !pidgood[s] && dblock[i][j][3] && dblock[i][l][3])
					      {
						flag_tmp[3] = true;
					      }
					    if(pidgood[i] && pidgood[j] && pidgood[l] && pidgood[n] && !pidgood[s] && ( flag_acc[i] || flag_acc[j] || flag_acc[l] || flag_acc[n]))
					      {
						//4E
						flag_tmp[4] = true;
					      }
					    if(pidaccgood[i] && pidgood[j] && pidgood[l] && pidgood[n] && !pidgood[s] && dblock[i][j][3] && dblock[i][l][3] && dblock[i][n][3])
					      {
						flag_tmp[5] = true;
					      }
					    if(pidgood[i] && pidgood[j] && pidgood[l] && pidgood[n] && pidgood[s])
					      {
						//5E
						flag_tmp[6] = true;
					      }
					    if(pidaccgood[i] && pidgood[j] && pidgood[l] && pidgood[n] && pidgood[s] && dblock[i][j][3] && dblock[i][l][3] && dblock[i][n][3] && dblock[i][s][3])
					      {
						flag_tmp[7] = true;
					      }
					  }

				      }

				  }
			      }
			  }
		      }
		    if(flag_tmp[0])
		      {
			track_tmp[5][6]++;
			track[5][1]++;
		      }
		    if(flag_tmp[1])
		      {
			track_tmp[5][7]++;
		      }
		    if(flag_tmp[2])
		      {
			track_tmp[5][8]++;
			track[5][1]++;
		      }
		    if(flag_tmp[3])
		      {
			track_tmp[5][9]++;
		      }
		    if(flag_tmp[4])
		      {
			track_tmp[5][10]++;
			track[5][1]++;
		      } 
		    if(flag_tmp[5])
		      {
			track_tmp[5][11]++;
		      }
		    if(flag_tmp[6])
		      {
			track[5][1]++;
			track_tmp[5][12]++;
		      }
		    if(flag_tmp[7])
		      {
			track_tmp[5][13]++;
		      }
		  }

		if(L_tr_n==6)
		  {
		    if(!pidgood[0] && !pidgood[1] && !pidgood[2] && !pidgood[3] && !pidgood[4] && !pidgood[5] )
		      {
			track_tmp[6][0]++;
		      }

		    for(i = 0;i < L_tr_n; i++)
		      {
			for(j = 0;j < L_tr_n;j++)
			  {
			    for(l = 0;l < j;l++)
			      {
				for(n = 0;n < l;n++)
				  {
				    for(s = 0;s < n;s++)
				      {
					for(t = 0;t < s;t++)
					  {
					    if(i != j && i != l && i != n  && i != s  && i != t)
					      {
						if(pidaccgood[i] && !pidgood[j] && !pidgood[l] && !pidgood[n] && !pidgood[s] && !pidgood[t] )
						  {
						    track_tmp[6][i+1]++;
						    track[6][1]++;
						  }
					      }
					  }
				      }
				  }
			      }
			  }

		      }

		    for(i = 0; i<12 ;i++)
		      {
			flag_tmp[i] = false;
		      }
		    for(i = 0;i <L_tr_n; i++)
		      {
			for(j = 0;j <L_tr_n; j++)
			  {
			    for(l = 0;l <L_tr_n; l++)
			      {
				for(n = 0;n <L_tr_n; n++)
				  {
				    for(s = 0;s <L_tr_n; s++)
				      {
					for(t = 0;t <L_tr_n; t++)
					  {
					    if(i != j && i != l && i != n  && i != s  && i != t && j != l && j != n && j != s && j != t && l != n  && l != s && l != t && n != s && n != t && s != t)
					      {
						if(pidgood[i] && pidgood[j] && (flag_acc[i] || flag_acc[j]) && !pidgood[l] && !pidgood[n] && !pidgood[s] && !pidgood[t])
						  {
						    //2E
						    flag_tmp[0] = true;
						  }
						if(pidgood[i] && pidgood[j] && dblock[i][j][3] && !pidgood[l] && !pidgood[n] && !pidgood[s] && !pidgood[t])
						  {
						    flag_tmp[1] = true;
						  }
						if(pidgood[i] && pidgood[j] && pidgood[l] && (flag_acc[i] || flag_acc[j]|| flag_acc[l])  && !pidgood[n] && !pidgood[s] && !pidgood[t])
						  {
						    //3E
						    flag_tmp[2] = true;
						  }
						if(pidaccgood[i] && pidgood[j] && pidgood[l] && !pidgood[n] && !pidgood[s] && !pidgood[t] && dblock[i][j][3] && dblock[i][l][3])
						  {
						    flag_tmp[3] = true;
						  }
						if(pidgood[i] && pidgood[j] && pidgood[l] && pidgood[n] && !pidgood[s] && !pidgood[t] && ( flag_acc[i] || flag_acc[j] || flag_acc[l] || flag_acc[n]))
						  {
						    //4E
						    flag_tmp[4] = true;
						  }
						if(pidaccgood[i] && pidgood[j] && pidgood[l] && pidgood[n] && !pidgood[s]  && !pidgood[t] && dblock[i][j][3] && dblock[i][l][3] && dblock[i][n][3])
						  {
						    flag_tmp[5] = true;
						  }
						if(pidgood[i] && pidgood[j] && pidgood[l] && pidgood[n] && pidgood[s] && !pidgood[t] && ( flag_acc[i] || flag_acc[j] || flag_acc[l] || flag_acc[n] || flag_acc[s]))
						  {
						    //5E
						    flag_tmp[6] = true;
						  }
						if(pidaccgood[i] && pidgood[j] && pidgood[l] && pidgood[n] && pidgood[s] && !pidgood[t] && dblock[i][j][3] && dblock[i][l][3] && dblock[i][n][3] && dblock[i][s][3])
						  {
						    flag_tmp[7] = true;
						  }
						if(pidgood[i] && pidgood[j] && pidgood[l] && pidgood[n] && pidgood[s] && pidgood[t])
						  {
						    //6E
						    flag_tmp[8] = true;
						  }
						if(pidaccgood[i] && pidgood[j] && pidgood[l] && pidgood[n] && pidgood[s] && pidgood[t] && dblock[i][j][3] && dblock[i][l][3] && dblock[i][n][3] && dblock[i][s][3] && dblock[i][t][3])
						  {
						    flag_tmp[9] = true;
						  }

					      }

					  }

				      }

				  }
			      }
			  }
		      }
		  
		    if(flag_tmp[0])
		      {
			track_tmp[6][7]++;
			track[6][1]++;
		      }
		    if(flag_tmp[1])
		      {
			track_tmp[6][8]++;
		      }
		    if(flag_tmp[2])
		      {
			track_tmp[6][9]++;
			track[6][1]++;
		      }
		    if(flag_tmp[3])
		      {
			track_tmp[6][10]++;
		      }
		    if(flag_tmp[4])
		      {
			track_tmp[6][11]++;
			track[6][1]++;
		      } 
		    if(flag_tmp[5])
		      {
			track_tmp[6][12]++;
		      }
		    if(flag_tmp[6])
		      {
			track[6][1]++;
			track_tmp[6][13]++;
		      }
		    if(flag_tmp[7])
		      {
			track_tmp[6][14]++;
		      }
		    if(flag_tmp[8])
		      {
			track[6][1]++;
			track_tmp[6][15]++;
		      }
		    if(flag_tmp[9])
		      {
			track_tmp[6][16]++;
		      }
		  }

		if(L_tr_n==7)
		  {
		    if(!pidgood[0] && !pidgood[1] && !pidgood[2] && !pidgood[3] && !pidgood[4] && !pidgood[5] && !pidgood[6])
		      {
			track_tmp[7][0]++;
		      }

		    for(i = 0;i < L_tr_n; i++)
		      {
			for(j = 0;j < L_tr_n ;j++)
			  {
			    for(l = 0;l < j;l++)
			      {
				for(n = 0;n < l;n++)
				  {
				    for(s = 0;s < n;s++)
				      {
					for(t = 0;t < s;t++)
					  {
					    for(u = 0;u < t;u++)
					      {
						if(i != j && i != l && i != n  && i != s  && i != t && i != u)
						  {
						    if(pidaccgood[i] && !pidgood[j] && !pidgood[l] && !pidgood[n] && !pidgood[s] && !pidgood[t] && !pidgood[u])
						      {
							track_tmp[7][i+1]++;
							track[7][1]++;
						      }
						  }
					      }
					  }
				      }
				  }
			      }
			  }
		      }
		    

		    count_tmp = 0;
		    for(i=0; i<L_tr_n; i++)
		      {
			if(pidaccgood[i])  
			  {
			    count_tmp++;
			  }
		      }
		    if(count_tmp == 2)
		      {
			track[7][1]++;
			track_tmp[7][8]++;
		      }
		    if(count_tmp == 3)
		      {
			track[7][1]++;
			track_tmp[7][9]++;
		      }
		    if(count_tmp == 4)
		      {
			track[7][1]++;
			track_tmp[7][10]++;
		      }
		    if(count_tmp == 5)
		      {
			track[7][1]++;
			track_tmp[7][11]++;
		      }
		    if(count_tmp == 6)
		      {
			track[7][1]++;
			track_tmp[7][12]++;
		      }
		    if(count_tmp == 7)
		      {
			track[7][1]++;
			track_tmp[7][13]++;
		      }
		  }

	      }  //within acc
	  } //n track >1

      }  //sample pid
	
    } //event loop


  event_total = k;
  for(i = 0; i<9; i++)
    {
      for(j = 0; j<20; j++)
	{
	  if(track[i][0]!=0)
	    {
	      track_tmp[i][j] = track_tmp[i][j]/track[i][0];
	    }
	}
    }
  //corrected the zero track loss, reconstruct loss. this is a loose calculation, the exact inefficiency is within this limit.
  //for large acceptance cut, all zero-track included
  //for small acceptance cut, use the correction
  if(track[0][0]/sample_pid_total>0.002)
    {
      track_tmp[0][0] = sample_total /sample_pid_total;
      track[0][0] =  track_tmp[0][0]* track[0][0];
    }
  sample_total = sample_total + track[0][0]+ (1 - track[1][1]/track[1][0]) * ( 2 * track[2][0] + 3 * track[3][0]);
  track[2][0] = track[2][0] + (1 - track[1][1]/track[1][0]) *  2 * track[2][0];
  track[3][0] = track[3][0] + (1 - track[1][1]/track[1][0]) *  3 * track[3][0];
  
  
  for(i = 0; i< 9; i++ )
    {
      track[i][0] = track[i][0]/sample_total;
      track[i][1] = track[i][1]/sample_total;
      track[i][3] = track[i][0] - track[i][1]; 
    }
  track[0][1]= track[0][0]/2.0;
  track[0][2]= track[0][0]/2.0; //low
  track[0][3]= track[0][0]/2.0; //up

  track[8][1]= track[8][0]/2.0;
  track[8][2]= track[8][0]/2.0; //low
  track[8][3]= track[8][0]/2.0; //up


  track[2][2]= (track_tmp[2][3] - track_tmp[2][4] - track_tmp[2][5]) * track[2][0]; 
  track[3][2]= (track_tmp[3][4] - track_tmp[3][5] - track_tmp[3][6] - track_tmp[3][7] + track_tmp[3][8] - track_tmp[3][9] - track_tmp[3][10] - track_tmp[3][11]) * track[3][0]; 
  track[4][2]= (track_tmp[4][5] - track_tmp[4][6] - track_tmp[4][7] - track_tmp[4][8] + track_tmp[4][9] - track_tmp[4][10] + track_tmp[4][11] - track_tmp[4][12]) * track[4][0]; 
  track[5][2]= (track_tmp[5][6] - track_tmp[5][7] + track_tmp[5][8] + track_tmp[5][10] + track_tmp[5][12] ) * track[5][0]; 
  track[6][2]= (track_tmp[6][7] - track_tmp[6][8] + track_tmp[6][9] + track_tmp[6][11] + track_tmp[6][13] + track_tmp[6][15]  ) * track[6][0]; 
  track[7][2]= (track_tmp[7][8] + track_tmp[7][9] + track_tmp[7][10] + track_tmp[7][11] + track_tmp[7][12] + track_tmp[7][13] ) * track[7][0]; 

  vdc_eff = 0;
  vdc_eff_low = 0;
  vdc_eff_up = 0;
  for(i = 0; i< 9; i++)
    {
      vdc_eff += track[i][1];
      vdc_eff_low += track[i][2];
      vdc_eff_up += track[i][3];
    }
  vdc_eff_up =  1 - vdc_eff;  //loose limit, include some events not distinguished


  //output the eff for each case
  cout<<runnum<<"\t"<<momentum<<"\t"<<event_total<<"\t"<<sample_total<<"\t"<<endl;
  for(i = 0; i< 9; i++)
    {
      cout<<track[i][0]<<"\t"<<track[i][1]<<"\t";
    }
  cout<<endl;


  for(i = 0; i< 9; i++)
    {
      cout<<"\n"<<i<<"\t";
      for(j = 0; j<17; j++)
	{
	  cout<<track_tmp[i][j]<<"\t";
	}
    }
  cout<<endl;

  TString singlefilename = filename + "lhrs_vdc_track_eff_";
  singlefilename += eng;
  singlefilename += "_";
  singlefilename += runnum;
  singlefilename += ".txt";
  ofstream outfile(singlefilename.Data(),ios::app);
  cout << singlefilename << endl;
  outfile<<runnum<<"\t"<<momentum<<"\t"<<event_total<<"\t"<<sample_total<<"\t"<<endl;

  for(i = 0; i< 9; i++)
    {
      outfile<<track[i][0]<<"\t"<<track[i][1]<<"\t";
    }
  outfile<<endl;

  for(i = 0; i< 9; i++)
    {
      outfile<<"\n track "<<i<<"\t";
      for(j = 0; j<20; j++)
	{
	  outfile<<track_tmp[i][j]<<"\t";
	}
    }
  outfile.close();

  TString outputfile = filename + "lhrs_vdc_track_eff_total_lists";
  outputfile += eng;
  outputfile += ".txt";
  ofstream resultfile(outputfile.Data(),ios::app);
  cout << outputfile << endl;
  cout<<"runnum"<<"\t"<<"momentum"<<"\t"<<"total event"<<"\t"<<"sample event"<<"\t"<<"one track"<<"\t"<<"total good track"<<"\t"<<"lower uncertainty"<<"\t"<<"upper uncertainty"<<endl;
  cout<<runnum<<"\t"<<momentum<<"\t"<<event_total<<"\t"<<sample_total<<"\t"<<track[1][0]<<"\t"<<vdc_eff<<"\t"<<vdc_eff_low<<"\t"<<vdc_eff_up<<endl;
  resultfile<<runnum<<"\t"<<momentum<<"\t"<<event_total<<"\t"<<sample_total<<"\t"<<track[1][0]<<"\t"<<vdc_eff<<"\t"<<vdc_eff_low<<"\t"<<vdc_eff_up<<endl;
  resultfile.close();
}


void getdb(Char_t* detname, Float_t* detpara, TDatime date)
{
  FILE* dbfile = THaAnalysisObject::OpenFile(detname,date,"","r",1);
  //  FILE* dbfile = THaAnalysisObject::OpenFile((char*)detname,date,"","r",1);
  char buf[256];
  char dum1[100];
  char dum2[100];
  for (Int_t j=0;j<2;j++) fgets(buf,256,dbfile); // get two dummy lines from data file
  Int_t k = fscanf(dbfile,"%f %f",detpara,detpara+1);
  if(k!=2){
    cout<<"Incorrect database, check it.. "<<endl;
  }
  while(fgets(buf,256,dbfile)){
    if(buf[0] == '#' || buf[0] == '\n' )
      continue;
    sscanf(buf,"%s%s",dum1,dum2);
    if((strcmp(dum1,"X,Y,Z")==0||strcmp(dum1,"x,y,z")==0)&&strcmp(dum2,"coords")==0){
      for (Int_t j=0;j<3;j++){
	fscanf(dbfile,"%f",detpara+2);
      }
      break;
    }
  }

  while(fgets(buf,256,dbfile)){
    if(buf[0] == '#' || buf[0] == '\n' )
      continue;
    sscanf(buf,"%s%s",dum1,dum2);
    if((strcmp(dum1,"x,y")==0||strcmp(dum1,"X,Y")==0)&&strcmp(dum2,"position")==0){
      for (Int_t j=0;j<2;j++){
	fscanf(dbfile,"%f",detpara+3+j);
      }
      break;
    }
  }

  while(fgets(buf,256,dbfile)){
    if(buf[0] == '#' || buf[0] == '\n' )
      continue;
    sscanf(buf,"%s%s",dum1,dum2);
    if(strcmp(dum1,"dx")==0&&strcmp(dum2,"and")==0){
      for (Int_t j=0;j<3;j++){
	fscanf(dbfile,"%f %f",detpara+5, detpara+6);
      }
      break;
    }
  }
  cout<<" \n lead glass "<<detname<<endl;
  for(Int_t m =0; m<7; m++)
    {
      cout<<detpara[m]<<endl;
    }
}

