/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug 21 16:54:59 2015 by ROOT version 5.34/13
// from TTree T/Hall A Analyzer Output DST
// found on file: gdh_1805.root
//////////////////////////////////////////////////////////

#ifndef temp_h
#define temp_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class temp {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types

   Int_t           Ndata_R_tr_r_ph;
   Double_t        R_tr_r_ph[10];   //[Ndata.R.tr.r_ph]
   Int_t           Ndata_R_tr_r_th;
   Double_t        R_tr_r_th[10];   //[Ndata.R.tr.r_th]
   Int_t           Ndata_R_tr_r_x;
   Double_t        R_tr_r_x[10];   //[Ndata.R.tr.r_x]
   Int_t           Ndata_R_tr_r_y;
   Double_t        R_tr_r_y[10];   //[Ndata.R.tr.r_y]




   Double_t        R_ps_e;
   Double_t        R_sh_e;
   Double_t        R_gold_dp;
   Double_t        R_gold_y;
   Double_t        R_gold_th;
   Double_t        R_gold_ph;
   Double_t        R_cer_asum_c;
   Double_t        R_tr_n;
   Double_t        beam_H_helicity;
   Double_t        evtype;
 //THaEvent        *Event_Branch;
   UInt_t          fEvtHdr_fEvtNum;
   Int_t           fEvtHdr_fEvtType;
   Int_t           fEvtHdr_fEvtLen;
   Double_t        fEvtHdr_fEvtTime;
   Int_t           fEvtHdr_fHelicity;
   Int_t           fEvtHdr_fRun;

   // List of branches

   TBranch        *b_Ndata_R_tr_r_ph;   //!
   TBranch        *b_R_tr_r_ph;   //!
   TBranch        *b_Ndata_R_tr_r_th;   //!
   TBranch        *b_R_tr_r_th;   //!
   TBranch        *b_Ndata_R_tr_r_x;   //!
   TBranch        *b_R_tr_r_x;   //!
   TBranch        *b_Ndata_R_tr_r_y;   //!
   TBranch        *b_R_tr_r_y;   //!
  //!


   TBranch        *b_R_gold_dp;   //!
   TBranch        *b_R_gold_y;   //!
   TBranch        *b_R_gold_th;   //!
   TBranch        *b_R_gold_ph;   //!
   TBranch        *b_R_ps_e;   //!
   TBranch        *b_R_sh_e;   //!
   TBranch        *b_R_cer_asum_c;   //!
   TBranch        *b_R_tr_n;   //!
   TBranch        *b_beam_H_helicity;   //!
   TBranch        *b_evtype;   //!

   temp(TTree *tree=0);
   virtual ~temp();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef temp_cxx
temp::temp(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("gdh_2080.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("gdh_2080.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

temp::~temp()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t temp::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t temp::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void temp::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);


   fChain->SetBranchAddress("Ndata.R.tr.r_ph", &Ndata_R_tr_r_ph, &b_Ndata_R_tr_r_ph);
   fChain->SetBranchAddress("R.tr.r_ph", R_tr_r_ph, &b_R_tr_r_ph);
   fChain->SetBranchAddress("Ndata.R.tr.r_th", &Ndata_R_tr_r_th, &b_Ndata_R_tr_r_th);
   fChain->SetBranchAddress("R.tr.r_th", R_tr_r_th, &b_R_tr_r_th);
   fChain->SetBranchAddress("Ndata.R.tr.r_x", &Ndata_R_tr_r_x, &b_Ndata_R_tr_r_x);
   fChain->SetBranchAddress("R.tr.r_x", R_tr_r_x, &b_R_tr_r_x);
   fChain->SetBranchAddress("Ndata.R.tr.r_y", &Ndata_R_tr_r_y, &b_Ndata_R_tr_r_y);
   fChain->SetBranchAddress("R.tr.r_y", R_tr_r_y, &b_R_tr_r_y);



   fChain->SetBranchAddress("R.gold.dp", &R_gold_dp, &b_R_gold_dp);
   fChain->SetBranchAddress("R.gold.y", &R_gold_y, &b_R_gold_y);
   fChain->SetBranchAddress("R.gold.th", &R_gold_th, &b_R_gold_th);
   fChain->SetBranchAddress("R.gold.ph", &R_gold_ph, &b_R_gold_ph);
   fChain->SetBranchAddress("R.ps.e", &R_ps_e, &b_R_ps_e);
   fChain->SetBranchAddress("R.sh.e", &R_sh_e, &b_R_sh_e);
   fChain->SetBranchAddress("R.cer.asum_c", &R_cer_asum_c, &b_R_cer_asum_c);
   fChain->SetBranchAddress("R.tr.n", &R_tr_n, &b_R_tr_n);
   fChain->SetBranchAddress("beam.H.helicity", &beam_H_helicity, &b_beam_H_helicity);
   fChain->SetBranchAddress("DR.evtypebits", &evtype, &b_evtype);


   Notify();
}

Bool_t temp::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void temp::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t temp::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef temp_cxx
