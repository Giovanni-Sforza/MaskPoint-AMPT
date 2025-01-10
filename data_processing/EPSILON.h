//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan  7 14:23:01 2022 by ROOT version 6.22/02
// from TTree epsilon_tree/epsilon_tree
// found on file: epsilonPbP.root
//////////////////////////////////////////////////////////

#ifndef EPSILON_h
#define EPSILON_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class EPSILON {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         epsilon1;
   Float_t         epsilon2;
   Float_t         epsilon3;

   // List of branches
   TBranch        *b_epsilon1;   //!
   TBranch        *b_epsilon2;   //!
   TBranch        *b_epsilon3;   //!

   EPSILON(TTree *tree=0);
   virtual ~EPSILON();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EPSILON_cxx
EPSILON::EPSILON(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("epsilonPbP.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("epsilonPbP.root");
      }
      f->GetObject("epsilon_tree",tree);

   }
   Init(tree);
}

EPSILON::~EPSILON()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EPSILON::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EPSILON::LoadTree(Long64_t entry)
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

void EPSILON::Init(TTree *tree)
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

   fChain->SetBranchAddress("epsilon1", &epsilon1, &b_epsilon1);
   fChain->SetBranchAddress("epsilon2", &epsilon2, &b_epsilon2);
   fChain->SetBranchAddress("epsilon3", &epsilon3, &b_epsilon3);
   Notify();
}

Bool_t EPSILON::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EPSILON::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EPSILON::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EPSILON_cxx
