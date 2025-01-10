#ifndef analysis_h
#define analysis_h
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TProfile.h"
// #include "string.h"
#include <string>
#include <TH2.h>
#include <fstream>
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "CorrelationCalculator.h"
#include "TTree.h"

#include <iostream>

using namespace std;
using std::vector;
using std::string;

const int MAXTRACK    = 3000000;

const int NHAR = 4;//v1, v2, v3
const double PI = acos(-1.0);
const double ETACUT = 2.5;
const double RAPCUT = 0.5;

const double ETALOW  = 0;
const double ETAHIGH = 2.4;

double RAPHigh = 2.4;
double RAPLow = -2.4;

const double PTLOW  = 0.4;
const double PTHIGH = 12.;

struct info { 
vector<int> Idxevent;
// std::vector<double> mID;
// std::vector<double> mx;
// std::vector<double> my;
// std::vector<double> mz;
// std::vector<double> mt;
// std::vector<double> psin;
// std::vector<double> fChcn2;
// std::vector<double> fChcn4;
// std::vector<double> weight2;
// std::vector<double> weight4;
std::vector<vector<double> > mp;
// std::vector<vector<double> > mr;
// std::vector<double> Psi;
}; 


class AMPTdata {
public :
    AMPTdata(string filename, int sys);
    AMPTdata(string filename);
    int find_track_charge(int pid);
    char fevname[200],fnchname[200];
    char name[200];
    int sysinp;
    int PID[102][2];
    int TotalParticleSpieces;
    bool chargeTracks;
    bool etaCut;
    bool ZPCt0;
   Int_t evno, totpart;
   CorrelationCalculator correlator;
   double Qcos[20][20];
   double Qsin[20][20];
   Int_t TrackAfter_P = {};
   Int_t TrackAfter_M = {};

   // double fChcn2[4];
   // double fChcn4[4];
   // double fChc422;
   // double fChc532;
   // double fChsc3232;
   // double fChsc4242;
   vector<double>rtrack;
   vector<double>ptrack;
	vector<int> NchNevent;

   //  vector<vector<double> >eventtrack;
   string filename1;
   TFile *file;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain


   Int_t           itrack;
   Int_t           nevent;
   Int_t           Event_nevent;
   Int_t           Event_nrun;
   Int_t           Event_multi;
   Float_t         Event_impactpar;
   Int_t           Event_NpartP;
   Int_t           Event_NpartT;
   Int_t           Event_NELP;
   Int_t           Event_NINP;
   Int_t           Event_NELT;
   Int_t           Event_NINT;
   Int_t           Indx[MAXTRACK];   //[multi]
   Int_t           ID[MAXTRACK];   //[multi]
   Float_t         Px[MAXTRACK];   //[multi]
   Float_t         Py[MAXTRACK];   //[multi]
   Float_t         Pz[MAXTRACK];   //[multi]
   Float_t         Mass[MAXTRACK];   //[multi]
   Float_t         X[MAXTRACK];   //[multi]
   Float_t         Y[MAXTRACK];   //[multi]
   Float_t         Z[MAXTRACK];   //[multi]
   Float_t         Time[MAXTRACK];   //[multi]
   Int_t           ievent[MAXTRACK];
   // Float_t         vpsi1[MAXTRACK];   //[multi]
   // Float_t         vpsi2[MAXTRACK];   //[multi]
   // Float_t         vpsi3[MAXTRACK];   //[multi]
   // Float_t         vpsi4[MAXTRACK];   //[multi]

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_Indx;   //!
   TBranch        *b_ID;   //!
   TBranch        *b_Px;   //!
   TBranch        *b_Py;   //!
   TBranch        *b_Pz;   //!
   TBranch        *b_Mass;   //!
   TBranch        *b_X;   //!
   TBranch        *b_Y;   //!
   TBranch        *b_Z;   //!
   TBranch        *b_Time;   //!
   TBranch        *b_nevent;
   TBranch        *b_ievent;
   // TBranch        *b_vpsi1;
   // TBranch        *b_vpsi2;
   // TBranch        *b_vpsi3;
   // TBranch        *b_vpsi4;
   TChain* chain  = new TChain("AMPT","");
   TChain* chainzpct0  = new TChain("AMPT","");
   TTree *trk_tree = new TTree("trk_event","track_event");

   vector<int> GetNevent(int NchLOW, int NchHIGH);
   info read(int Nchc, int DIMENSION, double PTLOWINP,double PTHIGHINP, int NchNevent);
   info tracks;
   vector<int> nentry;
   ifstream fin1;
    
   AMPTdata(TTree *tree=0);
   virtual ~AMPTdata();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Inittrk(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

class analysis: public AMPTdata {
public :
   using AMPTdata::AMPTdata;
   //  AMPTdata(){};
   //  void mkfile();
   
    void InitHist();
    void Nchcendef();
   //  void Getevent();
    
    void Getevent();
    void TrkLoop();
    void AMPT2vn();
    void FillQVector(Float_t phi = -999., Float_t eta = -999.);
    ofstream outData;

    int trk_total;
    float Psi[NHAR];
   //  std::vector<int> cent_event[10];
    vector<int>trk_event[501];
    int index;
    double cb[10]={10, 20, 30, 40, 50, 60, 70, 80, 90, 95};
   double IMpcut1[10]={0.};
   double IMpcut2[11]={0.}; //0,1,2,3,4,5,6,7,8,9
   
    TFile *fout;
   virtual void     Loop();
   TH1D *Nch = new TH1D("Nch","",500,0.,500.);

    vector<double>vpsi[NHAR];
   // virtual ~AMPTdata();
};

// class readdata: public AMPTdata {
// public :

//     using AMPTdata::AMPTdata;
//     info read(int Nchc);
//     info tracks;
//     int itrack,ievent;
//     vector<int> nentry;
//     ifstream fin1;
//    // virtual ~readdata();
// };
#endif

#ifdef analysis_cxx
AMPTdata::AMPTdata(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ampt_afterART.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ampt_afterART.root");
      }
      f->GetObject("AMPT",tree);

   }
   Init(tree);
}

AMPTdata::~AMPTdata()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AMPTdata::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AMPTdata::LoadTree(Long64_t entry)
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

void AMPTdata::Init(TTree *tree)
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

   fChain->SetBranchAddress("Event", &Event_nevent, &b_Event);
   fChain->SetBranchAddress("Indx", Indx, &b_Indx);
   fChain->SetBranchAddress("ID", ID, &b_ID);
   fChain->SetBranchAddress("Px", Px, &b_Px);
   fChain->SetBranchAddress("Py", Py, &b_Py);
   fChain->SetBranchAddress("Pz", Pz, &b_Pz);
   fChain->SetBranchAddress("Mass", Mass, &b_Mass);
   fChain->SetBranchAddress("X", X, &b_X);
   fChain->SetBranchAddress("Y", Y, &b_Y);
   fChain->SetBranchAddress("Z", Z, &b_Z);
   fChain->SetBranchAddress("Time", Time, &b_Time);
   Notify();
}

void AMPTdata::Inittrk(TTree *tree)
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

   fChain->SetBranchAddress("nevent", &nevent, &b_nevent);
   fChain->SetBranchAddress("ievent", &ievent, &b_ievent);
   // fChain->SetBranchAddress("vpsi1", &vpsi1, &b_vpsi1);
   // fChain->SetBranchAddress("vpsi2", &vpsi2, &b_vpsi2);
   // fChain->SetBranchAddress("vpsi3", &vpsi3, &b_vpsi3);
   // fChain->SetBranchAddress("vpsi4", &vpsi4, &b_vpsi4);
   Notify();
}

Bool_t AMPTdata::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AMPTdata::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AMPTdata::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analysis_cxx
