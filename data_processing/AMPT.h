#ifndef AMPT_h
#define AMPT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TChain.h"
#include "TBranch.h"
#include "TString.h"
#include "TArrayI.h"
#include "TArrayL.h"
#include "TArrayF.h"

const Int_t kMaxtrack = 50000;
typedef struct 
{

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   Int_t           nevent;
   Int_t           nrun;
   Int_t           multi;
   Float_t         impactpar;
   Int_t           NpartP;
   Int_t           NpartT;
   Int_t           NELP;
   Int_t           NINP;
   Int_t           NELT;
   Int_t           NINT;
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

}Cell_t;
class AMPT {
public :
   //TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   //Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   // List of branches
//------------------

    Cell_t cell;
   Float_t ImPar() const{return cell.impactpar;}
   Int_t Ntrack() const{return cell.multi;}
   Int_t Ninp() const{return cell.NINP;}
   Int_t Nint() const{return cell.NINT;}
    Int_t Npartp() const{return cell.NpartP;}
    Int_t Npartt() const{return cell.NpartT;}  
	
   Int_t mIndx(Int_t i) const{return cell.Indx[i];}       
   Int_t mID(Int_t i) const{return cell.ID[i];}   
   Float_t mPx(Int_t i) const{return cell.fPx[i];}
   Float_t mPy(Int_t i) const{return cell.fPy[i];}
   Float_t mPz(Int_t i) const{return cell.fPz[i];}
   Float_t mX(Int_t i) const{return cell.fX[i];}
   Float_t mY(Int_t i) const{return cell.fY[i];}
   
   
   Float_t mMass(Int_t i) const{return cell.fMass[i];}
   Float_t TIME(Int_t i) const{return cell.fT[i];}
   
   void SetData(TChain *);

};

#endif