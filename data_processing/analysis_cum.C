#define analysis_cxx
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <string>
#include <fstream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TProfile.h"

#include "EPSILON.h"
#include "analysis.h"
#include "CorrelationCalculator.h"
using namespace std;
double eta_shifttolab=0.465;
// double eta_shifttolab=0.;
AMPTdata::AMPTdata(string filename,int sys):filename1(filename),sysinp(sys){
        
    // TChain* chain  = new TChain("AMPT","");
		if(filename.find("PbPb") != string::npos)
	{
		eta_shifttolab=0.;
	}
	else
	{
		eta_shifttolab=0.465;
	}
// double eta_shifttolab=0.465;
// double eta_shifttolab=0.;
    char fname[400];
    sprintf(fname,"%s/ampt_afterART.list",filename.c_str());
    int fileNumber = 0;
    char FileList[512];
    ifstream* inputStream = new ifstream;
    inputStream->open(fname);
    if (!(inputStream))
    {
	printf("can not open list file\n");
	return;
    }
    for (;inputStream->good();)
    {
	inputStream->getline(FileList,512);
	if  ( inputStream->good() )
	{
	    TFile *ftmp = new TFile(FileList);
	    if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys()))
	    {
		printf(" file %s error in opening!!!\n",FileList);
	    }
	    else
	    {
		printf(" read in file %s\n",FileList);
		chain->Add(FileList);
		fileNumber++;
	    }
	    delete ftmp;
	}
    }
    printf(" files read in %d\n",fileNumber);

    Init(chain);
    etaCut       = true;
    chargeTracks = true;
    ZPCt0        = false;
    
     memset(PID, 0, sizeof PID);
    //==============================================================================
    ifstream fin;
    fin.open("PID.txt");
    int i=0;
    while(!fin.eof())
    {
        fin>>PID[i][0]>>PID[i][1];
        i++;
    }
    fin.close();
    //==============================================================================
    sprintf(fevname,"%s/trk_event_corr1.root", filename.c_str());
}

AMPTdata::AMPTdata(string filename):filename1(filename){
    
    char fname[400];
    // sprintf(fname,"%s/ampt_zpct0.list",filename.c_str());
    sprintf(fname,"%s/ampt_afterART.list",filename.c_str());

    int fileNumber = 0;
    char FileList[512];
    ifstream* inputStream = new ifstream;
    inputStream->open(fname);
    if (!(inputStream))
    {
	return;
    }
    for (;inputStream->good();)
    {
	inputStream->getline(FileList,512);
	if  ( inputStream->good() )
	{
	    TFile *ftmp = new TFile(FileList);
	    if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys()))
	    {
	    }
	    else
	    {
		chain->Add(FileList);
		fileNumber++;
	    }
	    delete ftmp;
	}
    }
    Init(chain);
    etaCut       = true;
    chargeTracks = false;
    ZPCt0        = false;

     memset(PID, 0, sizeof PID);
    //==============================================================================
    ifstream fin;
    fin.open("PID.txt");
    int i=0;
    while(!fin.eof())
    {
        fin>>PID[i][0]>>PID[i][1];
        i++;
    }
    fin.close();
    //==============================================================================
// }

// void analysis::mkfile(){ 
//     // sprintf(fevname,"%s/trk_event.list",filename1.c_str());    
//     // // sprintf(fevname,"trk_event.list",filename.c_str());
//     // outData.open(fevname);
    
    sprintf(fevname,"%s/trk_event_corr1.root", filename.c_str());
    file = new TFile(fevname);
    trk_tree = (TTree*)file->Get("trk_event");
    Inittrk(trk_tree);
		if(filename.find("PbPb") != string::npos)
	{
		eta_shifttolab=0.;
	}
	else
	{
		eta_shifttolab=0.465;
	}
}

int AMPTdata::find_track_charge(int pid){  //-1, 0, 1
    for(int i=0;i<102;i++)
    {
        if(PID[i][0]==pid)  return PID[i][1];
    }
    return 999;
}
void analysis::TrkLoop()
{
    // info tracks;
    trk_total = 0;
    int TotalTracks = Event_multi;

    for(int iTrk=0; iTrk<TotalTracks; iTrk++)
    {
        
        int _id     = ID[iTrk];
        double _px   = Px[iTrk];
        double _py   = Py[iTrk];
        double _pz   = Pz[iTrk];
        double _mass = Mass[iTrk];
        
        double _x    = X[iTrk];
        double _y    = Y[iTrk];
        double _z    = Z[iTrk];
        double _time = Time[iTrk];
        
        //===============================================
        double _phi  =  atan2(_py,_px);
        if( _px == 0 && _py==0 ) continue; 
        
        TVector3 mP(_px,_py,_pz);
        double _pt     = mP.Perp();
        // if(_pt < 1e-7) continue;
        if(_pt <= PTLOW) continue;
        double _eta = mP.PseudoRapidity()+eta_shifttolab;
        // double _eta = mP.PseudoRapidity();
         
        if(ZPCt0) 
        {
        chargeTracks = false; 
        etaCut = false; 
        }
        
        if(etaCut){
            if(_eta >= ETAHIGH || _eta <= ETALOW) continue;
        }
        if(chargeTracks){//only consider charged tracks
            int trk_chOld  = find_track_charge(_id); // -1, 0, 1
            if( trk_chOld == 0 || trk_chOld==999 ) continue;
        }
        //=================================================

    //    if(_pt >= PTHIGH || _pt <= PTLOW) continue;
    //    if(_eta >= ETAHIGH || _eta <= ETALOW) continue;

        trk_total++;
        
    }//end of track loop
    if (trk_total<500)    trk_event[trk_total].push_back(index);
	Nch->Fill(trk_total);
    
}//
// void analysis::Nchcendef(){
    

//     // sprintf(fevname,"psi_from%d_to%d.list", from, to);
//     // TFile *f1=new TFile("./out.root");
    
    
//     // TH1D *h1;
//     // TH1D *h2;

//     // h1=(TH1D*)f1->Get("Nch");
//     // //h2=(TH1D*)f2->Get("RefImp");
//     double mindel[10];
//     for(int j=9;j>=0;j--){
//     mindel[j]=10;
//     }
//     // char c[10]={'90', '80', '70', '60', '50', '40', '30', '20', '10', '5'};
//     // memset(mindel, 0.1, sizeof mindel);
//     // memset(mindel1, 0.1, sizeof mindel1);
//     double sum1=0.,sum2=0.;
//     for(int i=1;i<501;i++){
//     sum1+=Nch->GetBinContent(i);
//     //cout<<Nch->GetBinCenter(i)<<endl;
//     }
//     double a1=0.,a2=0.;
//     for(int i=1;i<501;i++){
//     double b1=Nch->GetBinLowEdge(i);
//     a1+=Nch->GetBinContent(i);
//     double cent=a1/sum1*100.;
//     //if(fabs(a1/sum1-0.05)<0.005){cout<<i<<" "<<b1<<" "<<a1/sum1<<endl;}
//     if(fabs(a1/sum1-0.1)<0.02){cout<<"90 "<<b1<<" "<<a1/sum1<<endl;}
//     if(fabs(a1/sum1-0.2)<0.01){cout<<"80 "<<b1<<" "<<a1/sum1<<endl;}
//     if(fabs(a1/sum1-0.3)<0.01){cout<<"70 "<<b1<<" "<<a1/sum1<<endl;}
//     if(fabs(a1/sum1-0.4)<0.01){cout<<"60 "<<b1<<" "<<a1/sum1<<endl;}
//     if(fabs(a1/sum1-0.5)<0.01){cout<<"50 "<<b1<<" "<<a1/sum1<<endl;}
//     if(fabs(a1/sum1-0.6)<0.01){cout<<"40 "<<b1<<" "<<a1/sum1<<endl;}
//     if(fabs(a1/sum1-0.7)<0.01){cout<<"30 "<<b1<<" "<<a1/sum1<<endl;}
//     if(fabs(a1/sum1-0.8)<0.01){cout<<"20 "<<b1<<" "<<a1/sum1<<endl;}
//     if(fabs(a1/sum1-0.9)<0.01){cout<<"10 "<<b1<<" "<<a1/sum1<<endl;}
//     if(fabs(a1/sum1-0.95)<0.015){cout<<"5 "<<b1<<" "<<a1/sum1<<endl;}
//     for(int j=0;j<10;j++){
//     if(cent-cb[j]>0 && fabs(cent-cb[j])<mindel[j]){
//         mindel[j]=fabs(cent-cb[j]);
//         IMpcut2[9-j]=b1;}
//         }
//     }
//     for(int j=0;j<10;j++)
//     {
//         IMpcut1[j]=IMpcut2[j+1];

//             outData<<trk_event[i][j]<<" "<<i;
//             outData<<endl;
//     }
// }


void analysis::Getevent(){
    
    file = new TFile(fevname,"recreate");
    trk_tree->Branch("nevent",&nevent,"nevent/I");
    trk_tree->Branch("ievent",ievent,"ievent[nevent]/I");
    // trk_tree->Branch("vpsi1",vpsi1,"vpsi1[nevent]/F");
    // trk_tree->Branch("vpsi2",vpsi2,"vpsi2[nevent]/F");
    // trk_tree->Branch("vpsi3",vpsi3,"vpsi3[nevent]/F");
    // trk_tree->Branch("vpsi4",vpsi4,"vpsi4[nevent]/F");
    // Long64_t nentries = chain->GetEntries();

    // for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //     for (int i=0;i<10;i++) {
    //         if(trk_event[jentry]>=IMpcut1[i] && trk_event[jentry]<IMpcut2[i])
    //         {
    //            cent_event[i].push_back(jentry);

    //         }
    //     }
    
    // }//end event loop
        for (int i=0;i<500;i++)
            {
                nevent=trk_event[i].size();
                // outData<<i<<" "<<trk_event[i].size();
                cout<<i<<" "<<trk_event[i].size()<<endl;
                for (unsigned j=0;j<trk_event[i].size();j++) {
                    ievent[j]=trk_event[i][j];
                    // vpsi1[j]=vpsi[0][ievent[j]];
                    // vpsi2[j]=vpsi[1][ievent[j]];
                    // vpsi3[j]=vpsi[2][ievent[j]];
                    // vpsi4[j]=vpsi[3][ievent[j]];

                // outData<<" "<<trk_event[i][j];
                }
            // outData<<endl;
            trk_tree->Fill();
            }
    trk_tree->Write();
    // file->Write();
    // outData.close();
}

vector<int> AMPTdata::GetNevent(int NchLOW, int NchHIGH)
{
    vector<int>().swap(NchNevent);
    for(int Nchc=NchLOW; Nchc<NchHIGH;Nchc++)
    {

        trk_tree->GetEntry(Nchc);
		NchNevent.push_back(nevent);
	}
	return NchNevent;
}

info AMPTdata::read(int Nchc, int DIMENSION, double PTLOWINP,double PTHIGHINP, int NchNevent )
{      
    vector<int>().swap(tracks.Idxevent);
    // vector<double>().swap(tracks.fChcn2);
    // vector<double>().swap(tracks.fChcn4);
    // vector<double>().swap(tracks.weight2);
    // vector<double>().swap(tracks.weight4);
    vector<vector<double> >().swap(tracks.mp);
	// vector<double>().swap(tracks.Psi);
    // for(int Nchc=Nchc_low; Nchc<Nchc_up;Nchc++)
    // {
        //==============================================================================
        // sprintf(fevname,"trk_event.list");
        // sprintf(fevname,"%s/trk_event.root",filename1.c_str());
        // fin1.open(fevname);
        // for (int i=1;i<Nchc;i++)
        // {
        //     fin1>>itrack>>nevent;
        //     for (int j=0;j<nevent;j++)
        //     {
        //         fin1>>ievent;
        //     }
        // }
        // fin1>>itrack>>nevent;

    // for(int Nchc=NchLOW; Nchc<NchHIGH;Nchc++)
    // {
		trk_tree->GetEntry(Nchc);
        chain->GetEntries();
        for (Long64_t jentry=0; jentry<NchNevent;jentry++)
        {
        // fin1>>ievent;
        chain->GetEntry(ievent[jentry]);
        // double Psin=-vpsi2[jentry];
        //     if(isinf(Psin)) {cout<< -vpsi2[jentry]<<" "<<jentry<<" "<<Nchc<<" "<<filename1<<endl;};
        int TotalTracks = Event_multi;


    //=====================================
    //=======   for Each Event   ==========
        
        for (int i = 0; i < 20; i++)
        {
            for (int j = 0; j < 20; j++)
            {
                Qcos[i][j] = 0;
                Qsin[i][j] = 0;
            }
        }
        TrackAfter_M = 0;
        TrackAfter_P = 0;

        vector<double>().swap(ptrack);
        for(int iTrk=0; iTrk<TotalTracks; iTrk++)
        {
            
            int _id     = ID[iTrk];
            double _px   = Px[iTrk];
            double _py   = Py[iTrk];
            double _pz   = Pz[iTrk];
            double _mass = Mass[iTrk];
            
            double _x    = X[iTrk];
            double _y    = Y[iTrk];
            double _z    = Z[iTrk];
            double _time = Time[iTrk];
            
            //===============================================



            // double _phi  =  atan2(_py,_px);
            // if( _px == 0 && _py==0 ) continue; 
            
            TVector3 mP(_px,_py,_pz);
            double _pt     = mP.Perp();
            double _pmag     = mP.Mag();
            // if(_pt < 1e-7) continue;
            if(_pt >= PTHIGHINP || _pt <= PTLOWINP) continue;
            
            Float_t phi = (_pt != 0.) ? TMath::Pi() + TMath::ATan2(-_py, -_px) : -1;
            Float_t theta = (_pmag != 0.) ? TMath::ACos(_pz / _pmag) : -999.;
            Float_t eta = -TMath::Log(TMath::Tan(0.5 * theta))+eta_shifttolab;
            // eta = -TMath::Log(TMath::Tan(0.5 * theta)) +eta_shifttolab;
            // etaCut       = false;
            if(etaCut){
                if(eta >= ETAHIGH || eta <= ETALOW) continue;
                // if(eta >= ETAHIGH || eta <= 0) continue;
            }
            
            // if(ZPCt0) chargeTracks = false;
            
            if(chargeTracks){//only consider charged tracks
                int trk_chOld  = find_track_charge(_id); // -1, 0, 1
                if( trk_chOld == 0 || trk_chOld==999 ) continue;
            }
            // if (_id!=321) continue;
            //=================================================

            
            if (eta < 0) TrackAfter_M++;
            if (eta > 0) TrackAfter_P++;
            for (int iharm = 0; iharm < 20; iharm++)
                {
                for (int ipow = 0; ipow < 20; ipow++)
                    {
                    Qcos[iharm][ipow] += TMath::Power(1, ipow) * TMath::Cos(iharm * phi);
                    Qsin[iharm][ipow] += TMath::Power(1, ipow) * TMath::Sin(iharm * phi);
                    }
                }

            // tracks.mID.push_back(_id );
            // tracks.mx.push_back(_x   );
            // tracks.my.push_back(_y   );
            // tracks.mz.push_back(_z   );
            // tracks.mt.push_back(_time);
            // tracks.Nch=trk_total;
            // vector<double>().swap(ptrack);

            // event plane
            // double px_event=_px*cos(Psin)-_py*sin(Psin);
            // double py_event=_px*sin(Psin)+_py*cos(Psin);

            // if(filename1.find("PbP5p02TeV") != string::npos)
            // {
            // double phiran=rand()*1.0/RAND_MAX*2.0*PI;
            // double px_rand=_px*cos(phiran)-_py*sin(phiran);
            // double py_rand=_px*sin(phiran)+_py*cos(phiran);
            // }

            // double phiran=rand()*1.0/RAND_MAX*2.0*TMath::Pi();
            // double px_rand=_x*cos(phiran)-_y*sin(phiran);
            // double py_rand=_x*sin(phiran)+_y*cos(phiran);
            // ptrack.push_back(px_rand );
            // ptrack.push_back(py_rand );
            // double px_norm=_px/_pt;
            // double py_norm=_py/_pt;
            // double px_norm=_px/_pmag;
            // double py_norm=_py/_pmag;
            // double pz_norm=_pz/_pmag;
			// if(eta < 0) continue;
            ptrack.push_back(_px );
            ptrack.push_back(_py );
            if (DIMENSION>2) ptrack.push_back(_pz );
            if (DIMENSION>3) ptrack.push_back(sqrt(_px*_px+_py*_py+_pz*_pz+_mass*_mass));
            // if (isnan(_px) ||isnan(_px)) cout<<_x<<" "<<_y<<" "<<Psin<<endl;
            // eventtrack.push_back(ptrack );
            }


            
                correlator.FillQVector(correlator.Qvector, Qcos, Qsin);


		// if (TrackAfter_M > 1 && TrackAfter_P > 1)
        //     {
                tracks.mp.push_back(ptrack  );
                tracks.Idxevent.push_back(ievent[jentry] );
                // if (TrackAfter_M > 0 && TrackAfter_P > 0) {
                // tracks.fChcn2.push_back(correlator.Two(2, -2).Re() );
                // tracks.fChcn2.push_back(correlator.Two(3, -3).Re() );
                // tracks.fChcn2.push_back(correlator.Two(4, -4).Re() );
                // tracks.fChcn2.push_back(correlator.Two(5, -5).Re() );
                // tracks.weight2.push_back(correlator.Two(0, 0).Re() );
                // }
				// else
				// {
                // tracks.fChcn2.push_back(0. );
                // tracks.fChcn2.push_back(0. );
                // tracks.fChcn2.push_back(0. );
                // tracks.fChcn2.push_back(0. );
                // tracks.weight2.push_back(0. );
				// }
                // if (TrackAfter_M > 1 && TrackAfter_P > 1){
                // tracks.fChcn4.push_back(correlator.Four(2, 2, -2, -2).Re() );
                // tracks.fChcn4.push_back(correlator.Four(3, 3, -3, -3).Re() );
                // tracks.fChcn4.push_back(correlator.Four(4, 4, -4, -4).Re() );
                // tracks.fChcn4.push_back(correlator.Four(5, 5, -5, -5).Re() );
                // tracks.weight4.push_back(correlator.Four(0, 0, 0, 0).Re() );
				// tracks.Psi.push_back((atan2(Qsin[2][2], Qcos[2][2]))/2.) ;
				// }
				// else
				// {
                // tracks.fChcn4.push_back(0. );
                // tracks.fChcn4.push_back(0. );
                // tracks.fChcn4.push_back(0. );
                // tracks.fChcn4.push_back(0. );
                // tracks.weight4.push_back(0. );
				// }
            // }



        }//end event loop
	// }

    
//    }
    // fin1.close();
    return tracks;
}

void analysis::Loop()
{
    //   In a ROOT session, you can do:
    //      root> .L analysis.C
    //      root> analysis t
    //      root> t.GetEntry(12); // Fill t data members with entry number 12
    //      root> t.Show();       // Show values of entry 12
    //      root> t.Show(16);     // Read and show values of entry 16
    //      root> t.Loop();       // Loop on all entries
    //
    
    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    chain->SetBranchStatus("*",0);  // disable all branches
    //    chain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    chain->GetEntry(jentry);       //read all branches
    //by  b_branchname->GetEntry(ientry); //read only this branch
    // if (chain == 0) return;
    unsigned nentries = int(chain->GetEntries());
    cout<<"We have "<<nentries<<" events"<<endl;
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        index=jentry;
        nb = chain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        if(jentry%10000==0) cout<<"working on "<<jentry<<endl;
        
        TrkLoop();
        
    }//end event loop
    
    // TFile *saveFile = new TFile("parton.root", "RECREATE");
    sprintf(fnchname,"%s.root", filename1.c_str());
    fout = new TFile(fnchname,"recreate");
    Nch->Write();
    fout -> Write();
    
    // writeHistograms(FileOutput);
    // saveFile->Write();
    delete fout;
    delete Nch;
    // outData.close();
}

//_____________________________________________________________________
TComplex CorrelationCalculator::Q(int n, int p)
{

	if(n>=0) return Qvector[n][p];
	else return TComplex::Conjugate(Qvector[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QsubLeft(int n, int p)
{

	if(n>=0) return QvectorSubLeft[n][p];
	else return TComplex::Conjugate(QvectorSubLeft[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QsubRight(int n, int p)
{

	if(n>=0) return QvectorSubRight[n][p];
	else return TComplex::Conjugate(QvectorSubRight[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QsubMiddle(int n, int p)
{

	if(n>=0) return QvectorSubMiddle[n][p];
	else return TComplex::Conjugate(QvectorSubMiddle[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::p(int n, int p)
{

	if(n>=0) return pvector[n][p];
	else return TComplex::Conjugate(pvector[n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::q(int n, int p)
{

	if(n>=0) return qvector[n][p];
	else return TComplex::Conjugate(qvector[n][p]);

}
//____________________________________________________________________
void CorrelationCalculator::ResetQ(const int nMaxHarm, const int nMaxPow)
{
	for(int i=0; i<nMaxHarm; i++)
	{
		for(int j=0; j<nMaxPow; j++)
		{
			Qvector[i][j] = TComplex(0.,0.);
		}
	}
}
//____________________________________________________________________
TComplex CorrelationCalculator::Two(int n1, int n2)
{

	TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Two_3SubLM(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubMiddle(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Two_3SubRM(int n1, int n2)
{

	TComplex formula = QsubMiddle(n1,1)*QsubRight(n2,1);
	return formula;

}
//
//____________________________________________________________________
TComplex CorrelationCalculator::Two_3SubLR(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubRight(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Three(int n1, int n2, int n3)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
		- Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Three_3Sub(int n1, int n2, int n3)
{

	TComplex formula = QsubLeft(n1,1)*QsubMiddle(n2,1)*QsubRight(n3,1);
	return formula;
}

//____________________________________________________________________
TComplex CorrelationCalculator::Four(int n1, int n2, int n3, int n4)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
		- Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
		+ Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
		+ 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Four_3SubMMLR(int n1, int n2, int n3, int n4)
{
	TComplex formula = QsubMiddle(n1,1)*QsubMiddle(n2,1)*QsubLeft(n3,1)*QsubRight(n4,1)-QsubMiddle(n1+n2,2)*QsubLeft(n3,1)*QsubRight(n4,1);
	return formula;
}

//____________________________________________________________________
TComplex CorrelationCalculator::Four_3SubLLMR(int n1, int n2, int n3, int n4)
{
	TComplex formula = QsubLeft(n1,1)*QsubLeft(n2,1)*QsubMiddle(n3,1)*QsubRight(n4,1)-QsubLeft(n1+n2,2)*QsubMiddle(n3,1)*QsubRight(n4,1);
	return formula;
}

//____________________________________________________________________
TComplex CorrelationCalculator::Four_3SubRRML(int n1, int n2, int n3, int n4)
{
	TComplex formula = QsubRight(n1,1)*QsubRight(n2,1)*QsubMiddle(n3,1)*QsubLeft(n4,1)-QsubRight(n1+n2,2)*QsubMiddle(n3,1)*QsubLeft(n4,1);
	return formula;
}

//____________________________________________________________________
TComplex CorrelationCalculator::Five(int n1, int n2, int n3, int n4, int n5)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
		- Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
		+ 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
		+ Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
		+ Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
		- Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
		- 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
		+ Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
		+ Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
		- Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
		+ Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
		- 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
		- 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
		+ Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
		+ Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
		+ 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
		+ 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
		- 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
		+ Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
		+ Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
		+ 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
		+ 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
		- 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)
		- 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
		- 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5);
	return formula;

}
//___________________________________________________________________
TComplex CorrelationCalculator::Six(int n1, int n2, int n3, int n4, int n5, int n6)
{


	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
		- Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
		+ 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
		+ Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
		+ Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
		- Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
		- 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
		+ Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
		+ Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
		- Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
		+ Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
		- 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
		- 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
		+ Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
		+ Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
		+ 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
		+ 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
		- 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
		+ Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
		+ Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
		+ 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
		+ 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
		- 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
		- 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
		- 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
		- Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
		+ Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
		- 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
		- Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
		- Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
		+ Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
		- 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
		+ 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
		+ Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
		+ Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
		+ Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
		+ Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
		- 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
		- Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
		- 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
		+ 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
		- 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
		+ 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
		+ Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
		+ Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
		+ Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
		+ Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
		- 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
		- Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
		- 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
		+ 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
		- 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
		+ 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
		- 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
		- 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
		- 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
		- Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
		+ Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
		- 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
		- Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
		- Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
		+ Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
		- 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
		+ 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
		- 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
		- 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
		+ 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
		- 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
		+ 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
		+ 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
		- 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
		- 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
		- 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
		- 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
		+ 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
		+ Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
		+ Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
		+ Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
		+ Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
		- 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
		- Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
		- 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
		+ 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
		- 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
		+ 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
		- 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
		- 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
		- 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
		+ 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
		- 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
		+ 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
		+ 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
		+ 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
		+ 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
		- 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
		+ 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
		+ 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
		+ 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
		- 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
		+ 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
		- 120.*Q(n1+n2+n3+n4+n5+n6,6);
	return formula;

}
//_________________________________________________________________________________
TComplex CorrelationCalculator::Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6};

	for(int k=7; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[6] = {0,1,2,3,4,5};
		int iPerm = 0;
		int argument = 0;
		int count = 0;

		// k==6: there is just one combination, we can add it manually
		if(k==6){
			Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
				Six(n1, n2, n3, n4, n5, n6)*Q(n7, 7-k);
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
						Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
								Narray[int(array[3])], Narray[int(array[4])])*
						Q(Narray[int(array[5])]+n7, 7-k);
				}
			}while(std::next_permutation(array, array+6));
		}// k==5

		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
									Narray[int(array[3])])*
							Q(Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==4

		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
							Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==3

		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Two(Narray[int(array[0])], Narray[int(array[1])])*
							Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
									+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==2

		else if(k == 1){
			Correlation = Correlation
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7, 7-k);
		}// k==1

		else if(k == 0){
			Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1+n2+n3+n4+n5+n6+n7, 7-k);
		}// k==0

		else{
			std::cout<<"invalid range of k"<<std::endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6, n7};

	for(int k=8; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[7] = {0,1,2,3,4,5,6};
		int iPerm = 0;
		int argument = 0;
		int count = 0;

		// k==7: there is just one combination, we can add it manually
		if(k==7){
			Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
				Seven(n1, n2, n3, n4, n5, n6, n7)*Q(n8, 8-k);
		}// k==7

		else if(k==6){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4] && array[4] < array[5]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
						Six(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
								Narray[int(array[3])], Narray[int(array[4])], Narray[int(array[5])])*
						Q(Narray[int(array[6])]+n8, 8-k);
				}
			}while(std::next_permutation(array, array+7));
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
									Narray[int(array[3])], Narray[int(array[4])])*
							Q(Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==5

		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])], Narray[int(array[3])])*
							Q(Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==4

		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
							Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==3

		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%120 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Two(Narray[int(array[0])], Narray[int(array[1])])*
							Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
									+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==2

		else if(k == 1){
			Correlation = Correlation
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n7, 1)*Q(n1+n2+n3+n4+n5+n6+n8, 8-k);
		}// k==1

		else if(k == 0){
			Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1+n2+n3+n4+n5+n6+n7+n8, 8-k);
		}// k==0

		else{
			std::cout<<"invalid range of k"<<std::endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;

}

//_____________________________________________________________________________

void CorrelationCalculator::FillQVector(TComplex _Qvector[20][20], double _Qcos[20][20], double _Qsin[20][20]) {
	for(int iharm=0; iharm<20; iharm++)
	{
		for(int ipow=0; ipow<20; ipow++)
		{
			_Qvector[iharm][ipow] = TComplex(_Qcos[iharm][ipow], _Qsin[iharm][ipow]);
		}
	}
}

