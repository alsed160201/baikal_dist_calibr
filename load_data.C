#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h" 
#include "TStyle.h" 
#include "TCanvas.h"
#include "TLegend.h"

#include "BMCEvent.h"
#include "BSource.h"
#include "BSecInteraction.h"
#include "BEventMask.h"
#include "BHelperFunctions.h"
#include "BMuonNamespace.h"
#include "BRecoMuon.h"
#include "BGeomTel.h"

//int NCLUSTER = 2;
TCanvas* cnmuon; // number of muons in event
TCanvas* ctt;    // true first muon time in ns
TCanvas* ctr;    // reco time at ref track point in ns -----------
TCanvas* cpa;    // polar angle of reconstructed muon
TCanvas* cly;    // reco OM LY in p.e.
TCanvas* ct;     // reco OM time in ns
TCanvas* cdt;    // OM time difference to ref time in ns
TCanvas* c2cd;   // reco LY vs distance no selections
TCanvas* c2cdns; //  reco LY vs distance, number of muon selection
TCanvas* c2cdas; //  reco LY vs distance, track angle selection
TCanvas* c2cdts; //  reco LY vs distance, time dif selection
TCanvas* c2cdls; //  reco LY vs distance, ly selection
//TCanvas* c2td;   // reco time vs distance

void load_data(std::string _filelist
		, int save_event = 1
		, int save_pulse = 1)
{
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);

  //----------------output file configuration-------------
  TString fout = "./output/";
  TFile* outputFile_event = new TFile(fout + "/mc_data_event.root","recreate");
  TFile* outputFile_pulse = new TFile(fout + "/mc_data_pulse.root","recreate");
  TTree *trOut_event = new TTree("eDataTree", "Postprocessed event mc data");
  TTree *trOut_pulse = new TTree("pDataTree", "Postprocessed pulse mc data");
  
  Int_t  eventId; 
  Int_t  nMuons;
  Float_t trueTime, theta, phi, refTime;

  trOut_event->Branch("eventId", &eventId);
  trOut_event->Branch("nMuons", &nMuons);
  trOut_event->Branch("trueTime", &trueTime);
  trOut_event->Branch("theta", &theta);
  trOut_event->Branch("phi", &phi);
  trOut_event->Branch("refTime", &refTime);

  Int_t  chanID;
  Float_t pulseLY, pulseTime, dTime, distToPoint_BH;
  trOut_pulse->Branch("eventId", &eventId); 
  trOut_pulse->Branch("chanID", &chanID); 
  trOut_pulse->Branch("pulseLY", &pulseLY); 
  trOut_pulse->Branch("pulseTime", &pulseTime); 
  trOut_pulse->Branch("dTime", &dTime); 
  trOut_pulse->Branch("distToPoint_BH", &distToPoint_BH); 

  //----------------read file------------------------------- 
  int ifile=0;
  char tmp[100];
  ifstream flist_mc;   
  flist_mc.open(_filelist.data());
  while (!flist_mc.eof()){
    ifile++;
    std::string fname_mc;
    flist_mc>>fname_mc; 
    flist_mc.getline(tmp,100,'\n');
    if (fname_mc.empty()) continue;
    TFile fmc(fname_mc.c_str());
    
    std::cout<<"processing "<<fname_mc<<std::endl;
    
    bool skip=false;
    if (fmc.GetSize()<2500000) skip=true;
    if (skip||fmc.IsZombie()) {
      std::cout<<"corrupted file"<<std::endl;
      continue;
    }
    
    TTree* trMC=(TTree*)fmc.Get("Events");
    if (trMC==NULL) continue;
    
    BEvent* bevt=0;
    trMC->SetBranchAddress("BEvent.",&bevt);   
    BGeomTel* bgeomtel=0;
    trMC->SetBranchAddress("BGeomTel.",&bgeomtel);
    BMCEvent* bmcev=0;
    trMC->SetBranchAddress("BMCEvent.",&bmcev);
    BRecoMuon* breco=0;
    trMC->SetBranchAddress("BRecoMuon.",&breco);

    std::cout<<"Events: "<<trMC->GetEntries()<<std::endl;
    
    for (int i=0; i<trMC->GetEntries(); i++){    

      trMC->GetEntry(i);
      eventId += 1;
      
      Double_t eventWeight=bmcev->GetEventWeight();
      if ( eventWeight == 0 ) continue;      
      if (bmcev->GetTrack(0) == NULL) continue;

      //---- mc info ----------------------
      //number of muons, which produced response in the detector 
      nMuons = bmcev->GetResponseMuonsN();

      // true time of first muon
      trueTime = bmcev->GetFirstMuonTime();

      //---- reco info ---------------------
      //Calculate distance from reconstructed track to channels
      theta = breco->GetThetaRec(); //_clHM();
      phi = breco->GetPhiRec(); //_clHM();
      refTime = breco->GetTimeXYZRec();

      Float_t trackThetaRad = TMath::Pi()*(theta)/180;
      Float_t trackPhiRad = TMath::Pi()*(phi)/180;

      TVector3 recoVec(sin(trackThetaRad)*cos(trackPhiRad),
		       sin(trackThetaRad)*sin(trackPhiRad),
		       cos(trackThetaRad));

      //reference point at the muon track and its time:
      TVector3 refPoint = breco->GetXYZRec();
     
      //loop over bevent pulses (fired OM's)
      for (int ipulse = 0; ipulse < bevt->NHits(); ipulse++){
	pulseLY = bevt->Q(ipulse);
	chanID = bevt->HitChannel(ipulse);
	pulseTime = bevt->GetImpulse(ipulse)->GetTime();         //ns
	dTime = pulseTime-refTime;

	TVector3 chanPos = TVector3(bgeomtel->At(chanID)->GetX(),
                                    bgeomtel->At(chanID)->GetY(),
                                    bgeomtel->At(chanID)->GetZ());
          
	distToPoint_BH = BHelperFunctions::GetTrackDistanceToPoint(refPoint, recoVec, chanPos);
        trOut_pulse->Fill();
      }
       trOut_event->Fill();
    }
  }

  outputFile_event->cd();
  trOut_event->Write();
  outputFile_event->Close();

  outputFile_pulse->cd();
  trOut_pulse->Write();
  outputFile_pulse->Close();

}


