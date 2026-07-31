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

#include "./helpers/help_functions.C"

void load_data(std::string _filelist,
		  int save_event = 1,
		  int save_pulse = 1,
		  int save_gentrk = 1,
	          float _minTD = -400, // in ns
	          float _maxTD = 400,
		  float _minLY = 5     // in p.e.
)
{
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);

  //----------------output file configuration-------------

  //----------------configuration for event data file-----------------
  TString fout = "./output/";
  TFile* outputFile_event = TFile::Open(fout + "/mc_reco_event.root","recreate");
  TTree *trOut_event = new TTree("eRecoTree", "Postprocessed event mc reco data");
  
  Int_t  eventId, clusterID, number = 0; 
  Int_t  nRespMuons, nTrueMuons;
  Float_t eventWeight, trueTime, trueFirstMuonEnergy, trueFirstMuonTheta, trueFirstMuonPhi, truePrimaryTheta, truePrimaryPhi, truePrimaryEnergy, trueBundleEnergy;
  Float_t RecoTheta, RecoPhi, RecoRefTime;

  trOut_event->Branch("eventId", &eventId);
  trOut_event->Branch("number", &number);
  trOut_event->Branch("clusterID", &clusterID);
  trOut_event->Branch("nRespMuons", &nRespMuons);
  trOut_event->Branch("nTrueMuons", &nTrueMuons);
  trOut_event->Branch("eventWeight", &eventWeight);

  trOut_event->Branch("trueTime", &trueTime);
  trOut_event->Branch("trueFirstMuonEnergy", &trueFirstMuonEnergy);
  trOut_event->Branch("trueFirstMuonTheta", &trueFirstMuonTheta);
  trOut_event->Branch("trueFirstMuonPhi", &trueFirstMuonPhi);

  trOut_event->Branch("truePrimaryTheta", &truePrimaryTheta);
  trOut_event->Branch("truePrimaryPhi", &truePrimaryPhi);
  trOut_event->Branch("truePrimaryEnergy", &truePrimaryEnergy);
  trOut_event->Branch("trueBundleEnergy", &trueBundleEnergy);
  

  trOut_event->Branch("RecoTheta", &RecoTheta);
  trOut_event->Branch("RecoPhi", &RecoPhi);
  trOut_event->Branch("RecoRefTime", &RecoRefTime);

  //----------------configuration for gen muons data file-----------------
  TFile* outputFile_gentrk = TFile::Open(fout + "/mc_gen_trk.root","recreate");
  TTree *trOut_gentrk = new TTree("tGenTree", "Postprocessed track mc gen data");

  Float_t trueMuonTheta, trueMuonPhi, trueMuonEnergy, trueMuonDelay;
  trOut_gentrk->Branch("eventId", &eventId);
  trOut_gentrk->Branch("number", &number);
  trOut_gentrk->Branch("trueMuonTheta", &trueMuonTheta);
  trOut_gentrk->Branch("trueMuonPhi", &trueMuonPhi);
  trOut_gentrk->Branch("trueMuonEnergy", &trueMuonEnergy);
  trOut_gentrk->Branch("trueMuonDelay", &trueMuonDelay);


  //----------------configuration for pulse data file-----------------

  TFile* outputFile_pulse = TFile::Open(fout + "/mc_reco_pulse.root","recreate");
  TTree *trOut_pulse = new TTree("pRecoTree", "Postprocessed pulse mc reco data");

  Int_t  chanID;
  Float_t pulseLY,  pulseTime, dTime, RecoDistToPoint_BH, trueFirstMuonDistToPoint_BH;
  trOut_pulse->Branch("eventId", &eventId); 
  trOut_pulse->Branch("number", &number);
  trOut_pulse->Branch("chanID", &chanID); 
  trOut_pulse->Branch("pulseLY", &pulseLY); 
  trOut_pulse->Branch("pulseTime", &pulseTime); 
  trOut_pulse->Branch("dTime", &dTime); 
  trOut_pulse->Branch("RecoDistToPoint_BH", &RecoDistToPoint_BH); 
  trOut_pulse->Branch("trueFirstMuonDistToPoint_BH", &trueFirstMuonDistToPoint_BH); 

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
      eventId =bmcev->GetEventN();
      number += 1; 
      eventWeight=bmcev->GetEventWeight();
      if ( eventWeight == 0 ) continue;      
      if (bmcev->GetTrack(0) == NULL) continue;

      //---------------- mc info ----------------------
      //number of muons, which produced response in the detector 
      nRespMuons = bmcev->GetResponseMuonsN();
      nTrueMuons = bmcev->GetMuonsN();
      // true time of first muon
      trueTime = bmcev->GetFirstMuonTime();

      // true first muon properties
      Int_t imuon=0;
      trueFirstMuonEnergy=bmcev->GetTrack(imuon)->GetMuonEnergy();
      trueFirstMuonTheta=bmcev->GetTrack(imuon)->GetTheta();
      trueFirstMuonPhi=bmcev->GetTrack(imuon)->GetPhi();

      Float_t trueFirstMuonThetaRad = TMath::Pi()*(trueFirstMuonTheta)/180;
      Float_t trueFirstMuonPhiRad = TMath::Pi()*(trueFirstMuonPhi)/180;

      TVector3 trueFirstMuonVec(sin(trueFirstMuonThetaRad)*cos(trueFirstMuonPhiRad),
                       sin(trueFirstMuonThetaRad)*sin(trueFirstMuonPhiRad),
                       cos(trueFirstMuonThetaRad));

      TVector3 trueFirstMuonPoint(bmcev->GetTrack(imuon)->GetX(),
                         bmcev->GetTrack(imuon)->GetY(),
                         bmcev->GetTrack(imuon)->GetZ());

      for (int itrk = 0; itrk < bmcev->GetResponseMuonsN(); itrk++){
	  trueMuonEnergy = bmcev->GetTrack(itrk)->GetMuonEnergy();
	  trueMuonTheta = bmcev->GetTrack(itrk)->GetTheta();
	  trueMuonPhi = bmcev->GetTrack(itrk)->GetPhi();
	  trueMuonDelay = bmcev->GetTrack(itrk)->GetDelay(); //with respect to first one (nanosec)

	  Float_t trueMuonThetaRad = TMath::Pi()*(trueMuonTheta)/180;
	  Float_t trueMuonPhiRad = TMath::Pi()*(trueMuonPhi)/180;

	  TVector3 trueMuonVec(sin(trueMuonThetaRad)*cos(trueMuonPhiRad),
		       sin(trueMuonThetaRad)*sin(trueMuonPhiRad),
		       cos(trueMuonThetaRad));

	  TVector3 trueMuonPoint(bmcev->GetTrack(itrk)->GetX(),
		         bmcev->GetTrack(itrk)->GetY(),
		         bmcev->GetTrack(itrk)->GetZ());


	  trOut_gentrk->Fill();
      }


      // Primary Particle Angles
      truePrimaryTheta = bmcev->GetPrimaryParticleTheta();
      truePrimaryPhi = bmcev->GetPrimaryParticlePhi();
      truePrimaryEnergy = bmcev->GetPrimaryParticleEnergy();
      trueBundleEnergy  = bmcev->GetSumEnergyBundleReg();
      // response cluster number
      clusterID = extractClusterId(fname_mc);



      //---------------- reco info ---------------------
      //Calculate distance from reconstructed track to channels
      RecoTheta = breco->GetThetaRec(); //_clHM();
      RecoPhi = breco->GetPhiRec(); //_clHM();
      RecoRefTime = breco->GetTimeXYZRec();

      Float_t trackThetaRad = TMath::Pi()*(RecoTheta)/180;
      Float_t trackPhiRad = TMath::Pi()*(RecoPhi)/180;

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
	dTime = pulseTime-RecoRefTime;

	TVector3 chanPos = TVector3(bgeomtel->At(chanID)->GetX(),
                                    bgeomtel->At(chanID)->GetY(),
                                    bgeomtel->At(chanID)->GetZ());
          
	RecoDistToPoint_BH = BHelperFunctions::GetTrackDistanceToPoint(refPoint, recoVec, chanPos);
	trueFirstMuonDistToPoint_BH = BHelperFunctions::GetTrackDistanceToPoint(trueFirstMuonPoint, trueFirstMuonVec, chanPos);
	
	if ( pulseLY < _minLY) continue;
	if ( dTime < _minTD || dTime > _maxTD ) continue;

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

  outputFile_gentrk->cd();
  trOut_gentrk->Write();
  outputFile_gentrk->Close();
}


