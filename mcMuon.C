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

void mcMuon(std::string _filelist,
	    int _minNM = 1,    // number of muons
	    int _maxNM = 1000,
	    float _minPA = 0,     // in degree
	    float _maxPA = 180,
	    float _minTD = -500, // in ns
	    float _maxTD = 500,
	    float _minLY = 5,     // in p.e.
	    int _nmuonnbin = 70,
	    float _nmuonmin = 0,
	    float _nmuonmax = 350,
	    int _thetanbin = 18,  
	    float _thetamin = 0,
	    float _thetamax = 180,// in degree
	    int _timenbin = 300,
	    float _timemin = 0, 
	    float _timemax = 6000,// in ns
	    int _dtimenbin = 200,
	    float _dtimemin = -4000, 
	    float _dtimemax = 4000,// in ns
	    int _lynbin = 100,
	    float _lymin = 10,
	    float _lymax = 110, // in p.e.
	    int _distnbin = 100,
	    float _distmin = 0, 
	    float _distmax = 1000 // in m ??????
)
{
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);

  char stmp[120];
  char sleg[120];
  char stit[120];
  snprintf(stit,sizeof stit,"2020 MC atmospheric muons");

  TFile* outputFile = new TFile("./output/outputFile.root","recreate");

  TH1F* hNmuon = new TH1F("hNmuon","Number of muons in event",_nmuonnbin,_nmuonmin,_nmuonmax);
  TH1F* htruetime = new TH1F("htruetime","MC time of first muon",_timenbin,_timemin,_timemax);
  TH1F* hreftime = new TH1F("hreftime","Time at reco track reference point",_timenbin,_timemin,_timemax);
  TH1F* htheta = new TH1F("htheta","Polar angle of reco muon track",_thetanbin,_thetamin,_thetamax);
  TH1F* hly = new TH1F("hly","Pulse LY",_lynbin,_lymin,_lymax);
  TH1F* htimes = new TH1F("htimes","Pulse times",_timenbin,_timemin,_timemax);
  TH1F* hdt = new TH1F("hdt","OM time - ref time",_dtimenbin,_dtimemin,_dtimemax);

  TH2F* hLYvsTrackDist = new TH2F("hLYvsTrackDist","LY vs dist to OM",
				  _distnbin,_distmin,_distmax,_lynbin,_lymin,_lymax);
  TH2F* hLYvsTrackDistNM = new TH2F("hLYvsTrackDistNM","LY vs dist to OM, cut on Nmuon",
				     _distnbin,_distmin,_distmax,_lynbin,_lymin,_lymax);
  TH2F* hLYvsTrackDistPA = new TH2F("hLYvsTrackDistPA","LY vs dist to OM, cut on polar angle",
				     _distnbin,_distmin,_distmax,_lynbin,_lymin,_lymax);
  TH2F* hLYvsTrackDistTD = new TH2F("hLYvsTrackDistTD","LY vs dist to OM, cut on time dif",
				     _distnbin,_distmin,_distmax,_lynbin,_lymin,_lymax);
  TH2F* hLYvsTrackDistLY = new TH2F("hLYvsTrackDistLY","LY vs dist to OM, cut on LY",
				     _distnbin,_distmin,_distmax,_lynbin,_lymin,_lymax);
  //TH2F* hLYvsDtime = new TH2F("hLYvsDtime","LY vs OM time dif to ref time",
  //			     _dtimenbin,_dtimemin,_dtimemax,_lynbin,_lymin,_lymax);

  float maxNmuon = 0;
  float maxTrueTime = 0;
  float maxRefTime = 0;
  float maxLY = 0;
  float maxRecoTime = 0;
  float maxDist = 0;
  float minDtime = 0;
  float maxDtime = 0;

  //---------- read file 
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
      
      Double_t eventWeight=bmcev->GetEventWeight();
      if ( eventWeight == 0 ) continue;      
      if (bmcev->GetTrack(0) == NULL) continue;

      //---- mc info ----------------------
      //number of muons, which produced response in the detector 
      Int_t nMuons = bmcev->GetResponseMuonsN();
      hNmuon->Fill(nMuons,eventWeight);
      if ( maxNmuon < nMuons ) maxNmuon = nMuons;

      // true time of first muon
      Double_t trueTime = bmcev->GetFirstMuonTime();
      htruetime->Fill(trueTime,eventWeight);
      if ( maxTrueTime < trueTime ) maxTrueTime = trueTime;

      //---- reco info ---------------------
      //Calculate distance from reconstructed track to channels
      Double_t theta = breco->GetThetaRec(); //_clHM();
      Double_t phi = breco->GetPhiRec(); //_clHM();
      double refTime = breco->GetTimeXYZRec();
      hreftime->Fill(refTime,eventWeight);
      if ( maxRefTime < refTime ) maxRefTime = refTime;

      Float_t trackThetaRad = TMath::Pi()*(theta)/180;
      Float_t trackThetaGrad = theta;
      Float_t trackPhiRad = TMath::Pi()*(phi)/180;
      htheta->Fill(trackThetaGrad,eventWeight);
      
      TVector3 recoVec(sin(trackThetaRad)*cos(trackPhiRad),
		       sin(trackThetaRad)*sin(trackPhiRad),
		       cos(trackThetaRad));

      //reference point at the muon track and its time:
      TVector3 refPoint = breco->GetXYZRec();
     
      //loop over bevent pulses (fired OM's)
      for (int ipulse = 0; ipulse < bevt->NHits(); ipulse++){
	float pulseLY = bevt->Q(ipulse);
	int chanID = bevt->HitChannel(ipulse);
	Float_t pulseTime = bevt->GetImpulse(ipulse)->GetTime();         //ns
	float dTime = pulseTime-refTime;

	hly->Fill(pulseLY,eventWeight);
	htimes->Fill(pulseTime,eventWeight);
	hdt->Fill(dTime,eventWeight);

	if ( maxLY < pulseLY ) maxLY = pulseLY;
	if ( maxRecoTime < pulseTime ) maxRecoTime = pulseTime;
	if ( minDtime > dTime ) minDtime = dTime;
	if ( maxDtime < dTime ) maxDtime = dTime;

	TVector3 chanPos = TVector3(bgeomtel->At(chanID)->GetX(),
                                    bgeomtel->At(chanID)->GetY(),
                                    bgeomtel->At(chanID)->GetZ());
            
	//Double_t fDistToPoint=BMuonNamespace::TrackDistanceToPoint(&trueTrack, &somePoint);
	Double_t distToPoint_BH = BHelperFunctions::GetTrackDistanceToPoint(refPoint, recoVec, chanPos);
	if ( maxDist < distToPoint_BH ) maxDist = distToPoint_BH;
      
	hLYvsTrackDist->Fill(distToPoint_BH,pulseLY,eventWeight);
	//hLYvsDtime->Fill(dTime,pulseLY,eventWeight);
	if ( nMuons >= _minNM && nMuons <= _maxNM )
	  hLYvsTrackDistNM->Fill(distToPoint_BH,pulseLY,eventWeight);
	if ( theta >= _minPA && theta <= _maxPA )
	  hLYvsTrackDistPA->Fill(distToPoint_BH,pulseLY,eventWeight);
	if ( dTime >= _minTD && dTime <= _maxTD )
	  hLYvsTrackDistTD->Fill(distToPoint_BH,pulseLY,eventWeight);
	if ( pulseLY >= _minLY )
	  hLYvsTrackDistLY->Fill(distToPoint_BH,pulseLY,eventWeight);
      }
    }
  }

  std::cout << "Max number of muons: " << maxNmuon << std::endl;
  std::cout << "Max true first muon time: " << maxTrueTime << " ns" << std::endl;
  std::cout << "Max reco time at ref point: " << maxRefTime << " ns" << std::endl;
  std::cout << "Max reco pulse LY: " << maxLY << " p.e." << std::endl;
  std::cout << "Max reco pulse time: " << maxRecoTime << " ns" << std::endl;
  std::cout << minDtime << " < time dif < " << maxDtime << " ns" << std::endl;
  std::cout << "Max distance to reco track: " << maxDist << " m" << std::endl;
 


  //------------- plot to canvas
  if ( gROOT->GetListOfCanvases()->FindObject("cnmuon") == NULL )
    cnmuon = new TCanvas("cnmuon","N muons", 10, 10, 400, 400);
  cnmuon->cd();
  snprintf(stmp,sizeof stmp,"%s, Nevt = %d",stit,int(hNmuon->GetEntries()));  
  hNmuon->SetTitle(stmp);
  hNmuon->GetXaxis()->SetTitle("Number of muons");
  hNmuon->GetYaxis()->SetTitle("entries");
  hNmuon->DrawCopy();

  if ( gROOT->GetListOfCanvases()->FindObject("cly") == NULL )
    cly = new TCanvas("cly","LY", 10, 510, 400, 400);
  cly->cd();
  cly->SetLogy(1);
  snprintf(stmp,sizeof stmp,"%s, N = %d",stit,int(hly->GetEntries()));  
  hly->SetTitle(stmp);
  hly->GetXaxis()->SetTitle("LY in OMs [p.e.]");
  hly->GetYaxis()->SetTitle("entries");
  hly->DrawCopy();

  if ( gROOT->GetListOfCanvases()->FindObject("ctt") == NULL )
    ctt = new TCanvas("ctt","True time", 510, 10, 400, 400);
  ctt->cd();
  TLegend* legctt = new TLegend(0.5,0.7,0.85,0.85);
  legctt->SetTextSize(0.045);
  snprintf(stmp,sizeof stmp,"%s, Nevt = %d",stit,int(htruetime->GetEntries()));  
  htruetime->SetTitle(stmp);
  htruetime->GetXaxis()->SetTitle("Time [ns]");
  htruetime->GetYaxis()->SetTitle("entries");
  htruetime->DrawCopy();
  snprintf(sleg, sizeof sleg,"mc first muon");
  legctt->AddEntry(htruetime,sleg,"l");
  hreftime->SetLineColor(2);
  hreftime->DrawCopy("same");
  snprintf(sleg, sizeof sleg,"reco at ref point");
  legctt->AddEntry(hreftime,sleg,"l");
  legctt->Draw("same");
  ctt->Update();

  if ( gROOT->GetListOfCanvases()->FindObject("ct") == NULL )
    ct = new TCanvas("ct","Time", 510, 310, 400, 400);
  ct->cd();
  snprintf(stmp,sizeof stmp,"%s, N = %d",stit,int(htimes->GetEntries()));  
  htimes->SetTitle(stmp);
  htimes->GetXaxis()->SetTitle("OM time [ns]");
  htimes->GetYaxis()->SetTitle("entries");
  htimes->DrawCopy();

  if ( gROOT->GetListOfCanvases()->FindObject("cdt") == NULL )
    cdt = new TCanvas("cdt","Time dif", 510, 610, 400, 400);
  cdt->cd();
  snprintf(stmp,sizeof stmp,"%s, N = %d",stit,int(hdt->GetEntries()));  
  hdt->SetTitle(stmp);
  hdt->GetXaxis()->SetTitle("OM time - Ref time [ns]");
  hdt->GetYaxis()->SetTitle("entries");
  hdt->DrawCopy();


  if ( gROOT->GetListOfCanvases()->FindObject("c2cd") == NULL )
    c2cd = new TCanvas("c2cd","LY vs Dist", 1010, 10, 600, 400);
  c2cd->cd();
  snprintf(stmp,sizeof stmp,"%s, no sel, N = %d",stit,int(hLYvsTrackDist->GetEntries()));  
  hLYvsTrackDist->SetTitle(stmp);
  hLYvsTrackDist->GetXaxis()->SetTitle("Distance from reco track to OM [m]");
  hLYvsTrackDist->GetYaxis()->SetTitle("LY [p.e.]");
  hLYvsTrackDist->DrawCopy("colz");

  if ( gROOT->GetListOfCanvases()->FindObject("c2cdns") == NULL )
    c2cdns = new TCanvas("c2cdns","LY vs Dist, Nmuon sel", 1010, 210, 600, 400);
  c2cdns->cd();
  snprintf(stmp,sizeof stmp,"%s, %d #leq Nmuon #leq %d, N = %d (%5.3f)",
	   stit,_minNM, _maxNM, int(hLYvsTrackDistNM->GetEntries()),
	   float(hLYvsTrackDistNM->GetEntries())/float(hLYvsTrackDist->GetEntries()));  
  hLYvsTrackDistNM->SetTitle(stmp);
  hLYvsTrackDistNM->GetXaxis()->SetTitle("Distance from reco track to OM [m]");
  hLYvsTrackDistNM->GetYaxis()->SetTitle("LY [p.e.]");
  hLYvsTrackDistNM->DrawCopy("colz");

  if ( gROOT->GetListOfCanvases()->FindObject("c2cdas") == NULL )
    c2cdas = new TCanvas("c2cdas","LY vs Dist, Polar ang sel", 1010, 410, 600, 400);
  c2cdas->cd();
  snprintf(stmp,sizeof stmp,"%s, %3.0f #leq #theta #leq %3.0f, N = %d (%5.3f)",
	   stit,_minPA, _maxPA, int(hLYvsTrackDistPA->GetEntries()),
	   float(hLYvsTrackDistPA->GetEntries())/float(hLYvsTrackDist->GetEntries()));  
  hLYvsTrackDistPA->SetTitle(stmp);
  hLYvsTrackDistPA->GetXaxis()->SetTitle("Distance from reco track to OM [m]");
  hLYvsTrackDistPA->GetYaxis()->SetTitle("LY [p.e.]");
  hLYvsTrackDistPA->DrawCopy("colz");

  if ( gROOT->GetListOfCanvases()->FindObject("c2cdts") == NULL )
    c2cdts = new TCanvas("c2cdts","LY vs Dist, time dif sel", 1010, 610, 600, 400);
  c2cdts->cd();
  snprintf(stmp,sizeof stmp,"%s, %4.0f #leq dt #leq %4.0f ns, N = %d (%5.3f)",
	   stit,_minTD, _maxTD, int(hLYvsTrackDistTD->GetEntries()),
	   float(hLYvsTrackDistTD->GetEntries())/float(hLYvsTrackDist->GetEntries()));  
  hLYvsTrackDistTD->SetTitle(stmp);
  hLYvsTrackDistTD->GetXaxis()->SetTitle("Distance from reco track to OM [m]");
  hLYvsTrackDistTD->GetYaxis()->SetTitle("LY [p.e.]");
  hLYvsTrackDistTD->DrawCopy("colz");


  outputFile->cd();
  hNmuon->Write();
  htruetime->Write();
  htheta->Write();
  hreftime->Write();
  hly->Write();
  htimes->Write();
  hdt->Write();
  hLYvsTrackDist->Write();
  hLYvsTrackDistNM->Write();
  hLYvsTrackDistPA->Write();
  hLYvsTrackDistTD->Write();
  hLYvsTrackDistLY->Write();
  outputFile->Close();
}


