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
TCanvas* cdtevsr;    // Reco minus Expected OM time difference in ns
TCanvas* cnoms;    // Number of fired OM vs dist from track
TCanvas* cnomt;    // Number of total OM vs dist from track
TCanvas* c2cd;   // reco LY vs distance no selections
TCanvas* c2cdns; //  reco LY vs distance, number of muon selection
TCanvas* c2cdas; //  reco LY vs distance, track angle selection
TCanvas* c2cdts; //  reco LY vs distance, time dif selection
TCanvas* c2cdtsevsr; //  reco LY vs distance, Expected minus Reco time dif selection
TCanvas* c2cdls; //  reco LY vs distance, ly selection
//TCanvas* c2td;   // reco time vs distance

void mcMuon(std::string _filelist,
	    int _minNM = 1,    // number of muons
	    int _maxNM = 1000,
	    float _minPA = 0,     // in degree
	    float _maxPA = 180,
	    float _minTD = -500, // in ns
	    float _maxTD = 500,
	    float _minTDEvsR = -20, 
	    float _maxTDEvsR = 40,// in ns
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
	    float _dtimeminEvsR = -100, 
	    float _dtimemaxEvsR = 200,// in ns
	    int _lynbin = 100,
	    float _lymin = 0,
	    float _lymax = 110, // in p.e.
	    int _distnbin = 100,
	    float _distmin = 0, 
	    float _distmax = 500, // in m ??????
	    float _rmin = 5,
	    float _rmax = 30,      
      float _nsteps = 10   // in m
)
{
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);

  char stmp[120];
  char sleg[120];
  char stit[120];
  snprintf(stit,sizeof stit,"2020 MC atmospheric muons");

  TString fout = "./output/";

  TFile* outputFile = new TFile(fout + "data/mcMuonFile.root","recreate");

  TH1F* hNmuon = new TH1F("hNmuon","Number of muons in event",_nmuonnbin,_nmuonmin,_nmuonmax);
  TH1F* htruetime = new TH1F("htruetime","MC time of first muon",_timenbin,_timemin,_timemax);
  TH1F* hreftime = new TH1F("hreftime","Time at reco track reference point",_timenbin,_timemin,_timemax);
  TH1F* htheta = new TH1F("htheta","Polar angle of reco muon track",_thetanbin,_thetamin,_thetamax);
  TH1F* hly = new TH1F("hly","Pulse LY",_lynbin,_lymin,_lymax);
  TH1F* htimes = new TH1F("htimes","Pulse times",_timenbin,_timemin,_timemax);
  TH1F* hdt = new TH1F("hdt","OM time - ref time",_dtimenbin,_dtimemin,_dtimemax);
  TH1F* hdtEvsR = new TH1F("hdt","OM time - ref time",_dtimenbin,_dtimeminEvsR,_dtimemaxEvsR);
  TH1F* hTrackDistSigOM = new TH1F("TrackDistSigOM ", "number of fired OM vs dist from track", _nsteps, _rmin, _rmax);
  TH1F* hTrackDistTotalOM = new TH1F("TrackDistTotalOM ", "number of total OM vs dist from track", _nsteps, _rmin, _rmax);

  TH2F* hLYvsTrackDist = new TH2F("hLYvsTrackDist","LY vs dist to OM",
				  _distnbin,_distmin,_distmax,_lynbin,_lymin,_lymax);
  TH2F* hLYvsTrackDistNM = new TH2F("hLYvsTrackDistNM","LY vs dist to OM, cut on Nmuon",
				     _distnbin,_distmin,_distmax,_lynbin,_lymin,_lymax);
  TH2F* hLYvsTrackDistPA = new TH2F("hLYvsTrackDistPA","LY vs dist to OM, cut on polar angle",
				     _distnbin,_distmin,_distmax,_lynbin,_lymin,_lymax);
  TH2F* hLYvsTrackDistTD = new TH2F("hLYvsTrackDistTD","LY vs dist to OM, cut on time dif",
				     _distnbin,_distmin,_distmax,_lynbin,_lymin,_lymax);
  TH2F* hLYvsTrackDistTDEvsR = new TH2F("hLYvsTrackDistTDEvsR","LY vs dist to OM, cut on Expected vs Reco time dif",
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

  /*Double_t r_x[_nsteps];
  Double_t r_y[_nsteps];
  Double_t r_ex[_nsteps];
  Double_t r_ey[_nsteps];

  for (int i = 0; i < _nsteps; ++i) {
	r_x[i] = _rmin + (i + 0.5) * _step;
        r_ex[i] = 0.5 * _rstep;
  }*/

  //---------- read file 
  int ifile=0;
  char tmp[100];
  int nRecoEvents= 0;
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

      //====================selection cuts===================================
      if (breco->GetNHits() < 8) continue;
      if (breco->GetNStrings() < 2) continue;
      if (breco->GetCovMatrixStatus() != 3) continue;
      if (breco->GetThetaRec() <= 100) continue;
      if (breco->GetZDist() < 200) continue;

      // if (breco->GetDEDX_energy() >= 3) continue;
      // if (bmcev->GetMuonsN() > 3) continue;
      //if (breco->GetClassBDT() > 0.25) continue;
      //if (breco->GetClassBDTLowE() > 0.25) continue;
      nRecoEvents++;
      //=====================================================================
      

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
        // if ( pulseLY <= 5) continue;

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

        Double_t distToPoint_BH = BHelperFunctions::GetTrackDistanceToPoint(refPoint, recoVec, chanPos);
        Double_t ExpectedOMTime = BHelperFunctions::GetPropagationTime(refPoint, recoVec, chanPos);
        Double_t dTimeExpVsRec = dTime - ExpectedOMTime;

        // if ( dTimeExpVsRec < -20 || dTimeExpVsRec > 40) continue;

        hdtEvsR->Fill(dTimeExpVsRec,eventWeight);

        if ( dTimeExpVsRec >= _minTDEvsR && dTimeExpVsRec <= _maxTDEvsR && pulseLY >= _minLY) {
          hTrackDistSigOM->Fill(distToPoint_BH,eventWeight);
        }

        if ( maxDist < distToPoint_BH ) maxDist = distToPoint_BH;
            
        hLYvsTrackDist->Fill(distToPoint_BH,pulseLY,eventWeight);
        if ( nMuons >= _minNM && nMuons <= _maxNM )
          hLYvsTrackDistNM->Fill(distToPoint_BH,pulseLY,eventWeight);
        if ( theta >= _minPA && theta <= _maxPA )
          hLYvsTrackDistPA->Fill(distToPoint_BH,pulseLY,eventWeight);
        if ( dTime >= _minTD && dTime <= _maxTD )
          hLYvsTrackDistTD->Fill(distToPoint_BH,pulseLY,eventWeight);
        if ( dTimeExpVsRec >= _minTDEvsR && dTimeExpVsRec <= _maxTDEvsR)
          hLYvsTrackDistTDEvsR->Fill(distToPoint_BH,pulseLY,eventWeight);
        if ( pulseLY >= _minLY )
          hLYvsTrackDistLY->Fill(distToPoint_BH,pulseLY,eventWeight);
      }

      for (int channel = 0; channel < bgeomtel->GetNumOMs(); channel++){

        TVector3 chanPos = TVector3(bgeomtel->At(channel)->GetX(),
		                    bgeomtel->At(channel)->GetY(),
		                    bgeomtel->At(channel)->GetZ()); 
        Double_t distToPoint_BH = BHelperFunctions::GetTrackDistanceToPoint(refPoint, recoVec, chanPos);
        hTrackDistTotalOM->Fill(distToPoint_BH,eventWeight);
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
  cnmuon->SaveAs(fout + "figures/" + "nmuons.pdf");

  if ( gROOT->GetListOfCanvases()->FindObject("cly") == NULL )
    cly = new TCanvas("cly","LY", 10, 510, 400, 400);
  cly->cd();
  cly->SetLogy(1);
  snprintf(stmp,sizeof stmp,"%s, N = %d",stit,int(hly->GetEntries()));  
  hly->SetTitle(stmp);
  hly->GetXaxis()->SetTitle("LY in OMs [p.e.]");
  hly->GetYaxis()->SetTitle("entries");
  hly->DrawCopy();
  cly->SaveAs(fout + "figures/"+ "LY.pdf");

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
  ctt->SaveAs(fout + "figures/" + "true_time.pdf");

  if ( gROOT->GetListOfCanvases()->FindObject("ct") == NULL )
    ct = new TCanvas("ct","Time", 510, 310, 400, 400);
  ct->cd();
  snprintf(stmp,sizeof stmp,"%s, N = %d",stit,int(htimes->GetEntries()));  
  htimes->SetTitle(stmp);
  htimes->GetXaxis()->SetTitle("OM time [ns]");
  htimes->GetYaxis()->SetTitle("entries");
  htimes->DrawCopy();
  ct->SaveAs(fout + "figures/" + "puls_time.pdf");

  if ( gROOT->GetListOfCanvases()->FindObject("cdt") == NULL )
    cdt = new TCanvas("cdt","Time dif", 510, 610, 400, 400);
  cdt->cd();
  snprintf(stmp,sizeof stmp,"%s, N = %d",stit,int(hdt->GetEntries()));  
  hdt->SetTitle(stmp);
  hdt->GetXaxis()->SetTitle("OM time - Ref time [ns]");
  hdt->GetYaxis()->SetTitle("entries");
  hdt->DrawCopy();
  cdt->SaveAs(fout + "figures/" + "rel_time.pdf");

  if ( gROOT->GetListOfCanvases()->FindObject("cdtevsr") == NULL )
    cdtevsr = new TCanvas("cdtevsr","Reco vs Expected Time dif", 510, 610, 400, 400);
  cdtevsr->cd();
  snprintf(stmp,sizeof stmp,"%s, N = %d",stit,int(hdtEvsR->GetEntries()));  
  hdtEvsR->SetTitle(stmp);
  hdtEvsR->GetXaxis()->SetTitle("Expected dtime - Reco dtime [ns]");
  hdtEvsR->GetYaxis()->SetTitle("entries");
  hdtEvsR->DrawCopy();
  cdtevsr->SaveAs(fout + "figures/" +  "evsr_dtime.pdf");

  if ( gROOT->GetListOfCanvases()->FindObject("cnomt") == NULL )
    cnomt = new TCanvas("cnomt","Number of total OM vs dist from track", 510, 610, 400, 400);
  cnomt->cd();
  TLegend* legnom = new TLegend(0.5,0.7,0.85,0.85);
  legnom->SetTextSize(0.045);
  hTrackDistSigOM->Scale(1.0/nRecoEvents);
  hTrackDistTotalOM->Scale(1.0/nRecoEvents);
  TH1F *hTrackDistZeroOM = (TH1F*) hTrackDistTotalOM->Clone("hTrackDistZeroOM");
  hTrackDistZeroOM->SetTitle("number of zero OM vs dist from track");
  hTrackDistZeroOM->Add(hTrackDistSigOM, -1.0);
  snprintf(stmp,sizeof stmp,"%s, N = %d",stit,int(hTrackDistSigOM->GetEntries()));  
  hTrackDistTotalOM->SetTitle(stmp);
  //hTrackDistTotalOM->SetMinimum(0);
  hTrackDistTotalOM->GetXaxis()->SetTitle("OM dist, m");
  hTrackDistTotalOM->GetYaxis()->SetTitle("hits");
  hTrackDistTotalOM->DrawCopy();
  snprintf(sleg, sizeof sleg,"total OM");
  legnom->AddEntry(hTrackDistTotalOM,sleg,"l");
  hTrackDistZeroOM->SetLineColor(2);
  hTrackDistZeroOM->SetMarkerColor(2);
  hTrackDistZeroOM->DrawCopy("same");
  snprintf(sleg, sizeof sleg,"zero OM");
  legnom->AddEntry(hTrackDistZeroOM,sleg,"l");
  legnom->Draw("same");
  cnomt->Update();
  cnomt->SaveAs(fout + "figures/" + "totalnom_vs_dist.pdf");

  if ( gROOT->GetListOfCanvases()->FindObject("cnoms") == NULL )
    cnoms = new TCanvas("cnoms","Number of signal OM vs dist from track", 510, 610, 400, 400);
  cnoms->cd();
  snprintf(stmp,sizeof stmp,"%s, N = %d",stit,int(hTrackDistSigOM->GetEntries()));  
  hTrackDistSigOM->SetTitle(stmp);
  hTrackDistSigOM->GetXaxis()->SetTitle("OM dist, m");
  hTrackDistSigOM->GetYaxis()->SetTitle("hits");
  hTrackDistSigOM->DrawCopy();
  cnoms->SaveAs(fout + "figures/" + "signom_vs_dist.pdf");

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


  if ( gROOT->GetListOfCanvases()->FindObject("c2cdtsevsr") == NULL )
    c2cdtsevsr = new TCanvas("c2cdtsevsr","LY vs Dist, reco vs expected time dif sel", 1010, 610, 600, 400);
  c2cdtsevsr->cd();
  snprintf(stmp,sizeof stmp,"%s, %4.0f #leq Exp vs Reco dt #leq %4.0f ns, N = %d (%5.3f)",
	   stit,_minTDEvsR, _maxTDEvsR, int(hLYvsTrackDistTDEvsR->GetEntries()),
	   float(hLYvsTrackDistTDEvsR->GetEntries())/float(hLYvsTrackDistTDEvsR->GetEntries()));  
  hLYvsTrackDistTDEvsR->SetTitle(stmp);
  hLYvsTrackDistTDEvsR->GetXaxis()->SetTitle("Distance from reco track to OM [m]");
  hLYvsTrackDistTDEvsR->GetYaxis()->SetTitle("LY [p.e.]");
  hLYvsTrackDistTDEvsR->DrawCopy("colz");

  outputFile->cd();
  hNmuon->Write();
  htruetime->Write();
  htheta->Write();
  hreftime->Write();
  hly->Write();
  htimes->Write();
  hdt->Write();
  hdtEvsR->Write();
  hLYvsTrackDist->Write();
  hLYvsTrackDistNM->Write();
  hLYvsTrackDistPA->Write();
  hLYvsTrackDistTD->Write();
  hLYvsTrackDistTDEvsR->Write();
  hLYvsTrackDistLY->Write();
  outputFile->Close();
}


