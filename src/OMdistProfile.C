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

#include "../helpers/help_functions.C"

//int NCLUSTER = 2;
TCanvas* cnoms;    // Number of fired OM vs dist from track
TCanvas* cnomt;    // Number of total OM vs dist from track
TCanvas* cphe;    // Number of ph. e.
TCanvas* cpheGroups;    // Number of ph. e. split by groups
TCanvas* cpheCuts;    // Number of ph. e. split by pulseLY cut level

void OMdistProfile(std::string _filelist,
	    float _rmin = 2.5,
	    float _rmax = 40,
      float _nsteps = 15,   // in m
      bool _doGroupPlots = true,   // build hTrackDistSigOM_g/hTrackDistTotalOM_g and the per-group ph.e. plot
      bool _doCutPlots = true      // build hTrackDistSigOM_c and the per-pulseLY-cut ph.e. plot
)
{
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);

  char stmp[120];
  char sleg[120];
  char stit[120];
  snprintf(stit,sizeof stit,"2020 MC atmospheric muons");

  TString fout = "./output/figures/";

  TFile* outputFile = new TFile("./output/data/OMdistProfileFile.root","recreate");

  TH1F* hTrackDistSigOM = new TH1F("TrackDistSigOM ", "number of fired OM vs dist from track", _nsteps, _rmin, _rmax);
  hTrackDistSigOM->GetXaxis()->SetTitle("OM dist, m");
  TH1F* hTrackDistTotalOM = new TH1F("TrackDistTotalOM ", "number of total OM vs dist from track", _nsteps, _rmin, _rmax);
  hTrackDistTotalOM->GetXaxis()->SetTitle("OM dist, m");

  //=========== grouping configuration ===========
  // holds every per-event object grouping variables may read from, so
  // GetGroupingVariable always has the same signature no matter which
  // object the chosen variable actually comes from.
  struct EventRefs {
    BMCEvent* bmcev;
    BRecoMuon* breco;
    BGeomTel* bgeomtel;
  };

  // variable to split events by, and its bin edges (like TH1 bin edges):
  // edges = {e0, e1, ..., eN} define N groups [e0,e1), [e1,e2), ..., [eN-1,eN)

  //========================= example 1: split by number of muons =========================
  // TString groupVarName = "N_{mc #mu}";
  // std::vector<Double_t> groupEdges = {0, 3, 10, 50, 100};
  // auto GetGroupingVariable = [](const EventRefs& ev) -> Double_t {
  //   return ev.bmcev->GetMuonsN();
  // };

  //========================= example 2: split by number of hitted strings =========================
  // TString groupVarName = "N_{strings}";
  // std::vector<Double_t> groupEdges = {2, 3, 4, 5, 9};
  // auto GetGroupingVariable = [](const EventRefs& ev) -> Double_t {
  //   return ev.breco->GetNStrings();
  // };

  //========================= example 3: split by reconstructed reconstructed theta angle =========================
  // TString groupVarName = "#theta rec, deg";
  // std::vector<Double_t> groupEdges = {120, 140, 160, 181};
  // auto GetGroupingVariable = [](const EventRefs& ev) -> Double_t {
  //   return ev.breco->GetThetaRec();
  // };

    //========================= example 4: split by reconstructed reconstructed source zenith angle =========================
  TString groupVarName = "zenith angle";
  std::vector<Double_t> groupEdges = {0, 20, 40, 81};
  auto GetGroupingVariable = [](const EventRefs& ev) -> Double_t {
    return ev.breco->GetSourceZenithRec();
  };


  //================================================

  const int nGroups = (int)groupEdges.size() - 1;
  std::vector<Int_t> groupColors = {kBlack, kRed, kBlue, kGreen+2, kMagenta+1, kOrange+7, kCyan+2};

  std::vector<TString> groupLabels(nGroups);
  for (int g=0; g<nGroups; g++){
    groupLabels[g] = TString::Format("%s = [%g,%g)", groupVarName.Data(), groupEdges[g], groupEdges[g+1]);
  }

  // returns the group index for val, or -1 if val falls outside [edges.front(), edges.back())
  auto GetGroupIndex = [&](Double_t val) -> int {
    if (val < groupEdges.front() || val >= groupEdges.back()) return -1;
    for (int g=0; g<nGroups; g++){
      if (val >= groupEdges[g] && val < groupEdges[g+1]) return g;
    }
    return -1;
  };

  std::vector<TH1F*> hTrackDistSigOM_g(nGroups);
  std::vector<TH1F*> hTrackDistTotalOM_g(nGroups);
  std::vector<int> nRecoEvents_g(nGroups, 0);

  if (_doGroupPlots){
    for (int g=0; g<nGroups; g++){
      hTrackDistSigOM_g[g] = new TH1F(Form("TrackDistSigOM_g%d",g), "number of fired OM vs dist from track", _nsteps, _rmin, _rmax);
      hTrackDistSigOM_g[g]->GetXaxis()->SetTitle("OM dist, m");
      hTrackDistTotalOM_g[g] = new TH1F(Form("TrackDistTotalOM_g%d",g), "number of total OM vs dist from track", _nsteps, _rmin, _rmax);
      hTrackDistTotalOM_g[g]->GetXaxis()->SetTitle("OM dist, m");
    }
  }

  //=========== pulseLY cut-level configuration ===========
  // pulseLY is a per-pulse quantity (not per-event), so it can't go through
  // GetGroupingVariable/EventRefs; instead, each threshold gets its own
  // "signal OM" histogram, filled with hits passing pulseLY > threshold.
  // hTrackDistTotalOM (denominator) doesn't depend on pulseLY, so it's shared.
  TString cutVarName = "Q_{cut}, p.e.";
  std::vector<Double_t> pulseLYCuts = {0, 3, 5, 7, 10};

  const int nCuts = (int)pulseLYCuts.size();
  std::vector<Int_t> cutColors = {kBlack, kRed, kBlue, kGreen+2, kMagenta+1, kOrange+7, kCyan+2};

  std::vector<TString> cutLabels(nCuts);
  for (int c=0; c<nCuts; c++){
    cutLabels[c] = TString::Format("%s > %g", cutVarName.Data(), pulseLYCuts[c]);
  }

  std::vector<TH1F*> hTrackDistSigOM_c(nCuts);
  if (_doCutPlots){
    for (int c=0; c<nCuts; c++){
      hTrackDistSigOM_c[c] = new TH1F(Form("TrackDistSigOM_c%d",c), "number of fired OM vs dist from track", _nsteps, _rmin, _rmax);
      hTrackDistSigOM_c[c]->GetXaxis()->SetTitle("OM dist, m");
    }
  }

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

      int groupIdx = GetGroupIndex(GetGroupingVariable(EventRefs{bmcev, breco, bgeomtel}));
      if (groupIdx >= 0) nRecoEvents_g[groupIdx]++;

      //---- mc info ----------------------
      //number of muons, which produced response in the detector 
      Int_t nMuons = bmcev->GetResponseMuonsN();

      // true time of first muon
      Double_t trueTime = bmcev->GetFirstMuonTime();

      //---- reco info ---------------------
      //Calculate distance from reconstructed track to channels
      Double_t theta = breco->GetThetaRec(); //_clHM();
      Double_t phi = breco->GetPhiRec(); //_clHM();
      double refTime = breco->GetTimeXYZRec();

      Float_t trackThetaRad = TMath::Pi()*(theta)/180;
      Float_t trackThetaGrad = theta;
      Float_t trackPhiRad = TMath::Pi()*(phi)/180;
      
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
        TVector3 chanPos = TVector3(bgeomtel->At(chanID)->GetX(),
                                          bgeomtel->At(chanID)->GetY(),
                                          bgeomtel->At(chanID)->GetZ());
        Double_t distToPoint_BH = BHelperFunctions::GetTrackDistanceToPoint(refPoint, recoVec, chanPos);
        Double_t ExpectedOMTime = BHelperFunctions::GetPropagationTime(refPoint, recoVec, chanPos);
        Double_t dTimeExpVsRec = dTime - ExpectedOMTime;

        if ( dTimeExpVsRec < -20 || dTimeExpVsRec > 40) continue;

        hTrackDistSigOM->Fill(distToPoint_BH,eventWeight);
        if (_doGroupPlots && groupIdx >= 0) hTrackDistSigOM_g[groupIdx]->Fill(distToPoint_BH,eventWeight);

        if (_doCutPlots){
          for (int c=0; c<nCuts; c++){
            if (pulseLY > pulseLYCuts[c]) hTrackDistSigOM_c[c]->Fill(distToPoint_BH,eventWeight);
          }
        }

      }

      for (int channel = 0; channel < bgeomtel->GetNumOMs(); channel++){

        TVector3 chanPos = TVector3(bgeomtel->At(channel)->GetX(),
		                    bgeomtel->At(channel)->GetY(),
		                    bgeomtel->At(channel)->GetZ());
        Double_t distToPoint_BH = BHelperFunctions::GetTrackDistanceToPoint(refPoint, recoVec, chanPos);
        hTrackDistTotalOM->Fill(distToPoint_BH,eventWeight);
        if (_doGroupPlots && groupIdx >= 0) hTrackDistTotalOM_g[groupIdx]->Fill(distToPoint_BH,eventWeight);
      }
  
    }
  }

  //-------------OM hist magic------------------------
  hTrackDistSigOM->Scale(1.0/nRecoEvents);
  hTrackDistTotalOM->Scale(1.0/nRecoEvents);

  TH1F *hTrackDistZeroOM = (TH1F*) hTrackDistTotalOM->Clone("hTrackDistZeroOM");
  hTrackDistZeroOM->SetTitle("number of zero OM vs dist from track");
  hTrackDistZeroOM->Add(hTrackDistSigOM, -1.0);

  TH1F *hTrackDistZeroOMfrac = (TH1F*) hTrackDistZeroOM->Clone("hTrackDistZeroOMfrac");
  hTrackDistZeroOMfrac->SetTitle("fraction of zero OM vs dist from track");
  hTrackDistZeroOMfrac->Divide(hTrackDistTotalOM);

  TGraphErrors *pheGraph = MakeLogGraph(hTrackDistZeroOMfrac);

  //-------------OM hist magic per group------------------------
  std::vector<TH1F*> hTrackDistZeroOM_g(nGroups);
  std::vector<TH1F*> hTrackDistZeroOMfrac_g(nGroups);
  std::vector<TGraphErrors*> pheGraph_g(nGroups);

  if (_doGroupPlots){
    for (int g=0; g<nGroups; g++){
      if (nRecoEvents_g[g] > 0){
        hTrackDistSigOM_g[g]->Scale(1.0/nRecoEvents_g[g]);
        hTrackDistTotalOM_g[g]->Scale(1.0/nRecoEvents_g[g]);
      }

      hTrackDistZeroOM_g[g] = (TH1F*) hTrackDistTotalOM_g[g]->Clone(Form("hTrackDistZeroOM_g%d",g));
      hTrackDistZeroOM_g[g]->SetTitle("number of zero OM vs dist from track");
      hTrackDistZeroOM_g[g]->Add(hTrackDistSigOM_g[g], -1.0);

      hTrackDistZeroOMfrac_g[g] = (TH1F*) hTrackDistZeroOM_g[g]->Clone(Form("hTrackDistZeroOMfrac_g%d",g));
      hTrackDistZeroOMfrac_g[g]->SetTitle("fraction of zero OM vs dist from track");
      hTrackDistZeroOMfrac_g[g]->Divide(hTrackDistTotalOM_g[g]);

      pheGraph_g[g] = MakeLogGraph(hTrackDistZeroOMfrac_g[g]);
      pheGraph_g[g]->SetName(Form("g_ln_hTrackDistZeroOMfrac_g%d",g));
      pheGraph_g[g]->SetTitle(groupLabels[g]+";OM dist, m;proxy ph. e.");
      pheGraph_g[g]->SetLineColor(groupColors[g % groupColors.size()]);
      pheGraph_g[g]->SetMarkerColor(groupColors[g % groupColors.size()]);
      pheGraph_g[g]->SetMarkerStyle(20+g);
    }
  }

  //-------------OM hist magic per pulseLY cut level------------------------
  std::vector<TH1F*> hTrackDistZeroOM_c(nCuts);
  std::vector<TH1F*> hTrackDistZeroOMfrac_c(nCuts);
  std::vector<TGraphErrors*> pheGraph_c(nCuts);

  if (_doCutPlots){
    for (int c=0; c<nCuts; c++){
      hTrackDistSigOM_c[c]->Scale(1.0/nRecoEvents);

      hTrackDistZeroOM_c[c] = (TH1F*) hTrackDistTotalOM->Clone(Form("hTrackDistZeroOM_c%d",c));
      hTrackDistZeroOM_c[c]->SetTitle("number of zero OM vs dist from track");
      hTrackDistZeroOM_c[c]->Add(hTrackDistSigOM_c[c], -1.0);

      hTrackDistZeroOMfrac_c[c] = (TH1F*) hTrackDistZeroOM_c[c]->Clone(Form("hTrackDistZeroOMfrac_c%d",c));
      hTrackDistZeroOMfrac_c[c]->SetTitle("fraction of zero OM vs dist from track");
      hTrackDistZeroOMfrac_c[c]->Divide(hTrackDistTotalOM);

      pheGraph_c[c] = MakeLogGraph(hTrackDistZeroOMfrac_c[c]);
      pheGraph_c[c]->SetName(Form("g_ln_hTrackDistZeroOMfrac_c%d",c));
      pheGraph_c[c]->SetTitle(cutLabels[c]+";OM dist, m;proxy ph. e.");
      pheGraph_c[c]->SetLineColor(cutColors[c % cutColors.size()]);
      pheGraph_c[c]->SetMarkerColor(cutColors[c % cutColors.size()]);
      pheGraph_c[c]->SetMarkerStyle(20+c);
    }
  }

  //------------- plot to canvas----------------------
  if ( gROOT->GetListOfCanvases()->FindObject("cnomt") == NULL )
    cnomt = new TCanvas("cnomt","Number of total OM vs dist from track", 510, 610, 400, 400);
  cnomt->cd();
  TLegend* legnom = new TLegend(0.2,0.7,0.55,0.85);
  legnom->SetTextSize(0.045);

  snprintf(stmp,sizeof stmp,"%s, N = %d",stit,int(hTrackDistSigOM->GetEntries()));  
  hTrackDistTotalOM->SetTitle(stmp);
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
  cnomt->SaveAs(fout+"totalnom_vs_dist.pdf");


  if ( gROOT->GetListOfCanvases()->FindObject("cnoms") == NULL )
    cnoms = new TCanvas("cnoms","Number of signal OM vs dist from track", 510, 610, 400, 400);
  cnoms->cd();
  snprintf(stmp,sizeof stmp,"%s, N = %d",stit,int(hTrackDistSigOM->GetEntries()));  
  hTrackDistSigOM->SetTitle(stmp);
  hTrackDistSigOM->GetXaxis()->SetTitle("OM dist, m");
  hTrackDistSigOM->GetYaxis()->SetTitle("hits");
  hTrackDistSigOM->DrawCopy();
  cnoms->SaveAs(fout+"signom_vs_dist.pdf");



  if ( gROOT->GetListOfCanvases()->FindObject("cphe") == NULL )
    cphe = new TCanvas("cphe","ph. e. estimation", 510, 610, 400, 400);
  cphe->cd();
  //cphe->SetLogy(1);
  pheGraph->Draw("AP");
  cphe->SaveAs(fout+"phe_estimation.pdf");



  if (_doGroupPlots){
    if ( gROOT->GetListOfCanvases()->FindObject("cpheGroups") == NULL )
      cpheGroups = new TCanvas("cpheGroups", TString("ph. e. estimation by "+groupVarName).Data(), 510, 610, 400, 400);
    cpheGroups->cd();
    //cpheGroups->SetLogy(1);
    TLegend* legphe = new TLegend(0.55,0.7,0.88,0.9);
    legphe->SetTextSize(0.035);
    // frame's y-range is taken only from the first ("AP") graph; SetRangeUser on
    // the already-created frame axis (unlike SetMaximum) reliably overrides it,
    // so points from the other groups ("P SAME") aren't clipped at the top
    pheGraph_g[0]->Draw("AP");
    pheGraph_g[0]->GetYaxis()->SetRangeUser(0.0, 1.0);
    legphe->AddEntry(pheGraph_g[0], groupLabels[0], "lp");
    for (int g=1; g<nGroups; g++){
      pheGraph_g[g]->Draw("P SAME");
      legphe->AddEntry(pheGraph_g[g], groupLabels[g], "lp");
    }
    if (pheGraph_g[0]->GetHistogram())
      pheGraph_g[0]->GetHistogram()->SetTitle(";OM dist, m;ph. e.");
    legphe->Draw("same");
    cpheGroups->Update();
    cpheGroups->SaveAs(fout+"phe_estimation_groups.pdf");
  }


  if (_doCutPlots){
    if ( gROOT->GetListOfCanvases()->FindObject("cpheCuts") == NULL )
      cpheCuts = new TCanvas("cpheCuts", TString("ph. e. estimation by "+cutVarName).Data(), 510, 610, 400, 400);
    cpheCuts->cd();
    //cpheCuts->SetLogy(1);
    TLegend* legpheCuts = new TLegend(0.55,0.7,0.88,0.9);
    legpheCuts->SetTextSize(0.035);
    pheGraph_c[0]->Draw("AP");
    pheGraph_c[0]->GetYaxis()->SetRangeUser(0.0, 1.0);
    legpheCuts->AddEntry(pheGraph_c[0], cutLabels[0], "lp");
    for (int c=1; c<nCuts; c++){
      pheGraph_c[c]->Draw("P SAME");
      legpheCuts->AddEntry(pheGraph_c[c], cutLabels[c], "lp");
    }
    if (pheGraph_c[0]->GetHistogram())
      pheGraph_c[0]->GetHistogram()->SetTitle(";OM dist, m;ph. e.");
    legpheCuts->Draw("same");
    cpheCuts->Update();
    cpheCuts->SaveAs(fout+"phe_estimation_pulseLY_cuts.pdf");
  }


  //------------write into file------------------
  outputFile->cd();
  hTrackDistTotalOM->Write();
  hTrackDistSigOM->Write();
  hTrackDistZeroOM->Write();
  hTrackDistZeroOMfrac->Write();
  pheGraph->Write();
  if (_doGroupPlots){
    for (int g=0; g<nGroups; g++){
      hTrackDistSigOM_g[g]->Write();
      hTrackDistTotalOM_g[g]->Write();
      hTrackDistZeroOM_g[g]->Write();
      hTrackDistZeroOMfrac_g[g]->Write();
      pheGraph_g[g]->Write();
    }
  }
  if (_doCutPlots){
    for (int c=0; c<nCuts; c++){
      hTrackDistSigOM_c[c]->Write();
      hTrackDistZeroOM_c[c]->Write();
      hTrackDistZeroOMfrac_c[c]->Write();
      pheGraph_c[c]->Write();
    }
  }
  outputFile->Close();
}


