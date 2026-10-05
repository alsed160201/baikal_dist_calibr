#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h" 
#include "TF1.h"  
#include "TROOT.h"
#include "TVector3.h"
#include "TMath.h"

#include "BMCEvent.h"
#include "BEvent.h"
#include "BRecoMuon.h"
#include "BGeomTel.h"
#include "BHelperFunctions.h"

#include "../helpers/help_functions.C"

//int NCLUSTER = 2;
TCanvas* cphe;    // Number of ph. e.
TCanvas* cprof;    // Number of ph. e. via signal profile

void OMdistFiter(std::string _filelist,
	    float _rmin = 5,
	    float _rmax = 30,
      float _nsteps = 10,   // in m
      int _fit_scenario = 2
)
{
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);

  char stmp[120];
  char sleg[120];
  char stit[120];
  snprintf(stit,sizeof stit,"2020 MC atmospheric muons");

  TString figures_out = "./output/figures/fits/";
  TString data_out = "./output/data/";

  if (!EnsureDirectoryExists(figures_out)) {
        return; // Exits macro and stops the program
    }

  if (!EnsureDirectoryExists(data_out)) {
      return; // Exits macro and stops the program
  }

  TFile* outputFile = new TFile(data_out + "OMdistProfileFile.root","recreate");

  TH1F* hTrackDistSigOM = new TH1F("TrackDistSigOM ", "number of fired OM vs dist from track", _nsteps, _rmin, _rmax);
  hTrackDistSigOM->GetXaxis()->SetTitle("OM dist, m");
  TH1F* hTrackDistTotalOM = new TH1F("TrackDistTotalOM ", "number of total OM vs dist from track", _nsteps, _rmin, _rmax);
  hTrackDistTotalOM->GetXaxis()->SetTitle("OM dist, m");
  TProfile* hprof  = new TProfile("hprof","Profile of ph.e. signal versus track dist", _nsteps, _rmin, _rmax);
  hprof->GetXaxis()->SetTitle("OM dist, m");

  //====================================Fit function definition===================================

  Double_t WR = BHelperFunctions::GetWRefraction();
  Double_t cos_c = 1/WR;
  Double_t sin_c = sqrt(1 - cos_c*cos_c);

  // declared outside the branches so it stays visible for the fits below;
  // each scenario sets up its own formula, parameter names, start values and fixed parameters
  TF1 *fitfunction = nullptr;
  if (_fit_scenario == 1){
    fitfunction = new TF1("fitfunction", "([0]/x)*exp(-x/([1]*[2]))", _rmin, _rmax);
    fitfunction->SetParNames("A", "Lambda", "sin_c");
    fitfunction->SetParameters(8, 20);
    fitfunction->FixParameter(2, sin_c);
    fitfunction->SetParLimits(1, 0.1, 200);   // absorption length must stay positive
  }
  else if (_fit_scenario == 2){
    fitfunction = new TF1("fitfunction", "([0])*exp(-x/([1]*[2]))", _rmin, _rmax);
    fitfunction->SetParNames("A", "Lambda", "sin_c");
    fitfunction->SetParameters(8, 20);
    fitfunction->FixParameter(2, sin_c);
    fitfunction->SetParLimits(1, 0.1, 200);   // absorption length must stay positive
  }
  else if (_fit_scenario == 3){
    // scenario 1 plus a constant pedestal B (e.g. noise hits at large distances)
    fitfunction = new TF1("fitfunction", "([0]/x)*exp(-x/([1]*[2])) + [3]", _rmin, _rmax);
    fitfunction->SetParNames("A", "Lambda", "sin_c", "B");
    fitfunction->SetParameters(8, 20, 0, 0.01);
    fitfunction->FixParameter(2, sin_c);
    fitfunction->SetParLimits(1, 0.1, 200);   // absorption length must stay positive
    fitfunction->SetParLimits(3, 0, 1);   // pedestal can't be negative
  }
  else if (_fit_scenario == 4){
    // scenario 2 plus a constant pedestal B (e.g. noise hits at large distances)
    fitfunction = new TF1("fitfunction", "[0]*exp(-x/([1]*[2])) + [3]", _rmin, _rmax);
    fitfunction->SetParNames("A", "Lambda", "sin_c", "B");
    fitfunction->SetParameters(8, 20, 0, 0.01);
    fitfunction->FixParameter(2, sin_c);
    fitfunction->SetParLimits(1, 0.1, 200);   // absorption length must stay positive
    fitfunction->SetParLimits(3, 0, 1);   // pedestal can't be negative
  }

  else{
    std::cout << "[ROOT] Fit function is not defined" << std::endl;
    return; // Exits macro and stops the program
  }

  // start values, to reset the function before each fit
  std::vector<Double_t> initPars(fitfunction->GetParameters(),
                                 fitfunction->GetParameters() + fitfunction->GetNpar());


  //==============================================================================================


  //---------- read file---------------
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
      if (breco->GetNHits() > 10) continue;
      if (breco->GetNStrings() > 2) continue;
      if (breco->GetCovMatrixStatus() != 3) continue;
      if (breco->GetThetaRec() <= 160) continue;
      if (breco->GetZDist() < 400) continue;

      if (breco->GetDEDX_energy() >= 3) continue;
      // if (bmcev->GetMuonsN() != 1) continue;
      //if (breco->GetClassBDT() > 0.25) continue;
      //if (breco->GetClassBDTLowE() > 0.25) continue;
      nRecoEvents++;
      //=====================================================================

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

      // signal amplitude (ph.e.) per channel, for channels with a pulse that
      // passes the time cut below; looked up in the all-OM loop further down
      // so hprof averages over every OM (0 for OMs that didn't fire), not just
      // the fired ones
      std::unordered_map<int, float> chanSignal;

      //loop over bevent pulses (fired OM's)
      for (int ipulse = 0; ipulse < bevt->NHits(); ipulse++){

        float pulseLY = bevt->Q(ipulse);
        //  if ( pulseLY <= 3) continue;

        int chanID = bevt->HitChannel(ipulse);
        Float_t pulseTime = bevt->GetImpulse(ipulse)->GetTime();         //ns
        float dTime = pulseTime-refTime;
        TVector3 chanPos = TVector3(bgeomtel->At(chanID)->GetX(),
                                          bgeomtel->At(chanID)->GetY(),
                                          bgeomtel->At(chanID)->GetZ());
        Double_t distToPoint_BH = BHelperFunctions::GetTrackDistanceToPoint(refPoint, recoVec, chanPos);
        Double_t ExpectedOMTime = BHelperFunctions::GetPropagationTime(refPoint, recoVec, chanPos);
        Double_t OMlightAngle = 180*BHelperFunctions::GetOMlightAngle(refPoint, recoVec, chanPos)/TMath::Pi();
        Double_t dTimeExpVsRec = dTime - ExpectedOMTime;

        if ( dTimeExpVsRec < -20 || dTimeExpVsRec > 40) continue;


        if (chanSignal.find(chanID) == chanSignal.end()){
          hTrackDistSigOM->Fill(distToPoint_BH,eventWeight);
          chanSignal[chanID] = pulseLY;
        }
        else{
          chanSignal[chanID] += pulseLY;
        }
        
      }

      for (int channel = 0; channel < bgeomtel->GetNumOMs(); channel++){

        TVector3 chanPos = TVector3(bgeomtel->At(channel)->GetX(),
		                    bgeomtel->At(channel)->GetY(),
		                    bgeomtel->At(channel)->GetZ());
        Double_t distToPoint_BH = BHelperFunctions::GetTrackDistanceToPoint(refPoint, recoVec, chanPos);
        hTrackDistTotalOM->Fill(distToPoint_BH,eventWeight);

        // 0 ph.e. for OMs that didn't fire (or whose pulse failed the time cut above),
        // so hprof averages over every OM at this distance, not only the fired ones
        auto chanSignalIt = chanSignal.find(channel);
        float channelSignal = (chanSignalIt != chanSignal.end()) ? chanSignalIt->second : 0.0f;
        hprof->Fill(distToPoint_BH, channelSignal, eventWeight);
      }

    }
  }

  //-------------OM hist magic------------------------
  hTrackDistSigOM->Scale(1.0/nRecoEvents);
  hTrackDistTotalOM->Scale(1.0/nRecoEvents);

  TH1F *hTrackDistZeroOM = (TH1F*) hTrackDistTotalOM->Clone("hTrackDistZeroOM");
  hTrackDistZeroOM->SetTitle("number of zero OM vs dist from track");
  hTrackDistZeroOM->Add(hTrackDistSigOM, -1.0);
  SetPoissonErrors(hTrackDistZeroOM, nRecoEvents);

  // old variant: Zero/Total without "B" (treats them as independent, overestimates the errors)
  // TH1F *hTrackDistZeroOMfrac = (TH1F*) hTrackDistZeroOM->Clone("hTrackDistZeroOMfrac");
  // hTrackDistZeroOMfrac->SetTitle("fraction of zero OM vs dist from track");
  // hTrackDistZeroOMfrac->Divide(hTrackDistTotalOM);
  TH1F *hTrackDistZeroOMfrac = MakeZeroFraction(hTrackDistSigOM, hTrackDistTotalOM, "hTrackDistZeroOMfrac");

  TGraphErrors *pheGraph = MakeLogGraph(hTrackDistZeroOMfrac);


  //------------- fit and save ----------------------

  if ( gROOT->GetListOfCanvases()->FindObject("cphe") == NULL )
    cphe = new TCanvas("cphe","ph. e. estimation", 510, 610, 400, 400);
  cphe->cd();
  // cphe->SetLogy(1);
  pheGraph->Draw("AP");
  fitfunction->SetParameters(initPars.data());
  pheGraph->Fit(fitfunction);
  // pheGraph->Fit(fitfunction, "EX0");
  cphe->SaveAs(figures_out+"fit_estimation.pdf");


  if ( gROOT->GetListOfCanvases()->FindObject("cprof") == NULL )
    cprof = new TCanvas("cprof","ph. e. estimation (via signal profile)", 510, 610, 400, 400);
  cprof->cd();
  // cprof->SetLogy(1);
  hprof->Draw("PE");
  fitfunction->SetParameters(initPars.data());
  hprof->Fit(fitfunction);
  cprof->SaveAs(figures_out+"fit_profile_estimation.pdf");

  //------------write into file------------------
  outputFile->cd();
  hTrackDistTotalOM->Write();
  hTrackDistSigOM->Write();
  hTrackDistZeroOM->Write();
  hTrackDistZeroOMfrac->Write();
  pheGraph->Write();
  hprof->Write();
  outputFile->Close();
}
