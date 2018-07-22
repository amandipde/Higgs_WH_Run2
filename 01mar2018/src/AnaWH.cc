#include "configana.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "AnaWH.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::abs;
using std::max;
using std::sqrt;
using std::sort;

using namespace vhtm;

#define NEL(x) (sizeof((x))/sizeof((x)[0]))


// -----------
// Constructor
// -----------
AnaWH::AnaWH()
  : AnaBase()
{}
// ----------
// Destructor
// ----------
AnaWH::~AnaWH() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool AnaWH::beginJob() 
{ 
  AnaBase::beginJob();
  count_sig_event = 0;
  for (unsigned int i = 0; i < NEL(nEvtMuMuTau); i++) {
    nEvtMuMuTau[i] = 0;
  }
  for (unsigned int i = 0; i < NEL(nEvtMuTauTau); i++) {
    nEvtMuTauTau[i] = 0;
  }
  for (unsigned int i = 0; i < NEL(nSig); i++) {
    nSig[i] = 0;
  }
  for (unsigned int i = 0; i < NEL(nEvtMuMuTauTau); i++) {
    nEvtMuMuTauTau[i]= 0;
  }
  for (unsigned int i = 0; i < NEL(nEvtMuMuMuTau); i++) {
    nEvtMuMuMuTau[i]= 0;
  }
  // Open the output ROOT file
  _histf->cd();
  bookHistograms();

  return true;
}
// ---------------
// Book histograms
// ---------------
void AnaWH::bookHistograms() 
{
  static const double pi = TMath::Pi();
  new TH1F("nMuon", "Number of Muons", 6, -0.5, 5.5);
  new TH1F("muonPt1", "Highest Pt Muon Pt Distribution", 100, -0.5, 99.5);
  new TH1F("muonPt2", "Second highest Pt Muon Pt Distribution", 100, -0.5, 99.5);
  new TH1F("nTau", "Number of Taus", 4, -0.5, 3.5);
  new TH1F("tauPt1", "Highest Pt Tau Pt Distribution", 100, -0.5, 99.5);
  new TH1F("mumuInvMass", "Invariant mass of the two muons", 100, 0, 140);
  //new TH1F("deltaPt1","minimum Pt difference  between detector and generator level muon1", 100, 0, 100);
  //new TH1F("deltaPt2","minimum Pt difference  between detector and generator level muon2", 100, 0, 100);
  //new TH1F("deltaEta1","minimum Eta difference  between detector and generator level muon1", 100, 0, 100);
  //new TH1F("deltaEta2","minimum Eta difference  between detector and generator level muon2", 100, 0, 100);
  //new TH1F("deltaPhi1","minimum Phi difference  between detector and generator level muon1", 100, 0, 100);
  //new TH1F("deltaPhi2","minimum Phi difference  between detector and generator level muon2", 100, 0, 100);
  //new TH1F("deltaDR","minimum DR difference  between detector and generator level muon", 100, 0, 100);
  new TH1F("genInvMass", "Invariant mass of two muons in generator level irrespective their mother", 100, 0, 140);
  new TH1F("oddInvMass", "Inv Mass of the two muons in generator level requiring that the mothers of two muons are either W or Tau", 100, 0, 140);
  
  new TH1F("PhiMu1","phi angle of heighst pt muon", 100, -7, 7);
  new TH1F("PhiMu2","phi angle of 2nd heighst pt muon", 100, -7, 7);
  new TH1F("PhiTau","phi angle of heighst pt Tau", 100, -7, 7);

  //  new TH1F("DRMuMu","DR difference  between two highest pt muon", 100, 0, 10);
  new TH1F("DRMuTau1","DR difference  between 1st muon and tau", 100, 0, 10);
  new TH1F("DRMuTau2","DR difference  between 2nd muon and tau", 100, 0, 10);
  new TH1F("DEtaMuMu","Eta difference  between two highest pt muon", 100, -10, 10);
  new TH1F("DEtaMuTau1","Eta difference  between 1st muon and tau", 100, -10, 10);
  new TH1F("DEtaMuTau2","Eta difference  between 2nd muon and tau", 100, -10, 10);
  new TH1F("DPhiMuMu","Phi difference  between two highest pt muon", 100, -10, 10);
  new TH1F("DPhiMuTau1","Phi difference  between 1st muon and tau", 100, -10, 10);
  new TH1F("DPhiMuTau2","Phi difference  between 2nd muon and tau", 100, -10, 10);
  new TH1F("mva1", "new variable muon pt(mu1+mu2)/pt1+pt2", 100, 0, 2.0);
  new TH1F("mva_eta", "eta(mu1 +mu2)", 60, -6.0, 6.0);
  new TH1F("DRMuMu", "DR diff between two muons when mva1>0.98 for signal event", 100, -10, 10);
  new TH1F("mva_gen", "new variable muon pt(mu1+mu2)/pt1+pt2", 100, 0, 2.0);
  new TH2F("muon_correlation","correlation plot", 100, 0, 100, 100, 0, 100);
  new TH1F("deleta1", "eta diff between z and highest pt  muon for mva1>0.98", 100, -10, 10);
  new TH1F("deleta2", "eta diff between z and 2nd highest pt  muon for mva1>0.98", 100, -10, 10);
  new TH1F("Deta1", "eta diff between z and highest pt muon for mva1<0.5", 100, -10, 10);
  new TH1F("Deta2", "eta diff between z and 2nd highest pt muon for mva1<0.5", 100, -10, 10);
  new TH2F("dmuon_correlation", "correlation plot at the detector level", 100, 0, 100, 100, 0, 100);
  new TH2F("pt covariance", "plot of muon pt sum vs difference", 100, 0, 100, 100, 0, 100);
  new TH2F("eta_correlation", "eta correlation of two highest pt muons", 100, -5, 5, 100, -5, 5);
  new TH2F("P_correlation", "energy correlation of two highest pt muons", 100, 0, 100, 100, 0, 100);
  new TH1F("muoa charge", "charge of the highest pt muon", 100, -10, 10);
  new TH1F("muob charge", "charge of the 2nd highest pt muon", 100, -10, 10);
  new TH1F("tau charge", "charge of the tau", 100, -10, 10);
  new TH1F("WH_DR", "DR difference between W and Higgs at the generator level", 100, 0, 10);
  new TH1F("Ht1_DR","DR diff between the Higgs and the Hadronic Tau at the generator level", 100, 0, 10);
  new TH1F("HM2_DR","DR diff betwen the Higgs and tau decayed Muon at the generator level", 100, 0, 10);
  new TH1F("WM1_DR","DR diff between the W and first muon at the generator level", 100, 0, 10);
  new TH1F("t1M2_DR","DR diff between the Hadronic Tau and Tau decayed Muon at the generator level", 100, 0, 10);
  new TProfile("WH_Wpt","DR diff between W and H as a function of W's pt", 100, 0, 100, 0, 10);
  new TProfile("t1M2_Hpt","DR diff between Hadronic Tau and Tau decayed Muon as a function of Higgs's pt",100, 0, 100, 0, 10);
  new TProfile("WM1_Wpt", "DR diff between W boson and daughter Muon as a function of W's pt", 100, 0, 100, 0,10);
  new TH1F("DPhi_WH", "Phi diff between W and H", 100, -10, 10);
  new TH1F("DEta_WH", "Eta diff between W and H", 100, -10, 10);
  new TH1F("DPhi_t1M2","Phi diff between Hadronic Tau and Tau decayed Muon", 100, -10, 10);
  new TH1F("DEta_t1M2","Eta diff between Hadronic Tau and Tau decayed Muon", 100, -10, 10);
  new TH1F("mumuInvMass_ZH", "inv mass of two muons coming from Z decay for ZH process", 100, 0, 140);

 
}
// -------------------
// The main event loop
// -------------------
void AnaWH::clearLists() {
  vtxList.clear();
  muoList.clear();
  eleList.clear();
  tauList.clear();
  bjetList.clear();
  genMuonList.clear();
  genTauList.clear();
  GEN_MUON_LIST.clear();
  genMuoaList.clear();
  genMuobList.clear();
  genHList.clear();
  genWList.clear();
  genTau_H_List.clear();
  genMuon_Z_List.clear();
}
void AnaWH::eventLoop() 
{
  // Initialize analysis
  if (!beginJob()) return;
  

  int nPrint = max(10000, nEvents/1000);

  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  int nEvt = 0, nEvtVtx = 0;
  for (int ev = 0; ev < nEvents; ev++) {
    bool select = false;

    clearEvent();
    int lflag = _chain->LoadTree(ev); 
    int nentries = getEntry(lflag);    // returns bytes read
    ++nEvt;
    
    string currentFile(gSystem->BaseName(_chain->GetCurrentFile()->GetName())); 

    const Event* evt = dynamic_cast<Event*>(eventA->At(0));
    assert(evt);
    int run   = evt->run;
    int event = evt->event;

    // THIS PART HAS TO BE DELETED, ITS FOR CHECKING 
    _fLog << "run = " << run
          << ", event = " << event 
          << ", nMuon = " << n_muon
          << ", nElectron = " << n_electron
          << ", nTau = " << n_tau
          << ", nJet = " << n_jet
          << std::endl;
    _fLog << "muon block" << std::endl;
    for (int indx = 0; indx < n_muon; ++indx) {
      const Muon* muon = dynamic_cast<Muon*>(muonA->At(indx));
      _fLog << "i = " << indx
            << ", pt = " << muon->pt 
            << ", eta = " << muon->eta 
            << ", phi = " << muon->phi
            << ", trkD0 = "<< muon->trkD0
            << ", vtxDistZ"<< muon->vtxDistZ
            << std::endl;
    }
    _fLog << "electron block" << std::endl;
    for (int indx = 0; indx < n_electron; ++indx) {
      const Electron* elec = dynamic_cast<Electron*>(electronA->At(indx));
      _fLog << "i = " << indx
            << ", pt = " << elec->pt 
            << ", eta = " << elec->eta 
            << ", phi = " << elec->phi
            << std::endl;
    }
    _fLog << "tau block" << std::endl;
    for (int indx = 0; indx < n_tau; ++indx) {
      const Tau* tau = dynamic_cast<Tau*>(tauA->At(indx));
      _fLog << "i = " << indx
            << ", pt = " << tau->pt 
            << ", eta = " << tau->eta 
            << ", phi = " << tau->phi
            << std::endl;
    }
    _fLog << "end of an event" << std::endl;
    //HAS TO BE DELETED

    // Show status of the run
    if (currentFile != lastFile) 
      cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
           << " ==> " << currentFile 
           << " <<< Run# " << run
           << " Event# " << setw(8) << event << " >>> " 
           << " Events proc. " << setw(8) << ev << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
       cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
            << " ==> " << _chain->GetCurrentFile()->GetName() 
            << " <<< Run# " << run 
            << " Event# " << setw(8) << event << " >>> " 
            << " Events proc. " << setw(8) << ev << endl;

    if (_logOption >> 0 & 0x1) 
    _fLog << "run: " << run
          << ", event: " << event
          << ", n_tau: "<< n_tau
          << ", n_muon: "<< n_muon
          << ", n_jet: " << n_jet
          << ", n_vertex: " << n_vertex
          << ", n_met: " << n_met
          << ", n_electron: " << n_electron << endl;

    clearLists();
 
    //findGenInfo();
    if (_logOption >> 6 & 0x1) {
      bool dumpAll = false; // change to true to dump all the events
      if (dumpAll || (genMuonList.size() > 1 && genTauList.size() > 0)) dumpGenInfo();
    }

    findVtxInfo(vtxList, _vtxCutMap);
    findTauInfo(tauList, _tauCutMap);
    findMuonInfo(muoList, _muonCutMap);
    findElectronInfo(eleList, _electronCutMap);
    findJetInfo(bjetList, _bjetCutMap);

    if (_logOption)
    _fLog << "run: " << run
          << ", event: " << event
          << ", n_vertex_good: "   << vtxList.size()
          << ", n_muon_selected: " << muoList.size()
          << ", n_electron: "      << eleList.size()
          << ", n_tau_selected: "  << tauList.size()
          << ", n_bjet: "          << bjetList.size()
          << endl;

    // Event Selection Starts here .....
    // presence of > 1 good vertex
    if (vtxList.size() < 1) continue;
    ++nEvtVtx;
    selectMuMuTauEvent(muoList, tauList, bjetList);
    selectMuTauTauEvent(muoList, tauList, bjetList);
    selectMuMuTauTauEvent(muoList, tauList, bjetList);
    selectMuMuMuTauEvent(muoList, tauList, bjetList);
    //computeSignificance(muoList, tauList, bjetList);
    computeDeltaR(muoList,tauList); 
  }  
  // Analysis is over
  nEvtMuMuTau[0] = nEvtMuTauTau[0] = nEvtMuMuTauTau[0] = nEvtMuMuMuTau[0] = nEvt;
  nEvtMuMuTau[1] = nEvtMuTauTau[1] = nEvtMuMuTauTau[1] = nEvtMuMuMuTau[1] = nEvtVtx;

  endJob();
}
void AnaWH::findGenInfo() {
      // Generator level information
    if (n_genparticle) {
      //if (_logOption >> 6 & 0x1) 
      //  _fLog << "indx    status    pdgId     eta      phi      pt     energy            mID                             dID"
      //        << endl;
      for (int indx = 0; indx < n_genparticle; ++indx) {
        const GenParticle* gp = dynamic_cast<GenParticle*>(genParticleA->At(indx));
        if (!gp) continue;
        
	// cout << "indx = " << indx << ", id=" << gp->pdgId << ", status = " << gp->status << endl;
     
	//looking for a Hadronic Tau whose mother is Higgs
        if (gp->status == 2 && abs(gp->pdgId)==15) {
           vector<int> m = gp->motherIndices;
           const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(m[0]));
           if (!mgp) continue;

           int pdgid = mgp->pdgId;
           if (abs(pdgid) == 15) {
             vector<int> m2 = mgp->motherIndices;
             const GenParticle* m2gp = dynamic_cast<GenParticle*>(genParticleA->At(m2[0]));
             if (!m2gp) continue;
             pdgid = m2gp->pdgId;
             if (abs(pdgid)== 25) {
	       genHList.push_back(*m2gp);
             }
           }
           
           if (abs(pdgid) != 25) continue;

           vector<int> d = gp->daughtIndices;
           for (size_t i = 0; i < d.size(); ++i) {
	     int di = d[i];
             if (di >= n_genparticle) continue;
             const GenParticle* dgp = dynamic_cast<GenParticle*>(genParticleA->At(di));
             if (!dgp) continue;
             if (abs(dgp->pdgId) == 16) continue;
             if (abs(dgp->pdgId)!=11 && abs(dgp->pdgId)!=13 && abs(dgp->pdgId) !=12 && abs(dgp->pdgId) !=14 && dgp->pdgId !=22) {
               //cout << "pId=" << gp->pdgId << ", dId=" << dgp->pdgId << endl;
  	       genTauList.push_back(*gp);
               break;
             }
           }
          
        }


        // Looking for Muon whose mother is either W boson or Tau 
        if (gp->status == 1 && abs(gp->pdgId)== 13) {
	  vector<int> m = gp->motherIndices;
          const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(m[0]));
          if (!mgp) continue;
          if (abs (mgp->pdgId)==15)
	    genMuobList.push_back(*gp);
          if (abs (mgp->pdgId)==24){
	    genMuoaList.push_back(*gp);
            genWList.push_back(*mgp);
          }

          int mmid= mgp->pdgId;
          if (abs(mmid)==13){
            vector<int> m2 = mgp->motherIndices;
            const GenParticle* m2gp = dynamic_cast<GenParticle*>(genParticleA->At(m2[0]));
            if (!m2gp) continue;
            mmid = m2gp->pdgId;
            if (abs(mmid)==15)
	      genMuobList.push_back(*gp);
            if (abs(mmid)==24){
              genWList.push_back(*m2gp);
              genMuoaList.push_back(*gp);
	    }
          }
          
          
	   if (abs(mmid)==24 || abs(mmid)==15)
             genMuonList.push_back(*gp);
       	}

      	//looking for a Hadronic Tau whose mother is H-Boson
        if (gp->status == 2 && abs(gp->pdgId)==15) {
           vector<int> m = gp->motherIndices;
           const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(m[0]));
           if (!mgp) continue;

           int pdgid = mgp->pdgId;
           if (abs(pdgid) == 15) {
             vector<int> m2 = mgp->motherIndices;
             const GenParticle* m2gp = dynamic_cast<GenParticle*>(genParticleA->At(m2[0]));
             if (!m2gp) continue;
             pdgid = m2gp->pdgId;
           }
           
           if (abs(pdgid) != 25) continue;

           vector<int> d = gp->daughtIndices;
           for (size_t i = 0; i < d.size(); ++i) {
	     int di = d[i];
             if (di >= n_genparticle) continue;
             const GenParticle* dgp = dynamic_cast<GenParticle*>(genParticleA->At(di));
             if (!dgp) continue;
             if (abs(dgp->pdgId) == 16) continue;
             if (abs(dgp->pdgId)!=11 && abs(dgp->pdgId)!=13 && abs(dgp->pdgId) !=12 && abs(dgp->pdgId) !=14 && dgp->pdgId !=22) {
               //cout << "pId=" << gp->pdgId << ", dId=" << dgp->pdgId << endl;
  	       genTau_H_List.push_back(*gp);
               break;
             }
           }
          
        }
        // Looking for 2 Muon whose mother is Z 
        if (gp->status == 1 && abs(gp->pdgId)== 13) {
	  vector<int> m = gp->motherIndices;
          const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(m[0]));
          if (!mgp) continue;

          int mmid= mgp->pdgId;
          if (abs(mmid)==13){
            vector<int> m2 = mgp->motherIndices;
            const GenParticle* m2gp = dynamic_cast<GenParticle*>(genParticleA->At(m2[0]));
            if (!m2gp) continue;
            mmid = m2gp->pdgId;
            if (abs (mmid)== 13) {
              vector<int>m3 = m2gp->motherIndices;
              const GenParticle* m3gp = dynamic_cast<GenParticle*>(genParticleA->At(m3[0]));
              if (!m3gp) continue;
              mmid = m3gp->pdgId;
	    }
          }
 
	  if (abs(mmid)==23)
            genMuon_Z_List.push_back(*gp);
       	}



	//calculation of invariant mass in the generator level
	if(gp->status==1 && abs(gp->pdgId)== 13)
	  GEN_MUON_LIST.push_back(*gp); 
     
      }
  //end of genParticle loop
      if (GEN_MUON_LIST.size()>1){
        sort (GEN_MUON_LIST.begin(), GEN_MUON_LIST.end(), PtComparator<GenParticle>());
        const GenParticle& mua1=GEN_MUON_LIST[0];
        const GenParticle& mua2=GEN_MUON_LIST[1];
        TLorentzVector x1, x2;
        x1.SetPtEtaPhiE (mua1.pt, mua1.eta, mua1.phi, mua1.energy);
        x2.SetPtEtaPhiE (mua2.pt, mua2.eta, mua2.phi, mua2.energy);
        TLorentzVector sum= x1+x2 ;
        double MASS= sum.M();
        fillHist1D ("genInvMass", MASS);
      }
      
      //fillHist1D ("muon_count", genMuonList.size(), 1.0);  
      if (genMuonList.size() > 1) 
        sort(genMuonList.begin(), genMuonList.end(), PtComparator<GenParticle>());
      //      if (genMuonList.size() > 0) {
      //   double pt = genMuonList[0].pt;
      //  fillHist1D ("genMuon_pt1", pt, 1.0);
      // }
      if (genTauList.size() > 1)
        sort(genTauList.begin(), genTauList.end(), PtComparator<GenParticle>());
      if (genMuoaList.size() > 1)
        sort(genMuoaList.begin(), genMuoaList.end(), PtComparator<GenParticle>());
      if (genMuobList.size() > 1)
        sort(genMuobList.begin(), genMuobList.end(), PtComparator<GenParticle>());
      if (genHList.size() > 1)
        sort(genHList.begin(), genHList.end(), PtComparator<GenParticle>());
      if (genWList.size() > 1)
        sort(genWList.begin(), genWList.end(), PtComparator<GenParticle>());
 
     if (genMuonList.size()>1 && genTauList.size()>0){         
         const GenParticle& muona=genMuonList[0];
         const GenParticle& muonb=genMuonList[1];
         const GenParticle& taua =genTauList[0];
         if (muona.pt<15) return;
         if (muonb.pt<10) return;
         TLorentzVector v1,v2,t1;
         v1.SetPtEtaPhiE (muona.pt, muona.eta, muona.phi, muona.energy);
         v2.SetPtEtaPhiE (muonb.pt, muonb.eta, muonb.phi, muonb.energy);
         t1.SetPtEtaPhiE (taua.pt , taua.eta , taua.phi , taua.energy );
         TLorentzVector z= v1+v2;
         double mass= z.M();
         if (mass < 20) return;
         double mva_gen= z.Pt()/(v1.Pt()+v2.Pt());
  
         fillHist1D("mva_gen", mva_gen, 1.0);
         fillHist2D("muon_correlation", v1.Pt(), v2.Pt(), 1.0);
	 fillHist1D("oddInvMass",mass, 1.0);
      }
      if (genMuoaList.size()>0 && genMuobList.size()>0 && genTauList.size()>0) {
      
       	TLorentzVector M1,M2,W1,H1,t1;
        M1.SetPtEtaPhiE (genMuoaList[0].pt, genMuoaList[0].eta, genMuoaList[0].phi, genMuoaList[0].energy);
        M2.SetPtEtaPhiE (genMuobList[0].pt, genMuobList[0].eta, genMuobList[0].phi, genMuobList[0].energy);
        W1.SetPtEtaPhiE (genWList[0].pt, genWList[0].eta, genWList[0].phi, genWList[0].energy);
        H1.SetPtEtaPhiE (genHList[0].pt, genHList[0].eta, genHList[0].phi, genHList[0].energy);
        t1.SetPtEtaPhiE (genTauList[0].pt, genTauList[0].eta, genTauList[0].phi, genTauList[0].energy);
        fillHist1D("WH_DR", AnaUtil::deltaR(W1, H1), 1.0);
        fillHist1D("Ht1_DR",AnaUtil::deltaR(t1, H1), 1.0);
        fillHist1D("HM2_DR",AnaUtil::deltaR(M2, H1), 1.0);
        fillHist1D("WM1_DR",AnaUtil::deltaR(M1, W1), 1.0);
        fillHist1D("t1M2_DR",AnaUtil::deltaR(t1,M2), 1.0);
        fillProfile("WH_Wpt",W1.Pt(),AnaUtil::deltaR(W1, H1), 1.0);
        fillProfile("t1M2_Hpt",H1.Pt(),AnaUtil::deltaR(t1,M2), 1.0);
        fillProfile("WM1_Wpt", W1.Pt(),AnaUtil::deltaR(W1,M1), 1.0);
        fillHist1D("DPhi_WH",AnaUtil::deltaPhi(W1,H1), 1.0);
        fillHist1D("DEta_WH",W1.Eta()-H1.Eta(), 1.0);
        fillHist1D("DPhi_t1M2",AnaUtil::deltaPhi(t1,M2), 1.0);
        fillHist1D("DEta_t1M2",t1.Eta()-M2.Eta(), 1.0);
      }

      
    }
  //end of generator level   
    /* if (genTau_H_List.size()>0 && genMuon_Z_List.size()>1){
      signal_ZH = signal_ZH + 1;
    }

    if (genMuonList.size()>1 && genTauList.size()>0){
      signal= signal+1;
    }*/
}
void AnaWH::dumpGenInfo() {
  if (!n_genparticle) return;
  _fLog << setprecision(2);
  _fLog << "indx    status    pdgId     eta      phi      pt     energy             mID                             dID"
        << endl;
  for (int indx = 0; indx < n_genparticle; ++indx) {
     const GenParticle* gp = dynamic_cast<GenParticle*>(genParticleA->At(indx));
     if (!gp) continue;

     ostringstream mID;
     vector<int> m = gp->motherIndices;
     for (size_t i = 0; i < m.size(); ++i) {
       int mi = m[i];
       if (mi >= n_genparticle) continue;
       const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(mi));
       if (!mgp) continue;
       mID << " " << mgp->pdgId; 
     }
     string ms = mID.str();
     if (!ms.length()) ms = " -";
     
     ostringstream dID;
     vector<int> d = gp->daughtIndices;
     for (size_t i = 0; i < d.size(); ++i) {
     	  int di = d[i];
       if (di >= n_genparticle) continue;
       const GenParticle* dgp = dynamic_cast<GenParticle*>(genParticleA->At(di));
       if (!dgp) continue;
       double energy = dgp->energy;
       int pdgid = dgp->pdgId;
       if (abs(pdgid) == 21 && energy <= 1) continue;
       dID << " " << dgp->pdgId; 
     }
     string ds = dID.str();
     if (!ds.length()) ds = " -";

     _fLog << setw(4)  << indx
           << setw(8)  << gp->status
           << setw(10) << gp->pdgId
           << setw(10) << gp->eta
           << setw(9)  << gp->phi
           << setw(9)  << gp->pt
           << setw(9)  << gp->energy
           << setw(16) << ms 
           << ds
           << endl;
#if 0 
     vector<int> m = gp->motherIndices;
     _fLog << "# of mother: "  << m.size() << endl;
     for (size_t i = 0; i < m.size(); ++i) {
       int mi = m[i];
       _fLog << "index: " << mi << endl;
       if (mi >= n_genparticle) continue;
       const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(mi));
       _fLog << "\tpdgId: " << mgp->pdgId << ", status:" << mgp->status << endl; 
     }

     vector<int> d = gp->daughtIndices;
     _fLog << "# of daughter: " << d.size() << endl;
     for (size_t i = 0; i < d.size(); ++i) {
     	  int di = d[i];
       _fLog << "index: " << di << endl;
       if (di >= n_genparticle) continue;
       const GenParticle* dgp = dynamic_cast<GenParticle*>(genParticleA->At(di));
       _fLog << "\tpdgId: " << dgp->pdgId << ", status: " << dgp->status << endl; 
     }
#endif
  }
}


void AnaWH::computeDeltaR(const vector<Muon>& muoList, const vector<Tau>& tauList ){
//findTauInfo();
//findMuoInfo();
 if(muoList.size()>1 && tauList.size()>0){
  double deltaDRMuMu = sqrt(pow((muoList[0].eta - muoList[1].eta),2) + pow((muoList[0].phi - muoList[1].phi),2));
  double deltaDRMuTau1 = sqrt(pow((muoList[0].eta - tauList[0].eta),2) + pow((muoList[0].phi - tauList[0].phi),2));
  double deltaDRMuTau2 = sqrt(pow((muoList[1].eta - tauList[0].eta),2) + pow((muoList[1].phi - tauList[0].phi),2));
  double deltaEtaMuMu = (muoList[0].eta - muoList[1].eta);
  double deltaEtaMuTau1 = (muoList[0].eta - tauList[0].eta);
  double deltaEtaMuTau2 = (muoList[1].eta - tauList[0].eta);
  double deltaPhiMuMu = (muoList[0].phi - muoList[1].phi);
  double deltaPhiMuTau1 = (muoList[0].phi - tauList[0].phi);
  double deltaPhiMuTau2 = (muoList[1].phi - tauList[0].phi);

  const Muon& muoa = muoList[0];
  const Muon& muob = muoList[1];

  if (muoa.pt<15) return;
  if (muob.pt<10) return;

  double pt_plus= muoa.pt+muob.pt;
  double pt_minus= muoa.pt-muob.pt;
  fillHist2D("pt covariance", pt_plus, pt_minus, 1.0); 

  TLorentzVector m1, m2;
  m1.SetPtEtaPhiE(muoa.pt, muoa.eta, muoa.phi, muoa.energy); 
  m2.SetPtEtaPhiE(muob.pt, muob.eta, muob.phi, muob.energy); 
  TLorentzVector z = m1 + m2;
  double mass = z.M();
  if(mass < 20) return;

  //just to check the charges of the muons 
  fillHist1D ("muoa charge", muoa.charge, 1.0);
  fillHist1D ("muob charge", muob.charge, 1.0);
  fillHist1D ("tau charge", tauList[0].charge, 1.0);

  if ( (muoa.charge + muob.charge) != 0) return;
  if ( (muob.charge + tauList[0].charge)  != 0) return;
  
  fillHist2D ("dmuon_correlation", m1.Pt(), m2.Pt(), 1.0);
  fillHist2D ("eta_correlation", m1.Eta(), m2.Eta(), 1.0);
  fillHist2D ("P_correlation", m1.P(), m2.P(), 1.0);
  

  double mva_eta= z.Eta();
  double mva1= z.Pt()/(m1.Pt()+m2.Pt());
  if(mva1>0.98) {
    double DRMuMu = AnaUtil::deltaR(m1, m2);
    double deleta1= z.Eta()-m1.Eta();
    double deleta2= z.Eta()-m2.Eta();

    fillHist1D ("DRMuMu", DRMuMu, 1.0);
    fillHist1D ("deleta1", deleta1, 1.0);
    fillHist1D ("deleta2", deleta2, 1.0);
  }
  if (mva1<0.5) {
    double Deta1= z.Eta()-m1.Eta();
    double Deta2= z.Eta()-m2.Eta();

    fillHist1D("Deta1", Deta1, 1.0);
    fillHist1D("Deta2", Deta2, 1.0);
  }

  fillHist1D("PhiMu1",muoList[0].phi,1.0);
  fillHist1D("PhiMu2",muoList[1].phi,1.0);
  fillHist1D("PhiTau",tauList[0].phi,1.0);
 
  //  fillHist1D("DRMuMu",deltaDRMuMu,1.0);
  fillHist1D("DRMuTau1",deltaDRMuTau1,1.0);
  fillHist1D("DRMuTau2",deltaDRMuTau2,1.0);

  fillHist1D("DEtaMuMu",deltaEtaMuMu,1.0);
  fillHist1D("DEtaMuTau1",deltaEtaMuTau1,1.0);
  fillHist1D("DEtaMuTau2",deltaEtaMuTau2,1.0);
  fillHist1D("DPhiMuMu",deltaPhiMuMu,1.0);
  fillHist1D("DPhiMuTau1",deltaPhiMuTau1,1.0);
  fillHist1D("DPhiMuTau2",deltaPhiMuTau2,1.0);
    
  fillHist1D("mva1", mva1, 1.0);
  fillHist1D("mva_eta", mva_eta, 1.0);

 }
}
/*void AnaBase::GenDetLevelMatch (const vector<Muon>& muoList,
                                const vector<GenParticle>& genMuonList){
for (int i = 0; i < 2; i++){
    double minDR = 990;
    int index = -1;
    for(int j = 0; j < genMuonList.size();j++){
       double deltaDR = sqrt(pow((genMuonList[j].eta - muolist[i].eta),2) + pow((genMuonList[j].phi - muolist[i].phi),2));
       if(deltaDR<minDr){
          minDR = deltaDR;
          index = j; //gen label index
          }
       }
    if(index < 0) continue;
    
    if(i = 0){
       fillHist1D("deltaPt1",muolist[0].pt - genMuonList[index].pt);
       fillHist1D("deltaPt1",muolist[0].eta - genMuonList[index].eta);
       fillHist1D("deltaPt1",muolist[0].phi - genMuonList[index].phi);
    }
    else{
       fillHist1D("deltaDR1",minDR);
       fillHist1D("deltaPt1",muolist[1].pt - genMuonList[index].pt);
       fillHist1D("deltaPt1",muolist[1].eta - genMuonList[index].eta);
       fillHist1D("deltaPt1",muolist[1].phi - genMuonList[index].phi);
       }  
    }
}
*/
void AnaWH::selectMuMuTauEvent(const vector<Muon>& muoList, 
                                 const vector<Tau>& tauList, 
                                 const vector<Jet>& bjetList) 
{
  fillHist1D("nMuon", muoList.size(), 1.0);
  // 2 muons
  if (muoList.size() < 2) return;
  ++nEvtMuMuTau[2];

  // muon Pt
  const Muon& muoa = muoList[0];
  const Muon& muob = muoList[1];
  fillHist1D("muonPt1", muoa.pt, 1.0);
  fillHist1D("muonPt2", muob.pt, 1.0);
  if (muoa.pt <= 15) return;
  ++nEvtMuMuTau[3];

  if (muob.pt <= 10) return;
  ++nEvtMuMuTau[4];

  TLorentzVector m1, m2;
  m1.SetPtEtaPhiE(muoa.pt, muoa.eta, muoa.phi, muoa.energy); 
  m2.SetPtEtaPhiE(muob.pt, muob.eta, muob.phi, muob.energy); 
  TLorentzVector z = m1 + m2;
  double mass = z.M();
  //fillHist1D("mumuInvMass", mass, 1.0);
  if (mass <= 20) return;
  ++nEvtMuMuTau[5];

  // 1 hadronic tau
  fillHist1D("nTau", tauList.size(), 1.0);
  if (tauList.size() < 1) return;
  ++nEvtMuMuTau[6];

  const Tau& tau = tauList[0];
  fillHist1D("tauPt1", tau.pt, 1.0);

  // Tau Pt
  if (tau.pt <= 15) return;
  ++nEvtMuMuTau[7];

  // no b-tagged jet
  if (bjetList.size()) return;
  ++nEvtMuMuTau[8];

  // Opposite charge for mu_2 & tau
  if ( (muob.charge + tau.charge) != 0) return;
  ++nEvtMuMuTau[9];

  // Same charge for muon_1 & muon_2 (to remove DY)
  if ( (muoa.charge + muob.charge) == 0) return;
  ++nEvtMuMuTau[10];
}
void AnaWH::selectMuTauTauEvent(const vector<Muon>& muoList, 
                                  const vector<Tau>& tauList, 
                                  const vector<Jet>& bjetList) 
{
  // 1 muons
  if (muoList.size() < 1) return;
  ++nEvtMuTauTau[2];

  // muon Pt
  const Muon& muo = muoList[0];
  if (muo.pt <= 17) return;
  ++nEvtMuTauTau[3];

  // 2 hadronic tau
  if (tauList.size() < 2) return;
  ++nEvtMuTauTau[4];

  const Tau& taua = tauList[0];
  const Tau& taub = tauList[1];

  // Tau Pt
  if (taua.pt <= 22) return;
  ++nEvtMuTauTau[5];

  if (taub.pt <= 15) return;
  ++nEvtMuTauTau[6];

  // no b-tagged jet
  if (bjetList.size()) return;
  ++nEvtMuTauTau[7];

  // Opposite charge for tau_1 & tau_2
  if ( (taua.charge + taub.charge) != 0) return;
  ++nEvtMuTauTau[8];
}
void AnaWH::selectMuMuTauTauEvent(const vector<Muon>& muoList, 
                                  const vector<Tau>& tauList, 
                                  const vector<Jet>& bjetList) 
{
  // 2 muons
  if (muoList.size() < 2) return;
  ++nEvtMuMuTauTau[2];

  // muon Pt
  const Muon& muoa = muoList[0];
  const Muon& muob = muoList[1];   
  if (muoa.pt <= 15) return;
  ++nEvtMuMuTauTau[3];
  if (muob.pt <= 10) return;
  ++nEvtMuMuTauTau[4];
  
  TLorentzVector m1, m2;
  m1.SetPtEtaPhiE(muoa.pt, muoa.eta, muoa.phi, muoa.energy); 
  m2.SetPtEtaPhiE(muob.pt, muob.eta, muob.phi, muob.energy); 
  TLorentzVector z = m1 + m2;
  double mass = z.M();
  fillHist1D("mumuInvMass_ZH", mass, 1.0);
  if (mass <= 20) return;
  ++nEvtMuMuTauTau[5];

  // 2 hadronic tau
  if (tauList.size() < 2) return;
  ++nEvtMuMuTauTau[6];

  const Tau& taua = tauList[0];
  const Tau& taub = tauList[1];

  // Tau Pt
  if (taua.pt <= 15) return;
  ++nEvtMuMuTauTau[7];

  if (taub.pt <= 10) return;
  ++nEvtMuMuTauTau[8];

  // no b-tagged jet
  if (bjetList.size()) return;
  ++nEvtMuMuTauTau[9];

  // Opposite charge for tau_1 & tau_2
  if ( (taua.charge + taub.charge) != 0) return;
  ++nEvtMuMuTauTau[10];
}
void AnaWH::selectMuMuMuTauEvent(const vector<Muon>& muoList, 
                                  const vector<Tau>& tauList, 
                                  const vector<Jet>& bjetList) 
{
  // 3 muons
  if (muoList.size() < 3) return;
  ++nEvtMuMuMuTau[2];

  // muon Pt
  const Muon& muoa = muoList[0];
  const Muon& muob = muoList[1];
  const Muon& muoc = muoList[2];   
  if (muoa.pt <= 15) return;
  ++nEvtMuMuMuTau[3];
  if (muob.pt <= 10) return;
  ++nEvtMuMuMuTau[4];
  if (muoc.pt <= 10) return;
  ++nEvtMuMuMuTau[5];
  
  TLorentzVector m1, m2;
  m1.SetPtEtaPhiE(muoa.pt, muoa.eta, muoa.phi, muoa.energy); 
  m2.SetPtEtaPhiE(muob.pt, muob.eta, muob.phi, muob.energy); 
  TLorentzVector z = m1 + m2;
  double mass = z.M();
  fillHist1D("mumuInvMass_ZH", mass, 1.0);
  if (mass <= 20) return;
  ++nEvtMuMuMuTau[6];

  // 1 hadronic tau
  if (tauList.size() < 1) return;
  ++nEvtMuMuMuTau[7];

  const Tau& taua = tauList[0];

  // Tau Pt
  if (taua.pt <= 15) return;
  ++nEvtMuMuMuTau[8];

  // no b-tagged jet
  if (bjetList.size()) return;
  ++nEvtMuMuMuTau[9];

}

void AnaWH::computeSignificance(const vector<Muon>& muoList, 
                                  const vector<Tau>& tauList, 
                                  const vector<Jet>& bjetList) 
{
  // 2 muons
  if (muoList.size() < 2) return;

  // muon Pt
  const Muon& muoa = muoList[0];
  const Muon& muob = muoList[1];

  if (muoa.pt <= 10) return;

  if (muob.pt <= 8) return;

  TLorentzVector m1, m2;
  m1.SetPtEtaPhiE(muoa.pt, muoa.eta, muoa.phi, muoa.energy); 
  m2.SetPtEtaPhiE(muob.pt, muob.eta, muob.phi, muob.energy); 
  TLorentzVector z = m1 + m2;
  double mass = z.M();
  if (mass <= 20) return;

  // 1 hadronic tau
  if (tauList.size() < 1) return;

  const Tau& tau = tauList[0];

  // Tau Pt
  if (tau.pt <= 12) return;

  // no b-tagged jet
  if (bjetList.size()) return;

  // Opposite charge for mu_2 & tau 
  if ( (muob.charge + tau.charge) != 0) return;

  // Same charge for muon_1 & muon_2 (to remove DY)
  if ( (muoa.charge + muob.charge) == 0) return;
 
  // Now study the effect of highest-pt muon pt
  ++nSig[0];
  double pt= muoa.pt;
  double vi = 5, step = 5, vf = 100;
  double value = vi;
  int i = 0;
  while (value <= vf) {
    if (pt > value) ++nSig[++i];
    value += step;
  }
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void AnaWH::endJob() 
{  
  const string tagsMuMuTau[] = 
  {
    "Total events processed",
    ">= 1 Good event vertex",
    ">= 2 selected muons",
    "First Muon Pt> 15 GeV cut",
    "Second Muon Pt>10 GeV cut",
    "Di-Muon invariant mass >= 20 GeV Cut",
    ">= 1 selected tau",
    "First Tau Pt>15 GeV cut",
    "no tagged b-jet",
    "mu+tau charge == 0",
    "mu+mu charge != 0"
  };

  _fLog << "Statistics for W->munu, H->tau+tau-; tau->mu,tau->hadron" << endl;
  for (unsigned int i = 0; i < NEL(nEvtMuMuTau); i++) {
    _fLog << setw(64) << tagsMuMuTau[i] << setw(10) << nEvtMuMuTau[i] << endl;
  }
  const string tagsMuTauTau[] = 
  {
    "Total events processed",
    ">= 1 Good event vertex",
    ">= 1 selected muons",
    "Muon Pt> 17 GeV cut",
    ">= 2 selected taus",
    "First Tau Pt>22 GeV cut",
    "First Tau Pt>15 GeV cut",
    "no tagged b-jet",
    "tau+tau charge == 0"
  };
  _fLog << "Statistics for W->munu, H->tau+tau-; tau->hadron,tau->hadron" << endl;
  for (unsigned int i = 0; i < NEL(nEvtMuTauTau); i++) {
    _fLog << setw(64) << tagsMuTauTau[i] << setw(10) << nEvtMuTauTau[i] << endl;
  }
  const string tagsMuMuTauTau[] =
  {
    "Total events processed",
    ">= 1 Good event vertex",
    ">= 2 selected Muons",
    "1st Muon Pt > 15 GeV Cut",
    "2nd Muon Pt > 10 Gev Cut",
    "di-muon inv mass> 20 GeV cut",
    ">= 2 selected taus",
    "1st Tau Pt >15 GeV cut",
    "2nd Tau PT >10 Gev cut",
    "no taged b-jet", 
    "tau+tau charge == 0"
  };
  _fLog << "Statistics for ZH process; Z to Mu Mu and H to Tau Tau"<< endl;
  for (unsigned int i=0; i< NEL(nEvtMuMuTauTau); i++) {
    _fLog << setw(64) << tagsMuMuTauTau[i] << setw(10) << nEvtMuMuTauTau[i]<< endl;
  }
  const string tagsMuMuMuTau[] =
  {
    "Total events processed",
    ">= 1 Good event vertex",
    ">= 3 selected Muons",
    "1st Muon Pt > 15 GeV Cut",
    "2nd Muon Pt > 10 Gev Cut",
    "3rd Muon Pt > 10 GeV Cut",
    "di-muon inv mass> 20 GeV cut",
    ">= 1 selected taus",
    "1st Tau Pt >15 GeV cut",
    "no taged b-jet"
  };
  _fLog << "Statistics for ZH process Z to Mu Mu and H to Mu Tau"<< endl;
  for (unsigned int i=0; i< NEL(nEvtMuMuMuTau); i++) {
    _fLog << setw(64) << tagsMuMuMuTau[i] << setw(10) << nEvtMuMuMuTau[i]<< endl;
  }

  // Print significance  numbers
  _fLog << "Events selected with N-1 cuts: " << nSig[0] << endl;
  _fLog << "Events selected with Pt cuts: " << endl;
  double vi = 10, step = 5;
  _fLog << setprecision(0);
  for (unsigned int i = 1; i < NEL(nSig); i++) {
    _fLog << setw(4) << vi << setw(10) << nSig[i] << endl;
    vi += step;
  }  
  _fLog << "Signal like events: " << count_sig_event << endl;
  _fLog << resetiosflags(ios::fixed);

  closeFiles();

  _histf->cd();
  _histf->Write();
  _histf->Close();
  delete _histf;
}
// ----------------------------------------------------------
// Perform event selection, For selection of Z -> e+e- events
// we need,
//   - > 0 Tight electron
//   - event within e+e- invariant mass window
// ----------------------------------------------------------
bool AnaWH::selectEvent() 
{
  return true;
}
//void AnaBase::fillTriggerHistograms(int run) 
//{
//  for (unsigned int ipath = 0; ipath  < _hltpaths->size(); ipath++) {
//    string pname = (*_hltpaths)[ipath];
//    for (vector<TProfile*>::iterator it = _triggerHistos.begin(); it != _triggerHistos.end(); it++) {
//      string hname = (*it)->GetName();
//      if (pname.find(hname) != string::npos) (*it)->Fill(run, (*_hltprescales)[ipath]);   
//    }
//  }
//}
// -------------------------------------------------------------------------------
// Poor man's way of a datacard. Each line between the 'START' and 'END' tags
// is read in turn, split into words, where the first element is the 'key' and
// the rest the value(s). If more than one values are present they are 
// stored in a vector. No safety mechanism is in place. Any line with an unknown 
// key is skipped. Comments lines should start with either '#' or '//', preferably
// in the first column. Empty lines are skipped. The file containing the datacards 
// is passed as the only argument of the program, there is no default datacard
// -------------------------------------------------------------------------------
bool AnaWH::readJob(const string& jobFile, int& nFiles)
{
  AnaBase::readJob(jobFile, nFiles);

  static const int BUF_SIZE = 256;

  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "Input File: " << jobFile << " could not be opened!" << endl;
    return false;
  }

  // note that you must use a pointer (reference!) to the cut map
  // in order to avoid scope related issues
  map<string, map<string, double>* > hmap;
  hmap.insert(pair<string, map<string, double>* >("vtxCutList", &_vtxCutMap));
  hmap.insert(pair<string, map<string, double>* >("electronCutList", &_electronCutMap));
  hmap.insert(pair<string, map<string, double>* >("muonCutList", &_muonCutMap));
  hmap.insert(pair<string, map<string, double>* >("tauCutList", &_tauCutMap));
  hmap.insert(pair<string, map<string, double>* >("bjetCutList", &_bjetCutMap));

  char buf[BUF_SIZE];
  vector<string> tokens;
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;   

    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;

    // Split the line into words
    AnaUtil::tokenize(line, tokens);
    int size = tokens.size();
    string key = tokens[0];

    map<string, map<string, double>* >::const_iterator pos = hmap.find(key);
    if (pos != hmap.end()) {
      map<string, double>* m = pos->second;        
      m->clear();
      for (int i = 1; i < size; ++i) {
        vector<string> cutstr;
        // Split the line into words
	AnaUtil::tokenize(tokens[i], cutstr, "=");
        m->insert( pair<string,double>(cutstr[0], atof(cutstr[1].c_str())));
      }
    }    
    tokens.clear();
  }
  // Close the file
  fin.close();

  printJob();
  return true;
}
void AnaWH::printJob(ostream& os)
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("vtxCutList", _vtxCutMap));
  hmap.insert(pair<string, map<string, double> >("electronCutList", _electronCutMap));
  hmap.insert(pair<string, map<string, double> >("muonCutList", _muonCutMap));
  hmap.insert(pair<string, map<string, double> >("tauCutList", _tauCutMap));
  hmap.insert(pair<string, map<string, double> >("bjetCutList", _bjetCutMap));

  for (map<string, map<string, double> >::const_iterator it  = hmap.begin(); 
                                                         it != hmap.end(); ++it)  
  {
    os << "=>>> " << it->first << endl; 
    map<string, double> m = it->second;
    os << setprecision(2);
    for (map<string,double>::const_iterator jt  = m.begin(); 
                                            jt != m.end(); ++jt)  
      os << jt->first << ": " << setw(7) << jt->second << endl;
    os << endl; 
  }
}
// ------------------------------------------------------------------------
// Convenience routine for filling 1D histograms. We rely on root to keep 
// track of all the histograms that are booked all over so that we do not 
// have to use any global variables to save the histogram pointers. Instead, 
// we use the name of the histograms and gROOT to retrieve them from the 
// Root object pool whenever necessary. This is the closest one can go to 
// hbook and ID based histogramming
// -------------------------------------------------------------------------
TH1* AnaWH::getHist1D(const string& hname) {
  TObject *obj = gDirectory->GetList()->FindObject(hname.c_str()); 
  if (!obj) {
    cerr << "**** Histogram for <<" << hname << ">> not found!" << endl;
    return 0;
  }
  TH1 *h = 0;
  if (obj->InheritsFrom("TH1D"))
    h = dynamic_cast<TH1D*>(obj);
  else if (obj->InheritsFrom("TH1C"))
    h = dynamic_cast<TH1C*>(obj);
  else if (obj->InheritsFrom("TH1K"))
    h = dynamic_cast<TH1K*>(obj);
  else if (obj->InheritsFrom("TH1S"))
    h = dynamic_cast<TH1S*>(obj);
  else if (obj->InheritsFrom("TH1I"))
    h = dynamic_cast<TH1I*>(obj);
  else
    h = dynamic_cast<TH1F*>(obj);

  if (!h) {
    cerr << "**** 1D Histogram <<" << hname << ">> not found" << endl;
    return 0;
  }
  return h;
}
template <class T>
bool AnaWH::fillHist1D(const string& hname, T value, double w) 
{
  TH1* h = getHist1D(hname);
  if (!h) return kFALSE;
  h->Fill(value, w);
  return kTRUE;
}
// ---------------------------------------------
// Convenience routine for filling 2D histograms
// ---------------------------------------------
template <class T1, class T2>
bool AnaWH::fillHist2D(const string& hname, T1 xvalue, T2 yvalue, double w) 
{
  TObject *obj = gDirectory->GetList()->FindObject(hname.c_str()); 
  TH2 *h = 0;

  if (obj->InheritsFrom("TH2D"))
    h = dynamic_cast<TH2D*>(obj);
  else if (obj->InheritsFrom("TH2C"))
    h = dynamic_cast<TH2C*>(obj);
  else if (obj->InheritsFrom("TH2S"))
    h = dynamic_cast<TH2S*>(obj);
  else if (obj->InheritsFrom("TH2I"))
    h = dynamic_cast<TH2I*>(obj);
  else
    h = dynamic_cast<TH2F*>(obj);

  if (!h) {
    cerr << "**** 2D Histogram <<" << hname << ">> not found" << endl;
    return kFALSE;
  }
  h->Fill(xvalue, yvalue, w);
  return kTRUE;
}
// --------------------------------------------------
// Convenience routine for filling profile histograms
// --------------------------------------------------
TProfile* AnaWH::getProfile(const string& hname) 
{
  TProfile *h = dynamic_cast<TProfile*>(gDirectory->GetList()->FindObject(hname.c_str()));
  if (!h) {
    cerr << "**** Profile Histogram <<" << hname << ">> not found" << endl;
    return 0;
  }
  return h;
}
bool AnaWH::fillProfile(const string& hname, float xvalue, float yvalue, double w) 
{
  TProfile *h = getProfile(hname);
  if (!h) return false;

  h->Fill(xvalue, yvalue, w);
  return true;
}
