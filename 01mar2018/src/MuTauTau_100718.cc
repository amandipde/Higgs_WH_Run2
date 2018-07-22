#include "configana.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>

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

#include "MuTauTau.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"
#include "MVASkim.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::abs;
using std::max;
using std::sqrt;
using std::sort;
using std::setprecision;
using std::setw;

using namespace vhtm;

// -----------
// Constructor
// -----------
MuTauTau::MuTauTau()
  : AnaBase(),
    _createMVATree(false),
    _mvaInputFile(""),
    _skimObj(0)
{}
// ----------
// Destructor
// ----------
MuTauTau::~MuTauTau() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MuTauTau::beginJob() 
{ 
  if (!AnaBase::beginJob()) return false;

  histf()->cd();
  bookHistograms();

  // Optionally write selected events to another tree
  if (_createMVATree) _skimObj = new MVASkim(_mvaInputFile);
  return true;
}
// ---------------
// Book histograms
// ---------------
void MuTauTau::bookHistograms() 
{
  static const double pi = TMath::Pi();

  if (isMC() && usePUWt()) {
    new TH1D("npu", "nPileUp", 50, 0, 50);
    new TH1D("puwt", "PileUp weight factor", 100, 0, 2.0);
    new TH1D("trueNInt", "True N-Interaction", 50, 0, 50);
  } 
  new TH1D("evcounter", "Selected event counter", 15, -0.5, 14.5);
  new TH1D("nvertex", "No of good vertex", 30, -0.5, 29.5);

  // STAGE 1
  new TH1F("muonPt1_stage1", "Highest Pt Muon Pt  after vtx and >=1 Muon requirment", 100, -0.5, 99.5);
  new TH1F("muonEta1_stage1", "Highest Pt Muon Eta  after vtx and >=1 Muon requirment", 100, -3, 3);
  new TH1F("met_stage1", "met distribution of the event before applying Muon Pt cut", 200, 0, 200);
  new TH1F("taupt1_afterMuon", "Highest Pt Tau Pt distribution after Muon requirement", 100, 0, 100);

  // STAGE 2 
  new TH1F("tauPt1_stage2", "Highest Pt Tau Pt distribution just before Tau1 Pt cut", 100, 0, 100);
  new TH1F("tauEta1_stage2", "Highest Pt Tau Eta distribution just before Tau1 Pt cut", 100, -3, 3);
  new TH1F("met_stage2", "met distribution of the event just before Tau1 Pt cut", 200, 0, 200);
  new TH1F("tauMVA1", "First Tau MVA against electrons", 2, -0.5, 1.5);
  new TH1F("tauZVertex1", "First Tau Z vertex", 100, 0, 100);

  // STAGE 3
  new TH1F("muonPt1_stage3", "Highest Pt Muon Pt after Tau1 Pt cut", 100, -0.5, 99.5);
  new TH1F("muonEta1_stage3", "Highest Pt Muon Eta after Tau1 Pt cut", 100, -3, 3);
  new TH1F("tauPt2_stage3", "2nd Highest Pt Tau Pt Distribution before Tau2  Pt cut", 100, -0.5, 99.5);
  new TH1F("tauEta2_stage3", "2nd Highest Pt Tau Eta Distribution before Tau2 Pt cut", 100, -3, 3);
  new TH1F("deltaR_mutau1_stage3", "delta R difference between the selected muon and highest pt tau just after Tau1 Pt cut", 100, 0, 10);
  new TH1F("met_stage3", "met distribution of the event just after Tau1 Pt cut", 200, 0, 200);
  new TH1F("deltaR_mutau2_stage3", "deltaR angle between mu and tau2 just before Tau2 pt cut", 100, 0, 10);
  new TH1F("deltaR_tau1tau2_stage3", "deltaR angle between tau1 and tau2 just before Tau2 pt cut", 100, 0, 10);
  new TH1F("tautauInvmass_stage3", "Tau Tau Invariant mass/vissible Higgs mass just before Tau2 Pt cut", 200, 0, 200);
  new TH1F("LT_stage3", "Scalar eT sum of three leptons before Tau2 pT cut", 200, 0, 200);
  new TH1F("vLT_stage3", "Vector eT sum of three leptons before Tau2 pT cut", 200, 0, 200);
  new TH1F("diffPhi_Mu_Met_stage3", "delta Phi difference of Mu and Met just before the Tau2 Pt cut", 100, 0, 2*pi);
  new TH1F("diffPhi_Mu_Tau1_stage3", "delta Phi difference of Mu and Tau1 just before the Tau2 Pt cut", 100, 0, 2*pi);
  new TH1F("DiTauPt_stage3", "Pt of the Di-Tau just before the Tau2 Pt cut", 200, 0, 200);
  new TH1F("diffPhi_Mu_DiTau_stage3", "delta Phi difference of Mu and Di-Tau just before the Tau2 Pt cut", 100, 0, 2*pi);
  new TH1F("tauMVA2", "Second Tau MVA against electrons", 2, -0.5, 1.5);
  new TH1F("tauZVertex2", "Second Tau Z vertex", 100, 0, 100);
  new TH1F("eventZ", "Z vertex for the 3 Leptons", 60, 0, 60);

  // STAGE 4
  new TH1F("deltaR_mutau2_stage4", "deltaR angle between mu and tau2 after Tau2 pt cut", 100, 0, 10);
  new TH1F("tautauInvmass_stage4", "Tau Tau Invariant mass// vissible Higgs mass just after Tau2 Pt cut", 200, 0, 200);
  new TH1F("met_stage4", "met distribution of the event just after Tau2 Pt cut", 200, 0, 200);
  new TH1F("deltaR_tau1tau2_stage4", "deltaR angle between tau1 and tau2 after Tau2 pt cut", 100, 0, 10);
  new TH1F("diffPhi_Mu_Met_stage4", "delta Phi difference of Mu and Met just after the Tau2 Pt cut", 100, 0, 2*pi);
  new TH1F("diffPhi_Mu_Tau1_stage4", "delta Phi difference of Mu and Tau1 just after the Tau2 Pt cut", 100, 0, 2*pi);
  new TH1F("DiTauPt_stage4", "Pt  of the Di-Tau just after the Tau2 Pt cut", 200, 0, 200);
  new TH1F("diffPhi_Mu_DiTau_stage4", "delta Phi difference of Mu and Di-Tau just after the Tau2 Pt cut", 100, 0, 2*pi);
  new TH1F("LT_stage4", "sum of Transverse Energy of the three leptons just after Tau2 Pt cut", 200, 0, 200);
  new TH1F("vLT_stage4", "Vector eT sum of three leptons after Tau2 pT cut", 200, 0, 200);
  new TH1F("var1_stage4", " plot of  Pt(Tau1+Tau2)/(tauiPt + tau2Pt) just after the Tau2 Pt cut", 100, 0, 2);
  new TH1F("siptSum", "Tau-pair signed IP significance", 60, 0, 60);

  // allcut
  new TH1F("nMuon_allcut", "Number of Muons", 6, -0.5, 5.5);
  new TH1F("muonPt1_allcut", "Highest Pt Muon Pt  after all cut applied", 100, -0.5, 99.5);
  new TH1F("muonEta1_allcut", "Highest Pt Muon Eta  after all cut applied", 100, -3, 3);
  new TH1F("tauPt1_allcut", "Highest Pt Tau Pt  after all cut applied", 100, -0.5, 99.5);
  new TH1F("tauEta1_allcut", "Highest Pt Tau Eta  after all cut applied", 100, -3, 3);
  new TH1F("tauPt2_allcut", "2nd Highest Pt Tau Pt  after all cut applied", 100, -0.5, 99.5);
  new TH1F("tauEta2_allcut", "2nd Highest Pt Tau Eta  after all cut applied", 100, -3, 3);
  new TH1F("tautauInvmass_allcut", "Tau Tau Invariant mass/vissible Higgs mass after all cut", 200, 0, 200);
  new TH1F("met_allcut", "met distribution of the event", 200, 0, 200); 
  new TH1F("diffPhi_Mu_Met_allcut", "delta Phi difference of Mu and Met", 100, 0, 2*pi);
  new TH1F("diffPhi_Mu_Tau1_allcut", "delta Phi difference of Mu and Tau1", 100, 0, 2*pi);
  new TH1F("DiTauPt_allcut", "Pt of the Di-Tau ", 140, 0, 140);
  new TH1F("diffPhi_Mu_DiTau_allcut", "delta Phi difference of Mu and Di-Tau", 100, 0, 2*pi);
  new TH1F("LT_allcut", "sum of Transverse Energy of the three leptons after all cuts", 200, 0, 200);
  new TH1F("vLT_allcut", "Vector eT sum of three leptons after all cuts", 200, 0, 200);
  new TH1F("var1_allcut", " plot of  Pt(Tau1+Tau2)/(tauiPt + tau2Pt)", 100, 0, 2);
  new TH1F("MT_Mu_Met_allcut","transverse mass plot for Muon and MET", 140, 0, 140);
  new TH1F("MT_Tau1_Met_allcut","transverse mass plot for Tau1 and MET", 140, 0, 140);
  new TH1F("MT_Tau2_Met_allcut","transverse mass plot for Tau2 and MET", 140, 0, 140);
}
// -------------------
// The main event loop
// -------------------
void MuTauTau::clearLists() {
  vtxList.clear();
  muoList.clear();
  eleList.clear();
  tauList.clear();
  bjetList.clear();
}
void MuTauTau::eventLoop() 
{
  // Initialize analysis
  if (!beginJob()) return;
  
  int nPrint = max(10000, nEvents()/1000);

  Options op;
  op.verbose = false;
  op.usesbit = true;  // Crucial
  op.printselected = false;

  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  ulong tbytes = 0;
  for (int ev = 0; ev < nEvents(); ++ev) {
    clearEvent();
    int lflag = chain()->LoadTree(ev); 
    tbytes += getEntry(lflag);    // returns total bytes read

    const Event& evt = event();

    histf()->cd();

    // PileUP weight
    puevWt_ = 1; // for data
    if (isMC() && usePUWt()) {
      int npu = 0;
      puevWt_ = wtPileUp(npu);
      AnaUtil::fillHist1D("npu", npu, 1.0);
      AnaUtil::fillHist1D("puwt", puevWt_, 1.0);

      vector<int> list = evt.trueNInt;
      AnaUtil::fillHist1D("trueNInt", (list.size() ? list.at(0) : 0), 1.0);
    }
    AnaUtil::fillHist1D("evcounter", 0, puevWt_);
    
    // Trigger selection
    if (useTrigger() && !isTriggered()) continue;
    AnaUtil::fillHist1D("evcounter", 1, puevWt_);
    
    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName())); 

    int run   = evt.run;
    int event = evt.event;
    int lumis = evt.lumis;
    // Show status of the run
    if (currentFile != lastFile) 
    cout << "Tree# " << setw(4) << chain()->GetTreeNumber()  
         << " ==> " << currentFile 
         << " <<< Run# " << run
         << " Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
    cout << "Tree# " << setw(4) << chain()->GetTreeNumber()  
         << " ==> " << chain()->GetCurrentFile()->GetName() 
         << " <<< Run# " << run 
         << " Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;

    if (logOption() > 0) 
    fLog() << "run: " << run
          << ", event: " << event
          << ", n_tau: "<< ntau()
	  << ", n_muon: "<< nmuon()
	  << ", n_jet: " << njet()
	  << ", n_vertex: " << nvertex()
	  << ", n_met: " << nmet()
	  << ", n_electron: " << nelectron()
          << endl;

    clearLists();
 
    op.verbose = (logOption() >> 1 & 0x1); 
    findVtxInfo(vtxList, op, fLog());
    int nvtx = vtxList.size();
    double vz = ( nvtx > 0 ) ? vtxList.at(0).z : -999;
    AnaUtil::fillHist1D("nvertex", nvtx, puevWt_);

    op.verbose = (logOption() >> 2 & 0x1); 
    findTauInfo(tauList, vz, op, fLog());

    op.verbose = (logOption() >> 3 & 0x1); 
    findMuonInfo(muoList, vz, op, fLog());

    op.verbose = (logOption() >> 4 & 0x1); 
    findElectronInfo(eleList, vz, op, fLog());

    op.verbose = (logOption() >> 5 & 0x1); 
    findJetInfo(bjetList, op, fLog());

    if (logOption() > 0)
    fLog() << "run: " << run
          << ", event: " << event
          << ", n_vertex_good: " << nvtx
          << ", n_muon_selected: " << muoList.size()
          << ", n_electron_selected: " << eleList.size()
          << ", n_tau_selected: " << tauList.size()
          << ", n_bjet_selected: " << bjetList.size()
          << endl;

    // Event Selection Starts here .....
    // presence of > 1 good vertex
    if (vtxList.size() < 1) continue;
    AnaUtil::fillHist1D("evcounter", 2, puevWt_);

    selectEvent();
  }  
  // Analysis is over
  endJob();
}
void MuTauTau::selectEvent()
{
  int nsmuon = muoList.size();
  // require atleast 1 muon 
  if (nsmuon < 1) return;
  AnaUtil::fillHist1D("evcounter", 3, puevWt_);

  // muon Pt
  const Muon& muo = muoList.at(0);

  // STAGE 1 : Just Before Using The Muon Pt Cut
  const MET& mt = metList()->at(0);

  AnaUtil::fillHist1D("muonPt1_stage1",  muo.pt, puevWt_);
  AnaUtil::fillHist1D("muonEta1_stage1", muo.eta, puevWt_);
  AnaUtil::fillHist1D("met_stage1", mt.met, puevWt_);

  if (muo.pt <= AnaUtil::cutValue(_evselCutMap, "ptMuon")) return;
  AnaUtil::fillHist1D("evcounter", 4, puevWt_);
  TLorentzVector M;
  M.SetPtEtaPhiE(muo.pt, muo.eta, muo.phi, muo.energy);

  // 1 hadronic tau (before mu-tau, tau-tau no-overlap requirement)
  int nstau = tauList.size();
  if (nstau < 1) return;
  AnaUtil::fillHist1D("evcounter", 5, puevWt_);
  AnaUtil::fillHist1D("taupt1_afterMuon", tauList.at(0).pt, puevWt_);

  // Find overlap with leading muon and between the tau's
  bool isGoodVtx; // no need to check 
  TVector3 vmu = findLeptonVtx(muo.vtxIndex, isGoodVtx);

  vector<Tau> fTauList;
  for (auto it  = tauList.begin(); it != tauList.end(); ++it) {
    const Tau& taua = (*it);

    TLorentzVector vta;
    vta.SetPtEtaPhiE(taua.pt, taua.eta, taua.phi, taua.energy);
    double drtm = AnaUtil::deltaR(vta, M); 
    if (drtm <= AnaUtil::cutValue(_evselCutMap, "drMuTau")) continue;

    double delz = abs(taua.zvertex - vmu.z()); 
    if (delz >= AnaUtil::cutValue(_evselCutMap, "dzMuTauVtx")) continue; 

    for (auto jt  = tauList.begin(); jt != tauList.end(); ++jt) {
      if (jt == it) continue;
      const Tau& taub = (*jt);

      TLorentzVector vtb;
      vtb.SetPtEtaPhiE(taub.pt, taub.eta, taub.phi, taub.energy);
      double drtt = AnaUtil::deltaR(vta, vtb); 
      if (drtt <= AnaUtil::cutValue(_evselCutMap, "drTauTau")) continue;

      double tau_zvdiff = abs(taua.zvertex - taub.zvertex);
      if (tau_zvdiff >= AnaUtil::cutValue(_evselCutMap, "dzTauTauVtx")) continue;

      fTauList.push_back(taua);
      break;
    }  
  }
  // Find the first hadronic tau (no mu-tau overlap)
  uint indx = 0;
  for (auto it  = fTauList.begin(); it != fTauList.end(); ++it,++indx) {
    const Tau& tau = (*it);
    if (tau.byTightCombinedIsolationDeltaBetaCorr > 0.5) break;
  }
  if (indx >= fTauList.size()) return; // the first tau not found
  AnaUtil::fillHist1D("evcounter", 6, puevWt_);

  const Tau& taua = fTauList.at(indx);
  // STAGE 2: Just Before Applying the Tau1 Pt Cut
  AnaUtil::fillHist1D("tauPt1_stage2", taua.pt, puevWt_);
  AnaUtil::fillHist1D("tauEta1_stage2", taua.eta, puevWt_);
  AnaUtil::fillHist1D("met_stage2", mt.met, puevWt_);
  AnaUtil::fillHist1D("tauMVA1", taua.againstElectronMVA, puevWt_);
  AnaUtil::fillHist1D("tauZVertex1", taua.zvertex, puevWt_);

  // Highest pT Tau
  if (taua.pt <= AnaUtil::cutValue(_evselCutMap, "ptTau1")) return;
  AnaUtil::fillHist1D("evcounter", 7, puevWt_);

  // STAGE 3: Just After Tau1 Pt Cut
  TLorentzVector T1;
  T1.SetPtEtaPhiE(taua.pt, taua.eta, taua.phi, taua.energy);
  double drM1T1 = AnaUtil::deltaR(M, T1);
  AnaUtil::fillHist1D("muonPt1_stage3", muo.pt, puevWt_);
  AnaUtil::fillHist1D("muonEta1_stage3", muo.eta, puevWt_);
  AnaUtil::fillHist1D("deltaR_mutau1_stage3", drM1T1, puevWt_);
  AnaUtil::fillHist1D("met_stage3", mt.met, puevWt_);

  // 2nd hadronic tau
  uint jndx = 0;
  for (auto it  = fTauList.begin(); it != fTauList.end(); ++it,++jndx) {
    if (jndx == indx) continue;

    const Tau& tau = (*it);
    if (tau.byMediumCombinedIsolationDeltaBetaCorr > 0.5 &&
        tau.againstElectronMedium > 0.5) break;
  }

  if (jndx >= fTauList.size()) return; // the second tau not found
  AnaUtil::fillHist1D("evcounter", 8, puevWt_);

  const Tau& taub = fTauList.at(jndx);

  TLorentzVector T2;
  T2.SetPtEtaPhiE(taub.pt, taub.eta, taub.phi, taub.energy);  
  double drM1T2 = AnaUtil::deltaR(M, T2);
  double drT1T2 = AnaUtil::deltaR(T1, T2);
  AnaUtil::fillHist1D("deltaR_mutau2_stage3", drM1T2, puevWt_);
  AnaUtil::fillHist1D("deltaR_tau1tau2_stage3", drT1T2, puevWt_);
  AnaUtil::fillHist1D("tauPt2_stage3", taub.pt, puevWt_);
  AnaUtil::fillHist1D("tauEta2_stage3", taub.eta, puevWt_);
  TLorentzVector t = T1 + T2;
  double invmass = t.M();
  AnaUtil::fillHist1D("tautauInvmass_stage3", invmass, puevWt_); 

  // Scalar sum
  double lt = T1.Et() + T2.Et() + M.Et();
  AnaUtil::fillHist1D("LT_stage3", lt, puevWt_);

  // Vector sum
  TLorentzVector V = T1 + T2 + M;
  AnaUtil::fillHist1D("vLT_stage3", V.Et(), puevWt_);

  double dphi = AnaUtil::deltaPhi(M.Phi(), mt.metphi); // acoplanarity
  double dphiM1T1 = AnaUtil::deltaPhi(M.Phi(), T1.Phi());
  AnaUtil::fillHist1D("diffPhi_Mu_Met_stage3", dphi, puevWt_);
  AnaUtil::fillHist1D("diffPhi_Mu_Tau1_stage3", dphiM1T1, puevWt_);
  AnaUtil::fillHist1D("DiTauPt_stage3", t.Pt(), puevWt_);
  AnaUtil::fillHist1D("diffPhi_Mu_DiTau_stage3", AnaUtil::deltaPhi(M, t), puevWt_);

  AnaUtil::fillHist1D("tauMVA2", taub.againstElectronMVA, puevWt_);
  AnaUtil::fillHist1D("tauZVertex2", taub.zvertex, puevWt_);

  // Common vertex
  double vz = vtxList.at(0).z;
  double evz = pow(abs(taua.zvertex - vz), 2)
             + pow(abs(taub.zvertex - vz), 2)
             + pow(abs(vmu.z() - vz), 2);
  AnaUtil::fillHist1D("eventZ", sqrt(evz), puevWt_);

  if (taub.pt <= AnaUtil::cutValue(_evselCutMap, "ptTau2")) return;
  AnaUtil::fillHist1D("evcounter", 9, puevWt_);

  // STAGE 4: Just After Tau2 Pt Cut
  AnaUtil::fillHist1D("deltaR_mutau2_stage4", drM1T2, puevWt_);
  AnaUtil::fillHist1D("tautauInvmass_stage4", invmass, puevWt_); 
  AnaUtil::fillHist1D("met_stage4", mt.met, puevWt_);
  AnaUtil::fillHist1D("deltaR_tau1tau2_stage4", drT1T2, puevWt_);

  AnaUtil::fillHist1D("diffPhi_Mu_Met_stage4", dphi, puevWt_);
  AnaUtil::fillHist1D("diffPhi_Mu_Tau1_stage4", dphiM1T1, puevWt_);
  AnaUtil::fillHist1D("DiTauPt_stage4", t.Pt(), puevWt_);
  AnaUtil::fillHist1D("diffPhi_Mu_DiTau_stage4", AnaUtil::deltaPhi(M, t), puevWt_);
  AnaUtil::fillHist1D("LT_stage4", lt, puevWt_);
  AnaUtil::fillHist1D("vLT_stage4", V.Et(), puevWt_);

  double var = t.Pt()/(T1.Pt() + T2.Pt());
  AnaUtil::fillHist1D("var1_stage4", var, puevWt_); 
  
  double sipt = pow(taua.ltsipt, 2) + pow(taub.ltsipt, 2); 
  AnaUtil::fillHist1D("siptSum", sqrt(sipt), puevWt_); 

  if ( (taua.charge + taub.charge) != 0) return;
  AnaUtil::fillHist1D("evcounter", 10, puevWt_);
  
  // muon veto
  int nm = vetoMuon(taua.zvertex, AnaUtil::cutValue(_evselCutMap, "muVetoPt"), 
		    AnaUtil::cutValue(_evselCutMap, "dzMuTauVtx"));
  if (nm != 1) return;
  AnaUtil::fillHist1D("evcounter", 11, puevWt_);

  // electron veto
  int nele = vetoElectron(taua.zvertex, AnaUtil::cutValue(electronCutMap(), "pt"),
			  AnaUtil::cutValue(_evselCutMap, "dzEleTauVtx"));
  if (nele > 0) return;
  AnaUtil::fillHist1D("evcounter", 12, puevWt_);

  // b-tagged jet veto
  if (bjetList.size() > 0) return;
  AnaUtil::fillHist1D("evcounter", 13, puevWt_);

  // mass and tau pait pT veto
  TLorentzVector a = M + T1;
  double muTauMass = a.M();
  double tauPairPt = t.Pt();
  if (muTauMass  < AnaUtil::cutValue(_evselCutMap, "muTauMass") && 
      tauPairPt < AnaUtil::cutValue(_evselCutMap, "tauPairPt")) return;
  AnaUtil::fillHist1D("evcounter", 14, puevWt_);

  AnaUtil::fillHist1D("nMuon_allcut", nsmuon, puevWt_);
  AnaUtil::fillHist1D("muonPt1_allcut", muo.pt, puevWt_);
  AnaUtil::fillHist1D("muonEta1_allcut", muo.eta, puevWt_);
  AnaUtil::fillHist1D("tauPt1_allcut", taua.pt, puevWt_);
  AnaUtil::fillHist1D("tauEta1_allcut", taua.eta, puevWt_);
  AnaUtil::fillHist1D("tauPt2_allcut", taub.pt, puevWt_);
  AnaUtil::fillHist1D("tauEta2_allcut", taub.eta, puevWt_);
  AnaUtil::fillHist1D("tautauInvmass_allcut", invmass, puevWt_); 
  AnaUtil::fillHist1D("met_allcut", mt.met, puevWt_);
  AnaUtil::fillHist1D("diffPhi_Mu_Met_allcut", dphi, puevWt_);
  AnaUtil::fillHist1D("diffPhi_Mu_Tau1_allcut", dphiM1T1, puevWt_);
  AnaUtil::fillHist1D("DiTauPt_allcut", t.Pt(), puevWt_);
  AnaUtil::fillHist1D("diffPhi_Mu_DiTau_allcut", AnaUtil::deltaPhi(M, t), puevWt_);
  AnaUtil::fillHist1D("LT_allcut", lt, puevWt_);
  AnaUtil::fillHist1D("vLT_allcut", V.Et(), puevWt_);
  AnaUtil::fillHist1D("var1_allcut", var, puevWt_); 
  double mass1 = sqrt(pow(M.Pt() + mt.met,2) - pow(M.Pt()*cos(M.Phi()) + mt.met * cos(mt.metphi),2) - 
                      pow(M.Pt()*sin(M.Phi()) + mt.met* sin(mt.metphi),2));   
  double mass2 = sqrt(pow(T1.Pt() + mt.met,2) - pow(T1.Pt()*cos(T1.Phi()) + mt.met * cos(mt.metphi),2) - 
                      pow(T1.Pt()*sin(T1.Phi()) + mt.met* sin(mt.metphi),2));   
  double mass3 = sqrt(pow(T2.Pt() + mt.met,2) - pow(T2.Pt()*cos(T2.Phi()) + mt.met * cos(mt.metphi),2) - 
                      pow(T2.Pt()*sin(T2.Phi()) + mt.met* sin(mt.metphi),2));   
  AnaUtil::fillHist1D("MT_Mu_Met_allcut", mass1, puevWt_);
  AnaUtil::fillHist1D("MT_Tau1_Met_allcut", mass2, puevWt_);
  AnaUtil::fillHist1D("MT_Tau2_Met_allcut", mass3, puevWt_);

  if (_skimObj) {
    TreeVariables varList;
    varList.muEta      = M.Eta();
    varList.muPt       = M.Pt();
    varList.tau1Eta    = T1.Eta();
    varList.tau1Pt     = T1.Pt();
    varList.tau2Eta    = T2.Eta();
    varList.tau2Pt     = T2.Pt();
    varList.diTauEta   = t.Eta();
    varList.diTauPt    = t.Pt();
    varList.dphiMuTau1 = AnaUtil::deltaPhi(M, T1);
    varList.met        = mt.met;
    TVector3 Mu_Alpha(M.Px(), M.Py(), M.Pz());
    TVector3 Z_Alpha(0, 0, 1.0);
    TVector3 d1(T1.Px(), T1.Py(), T2.Pz());
    TVector3 d2(T2.Px(), T2.Py(), T2.Pz());
    TVector3 ditauCrossProd = d1.Cross(d2);
    TVector3 muZ = Z_Alpha.Cross(Mu_Alpha);
    double angle = muZ.Angle(ditauCrossProd);
    varList.alpha       = angle; 
    varList.acop        = dphi;
    varList.dphiMuDiTau = AnaUtil::deltaPhi(M, t);
    _skimObj->fill(varList);
  }
  // Dump seleted event information
  // Event and object information
  dumpEvent("11111", evLog(), false);

  // event selection information
  evLog() << setprecision(3);
  evLog() << "nVertex = " << vtxList.size() << endl
         << "nMuon = " << nsmuon << endl
         << "pT Muon 1 = " << setw(8) << muo.pt << " GeV" << endl;
  evLog() << "nTau = " << nstau << endl
         << "pT Tau 1 = " << setw(8) << taua.pt << " GeV" << endl
         << "pT Tau 2 = " << setw(8) << taub.pt << " GeV" << endl
         << "dR(m,t1) = " << setw(7) << drM1T1 << endl
         << "dR(m,t2) = " << setw(7) << drM1T2 << endl
         << "dR(t1,t2) = " << setw(7) << drT1T2 << endl
         << "mass(t1,t2) = " << setw(8) << invmass << " GeV" << endl
         << "MET = " << setw(8) << mt.met << " GeV" << endl
         << "Mu-Tau1 mass = " << setw(8) << muTauMass  << " GeV" << endl
         << "Tau pair pT = " << setw(8) << tauPairPt  << " GeV"
         << endl;
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void MuTauTau::endJob() 
{  
  ostringstream ostra, ostrb, ostrc, ostrd, ostre, ostrf;
  ostra << "Leading Muon(Id, Iso, pT > " 
        <<  AnaUtil::cutValue(_evselCutMap, "ptMuon") << " GeV)";
  ostrb << ">= 1 Taus(Id, Iso, dR(M,T)>" 
	<< AnaUtil::cutValue(_evselCutMap, "drMuTau") 
        << ", dzv(M,T)<"
	<< AnaUtil::cutValue(_evselCutMap, "dzMuTauVtx") 
        << ", dR(T,T)>"
	<< AnaUtil::cutValue(_evselCutMap, "drTauTau") 
        << ", dzv(T,T)<"
	<< AnaUtil::cutValue(_evselCutMap, "dzTauTauVtx") 
      	<< ")";
  ostrc << "Highest pT Tau pT > " 
        << AnaUtil::cutValue(_evselCutMap, "ptTau1") << " GeV";
  ostre << "Second Highest pT Tau pT > " 
        << AnaUtil::cutValue(_evselCutMap, "ptTau2") << " GeV";
  ostrf << "!(Mu-Tau1 mass >= " 
        << AnaUtil::cutValue(_evselCutMap, "muTauMass") 
        << " GeV && Tau pair pT >= "
        << AnaUtil::cutValue(_evselCutMap, "tauPairPt")
        << " GeV)";
  string slist[] = 
  {
    "Total events processed",
    "After Trigger selection",
    ">= 1 good Event Vertex",
    "> 0 selected muons (Id, Iso)",
    ostra.str(),
    "> 0 Taus(Id, Iso)",
    ostrb.str(),
    ostrc.str(),
    ">= 2 Taus(Id, Iso, dR, dzv)", 
    ostrd.str(),
    "Tau pair charge = 0",
    "no extra Muon",
    "Electron veto",
    "b-jet veto",
    ostrf.str()
  };
  int prec = (isMC()) ? 1 : 0;
  fLog() << setprecision(prec);
  fLog() << "============================================================" << endl
        << "Statistics for W->munu, H->tau+tau-; tau->hadron,tau->hadron" << endl
        << "============================================================" << endl;

  TH1 *h = AnaUtil::getHist1D("evcounter");
  if (h) {
    for (int i = 1; i < h->GetNbinsX()+1; ++i)
      fLog() << setw(80) << slist[i-1] 
            << setw(12) << h->GetBinContent(i)
            << endl;
  }
  fLog() << resetiosflags(ios::fixed);

  closeFiles();

  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();

  if (_skimObj) _skimObj->close();
}
// -------------------------------------------------------------------------------
// Poor man's way of a datacard. Each line between the 'START' and 'END' tags
// is read in turn, split into words, where the first element is the 'key' and
// the rest the value(s). If more than one values are present they are 
// stored in a vector. No safety mechanism is in place. Any line with an unknown 
// key is skipped. Comments lines should start with either '#' or '//', preferably
// in the first column. Empty lines are skipped. The file containing the datacards 
// is passed as the only argument of the program, there is no default
// -------------------------------------------------------------------------------
bool MuTauTau::readJob(const string& jobFile, int& nFiles)
{
  if (!AnaBase::readJob(jobFile, nFiles)) return false;

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
  hmap.insert(pair<string, map<string, double>* >("evselCutList", &_evselCutMap));

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
    assert(tokens.size() > 1);
    string key = tokens[0];
    string value = tokens[1];
    if (key == "createMVATree")
      _createMVATree = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "mvaInputFile")
      _mvaInputFile = value;
    else if (key == "evselCutList")
      AnaUtil::storeCuts(tokens, hmap);

    tokens.clear();
  }
  // Close the file
  fin.close();

  printJob();

  return true;
}
void MuTauTau::printJob(ostream& os) const
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("evselCutList", _evselCutMap));
  AnaUtil::showCuts(hmap, os);
}
