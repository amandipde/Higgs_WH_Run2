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

#include "ElectronTriggerEfficiency.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"

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
ElectronTriggerEfficiency::ElectronTriggerEfficiency()
  : AnaBase(),
    _prbTrigPath(""),
    _dumpEvent(false),
    _matchTau(true),
    _iFlagL1(-1),
    _iFlagL2(-1)
{
  _tagTrigPathList.clear();  
}
// ----------
// Destructor
// ----------
ElectronTriggerEfficiency::~ElectronTriggerEfficiency() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool ElectronTriggerEfficiency::beginJob() 
{ 
  AnaBase::beginJob();

  // Open the output ROOT file
  histf()->cd();
  bookHistograms();

  return true;
}
// ---------------
// Book histograms
// ---------------
void ElectronTriggerEfficiency::bookHistograms() 
{
  new TH1D("evcounter", "Selected event counter", 11, -0.5, 10.5);
  new TH1F("invMass","Invariant Mass of Tag and Probe",100,30.0,130.0);
  new TH1F("invMassPass","Invariant Mass of Tag and Probe(Events Selected By Trigger)",100,30.0,130.0);
  new TH1F("invMassFail","Invariant Mass of Tag and Probe(Events Failed By Trigger)",100,30.0,130.0);
  new TH1F("invMassTrigPass","Invariant Mass of Tag and Probe TriggerObject",100,30.0,130.0);
  new TH1F("invMassTrigFail","Invariant Mass of Tag and Probe TriggerObject(failed events)",100,30.0,130.0);
  new TH1F("deltaZVtx", "Z-distance Vertices of Tag and Probe", 100, -0.01, 0.01);

  int nBin = 13;
  Double_t xbins[]= {0.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 65.0, 75.0, 85.0, 100.0};

  new TH1F("tagPt","Pt of Tag", 100, -0.5, 199.5);
  new TH1F("tagPtv","Pt of Tag (var. bin)", nBin, xbins);
  new TH1F("tagEta","Eta of Tag", 84, -2.1, 2.1);
  new TH1F("tagPhi", "Phi of Tag", 128, -3.2, 3.2);
  new TH1F("tagDeltaR", "DeltaR of Tag wrt Matched Trigger Object", 250, 0.0, 0.5);
  new TH1F("tagDeltaPt", "DeltaPt of Tag wrt Matched Trigger Object", 100, -25.0, 25.0);
  new TH1F("trigTagPt","Pt of Trigger Object matched with Tag", 100, -0.5, 199.5);
  new TH1F("trigTagEta","Eta of Trigger Object matched with Tag", 84, -2.1, 2.1);
  new TH1F("trigTagPhi","Phi of Trigger Object matched with Tag", 128, -3.2, 3.2);
  new TH1F("probePt","Pt of Probe ", 100, -0.5, 199.5);
  new TH1F("probePtv","Pt of Probe(var. bin)", nBin, xbins);
  new TH1F("probeEta","Eta of Probe ", 84, -2.1, 2.1);
  new TH1F("probePhi","Phi of Probe ", 128, -3.2, 3.2);
  new TH1F("probeDeltaR", "DeltaR of Probe wrt Matched Trigger Object", 250, 0.0, 0.5);
  new TH1F("probeDeltaPt", "DeltaPt of Probe wrt Matched Trigger Object", 100, -25.0, 25.0);
  new TH1F("matchedProbePt","Pt of Trigger Object matched with Probe",100,-0.5,199.5);
  new TH1F("matchedProbePtv","Pt of Trigger Object matched with Probe(var. bin)", nBin, xbins);
  new TH1F("matchedProbeEta","Eta of Trigger Object matched with Probe",84,-2.1,2.1);
  new TH1F("matchedProbePhi","Phi of Trigger Object matched with Probe", 128,-3.2,3.2);

  new TH2F("tagVsTrigPt","Pt Corr Tag Vs Trigger",100,-0.5,199.5,100,-0.5,199.5);
  new TH2F("tagVsTrigEta","Eta Corr Tag Vs Trigger",84,-2.1,2.1,84,-2.1,2.1);
  new TH2F("tagVsTrigPhi","Phi Corr Tag Vs Trigger",128,-3.2,3.2,128,-3.2,3.2);
  new TH2F("tagDRVsPt","Tag dR Vs Pt",50,0.0,100.0, 250, 0.0, 0.5);
  new TH2F("tagDRVsInvMass","Tag dR Vs Invariant Mass",100, 30.0, 130, 250, 0.0, 0.5);
  new TH2F("probeVsTrigPt","Pt Corr Probe Vs Trigger",100,-0.5,199.5,100,-0.5,199.5);
  new TH2F("probeVsTrigEta","Eta Corr Probe Vs Trigger",84,-2.1,2.1,84,-2.1,2.1);
  new TH2F("probeVsTrigPhi","Phi Corr Probe Vs Trigger",128,-3.2,3.2,128,-3.2,3.2);
  new TH2F("probeDRVsPt","ProbedR Vs Pt", 50, 0.0, 100.0, 250, 0.0, 0.5);
  new TH2F("probeDRVsInvMass","Probe dR Vs Invariant Mass",100, 30.0, 130, 250, 0.0, 0.5);
  new TH2F("probeVsTagDz", "Probe dZ Vs Tag dZ",100,-0.1,0.1,100,-0.1,0.1);
  new TH2F("probeVsTagDR", "Probe dR Vs Tag dR",250,0.0,0.5,250,0.0,0.5);
  new TH1F("probeDRFail", "Probe DR (failed events)", 400, 0.0,4.0); 
  new TH2F("probeVsTagDRFail", "Probe dR Vs Tag DR (failed events)", 400, 0.0,4.0, 400, 0.0,4.0);
  new TH1F("probeTrigDeltaPtFail", "deta Pt of Probe and TrigObject (failed events)", 100, -25.0, 25.0);
  new TH2F("probeVsTrigPtFail","Pt Corr Probe Vs Trigger (failed events)",100,-0.5,199.5,100,-0.5,199.5);    
  new TH2F("probeVsTrigEtaFail","Eta Corr Probe Vs Trigger (failed events)",84,-2.1,2.1,84,-2.1,2.1);
  new TH2F("probeVsTrigPhiFail","Phi Corr Probe Vs Trigger (failed events)",128,-3.2,3.2,128,-3.2,3.2);

  new TH2F("invMassVsPtPass","Invariant Mass of T&P vs pT(passed events)", nBin, xbins, 100,30.0,130.0);
  new TH2F("invMassVsPtFail","Invariant Mass of T&P Vs pT(failed events)", nBin, xbins, 100,30.0,130.0);
}

// -------------------
// The main event loop
// -------------------
void ElectronTriggerEfficiency::clearLists() {
  vtxList.clear();
  muoList.clear();
  eleList.clear();
  tauList.clear();
  bjetList.clear();
  trigObjList.clear();
}
void ElectronTriggerEfficiency::eventLoop() 
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
    tbytes += getEntry(lflag);    // returns bytes read
    
    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName())); 

    const Event& evt = event();

    int run   = evt.run;
    int event = evt.event;
    int lumis = evt.lumis;

    // Pile-up reweight
    _puevWt = 1; // for data
    if (_isMC && _usePUWt) {
      int npu = 0;
      _puevWt = wtPileUp(npu);
    }

    // Show status of the run
    if (currentFile != lastFile) 
    cout << "Tree# " << setw(4) << chain()->GetTreeNumber()  
         << " ==> " << currentFile 
         << " <<< Run# " << run
         << " <<< Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
    cout << "Tree# " << setw(4) << chain()->GetTreeNumber()  
         << " ==> " << chain()->GetCurrentFile()->GetName() 
         << " <<< Run# " << run 
         << " <<< Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;

    if (_logOption > 0) 
    fLog() << "run: " << run
          << ", event: " << event
          << ", n_tau: "<< n_tau
          << ", n_muon: "<< n_muon
          << ", n_jet: " << n_jet
          << ", n_vertex: " << n_vertex
          << ", n_met: " << n_met
          << ", n_electron: " << n_electron 
          << endl;

    AnaUtil::fillHist1D("evcounter", 0, _puevWt);

    // Trigger selection, do not check prescale
    if (_useTrigger && !isTriggered(false)) continue;
    AnaUtil::fillHist1D("evcounter", 1, _puevWt);

    clearLists();
    op.verbose = (_logOption >> 1 & 0x1); 
    findVtxInfo(vtxList, op, fLog());
    int nvtx = vtxList.size();
    double vz = ( nvtx > 0 ) ? vtxList.at(0).z : -999;

    op.verbose = (_logOption >> 4 & 0x1); 
    findElectronInfo(eleList, vz, op, fLog());

    op.verbose = (_logOption >> 2 & 0x1); 
    findTauInfo(tauList, vz, op, fLog());

    op.verbose = (_logOption >> 3 & 0x1); 
    findMuonInfo(muoList, vz, op, fLog());

    op.verbose = (_logOption >> 5 & 0x1); 
    findJetInfo(bjetList, op, fLog());

    findTriggerObjectInfo(trigObjList);

    if (_logOption > 0)
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
    AnaUtil::fillHist1D("evcounter", 2, _puevWt);

    fillTriggerFlags();
    selectEvent();
  }  
  // Analysis is over
  endJob();
}
// ----------------------------------------------------------
// Perform event selection, For selection of Z -> e+e- events
// we need,
//   - > 0 Tight electron
//   - event within e+e- invariant mass window
// ----------------------------------------------------------
void ElectronTriggerEfficiency::selectEvent() 
{
  // may also use nele >= 2
  if (eleList.size() != 2) return;
  AnaUtil::fillHist1D("evcounter", 3, _puevWt);

  if (_matchTau && tauList.size() < 1) return;
  AnaUtil::fillHist1D("evcounter", 4, _puevWt);

  if (muoList.size() || bjetList.size()) return;
  AnaUtil::fillHist1D("evcounter", 5, _puevWt);

  const MET& mt = metColl()[0];
  assert(mt);
  if (mt->met > AnaUtil::cutValue(_evselCutMap, "maxMET")) return;
  AnaUtil::fillHist1D("evcounter", 6, _puevWt);

  if (_dumpEvent) {
    dumpEvent("01010", fLog(), true);
    dumpTriggerPaths(fLog(), false);
    dumpTriggerObjectInfo(trigObjList, fLog());
  }

  int ntobj = trigObjList.size();
  double maxPtDiff = AnaUtil::cutValue(_evselCutMap, "maxPtDiff");

  // Find the Tag electron using offline information only
  // In case of TauPlusX match electron with a Tau (tag)
  // use an option _matchTau
  int nmatch[] = {0,0,0,0};
  random_shuffle (eleList.begin(), eleList.end());
  for (auto it  = eleList.begin(); it != eleList.end(); ++it) {
    const Electron& elec_t = (*it);
    if (elec_t.charge != AnaUtil::cutValue(_evselCutMap, "charge")) continue;
    ++nmatch[0];

    TLorentzVector taglv;
    taglv.SetPtEtaPhiE(elec_t.pt, elec_t.eta, elec_t.phi, elec_t.energy);

    // find if it is a tag candidate, we do not use it for DYToEE
    if (_matchTau) {
      double drmin = 999;
      double dpt = -1;
      for (auto jt = tauList.begin(); jt != tauList.end(); ++jt) {
        const Tau& tau = (*jt);
        TLorentzVector lvt;
        lvt.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
        double dr = AnaUtil::deltaR(taglv, lvt);
        if (dr < drmin) {
          drmin = dr;
          dpt = elec_t.pt - tau.pt;
        }
      }
      if (drmin >= AnaUtil::cutValue(_evselCutMap, "maxDrEleTau") ||
          abs(dpt) >= maxPtDiff) continue; // not a tag candidate
    }
    ++nmatch[1];

    // Ensure that the tag electron passes the event trigger path 
    // stop after the first trigger object match
    bool trigger_match_ev = false;
    for (vector<string>::const_iterator kt  = _trigPathList.begin();
	                                kt != _trigPathList.end(); ++kt) {
      const string& patt = (*kt);
      int tindx = -1;
      double drTag = matchTriggerObject(trigObjList, taglv, patt, -1, maxPtDiff, tindx);
      if (tindx < 0 || tindx >= ntobj) continue;
      if (drTag >= AnaUtil::cutValue(_evselCutMap, "maxDrTag")) continue;
      trigger_match_ev = true;
      break;
    }
    if (!trigger_match_ev) continue;
    ++nmatch[2];

    // Now consider (a) cross-trigger path(s) for the tag
    // stop after the first trigger object match
    bool trigger_match_x = false;
    double drTag = 999;
    int tindx = -1;
    for (auto kt = _tagTrigPathList.begin(); kt != _tagTrigPathList.end(); ++kt) {
      const string& patt = (*kt);
      drTag = matchTriggerObject(trigObjList, taglv, patt, -1, maxPtDiff, tindx);
      if (tindx < 0 || tindx >= ntobj) continue;
      if (drTag >= AnaUtil::cutValue(_evselCutMap, "maxDrTag")) continue;
      trigger_match_x = true;
      break;
    }
    if (!trigger_match_x) continue;
    ++nmatch[3];

    // Tag found
    AnaUtil::fillHist1D("tagDeltaR", drTag, _puevWt);
    AnaUtil::fillHist2D("tagDRVsPt", elec_t.pt, drTag, _puevWt);

    AnaUtil::fillHist1D("tagPt",  elec_t.pt, _puevWt);
    AnaUtil::fillHist1D("tagPtv", elec_t.pt, _puevWt);
    AnaUtil::fillHist1D("tagEta", elec_t.eta, _puevWt);
    AnaUtil::fillHist1D("tagPhi", elec_t.phi, _puevWt);
    
    const TriggerObject& tagTrigObj = trigObjList.at(tindx); 
    AnaUtil::fillHist1D("trigTagPt",    tagTrigObj.pt,  _puevWt);
    AnaUtil::fillHist1D("trigTagEta",   tagTrigObj.eta, _puevWt);
    AnaUtil::fillHist1D("trigTagPhi",   tagTrigObj.phi, _puevWt);
    AnaUtil::fillHist2D("tagVsTrigPt",  elec_t.pt,  tagTrigObj.pt,  _puevWt);
    AnaUtil::fillHist2D("tagVsTrigEta", elec_t.eta, tagTrigObj.eta, _puevWt);
    AnaUtil::fillHist2D("tagVsTrigPhi", elec_t.phi, tagTrigObj.phi, _puevWt);

    double drPtTag = elec_t.pt - tagTrigObj.pt;
    AnaUtil::fillHist1D("tagDeltaPt", drPtTag, _puevWt);

    TLorentzVector tagTrigTL;
    tagTrigTL.SetPtEtaPhiE(tagTrigObj.pt, tagTrigObj.eta, tagTrigObj.phi, tagTrigObj.energy);

    bool igtvtx; // status may not be checked, already done earlier
    TVector3 vtxtag = findLeptonVtx(elec_t.vtxIndex, igtvtx);

    // Now find the probe
    //    bool probe_found = false;
    for (auto jt  = eleList.begin(); jt != eleList.end(); ++jt) {
      if (jt == it) continue;
      const Electron& elec_p = (*jt);

      bool igpvtx;  // no need to check status, already done earlier
      TVector3 vtxprb = findLeptonVtx(elec_p.vtxIndex, igpvtx);

      double delz = vtxtag.z() - vtxprb.z();
      AnaUtil::fillHist1D("deltaZVtx", delz, _puevWt);
      if (abs(delz) > AnaUtil::cutValue(_evselCutMap, "maxDz")) continue;
      AnaUtil::fillHist2D("probeVsTagDz", elec_t.vtxDistZ, elec_p.vtxDistZ, _puevWt);

      TLorentzVector prblv;
      prblv.SetPtEtaPhiE(elec_p.pt, elec_p.eta, elec_p.phi, elec_p.energy);
      TLorentzVector lv = taglv + prblv;
      double mass = lv.M();
      AnaUtil::fillHist1D("invMass", mass, _puevWt);
      AnaUtil::fillHist2D("tagDRVsInvMass", mass, drTag, _puevWt);
      AnaUtil::fillHist2D("probeDRVsPt", elec_p.pt, drTag, _puevWt);

      // Charge 
      if (elec_t.charge + elec_p.charge != 0) continue;

      // Invariant mass
      double deltaM = abs(mass - 91.2);
      if (deltaM > AnaUtil::cutValue(_evselCutMap, "massWindowWide")) continue;

      // probe found
      if (deltaM < AnaUtil::cutValue(_evselCutMap, "massWindowNarrow")) {
        AnaUtil::fillHist1D("probePt",  elec_p.pt,  _puevWt);
        AnaUtil::fillHist1D("probePtv", elec_p.pt,  _puevWt);
        AnaUtil::fillHist1D("probeEta", elec_p.eta, _puevWt);
        AnaUtil::fillHist1D("probePhi", elec_p.phi, _puevWt);
      }
      if (_iFlagL2 == 1) 
        AnaUtil::fillHist1D("invMassPass", mass, _puevWt);
      else 
        AnaUtil::fillHist1D("invMassFail", mass, _puevWt);

      int pindx = -1;
      double drProbe = matchTriggerObject(trigObjList, prblv, _prbTrigPath, tindx, maxPtDiff, pindx);
      if (deltaM < AnaUtil::cutValue(_evselCutMap, "massWindowNarrow")) {
        AnaUtil::fillHist2D("probeDRVsInvMass", mass, drProbe, _puevWt); 
        AnaUtil::fillHist1D("probeDeltaR", drProbe, _puevWt);
        AnaUtil::fillHist2D("probeVsTagDR", drTag, drProbe, _puevWt);
      }
      if (pindx < 0 || pindx >= ntobj) {
        fLog() << "=> WARNING. No Trigger Object found for Probe! total # of TriggerObjects: " 
              << ntobj 
              << ", dr: " << drProbe 
              << ", pindx: " << pindx 
              << endl;
        continue;
      }
      assert(pindx != tindx);
      //probe_found = true; 

      const TriggerObject& probeTrigObj = trigObjList.at(pindx);
      TLorentzVector probeTrigTL;
      probeTrigTL.SetPtEtaPhiE(probeTrigObj.pt, probeTrigObj.eta, probeTrigObj.phi, probeTrigObj.energy);
      double massTrig = (tagTrigTL + probeTrigTL).M();
      double drPtProbe = elec_p.pt - probeTrigObj.pt;
      fLog() << "=> Summary:" 
            << setprecision(2)
            << " TagPt: "  << elec_t.pt
            << " TagTrigPt: " << tagTrigObj.pt 
            << " ProbePt: " << elec_p.pt 
            << " ProbeTrigPt: " << probeTrigObj.pt 
            << " Mass: " << mass 
            << " Mass Trig: " << massTrig 
            << setprecision(4)
            << " drTag: " << drTag 
            << " drProbe: " << drProbe 
            << endl;

      if (drProbe < AnaUtil::cutValue(_evselCutMap, "maxDrProbe")) {
        if (deltaM < AnaUtil::cutValue(_evselCutMap, "massWindowNarrow")) {
          AnaUtil::fillHist1D("matchedProbePt",  elec_p.pt,  _puevWt);
          AnaUtil::fillHist1D("matchedProbePtv", elec_p.pt,  _puevWt);
          AnaUtil::fillHist1D("matchedProbeEta", elec_p.eta, _puevWt);
          AnaUtil::fillHist1D("matchedProbePhi", elec_p.phi, _puevWt);
          AnaUtil::fillHist1D("probeDeltaPt",    drPtProbe,  _puevWt);
          AnaUtil::fillHist2D("probeVsTrigPt",   elec_p.pt,  probeTrigObj.pt, _puevWt);
          AnaUtil::fillHist2D("probeVsTrigEta",  elec_p.eta, probeTrigObj.eta, _puevWt);
          AnaUtil::fillHist2D("probeVsTrigPhi",  elec_p.phi, probeTrigObj.phi, _puevWt);
          AnaUtil::fillHist1D("invMassTrigPass", massTrig, _puevWt);
        }  
        AnaUtil::fillHist2D("invMassVsPtPass", elec_p.pt, mass, _puevWt);
        break;
      }
      else {
        if (deltaM < AnaUtil::cutValue(_evselCutMap, "massWindowNarrow")) {
          AnaUtil::fillHist1D("probeDRFail",          drProbe, _puevWt);
          AnaUtil::fillHist2D("probeVsTagDRFail",     drTag, drProbe, _puevWt);
          AnaUtil::fillHist1D("probeTrigDeltaPtFail", elec_p.pt - probeTrigObj.pt,  _puevWt);
          AnaUtil::fillHist2D("probeVsTrigPtFail",    elec_p.pt,  probeTrigObj.pt,  _puevWt);
          AnaUtil::fillHist2D("probeVsTrigEtaFail",   elec_p.eta, probeTrigObj.eta, _puevWt);
          AnaUtil::fillHist2D("probeVsTrigPhiFail",   elec_p.phi, probeTrigObj.phi, _puevWt);
          AnaUtil::fillHist1D("invMassTrigFail",      massTrig, _puevWt);
        }
        AnaUtil::fillHist2D("invMassVsPtFail", elec_p.pt, mass, _puevWt); 
      }
    }
    //if (probe_found) break;
  }
  for (uint i = 0; i < NEL(nmatch); ++i) {
    if (nmatch[i] > 0) AnaUtil::fillHist1D("evcounter", 7+i, _puevWt);
  }
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void ElectronTriggerEfficiency::endJob() {  
  fLog() << resetiosflags(ios::fixed);

  closeFiles();

  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();
}
void ElectronTriggerEfficiency::fillTriggerFlags(bool check_prescale, bool verbose) {
  int nFlag1 = 0,
      nFlag2 = 0;
  for (uint i = 0; i < hltpaths()->size(); ++i) {
    string path_name = (*hltpaths()).at(i);
    int prescale     = (check_prescale) ? (*hltprescales()).at(i) : 1;  
    int result       = (*hltresults()).at(i);  
    if (result != 1 || prescale != 1) continue;
    if (matchTriggerPath(_tagTrigPathList, path_name)) ++nFlag1;
    if (path_name.find(_prbTrigPath) != string::npos) ++nFlag2;
  }
  if (verbose && (nFlag1 > 1 || nFlag2 > 1)) 
    cout << "=> Multiple Match: nFlag1: " << nFlag1 
         << ", nFlag2: " << nFlag2 
         << endl;
  _iFlagL1 = (nFlag1 > 0) ? 1 : 0;
  _iFlagL2 = (nFlag2 > 0) ? 1 : 0;
}
bool ElectronTriggerEfficiency::readJob(const string& jobFile, int& nFiles)
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

    AnaUtil::tokenize(line, tokens);
    int vsize = tokens.size();
    assert(vsize > 1);
    string key = tokens.at(0);

    if (key == "tagTrigPathList")
      AnaUtil::buildList(tokens, _tagTrigPathList);
    else if (key == "probeTrigPath")
      _prbTrigPath = tokens[1];
    else if (key == "dumpEvent")
      _dumpEvent = atoi(tokens[1].c_str()) > 0 ? true : false;
    else if (key == "matchTau")
      _matchTau = atoi(tokens[1].c_str()) > 0 ? true : false;
    else if (key == "evselCutList")
      AnaUtil::storeCuts(tokens, hmap);

    tokens.clear();
  }
  // Close the file
  fin.close();

  printJob();

  return true;
}
void ElectronTriggerEfficiency::printJob(ostream& os) const
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("evselCutList", _evselCutMap));

  AnaUtil::showCuts(hmap, os);
}
