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

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "MuonTriggerEfficiency.h"
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

using vhtm::Event;
using vhtm::Vertex;
using vhtm::TriggerObject;
using vhtm::Muon;
using vhtm::Electron;
using vhtm::Tau;
using vhtm::Jet;

// -----------
// Constructor
// -----------
MuonTriggerEfficiency::MuonTriggerEfficiency()
  :  AnaBase(),
    _neventsMuTrig(0),
    _runSave(0),
    _lumiSave(0),
    _dumpEvent(false)
{
  _triggerPathLeg1.clear();  
  _triggerPathLeg2.clear();  

  initializeHistograms();
}
//--------------------------
// Initialize Histograms
//-------------------------
void MuonTriggerEfficiency::initializeHistograms() {
  _firstMuPt       = 0; 
  _secondMuPt      = 0; 

  _firstElPt       = 0; 
  _secondElPt      = 0; 

  _invMass         = 0;
  _invMassPass     = 0;
  _invMassFail     = 0;
  _invMassTrigPass = 0;
  _invMassTrigFail = 0;
  _deltaVtx        = 0;
  _invMassVsPtPass = 0;
  _invMassVsPtFail = 0;

  _tagPt           = 0;
  _tagPtv          = 0;
  _tagEta          = 0;
  _tagPhi          = 0;
  _tagDeltaR       = 0;
  _tagDeltaPt      = 0; 
  _trigTagPt       = 0;
  _trigTagEta      = 0;
  _trigTagPhi      = 0;
  _probePt         = 0;
  _probePtv        = 0;
  _probeEta        = 0;
  _probePhi        = 0;
  _probeDeltaR     = 0;
  _probeDeltaPt    = 0; 
  _matchedProbePt  = 0;
  _matchedProbePtv = 0;
  _matchedProbeEta = 0;
  _matchedProbePhi = 0;

  _probeDeltaRFail = 0;
  _probeTrigDeltaPtFail = 0;

  _tagVsTrigPt     = 0;
  _tagVsTrigEta    = 0;
  _tagVsTrigPhi    = 0;
  _tagDRVsPt       = 0; 
  _tagDRVsInvMass  = 0;
  _probeVsTrigPt   = 0;
  _probeVsTrigEta  = 0;
  _probeVsTrigPhi  = 0;
  _probeDRVsPt     = 0; 
  _probeDRVsInvMass= 0;
  _probeVsTagDz    = 0;
  _probeVsTagDR    = 0;
  _probeVsTagDRFail = 0;
  _probeVsTrigPtFail = 0;
  _probeVsTrigEtaFail = 0;
  _probeVsTrigPhiFail = 0;
}
//---------
// Destructor
// ----------
MuonTriggerEfficiency::~MuonTriggerEfficiency() 
{
  clearEvent();
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MuonTriggerEfficiency::beginJob() 
{
  if (!AnaBase::beginJob()) return false;

  _histf->cd();
  bookHistograms();

  return true;
}
// ---------------
// Book histograms
// ---------------
void MuonTriggerEfficiency::bookHistograms() 
{
  new TH1D("evcounter", "Selected event counter", 6, -0.5, 5.5);
  _firstMuPt       = new TH1F("firstMuPt","Pt of Most Energetic Muon",100,-0.5,199.5);
  _secondMuPt      = new TH1F("secondMuPt", "Pt of 2nd Most Energetic Muon", 100, -0.5, 199.5);
  _invMass         = new TH1F("InvMass","Invariant Mass of Tag and Probe",100,30.0,130.0);
  _invMassPass     = new TH1F("InvMassPass","Invariant Mass of Tag and Probe(Events Selected By Trigger)",100,30.0,130.0);
  _invMassFail     = new TH1F("InvMassFail","Invariant Mass of Tag and Probe(Events Failed By Trigger)",100,30.0,130.0);
  _invMassTrigPass = new TH1F("InvMassTrigPass","Invariant Mass of Tag and Probe TriggerObject",100,30.0,130.0);
  _invMassTrigFail = new TH1F("InvMassTrigFail","Invariant Mass of Tag and Probe TriggerObject(failed events)",100,30.0,130.0);
  _deltaVtx        = new TH1F("deltaZVtx", "Z-distance Vertices of Tag and Probe", 100, -0.01, 0.01);

  int nBin = 17;
  Double_t xmin = -1.5;
  Double_t xbins[nBin+1];
  for (int i = 0; i <nBin+1; i++) {
    if (i <= 11) xbins[i] = xmin + i*5.0;
    else  xbins[i] = xbins[i-1] + 10.0;
    std::cout << "Bin #" << i << "  Xval " << xbins[i] << std::endl;
  }
  _invMassVsPtPass = new TH2F("InvMassVsPtPass","Invariant Mass of Tag and Probe Vs Pt(passed events)",
			   nBin, xbins, 100,30.0,130.0);
  _invMassVsPtFail = new TH2F("InvMassVsPtFail","Invariant Mass of Tag and Probe Vs Pt(failed events)",
			   nBin, xbins, 100,30.0,130.0);
  _tagPt           = new TH1F("tagPt","Pt of Tag",100,-0.5,199.5);
  _tagPtv          = new TH1F("tagPtv","Pt of Tag (var. bin)",nBin, xbins);
  _tagEta          = new TH1F("tagEta","Eta of Tag",84,-2.1,2.1);
  _tagPhi          = new TH1F("tagPhi", "Phi of Tag",128,-3.2,3.2);
  _tagDeltaR       = new TH1F("tagDeltaR", "DeltaR of Tag wrt Matched Trigger Object", 250, 0.0, 0.5);
  _tagDeltaPt      = new TH1F("tagDeltaPt", "DeltaPt of Tag wrt Matched Trigger Object", 100, -25.0, 25.0);
  _trigTagPt       = new TH1F("trigTagPt","Pt of Trigger Object matched with Tag",100,-0.5,199.5);
  _trigTagEta      = new TH1F("trigTagEta","Eta of Trigger Object matched with Tag",84,-2.1,2.1);
  _trigTagPhi      = new TH1F("trigTagPhi","Phi of Trigger Object matched with Tag",128,-3.2,3.2);
  _probePt         = new TH1F("probePt","Pt of Probe ",100,-0.5,199.5);
  _probePtv        = new TH1F("probePtv","Pt of Probe(var. bin) ",nBin, xbins);
  _probeEta        = new TH1F("probeEta","Eta of Probe ",84,-2.1,2.1);
  _probePhi        = new TH1F("probePhi","Phi of Probe ", 128,-3.2,3.2);
  _probeDeltaR     = new TH1F("probeDeltaR", "DeltaR of Probe wrt Matched Trigger Object", 250, 0.0, 0.5);
  _probeDeltaPt    = new TH1F("probeDeltaPt", "DeltaPt of Probe wrt Matched Trigger Object", 100, -25.0, 25.0);
  _matchedProbePt  = new TH1F("matchedProbePt","Pt of Trigger Object matched with Probe",100,-0.5,199.5);
  _matchedProbePtv = new TH1F("matchedProbePtv","Pt of Trigger Object matched with Probe(var. bin)",nBin, xbins);
  _matchedProbeEta = new TH1F("matchedProbeEta","Eta of Trigger Object matched with Probe",84,-2.1,2.1);
  _matchedProbePhi = new TH1F("matchedProbePhi","Phi of Trigger Object matched with Probe", 128,-3.2,3.2);
    
  _tagVsTrigPt     = new TH2F("tagVsTrigPt","Pt Corr Tag Vs Trigger",100,-0.5,199.5,100,-0.5,199.5);
  _tagVsTrigEta    = new TH2F("tagVsTrigEta","Eta Corr Tag Vs Trigger",84,-2.1,2.1,84,-2.1,2.1);
  _tagVsTrigPhi    = new TH2F("tagVsTrigPhi","Phi Corr Tag Vs Trigger",128,-3.2,3.2,128,-3.2,3.2);
  _tagDRVsPt       = new TH2F("tagDRVsPt","Tag dR Vs Pt",50,0.0,100.0, 250, 0.0, 0.5);
  _tagDRVsInvMass  = new TH2F("tagDRVsInvMass","Tag dR Vs Invariant Mass",100, 30.0, 130, 250, 0.0, 0.5);
  _probeVsTrigPt   = new TH2F("probeVsTrigPt","Pt Corr Probe Vs Trigger",100,-0.5,199.5,100,-0.5,199.5);
  _probeVsTrigEta  = new TH2F("probeVsTrigEta","Eta Corr Probe Vs Trigger",84,-2.1,2.1,84,-2.1,2.1);
  _probeVsTrigPhi  = new TH2F("probeVsTrigPhi","Phi Corr Probe Vs Trigger",128,-3.2,3.2,128,-3.2,3.2);
  _probeDRVsPt     = new TH2F("probeDRVsPt","ProbedR Vs Pt",50,0.0,100.0, 250, 0.0, 0.5);
  _probeDRVsInvMass = new TH2F("probeDRVsInvMass","Probe dR Vs Invariant Mass",100, 30.0, 130, 250, 0.0, 0.5);
  _probeVsTagDz    = new TH2F("probeVsTagDz", "Probe dZ Vs Tag dZ",100,-0.1,0.1,100,-0.1,0.1);
  _probeVsTagDR    = new TH2F("probeVsTagDR", "Probe dR Vs Tag dR",250,0.0,0.5,250,0.0,0.5);
  _probeDeltaRFail = new TH1F("probeDRFailed", "Probe DR (failed events)", 400, 0.0,4.0); 
  _probeVsTagDRFail = new  TH2F("probeVsTagDRFailed", "Probe dR Vs Tag DR (failed events)", 400, 0.0,4.0, 400, 0.0,4.0);
  _probeTrigDeltaPtFail = new TH1F("probeTrigDeltaPtFailed", "deta Pt of Probe and TrigObject (failed events)", 100, -25.0, 25.0);
  _probeVsTrigPtFail = new TH2F("probeVsTrigPtFailed","Pt Corr Probe Vs Trigger (failed events)",100,-0.5,199.5,100,-0.5,199.5);    
  
  _probeVsTrigEtaFail  = new TH2F("probeVsTrigEtaFailed","Eta Corr Probe Vs Trigger (failed events)",84,-2.1,2.1,84,-2.1,2.1);
  _probeVsTrigPhiFail  = new TH2F("probeVsTrigPhiFailed","Phi Corr Probe Vs Trigger (failed events)",128,-3.2,3.2,128,-3.2,3.2);
}
// -------------------
// The main event loop
// -------------------
void MuonTriggerEfficiency::eventLoop() 
{
  // Initialize analysis
  if (!beginJob()) return;

  int nPrint = max(10000, nEvents/1000);

  Options op;
  op.verbose = false;
  op.usesbit = true;  // Crucial
  op.printselected = false;

  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  for (int ev = 0; ev < nEvents; ++ev) {
    clearEvent();
    int lflag = _chain->LoadTree(ev); 
    int nentries = getEntry(lflag);    // returns bytes read
    
    string currentFile(gSystem->BaseName(_chain->GetCurrentFile()->GetName())); 

    const Event* evt = dynamic_cast<Event*>(eventA->At(0));
    assert(evt);
    int run   = evt->run;
    int event = evt->event;
    int lumi  = evt->lumis; 

    // PileUP weight
    _puevWt = 1; // for data
    if (_isMC) {
      int npu = 0;
      _puevWt = wtPileUp(npu);
    }

    // Show status of the run
    if (currentFile != lastFile) 
      cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
           << " ==> " << currentFile 
           << " <<< Run# " << run
           << " Event# " << setw(12) << event << " >>> " 
           << " Events proc. " << setw(12) << ev << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
       cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
            << " ==> " << _chain->GetCurrentFile()->GetName() 
            << " <<< Run# " << run 
            << " Event# " << setw(12) << event << " >>> " 
            << " Events proc. " << setw(12) << ev << endl;

    if (_logOption >> 0 & 0x1) 
    _fLog << "run: " << run
          << ", event: " << event
          << ", n_tau: "<< n_tau
          << ", n_muon: "<< n_muon
          << ", n_jet: " << n_jet
          << ", n_vertex: " << n_vertex
          << ", n_triggerobj: " << n_triggerobj
          << ", n_met: " << n_met
          << ", n_electron: " << n_electron << endl;

    AnaUtil::fillHist1D("evcounter", 0, _puevWt);

    // Trigger selection, do not check prescale
    if (_useTrigger && !isTriggered(false)) continue;
    AnaUtil::fillHist1D("evcounter", 1, _puevWt);

    clearCollection();
    op.verbose = (_logOption >> 1 & 0x1); 
    findVtxInfo(vertexVec, op, _fLog);
    if (vertexVec.size() < 1) continue;
    AnaUtil::fillHist1D("evcounter", 2, _puevWt);
    double vz = vertexVec.at(0).z;

    op.verbose = (_logOption >> 4 & 0x1); 
    findElectronInfo(electronVec, vz, op, _fLog);

    op.verbose = (_logOption >> 2 & 0x1); 
    findTauInfo(tauVec, vz, op, _fLog);

    op.verbose = (_logOption >> 3 & 0x1); 
    findMuonInfo(muonVec, vz, op, _fLog);
    random_shuffle(muonVec.begin(), muonVec.end());

    op.verbose = (_logOption >> 5 & 0x1); 
    findJetInfo(jetVec, op, _fLog);

    findTriggerObjectInfo(trigObjVec);
    fillTriggerPrescaleHistograms(run,lumi);

    if (muonVec.size() > 0) _firstMuPt->Fill(muonVec.at(0).pt, _puevWt);
    if (muonVec.size() > 1) _secondMuPt->Fill(muonVec.at(1).pt, _puevWt);

    if (_logOption)
    _fLog << "run: " << run
          << ", event: " << event
          << ", n_muon_selected: " << muonVec.size()
          << endl;
    selectEvent();
  }  
  // Analysis is over
  endJob();
}
void MuonTriggerEfficiency::selectEvent() 
{
  int nMuonObj = muonVec.size();
  int nElecObj = electronVec.size();
  int nTauObj = tauVec.size();
  int nJetObj = jetVec.size();

  if (nMuonObj != 2) return;
  AnaUtil::fillHist1D("evcounter", 3, _puevWt);
  if (nTauObj) return;
  AnaUtil::fillHist1D("evcounter", 4, _puevWt);
  if (nElecObj || nJetObj) return;
  AnaUtil::fillHist1D("evcounter", 5, _puevWt);

  if (_dumpEvent) {
    dumpTriggerPaths(_fLog);
    dumpEvent("01010", _fLog, true);
    dumpTriggerObjectInfo(trigObjVec, _fLog);
  }
  fillTriggerTPHistograms();
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void MuonTriggerEfficiency::endJob() 
{
  map<TH1F*, TH1F*> histo_map;
  histo_map.insert(std::make_pair<TH1F*, TH1F*>(_probePtv, _matchedProbePtv));
  histo_map.insert(std::make_pair<TH1F*, TH1F*>(_probeEta, _matchedProbeEta));
  histo_map.insert(std::make_pair<TH1F*, TH1F*>(_probePhi, _matchedProbePhi));
  fillTriggerEfficiencies(histo_map);
  closeFiles();

  _histf->cd();
  _histf->Write();
  _histf->Close();
  delete _histf;
}
void MuonTriggerEfficiency::fillTriggerPrescaleHistograms(int run, int lumi) {
  _iFlagL1 = 0;
  _iFlagL2 = 0;
  _iPreScaleL1 = -1;
  _iPreScaleL2 = -1;
  int nFlag1 = 0;
  int nFlag2 = 0;

  for (uint i = 0; i < _hltpaths->size(); i++) {
    string path_name = (*_hltpaths).at(i);
    int prescale     = (*_hltprescales).at(i);  
    int result       = (*_hltresults).at(i);  
    if (matchTriggerPath(_triggerPathLeg1, path_name)) {
      if (prescale == 1) _iPreScaleL1 = 1;
      if (prescale != 1 && _runSave != run && _lumiSave != lumi) {
	_runSave  = run;
	_lumiSave = lumi;
      }
      if ( result == 1 )  {
	_iFlagL1 = 1;
	++nFlag1;
      }
    }
    if ( (matchTriggerPath(_triggerPathLeg1, path_name)) && ( result == 1 ) ) {
      if (prescale ==1) _iPreScaleL1 = 1;
      _iFlagL2 = 1;
      ++nFlag2;
    }
    for (vector<string>::iterator it  = _trigPathList.begin(); 
                                  it != _trigPathList.end(); it++) 
    {
      if (path_name.find((*it)) != string::npos) {
	string temp_str = path_name.substr(path_name.find("HLT_")+4);
	string trig_tag = temp_str.substr(0,temp_str.find_last_of("_v")-1);
	map<string, TriggerHistos>::iterator iPos = _triggerHistoMap.find(trig_tag);
	if (iPos == _triggerHistoMap.end()) bookTriggerPrescaleHistograms(trig_tag);
	iPos = _triggerHistoMap.find(trig_tag);
	if (iPos != _triggerHistoMap.end()) {
	  TProfile* tp = iPos->second.prescaleHist;
	  tp->Fill(run, prescale);
          if (prescale == 1 && result == 1) {
	    iPos->second.eventPassed++;
	    if (muonVec.size() > 0 )  iPos->second.firstMuPtAcc->Fill(muonVec.at(0).pt, _puevWt);
	    if (muonVec.size() > 1 )  iPos->second.secondMuPtAcc->Fill(muonVec.at(1).pt, _puevWt);
          }
	}      
      }
    }
  }
  if (nFlag1 > 1 || nFlag2 > 1) 
    std::cout << " Multiple Match :  nFlag1: " << nFlag1 << " nFlag2: " << nFlag2 << std::endl;
}
void MuonTriggerEfficiency::fillTriggerTPHistograms() {
  double maxPtDiff = AnaUtil::cutValue(_evselCutMap, "maxPtDiff");
  int ntobj = trigObjVec.size();
  for (uint iobj = 0; iobj < muonVec.size(); ++iobj) {
    const Muon& imuon = muonVec.at(iobj);     
    TLorentzVector tagv;
    tagv.SetPtEtaPhiE(imuon.pt, imuon.eta, imuon.phi, imuon.energy);

    // Find the Tag
    int tag_trig_indx = -1;
    uint trg_flg1 = 0; 
    double dRTag = matchTriggerObject(trigObjVec, tagv, _triggerPathLeg1, -1, 
                                      maxPtDiff,  tag_trig_indx, trg_flg1);

    if (tag_trig_indx < 0 || tag_trig_indx >= ntobj) continue;
    _tagDeltaR->Fill(dRTag, _puevWt);
    TriggerObject tagTrigObj = trigObjVec.at(tag_trig_indx); 
    double dRPtTag = tagv.Pt() - tagTrigObj.pt;
    _tagDRVsPt->Fill(imuon.pt, dRTag, _puevWt);
    if (dRTag <= AnaUtil::cutValue(_evselCutMap, "maxDr") && _iFlagL1 == 1) {
      _tagPt->Fill(tagv.Pt(), _puevWt);
      _tagPtv->Fill(tagv.Pt(), _puevWt);
      _tagEta->Fill(tagv.Eta(), _puevWt);
      _tagPhi->Fill(tagv.Phi(), _puevWt);
    
      TLorentzVector tagTrigTL;
      tagTrigTL.SetPtEtaPhiE(tagTrigObj.pt, tagTrigObj.eta, tagTrigObj.phi, tagTrigObj.energy);
      _trigTagPt->Fill(tagTrigObj.pt, _puevWt);
      _trigTagEta->Fill(tagTrigObj.eta, _puevWt);
      _trigTagPhi->Fill(tagTrigObj.phi, _puevWt);
      _tagVsTrigPt->Fill(tagv.Pt(), tagTrigObj.pt, _puevWt);
      _tagDeltaPt->Fill(dRPtTag, _puevWt);
      _tagVsTrigEta->Fill(tagv.Eta(), tagTrigObj.eta, _puevWt);
      _tagVsTrigPhi->Fill(tagv.Phi(), tagTrigObj.phi, _puevWt);

      bool isGoodVtx;
      TVector3 v1 = findLeptonVtx(imuon.vtxIndex, isGoodVtx);
      bool probe_found = false;
      // Find the Probe
      for (uint jobj = 0; jobj < muonVec.size(); ++jobj) {
	if (jobj == iobj) continue;
	const Muon& jmuon = muonVec.at(jobj);     
	TLorentzVector probev;
	probev.SetPtEtaPhiE(jmuon.pt, jmuon.eta, jmuon.phi, jmuon.energy);
        double mass = (tagv + probev).M();

	_invMass->Fill(mass, _puevWt);
        _tagDRVsInvMass->Fill(mass, dRTag, _puevWt);
	_probeDRVsPt->Fill(jmuon.pt, dRTag, _puevWt);
	if (_iFlagL2 == 1) 
          _invMassPass->Fill(mass, _puevWt);
	else 
          _invMassFail->Fill(mass, _puevWt);
        double deltaM = abs(mass - 91.2);
	if (deltaM < AnaUtil::cutValue(_evselCutMap, "massWindowWide")) {
	  probe_found = true;	  
	  int probe_trig_indx = -1;
          uint trg_flg2 = 0;
	  double dRProbe = matchTriggerObject(trigObjVec, probev, _triggerPathLeg2, tag_trig_indx, 
					      maxPtDiff, probe_trig_indx, trg_flg2);
          if (probe_trig_indx < 0 || probe_trig_indx >= ntobj) {
	    _fLog << " No Trigger Object found for Probe !!! total # of TriggerObjects " 
                  << ntobj
                  << ", drProbe: " << dRProbe << std::endl;
	    continue;
	  }
	  TriggerObject probeTrigObj = trigObjVec.at(probe_trig_indx);
	  TLorentzVector probeTrigTL;
	  probeTrigTL.SetPtEtaPhiE(probeTrigObj.pt, probeTrigObj.eta, probeTrigObj.phi, probeTrigObj.energy);
          if (deltaM < AnaUtil::cutValue(_evselCutMap, "massWindowNarrow")) {
	    _probePt->Fill(probev.Pt(), _puevWt);
	    _probePtv->Fill(probev.Pt(), _puevWt);
	    _probeEta->Fill(probev.Eta(), _puevWt);
	    _probePhi->Fill(probev.Phi(), _puevWt);
	    _probeDRVsInvMass->Fill(mass, dRProbe, _puevWt); 
	    _probeDeltaR->Fill(dRProbe, _puevWt);
	    _probeVsTagDR->Fill(dRTag, dRProbe, _puevWt);
	  }
	  double massTrig = (tagTrigTL + probeTrigTL).M();
	  double dRPtProbe = probev.Pt() - probeTrigObj.pt; 
	  _fLog << "=> Summary: " 
  	        << " TagPt: "  << tagv.Pt() 
	        << " TagTrigPt: " << tagTrigObj.pt 
	        << " drTag: "<< dRTag 
	        << " ProbePt: " << probev.Pt() 
	        << " ProbeTrigPt: " << probeTrigObj.pt 
	        << " drProbe: " << dRProbe 
	        << " Flag: " << trg_flg2 
	        << " Mass: " << mass 
	        << " MassTrig: " << massTrig 
	        << endl;
	  
	  if (dRProbe <= AnaUtil::cutValue(_evselCutMap, "maxDr") && _iFlagL2 == 1) {
            if (deltaM < AnaUtil::cutValue(_evselCutMap, "massWindowNarrow")) {
	      _matchedProbePt->Fill(probev.Pt(), _puevWt);
	      _matchedProbePtv->Fill(probev.Pt(), _puevWt);
	      _matchedProbeEta->Fill(probev.Eta(), _puevWt);
	      _matchedProbePhi->Fill(probev.Phi(), _puevWt);
	      _probeVsTrigPt->Fill(probev.Pt(), probeTrigObj.pt, _puevWt);
	      _probeDeltaPt->Fill(dRPtProbe, _puevWt);
	      _probeVsTrigEta->Fill(probev.Eta(), probeTrigObj.eta, _puevWt);
	      _probeVsTrigPhi->Fill(probev.Phi(), probeTrigObj.phi, _puevWt);
	      _invMassTrigPass->Fill(massTrig, _puevWt);
	    }
	    _invMassVsPtPass->Fill(probev.Pt(), mass, _puevWt); 
	    break;
	  } 
          else {
	    if (deltaM < AnaUtil::cutValue(_evselCutMap, "massWindowNarrow")) {
	      _probeDeltaRFail->Fill(dRProbe, _puevWt);
	      _probeVsTagDRFail->Fill(dRTag, dRProbe, _puevWt);
	      _probeTrigDeltaPtFail->Fill(probev.Pt() - probeTrigObj.pt, _puevWt);
	      _probeVsTrigPtFail->Fill(probev.Pt(), probeTrigObj.pt, _puevWt);
	      _probeVsTrigEtaFail->Fill(probev.Eta(), probeTrigObj.eta, _puevWt);
	      _probeVsTrigPhiFail->Fill(probev.Phi(), probeTrigObj.phi, _puevWt);
	      _invMassTrigFail->Fill(massTrig, _puevWt);
	    }
	    _invMassVsPtFail->Fill(probev.Pt(), mass, _puevWt); 
          } 
	}
      }
      if (probe_found) break;
    }
  }
}
void MuonTriggerEfficiency::fillTriggerEfficiencies(std::map<TH1F*, TH1F*> & histo_map) {
  int nBin1;
  nBin1 = _firstMuPt->GetNbinsX();
  for (map<string, TriggerHistos>::iterator it  = _triggerHistoMap.begin(); 
                                            it != _triggerHistoMap.end(); it++) {
    string tname = it->first;
    double muPt1, muPt2, ratio1, ratio2;  
    cout << " Trigger Path " << tname << " events " <<  it->second.eventPassed << endl;
    _fLog << " Trigger Path " << tname << " events " <<  it->second.eventPassed << endl;
    TriggerHistos local_histos = it->second;
    for (int i = 1; i < nBin1+1; i++) {
      muPt1 =  _firstMuPt->GetBinContent(i);
      muPt2 =  _secondMuPt->GetBinContent(i);

      ratio1 = 0.0;
      ratio2 = 0.0;     
      if (muPt1 > 0) ratio1 = local_histos.firstMuPtAcc->GetBinContent(i)/muPt1;    
      if (muPt2 > 0) ratio2 = local_histos.secondMuPtAcc->GetBinContent(i)/muPt2;
      local_histos.firstMuEff->SetBinContent(i,ratio1);  
      local_histos.secondMuEff->SetBinContent(i,ratio2);  
    }
  }
  cout << " Trigger Path 1 : " ;
  _fLog << " Trigger Path 1 : " ;
  for (vector<string>::iterator il1  = _triggerPathLeg1.begin(); 
                                il1 != _triggerPathLeg1.end(); ++il1){
    cout << (*il1) << " " ;
    _fLog << (*il1) << " " ;
  }  
  cout << endl;
  _fLog << endl;
  
  cout << " Trigger Path 2 : " ;
  _fLog << " Trigger Path 2 : " ;
  for (vector<string>::iterator il2  = _triggerPathLeg2.begin(); 
                                il2 != _triggerPathLeg2.end(); ++il2) {
    cout << (*il2) << " " ;
    _fLog << (*il2) << " " ;
  }  
  cout << endl;
  _fLog << endl;
  if (histo_map.size() == 0) return;
  for (std::map<TH1F*, TH1F*>::iterator it  = histo_map.begin(); 
                                        it != histo_map.end(); ++it) {
    TH1F* hist1 = it->first;
    TH1F* hist2 = it->second;
    int nBin = hist1->GetNbinsX();
    string hname1 = hist1->GetName();
    string hname2;
    if (hname1.find("Pt") != string::npos)  hname2 = "EfficiencyPt"; 
    else if (hname1.find("Eta") != string::npos)  hname2 = "EfficiencyEta"; 
    else if (hname1.find("Phi") != string::npos)  hname2 = "EfficiencyPhi"; 
    TH1F* heff = (TH1F*) gROOT->FindObject(hname2.c_str());
    if (heff) heff->Reset();
    else {
      if (hist1->GetXaxis()->IsVariableBinSize()) {
	std::cout << " Variable Bin " << hist1->GetNbinsX() << std::endl;
	heff = new TH1F(hname2.c_str(), hname2.c_str(), hist1->GetNbinsX(), hist1->GetXaxis()->GetXbins()->GetArray());
      }
      else  heff = new TH1F(hname2.c_str(), hname2.c_str(), nBin, hist1->GetXaxis()->GetXmin(), hist1->GetXaxis()->GetXmax());
    }
    double evt1 = hist1->GetEntries();
    double evt2 = hist2->GetEntries();
    double av_eff = 0.0;
    double av_eff_err = 0.0; 
    if (evt1 > 0.0 && hname1.find("Pt") != string::npos) {
      av_eff = evt2/evt1;
      av_eff_err = (1.0/evt1) * sqrt(evt2 * (1 - av_eff));
      cout << " Probe Events " << evt1 << "  Matched Probe Events " << evt2 << endl;
      _fLog << " Probe Events " << evt1 << "  Matched Probe Events " << evt2 << endl;
      cout << " Average Efficiency " << av_eff << " +- " << av_eff_err << endl;
      _fLog << " Average Efficiency " << av_eff << " +- " << av_eff_err << endl;
    }
    for (int j = 1; j < nBin+1; j++) {
      double val1 =  hist1->GetBinContent(j);
      double val2 =  hist2->GetBinContent(j);
      double ratio = 0.0;
      double error = 0.0;
      if (val1 > 0.0 ) {
	ratio = val2/val1; 
        error = (1.0 / val1) *  sqrt (val2 * (1-ratio));
	if (heff) {
	  heff->SetBinContent(j, ratio);
	  heff->SetBinError(j, error);
	}
      }
      if (hname1.find("Pt") != string::npos) {
	std::cout << j << " " << hist1->GetBinCenter(j) << " " << val1 << " " << val2 << " " << ratio << std::endl;
            _fLog << j << " " << hist1->GetBinCenter(j) << " " << val1 << " " << val2 << " " << ratio << std::endl;
      }
    }
  }
}
void MuonTriggerEfficiency::bookTriggerPrescaleHistograms(string tname) {

  TriggerHistos local_histos;

  TDirectory* topDir = gDirectory;
  TDirectory* newDir = topDir->mkdir(tname.c_str()); 
  if (newDir) newDir->cd();
  string hname = "Prescale_";
  hname += tname;
  string htitle = "Prescale Factor vs Run #(";
  htitle += tname +")";
  local_histos.prescaleHist = new TProfile(hname.c_str(),htitle.c_str(),1000,160400.0,180400.0,0.0,2000.0);  
  
  string hname1, hname2, htitle1, htitle2;

  hname1 = "ptMu1Acc";
  hname1 += tname;
  htitle1 = "First Muon Pt for ";
  htitle1 += tname ;
  local_histos.firstMuPtAcc = new TH1F(hname1.c_str(),htitle1.c_str(),100,-0.5, 199.5);
  
  hname2 = "ptMu2Acc";
  hname2 += tname;
  htitle2 = "Second Muon Pt for ";
  htitle2 += tname;
  local_histos.secondMuPtAcc = new TH1F(hname2.c_str(),htitle2.c_str(),100,-0.5, 199.5);
  
  hname1 = "effMu1";
  hname1 += tname;
  htitle1 = "Tigger Efficiency vs First Muon Pt for";
  htitle1 += tname;
  local_histos.firstMuEff = new TH1F(hname1.c_str(),htitle1.c_str(),100,-0.5, 199.5);
  
  hname2 = "effMu2";
  hname2 += tname;
  htitle2 = "Trigger Efficiency vs Second Muon Pt for";
  htitle2 += tname;
  local_histos.secondMuEff = new TH1F(hname2.c_str(),htitle2.c_str(),100,-0.5, 199.5);

  local_histos.eventPassed = 0;
  _triggerHistoMap.insert(std::make_pair(tname, local_histos)); 
  topDir->cd();
}
// -------------------------------------------------------------------------------
// Poor man's way of a datacard. Each line between the 'START' and 'END' tags
// is read in turn, split into words, where the first element is the 'key' and
// the rest the value(s). If more than one values are present they are 
// stored in a vector. No safety mechanism is in place. Any line with an unknown 
// key is skipped. Comments lines should start with either '#' or '//', preferably
// in the first column. Empty lines are skipped. The file containing the datacards 
// is passed as the only argument of the program, there is no default datacard
// -------------------------------------------------------------------------------
bool MuonTriggerEfficiency::readJob(const string& jobFile, int& nFiles)
{
  if (!AnaBase::readJob(jobFile, nFiles)) return false;

  static const int BUF_SIZE = 512;

  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "Input File: " << jobFile << " could not be opened!" << endl;
    return false;
  }

  // note that you must use a pointer (reference!) to the cut map
  // in order to avoid scope related issues
  map<string, map<string, double>* > hmap;
  hmap.insert (pair<string, map<string, double>* >("evselCutList", &_evselCutMap));

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
    string key = tokens.at(0);

    if (key == "triggerPathForLeg1")
      AnaUtil::buildList(tokens, _triggerPathLeg1);
    else if (key == "triggerPathForLeg2")
      AnaUtil::buildList(tokens, _triggerPathLeg2);
    else if (key == "dumpEvent")
      _dumpEvent = atoi(tokens.at(1).c_str()) > 0 ? true : false;
    else if (key == "evselCutList")
      AnaUtil::storeCuts(tokens, hmap);

    tokens.clear();
  }
  // Close the file
  fin.close();

  printJob();
  return true;
}
void MuonTriggerEfficiency::printJob(ostream& os) const
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("evselCutList", _evselCutMap));

  AnaUtil::showCuts(hmap, os);
}
void MuonTriggerEfficiency::clearCollection() {
  vertexVec.clear();
  muonVec.clear();
  electronVec.clear();
  tauVec.clear();
  jetVec.clear();
  trigObjVec.clear();
}
