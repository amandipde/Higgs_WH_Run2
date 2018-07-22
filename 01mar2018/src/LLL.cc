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
#include "LLL.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"
#include "MVASkim.h"
#include <cmath>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

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
LLL::LLL()
  : AnaBase(),
    _createMVATree(false),
    _readMVA(false),
    _readMVAFK(false),
    _mvaInputFile(""),
    _MVAdisFile(""),
    _MVAFKdisFile(""),
    _skimObj(0)
{}
// ----------
// Destructor
// ----------
LLL::~LLL() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool LLL::beginJob() 
{ 
  AnaBase::beginJob();

  histf()->cd();
  bookHistograms();

  return true;
}
// ---------------
// Book histograms
// ---------------
void LLL::bookHistograms() 
{
  new TH1D("counter_gen", "Selected event counter", 31, -0.5, 30.5);
  new TH1F("npu", "nPU distribution", 100, 0, 100);
  new TH1F("muPt", "pt distribution of the muon after all cut", 140, 0, 140);
  new TH1F("tau1Pt", "pt distribution of the tau1 after all cut", 140, 0, 140);
  new TH1F("mu2Pt", "pt distribution of the mu2 after all cut", 140, 0, 140);
  new TH1F("noMuon", "number of events with null muon collection", 3, -0.5, 2.5);
  new TH1F("muPt_Vtx", "pt distribution of the muons after vertex selection", 200, 0, 200);

  new TH1F("mueta", "eta distribution of the muon after all cut", 100, -2.5, 2.5);
  new TH1F("tau1eta", "eta distribution of the tau1 after all cut", 100, -2.5, 2.5);
  new TH1F("mu2eta", "eta distribution of the mu2 after all cut", 100, -2.5, 2.5);

  new TH1F("mvamet", "met distribution of the event after all cut", 350, 0, 350);
  new TH1F("pfmet", "met distribution of the event after all cut", 70, 0, 350);
  new TH1F("MuTauMass", "vissible mass distribution of the muon and leading tau after all cut", 200, 0, 200);
  
  new TH1F("filter", "WH normal vs WH contamination", 5, 0.5, 5.5);
  new TH1F("LT", "lepton pt sum", 200, 0, 200);
  new TH1F("WMuonPt_Gen", "Pt of the Muon coming from W Decay in Gen Level", 140, 0, 140);
  new TH1F("HMuonPt_Gen", "Pt of the Muon coming from H Decay in Gen Level", 140, 0, 140);
  new TH1F("Diff_MuPt_Gen", "Difference in Muon Pt coming from W and H in Gen Level", 100, -50, 50);
  new TH2D("WMuPt_vs_Diff_MuPt_Gen", "Difference in Muon Pt coming from W and H as a function of WMuPt in Gen Level", 140, 0, 140, 100, -50, 50);

}
// -------------------
// The main event loop
// -------------------
void LLL::clearLists() {
  vtxList.clear();
  muoList.clear();
  eleList.clear();
  tauList.clear();
  bjetList.clear();

  genHDecMuonList.clear();
  genWDecMuonList.clear();
  genHDecTauList.clear();
  genWDecTauList.clear();
  genHDecEleList.clear();
  genWDecEleList.clear();

  genMuonList.clear();
  genTauList.clear();
}
void LLL::eventLoop() 
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
  for (int ev = 0; ev < nEvents(); ++ev) {
    clearEvent();
    int lflag = chain()->LoadTree(ev); 
    int nbytes = getEntry(lflag);    // returns total bytes read

    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName())); 

    const Event& evt = eventColl()->at(0);

    histf()->cd();

    // PileUP weight
    puevWt_ = 1; // for data

    /*
    if (isMC()) {
      int npu = 0;
      puevWt_ = wtPileUp(npu);
    }
    */
    AnaUtil::fillHist1D("counter_gen", 0, puevWt_);

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

    if (eventIdMap().size()) {
      std::ostringstream mkey;
      mkey << event << "-" << lumis << "-" << run;
      if (eventIdMap().find(mkey.str()) == eventIdMap().end()) continue;
      evLog() << "Event " << event
            << " Lumis " << lumis
            << " Run " << run << endl;
      //dumpGenInfo(evLog()); 
      //dumpEvent("001111", evLog(), false);
      //cout << "Event = " << event << " ,nMuon = " << nmuon() << std::endl;
    }

    if (logOption() > 0) {
      fLog() << "Event " << event
            << " Lumis " << lumis
            << " Run " << run << endl;
      fLog() << "n_tau: "<< ntau()
	     << ", n_muon: "<< nmuon()
	     << ", n_jet: " << njet()
	     << ", n_vertex: " << nvertex()
	     << ", n_met: " << nmet()
             << ", n_electron: " << nelectron()
             << endl;
    }
    clearLists();
 
    //findGenInfo(genMuonList, genTauList);
    //if (genMuonList.size() != 1 || genTauList.size() != 2) continue;
    if (logOption() >> 1 & 0x1) dumpGenInfo(fLog()); 
  
    //Define a WmHmh Decay
    //if (!WmHmhFilter()) continue;                                        //Put this off if not required::ATTENTION
    findSigDecayGenInfo(genHDecMuonList, genWDecMuonList,
                        genHDecTauList,  genWDecTauList,
                        genHDecEleList,  genWDecEleList);
    std::cout << "a" << std::endl;
    if (!(genWDecMuonList.size() == 1 && genHDecMuonList.size() == 2 && genHDecTauList.size() == 0 &&
        genWDecTauList.size() == 0 && genWDecEleList.size() == 0 && genHDecEleList.size() == 0)) continue;  
    std::cout << "b" << std::endl;

    AnaUtil::fillHist1D("WMuonPt_Gen", genWDecMuonList[0].pt, puevWt_);
    AnaUtil::fillHist1D("HMuonPt_Gen", genHDecMuonList[0].pt, puevWt_);
    AnaUtil::fillHist1D("Diff_MuPt_Gen", genWDecMuonList[0].pt-genHDecMuonList[0].pt, puevWt_);
    AnaUtil::fillHist2D("WMuPt_vs_Diff_MuPt_Gen", genWDecMuonList[0].pt, genWDecMuonList[0].pt-genHDecMuonList[0].pt, puevWt_);

    AnaUtil::fillHist1D("counter_gen", 1, puevWt_);
    dumpEvent("000000", evLog(), false);
    AnaUtil::fillHist1D("npu", evt.nPU.at(0), puevWt_);

    // Trigger selection
    //if (useTrigger() && !isTriggered()) continue;
    AnaUtil::fillHist1D("counter_gen", 2, puevWt_);

    op.verbose = (logOption() >> 2 & 0x1); 
    findVtxInfo(vtxList, op, fLog());
    int nvtx = vtxList.size();
    double vz = ( nvtx > 0 ) ? vtxList.at(0).z : -999;

    if (logOption() > 0)
    fLog() << "Event " << event
          << " Lumis " << evt.lumis
          << " Run " << run 
          << " n_vertex_good " << nvtx
          << endl;

    op.verbose = (logOption() >> 5 & 0x1);
    findJetInfo(bjetList, op, fLog());

    // Event Selection Starts here .....
    // presence of > 1 good vertex
    //if (vtxList.size() < 1) continue;
    AnaUtil::fillHist1D("counter_gen", 3, puevWt_);

    std::cout << "c" << std::endl;
    selectEvent();
    std::cout << "d" << std::endl;

  }  
  // Analysis is over
  endJob();
}
void LLL::selectEvent() {
  study_eff();
  studyTauID();
}
void LLL::studyTauID(){
  findZtype(genMuonList_z, genEleList_z, genTauList_z);
  //ZToTauTau hadronic decay
  if (genMuonList_z.size() == 0 && genEleList_z.size() == 0 && genTauList_z.size() == 2){
    for (auto it = tauColl()->begin(); it != tauColl()->end(); ++it) {
      const Tau& tau = (*it);

      //if (fabs(tau.eta) < 2.3 && tau.pt > 20.0 && tau.decayModeFinding > 0.5)
    }
  }
  //ZToEE decay
  if (genMuonList_z.size() == 0 && genEleList_z.size() == 2 && genTauList_z.size() == 0){
  }
  //ZToMuMu decay
  if (genMuonList_z.size() == 2 && genEleList_z.size() == 0 && genTauList_z.size() == 0){
  }

}
void LLL::study_eff()
{
  double vz = vtxList.at(0).z;

  int mcounters[] = {0,0,0,0,0,0,0,0,0,0};
  vector<Muon> fMuoList;
  //std::cout << "nMuon " << muonColl()->size() << std::endl;
  if (muonColl()->size() == 0) AnaUtil::fillHist1D("noMuon", 1.0, puevWt_);

  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    AnaUtil::fillHist1D("muPt_Vtx", muon.pt, puevWt_);

    if (muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt") 
      || fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta")) continue;
    ++mcounters[0];

    if (!muon.isTrackerMuon || !muon.isGlobalMuonPromptTight) continue;
    ++mcounters[1];
    //if (mcounters[1] == 1) dumpEvent("000000", evLog(), false);

    if (muon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) continue;
    ++mcounters[2];

    if (muon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches") ||
        muon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) continue;
    ++mcounters[3];

    if (muon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) continue;
    ++mcounters[4];

    if (muon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) continue;
    ++mcounters[5];

    if (muon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) continue;
    ++mcounters[6];

    if (fabs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0")) continue;
    ++mcounters[7];

    if (abs(muon.vtxDistZ) >= AnaUtil::cutValue(muonCutMap(), "vtxDistZ")) continue;
    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;
    ++mcounters[8];

    //if (abs(muon.pfRelIso) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso")) continue;
    //if (abs(muon.relIso) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso")) continue;
    //if (fabs(muon.trkIso/muon.pt) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso")) continue;
    if (fabs(muon.pfChargedIsoR04/muon.pt) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso")) continue;
    ++mcounters[9];

    fMuoList.push_back(muon);
  }
  int ishift = 4;
  for (size_t i = 0; i < NEL(mcounters); ++i) {
    if (mcounters[i]) AnaUtil::fillHist1D("counter_gen", ishift + i, puevWt_);
  }
  
  // atleast 3 good muons
  if (fMuoList.size() < 3) return;
  AnaUtil::fillHist1D("counter_gen", 14, puevWt_);
  const Muon& muo = fMuoList.at(0);
  if (muo.pt <= 18.0) return;
  AnaUtil::fillHist1D("counter_gen", 15, puevWt_);
  const Muon& muo2 = fMuoList.at(1);

  //SS Muons
  //if ((muo.charge + muo2.charge) == 0) return;
  AnaUtil::fillHist1D("counter_gen", 16, puevWt_);

  TLorentzVector M, M2;
  M.SetPtEtaPhiE(muo.pt, muo.eta, muo.phi, muo.energy);
  M2.SetPtEtaPhiE(muo2.pt, muo2.eta, muo2.phi, muo2.energy);

  // Here comes MET
  //const MET& mvamet = metColl()->at(0);
  const MET& pfmet = metColl()->at(0);

  //---------------------------------------------
  //
  //              Muon Veto
  //
  //---------------------------------------------
  //if (vetoMuon(taua.zvertex, 10, 0.2) != 3) return;
  AnaUtil::fillHist1D("counter_gen", 27, puevWt_);
  
  //---------------------------------------------
  //
  //              Electron Veto
  //
  //---------------------------------------------
  //if (vetoElectron(taua.zvertex, 10, 0.14) != 0) return;
  AnaUtil::fillHist1D("counter_gen", 28, puevWt_);

  //---------------------------------------------
  //
  //              b-tagged Jets Veto
  //
  //---------------------------------------------
  //if (bjetList.size() > 0) return; 
  AnaUtil::fillHist1D("counter_gen", 29, puevWt_);


  //if (sumPt <= 80) return;                        
  AnaUtil::fillHist1D("counter_gen", 30, puevWt_);

  //After All cut
  AnaUtil::fillHist1D("muPt", muo.pt, puevWt_);
  AnaUtil::fillHist1D("mueta", muo.eta, puevWt_);
  AnaUtil::fillHist1D("mu2Pt", muo2.pt, puevWt_);
  AnaUtil::fillHist1D("mu2eta", muo2.eta, puevWt_);
  //AnaUtil::fillHist1D("mvamet", mvamet.met, puevWt_);
  AnaUtil::fillHist1D("pfmet", pfmet.met, puevWt_);

  //dumpEvent("000000", evLog(), false);

}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void LLL::endJob() 
{  
  histf()->cd();

  TH1 *h = AnaUtil::getHist1D("counter_gen");
  if (h) {
    fLog() << setprecision(3);
    fLog() << "Events Processed: " << int(h->GetBinContent(1)) << endl;
    for (int i = 2; i <= h->GetNbinsX(); ++i)
      fLog() << "Cut: " 
            << setw(3) << i 
            << setw(10) << int(h->GetBinContent(i))
            << setw(10) << ((h->GetBinContent(i-1)>0) ? h->GetBinContent(i)/h->GetBinContent(i-1) : 0.0)
            << setw(10) << ((h->GetBinContent(2)>0) ? h->GetBinContent(i)/h->GetBinContent(2) : 0.0)
            << endl;
  }
  histf()->Write();
  histf()->Close();
  delete histf();

  fLog() << resetiosflags(ios::fixed);

  if (_skimObj) _skimObj->close();

  closeFiles();
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
bool LLL::readJob(const string& jobFile, int& nFiles)
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
  hmap.insert(pair<string, map<string, double>* >("tau1CutList", &_tau1CutMap));
  hmap.insert(pair<string, map<string, double>* >("tau2CutList", &_tau2CutMap));

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
    string key = tokens.at(0);
    string value = tokens.at(1);

    if (key == "tau1CutList" || key == "tau2CutList")
      AnaUtil::storeCuts(tokens, hmap);
    else if (key == "createMVATree")
      _createMVATree = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "mvaInputFile")
      _mvaInputFile = value;
    else if (key == "readMVAFK")
      _readMVAFK = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "readMVA")
      _readMVA = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "MVAdisFile")
      _MVAdisFile = value;
    else if (key == "MVAFKdisFile")
      _MVAFKdisFile = value;


    tokens.clear();
  }
  // Close the file
  fin.close();

  printJob();

  return true;
}
void LLL::printJob(ostream& os) const
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("tau1CutList", _tau1CutMap));
  hmap.insert(pair<string, map<string, double> >("tau2CutList", _tau2CutMap));
  AnaUtil::showCuts(hmap, os);
}
