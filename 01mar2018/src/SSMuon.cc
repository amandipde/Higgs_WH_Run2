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

#include "SSMuon.h"
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
SSMuon::SSMuon()
  : AnaBase()
{}
// ----------
// Destructor
// ----------
SSMuon::~SSMuon() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool SSMuon::beginJob() 
{ 
  if (!AnaBase::beginJob()) return false;

  histf()->cd();
  bookHistograms();

  // Optionally write selected events to another tree
  return true;
}
// ---------------
// Book histograms
// ---------------
void SSMuon::bookHistograms() 
{
  static const double pi = TMath::Pi();

  if (isMC() && usePUWt()) {
    new TH1D("npu", "nPileUp", 50, 0, 50);
    new TH1D("puwt", "PileUp weight factor", 100, 0, 2.0);
    new TH1D("trueNInt", "True N-Interaction", 50, 0, 50);
  } 
  new TH1D("evcounter", "Selected event counter", 8, -0.5, 7.5);
  new TH1D("nvertex", "No of good vertex", 30, -0.5, 29.5);
  new TH1F("muPt1",  "Highest Pt Muon Pt  after vtx and >= 2 Muon requirment", 100, -0.5, 99.5);
  new TH1F("muEta1", "Highest Pt Muon Eta  after vtx and >= 2 Muon requirment", 100, -3, 3);
  new TH1F("muPt2",  "Highest Pt Muon Pt(vtx, >= 2 Muons, muPt1)", 100, -0.5, 99.5);
  new TH1F("muEta2", "Highest Pt Muon Eta(vtx,>= 2 Muons, muPt1)", 100, -3, 3);
  new TH1F("met", "Missing eT distribution", 100, 0, 100);
}
// -------------------
// The main event loop
// -------------------
void SSMuon::clearLists() {
  vtxList.clear();
  muoList.clear();
  eleList.clear();
  tauList.clear();
  bjetList.clear();
}
void SSMuon::eventLoop() 
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

    const Event& evt = eventColl()->at(0);

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
	     << ", ntau: "<< ntau()
	     << ", nmuon: "<< nmuon()
	     << ", njet: " << njet()
	     << ", nvertex: " << nvertex()
	     << ", nmet: " << nmet()
	     << ", nelectron: " << nelectron()
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
void SSMuon::selectEvent()
{ 
  int nsmuon = muoList.size();
  // require atleast 2 muon 
  if (nsmuon < 2) return;
  AnaUtil::fillHist1D("evcounter", 3, puevWt_);

  // muon Pt
  const Muon& muo1 = muoList.at(0);
  double pt1 = muo1.pt;
  AnaUtil::fillHist1D("muPt1", pt1, puevWt_);
  AnaUtil::fillHist1D("muEta1", muo1.eta, puevWt_);
  if (pt1 <= AnaUtil::cutValue(_evselCutMap, "ptMuon1")) return;
  AnaUtil::fillHist1D("evcounter", 4, puevWt_);

  const Muon& muo2 = muoList.at(1);
  double pt2 = muo2.pt;
  AnaUtil::fillHist1D("muPt2", pt2, puevWt_);
  AnaUtil::fillHist1D("muEta2", muo2.eta, puevWt_);
  if (pt2 <= AnaUtil::cutValue(_evselCutMap, "ptMuon2")) return;
  AnaUtil::fillHist1D("evcounter", 5, puevWt_);

  // missing_et cut
  const MET& mt = metColl()->at(0);

  int metv = mt.met;
  AnaUtil::fillHist1D("met", metv, puevWt_);
  if (metv > AnaUtil::cutValue(_evselCutMap, "maxMET")) return;
  AnaUtil::fillHist1D("evcounter", 6, puevWt_);

  // same signed muon
  if ((muo1.charge+muo2.charge) == 0) return;
  AnaUtil::fillHist1D("evcounter", 7, puevWt_);

  // Dump seleted event information
  // Event and object information
  dumpEvent("11111", evLog(), false);
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void SSMuon::endJob() 
{  
  string slist[] = 
    {
      "Total events processed",
      "After Trigger selection",
      ">= 1 good Event Vertex",
      "> 1 selected muons (Id, Iso)",
      "Highest pT Muon pT > 25 GeV",
      "2nd highest pT Muon pT > 20 GeV",
      "missing eT < 100 GeV",
      "Same-signed muons"
    };
  int prec = (isMC()) ? 1 : 0;
  fLog() << setprecision(prec);
  fLog() << "============================================================" << endl
        << "Statistics for 2 same-signed muon events                    " << endl
        << "============================================================" << endl;

  TH1 *h = AnaUtil::getHist1D("evcounter");
  if (h) {
    for (int i = 1; i <= h->GetNbinsX(); ++i)
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
bool SSMuon::readJob(const string& jobFile, int& nFiles)
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
    if (key == "evselCutList")
      AnaUtil::storeCuts(tokens, hmap);

    tokens.clear();
  }
  // Close the file
  fin.close();

  printJob();

  return true;
}
void SSMuon::printJob(ostream& os) const
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("evselCutList", _evselCutMap));
  AnaUtil::showCuts(hmap, os);
}
