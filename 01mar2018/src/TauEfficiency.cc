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

#include "TauEfficiency.h"
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
TauEfficiency::TauEfficiency(const string& filename)
  : AnaBase(filename)
{}
// ----------
// Destructor
// ----------
TauEfficiency::~TauEfficiency() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool TauEfficiency::beginJob() 
{ 
  AnaBase::beginJob();
  for (unsigned int i = 0; i < NEL(nSteps); i++) {
    nSteps[i] = 0;
  }
  // Open the output ROOT file
  _histf->cd();
  bookHistograms();

  return true;
}
// ---------------
// Book histograms
// ---------------
void TauEfficiency::bookHistograms() 
{
  static const double pi = TMath::Pi();
  new TH1F("lead_pt","pt of the lead particle in tau jet", 100, 0 , 100);
  new TH1F("iso", "iso cone wrt charged and photon", 100, 0, 10);
  new TH1F("met","event met ", 100, 0, 100);
  new TH1F("decay_mode","decay mode of hadronic tau",10, 0, 10);
  new TH1F("loose_iso","loose isolation cut-tau id discriminator", 10, 0,10);
  new TH1F("medium_iso","medium isolation cut-tau id discriminator", 10, 0,10);
  new TH1F("tight_iso","tight isolation cut-tau id discriminator", 10, 0,10);
  new TH1F("vis_mass","visible mass of the tau hadron and muon", 100, 0, 100);
  new TH1F("deltaR","deltaR angle between the muon and tau hadron",100, 0, 20);
  new TH1F("trans_mass","transverse mass of the muon and missing ET", 100, 0, 100); 
}
// -------------------
// The main event loop
// -------------------
void TauEfficiency::clearLists() {
  vtxList.clear();
  muoList.clear();
  tauList.clear();
  eleList.clear();
  bjetList.clear();
}
void TauEfficiency::eventLoop() 
{
  // Initialize analysis
  if (!beginJob()) return;
  

  int nPrint = max(10000, nEvents/1000);

  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  for (int ev = 0; ev < nEvents; ev++) {
    bool select = false;

    clearEvent();
    int lflag = _chain->LoadTree(ev); 
    int nentries = getEntry(lflag);    // returns bytes read
    
    string currentFile(gSystem->BaseName(_chain->GetCurrentFile()->GetName())); 

    const Event* evt = dynamic_cast<Event*>(eventA->At(0));
    assert(evt);
    int run   = evt->run;
    int event = evt->event;

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

    clearLists();
    findVtxInfo(vtxList, _vtxCutMap);
    findMuonInfo(muoList, _muonCutMap);
    findTauInfo(tauList, _tauCutMap);
    findJetInfo(bjetList, _bjetCutMap);
    findElectronInfo(eleList, _electronCutMap);
    // Event Selection Starts here .....
    // presence of > 1 good vertex
    ++nSteps[0];
    if (vtxList.size() < 1) continue;
    if (vtxList[0].ndf < 4.0 || vtxList[0].rho >= 2) continue;
    ++nSteps[1];

/*    //event selection by HLT path,,,specific to Tau ID analysis
    bool matched = false;
    for (int i = 0; i < n_triggerobj; ++i){
       const TriggerObject* tobj = dynamic_cast<TriggerObject*>(triggerobjA->At(i));
       std::map<std::string, unsigned int> pathlist = tobj->pathList;    
       for (std::map< std::string, unsigned int>::iterator it = pathlist.begin(); it != pathlist.end(); it++){
	  string pname = it->first;
          int flag = it->second;
          std::cout << "path name =" << pname << "flag = " << flag << std::endl;
          string str ("HLT_IsoMu15_L1ETM20");
          if (pname.find(str) == string::npos) continue;
          if (flag == 0) continue;
          matched = true;
          break;
       } 
       if (matched == true) break;
    }
    if (matched == false) continue;
*/  

    // no electron is required
    if (eleList.size()) continue;
    ++nSteps[2];

    // no bjet
    //if (bjetList.size()) continue;
    ++nSteps[3];

    if (muoList.size() < 1) continue;
    ++nSteps[4];
    const MET* taumet= dynamic_cast<MET*>(metA->At(0));
    fillHist1D("met", taumet->met, 1.0);

    if (taumet->met < 20) continue;
    ++nSteps[5];

    computeTauEff(muoList,tauList);
  }  
  // Analysis is over

  endJob();
}
void TauEfficiency::computeTauEff(const vector<Muon>& muoList, const vector<Tau>& tauList){

  const Muon& muon = muoList[0]; 
  if (muon.pt < 15) return;
  if (muoList.size() > 1 && muoList[1].pt >=  15) return;
  TLorentzVector M;
  M.SetPtEtaPhiE (muon.pt, muon.eta, muon.phi, muon.energy);
  for (unsigned int i = 0; i < tauList.size(); i++){
    const Tau& tau = tauList[i];
    fillHist1D("lead_pt",tau.leadParticlePt, 1.0);
    if (tau.leadParticlePt < 5) continue;
    fillHist1D("iso",tau.ptSumPFChargedHadronsIsoCone + tau.ptSumPhotonsIsoCone, 1.0);
    if ((tau.ptSumPFChargedHadronsIsoCone + tau.ptSumPhotonsIsoCone) >= 5) continue;
    fillHist1D("decay_mode", tau.decayModeFinding, 1.0);
    fillHist1D("loose_iso", tau.looseIsolation, 1.0);
    fillHist1D("medium_iso", tau.mediumIsolation, 1.0);
    fillHist1D("tight_iso", tau.tightIsolation, 1.0);
    TLorentzVector T,z;
    T.SetPtEtaPhiE (tau.pt, tau.eta, tau.phi, tau.energy);
    z = T+M;
    double vis_mass = z.M();
    fillHist1D("vis_mass", vis_mass, 1.0);
    double dR= AnaUtil::deltaR(M,T);
    fillHist1D("deltaR",dR, 1.0);
    //TVector3 v1,v2,v3;
    //v1.SetMagEtaPhi (1.0, muon.Eta(), muon.phi);
    //v2.SetMagEtaPhi (1.0, tau.Eta(), tau.phi);
        
  }
  const MET* taumet= dynamic_cast<MET*>(metA->At(0));
  double trans_mass = sqrt(pow(muon.pt+taumet->met,2)-pow(muon.pt*cos(muon.phi)+taumet->met*cos(taumet->metphi),2)-
                           pow(muon.pt*sin(muon.phi)+taumet->met*sin(taumet->metphi),2));
  fillHist1D("trans_mass", trans_mass, 1.0);
} 





// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void TauEfficiency::endJob() 
{
  const string steps[] =
  {
  "after initialization of job ",
  "after the vertex cuts ",
  "no electron requirment ",
  "no bjet requirment ",
  "after atleast 1 global muon requirment",
  "after met cut "
  }; 
  _fLog << "Tau ID efficiency, checking steps"<< endl;
  for (unsigned int i=0; i< NEL(steps); i++) {
    _fLog << setw(64) << steps[i] << setw(10) << nSteps[i]<< endl;
  }
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
bool TauEfficiency::selectEvent() 
{
  return true;
}
bool TauEfficiency::readJob(const string& jobFile, int& nFiles)
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
void TauEfficiency::printJob(ostream& os)
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
TH1* TauEfficiency::getHist1D(const string& hname) {
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
bool TauEfficiency::fillHist1D(const string& hname, T value, double w) 
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
bool TauEfficiency::fillHist2D(const string& hname, T1 xvalue, T2 yvalue, double w) 
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
TProfile* TauEfficiency::getProfile(const string& hname) 
{
  TProfile *h = dynamic_cast<TProfile*>(gDirectory->GetList()->FindObject(hname.c_str()));
  if (!h) {
    cerr << "**** Profile Histogram <<" << hname << ">> not found" << endl;
    return 0;
  }
  return h;
}
bool TauEfficiency::fillProfile(const string& hname, float xvalue, float yvalue, double w) 
{
  TProfile *h = getProfile(hname);
  if (!h) return false;

  h->Fill(xvalue, yvalue, w);
  return true;
}
