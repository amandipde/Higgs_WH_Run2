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
#include <cmath>
#include <sstream>
#include <tuple>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TProfile.h"

#include "AnaBase.h"
#include "AnaUtil.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::abs;
using std::max;
using std::sqrt;
using std::sort;
using std::setprecision;
using std::setw;
using std::ostream;
using std::ifstream;
using namespace vhtm;

// -----------
// Constructor
// -----------
AnaBase::AnaBase()
  : chain_(new TChain("treeCreator/vhtree")),
    histf_(0),
    eventList_(new vector<vhtm::Event>()),
    vertexList_(new vector<vhtm::Vertex>()),    
    genEventList_(new vector<vhtm::GenEvent>()),   
    tauList_(new vector<vhtm::Tau>()),    
    electronList_(new vector<vhtm::Electron>()),
    muonList_(new vector<vhtm::Muon>()),    
    jetList_(new vector<vhtm::Jet>()),    
    metList_(new vector<vhtm::MET>()),    
    mvametList_(new vector<vhtm::MET>()),    
    genParticleList_(new vector<vhtm::GenParticle>()),  
    genJetList_(new vector<vhtm::GenJet>()),   
    genMetList_(new vector<vhtm::GenMET>()),
    triggerObjList_(new vector<vhtm::TriggerObject>()),
    //trackList_(new vector<vhtm::Track>()),
    l1physbits_(new vector<int>()),
    l1techbits_(new vector<int>()),
    hltpaths_(new vector<string>()),
    hltresults_(new vector<int>()),
    hltprescales_(new vector<int>()),
    isMC_(false),
    isSignal_(false),
    readTrk_(false),
    readTrigObject_(true),
    logOption_(0),
    useTrigger_(true),
    usePUWt_(false),
    usenPV_(false),
    histFile_("default.root"),
    puHistFile_("./Reweight_PU_DY_2017.root"),
    useTrueNInt_(true),
    logFile_("default.out"),
    evFile_("default_event.out"),
    maxEvt_(0)
{
  cout << setiosflags(ios::fixed); 
  cout << "=== Start of Analysis === " << endl;
  fileList_.clear();
  brList_.clear();
  puWtList_.clear();
  trigPathList_.clear();
}
// ----------
// Destructor
// ----------
AnaBase::~AnaBase() 
{
  clearEvent();

  delete eventList_;
  delete vertexList_;
  delete genEventList_;
  delete tauList_;
  delete electronList_;
  delete muonList_;
  delete jetList_;
  delete metList_;
  delete mvametList_;
  delete genParticleList_;
  delete genJetList_;
  delete genMetList_;
  delete triggerObjList_;
  //delete trackList_;
}
// ------------------------
// Clear the clones arrays
// ------------------------
void AnaBase::clearEvent() 
{
  if (eventList_) eventList_->clear();
  if (genEventList_) genEventList_->clear();
  if (vertexList_) vertexList_->clear();
  if (tauList_) tauList_->clear();
  if (electronList_) electronList_->clear();
  if (muonList_) muonList_->clear();
  if (jetList_) jetList_->clear();
  if (metList_) metList_->clear();
  if (mvametList_) metList_->clear();
  if (genParticleList_) genParticleList_->clear();
  if (genJetList_) genJetList_->clear();
  if (genMetList_) genMetList_->clear();
  if (triggerObjList_) triggerObjList_->clear();
  //if (trackList_) trackList_->clear();

  if (l1physbits_)   l1physbits_->clear();
  if (l1techbits_)   l1techbits_->clear();
  if (hltpaths_)     hltpaths_->clear();
  if (hltresults_)   hltresults_->clear();
  if (hltprescales_) hltprescales_->clear();
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool AnaBase::beginJob() 
{ 
  if (isMC_ && usePUWt_ && !readPileUpHist()) return false;

  // Open the output ROOT file
  histf_ = TFile::Open(histFile_.c_str(), "RECREATE");

  setAddresses();
  nEvents_ = static_cast<int>(chain_->GetEntries()); 
  if (nEvents_ <= 0) {
    cerr << "******* nEvents = " << nEvents_ << ", returning!" << endl;
    return false;
  }
  if (maxEvt_ > 0) nEvents_ = std::min(nEvents_, maxEvt_);
  cout << " ===== # of events to analyse, nEvents = " << nEvents_ << endl;

  openFiles();

  return true;
}

// ------------------------------------
// Get Run number for the present event
// ------------------------------------
int AnaBase::getRunNumber() const 
{
  const Event& event = eventList_->at(0);
  return event.run;
}    
// ---------------------------------
// Add input Root files to the chain
// ---------------------------------
int AnaBase::setInputFile(const string& fname) 
{
  size_t found = fname.find("root:");
  if (found == string::npos && gSystem->AccessPathName(fname.c_str())) {
    cerr << ">>> Warning: File <<" << fname << ">> was not found!!" << endl;
    return static_cast<int>(chain_->GetEntries()); 
  }
  chain_->AddFile(fname.c_str(), -1);
  return static_cast<int>(chain_->GetEntries()); 
}
// ---------------------------------------
// Get total number of events in the chain
// --------------------------------------
int AnaBase::getEntries() const 
{
  return static_cast<int>(chain_->GetEntries());
}

// ------------------------------------------------------
// Open the output file with a global filehandle, C++ way
// ------------------------------------------------------
bool AnaBase::openFiles() 
{
  fLog_.open(logFile_.c_str(), ios::out);
  if (!fLog_) {
    cerr << "File: " << logFile_ << " could not be opened!" << endl;
    return false;
  }
  fLog_ << setiosflags(ios::fixed);

  evLog_.open(evFile_.c_str(), ios::out);
  if (!evLog_) {
    cerr << "File: " << evFile_ << " could not be opened!" << endl;
    return false;
  }
  evLog_ << setiosflags(ios::fixed);
  return true;
}
// ------------------------
// Close the output file
// ------------------------
void AnaBase::closeFiles() 
{
  if (fLog_) {
    fLog_ << resetiosflags(ios::fixed); 
    fLog_.close();
  }
  if (evLog_) {
    evLog_ << resetiosflags(ios::fixed); 
    evLog_.close();
  }
}
void AnaBase::setAddresses() 
{
  if (branchFound("Event")) chain_->SetBranchAddress("Event", &eventList_);
  if (branchFound("Vertex")) chain_->SetBranchAddress("Vertex", &vertexList_);
  if (branchFound("Tau"))  chain_->SetBranchAddress("Tau", &tauList_);
  if (branchFound("Electron")) chain_->SetBranchAddress("Electron", &electronList_);
  if (branchFound("Muon")) chain_->SetBranchAddress("Muon", &muonList_);
  if (branchFound("Jet")) chain_->SetBranchAddress("Jet", &jetList_);
  if (branchFound("MET")) chain_->SetBranchAddress("MET", &metList_);
  if (branchFound("mvaMET")) chain_->SetBranchAddress("mvaMET", &mvametList_);
  if (readTrigObject_ && branchFound("TriggerObject")) 
    chain_->SetBranchAddress("TriggerObject", &triggerObjList_);
  //if (readTrk_ && branchFound("Track")) 
  //  chain_->SetBranchAddress("Track", &trackList_);

  if (isMC_) {
    if (branchFound("GenEvent")) chain_->SetBranchAddress("GenEvent", &genEventList_);
    if (branchFound("GenParticle")) chain_->SetBranchAddress("GenParticle", &genParticleList_);
    if (branchFound("GenJet"))    chain_->SetBranchAddress("GenJet", &genJetList_);
  }
  // Now the trigger variables
  if (branchFound("l1physbits")) chain_->SetBranchAddress("l1physbits", &l1physbits_);
  if (branchFound("l1techbits")) chain_->SetBranchAddress("l1techbits", &l1techbits_);
  if (branchFound("hltresults")) chain_->SetBranchAddress("hltresults", &hltresults_);
  if (branchFound("hltprescales")) chain_->SetBranchAddress("hltprescales", &hltprescales_);
  if (branchFound("hltpaths")) chain_->SetBranchAddress("hltpaths", &hltpaths_);
}
bool AnaBase::branchFound(const string& b)
{
  TBranch* branch = chain_->GetBranch(b.c_str());  // Get branch pointer                                                                                  
  if (!branch) {
    cout << ">>> SetBranchAddress: <" << b << "> not found!" << endl;
    return false;
  }
  cout << ">>> SetBranchAddress: <" << b << "> found!" << endl;
  brList_.push_back(b);
  return true;
}
int AnaBase::getEntry(int lflag) const
{
  int nbytes = 0;
  for (auto it  = brList_.begin(); it != brList_.end(); ++it) {
    TBranch* branch = chain_->GetBranch((*it).c_str());
    if (!branch) {
      cout << ">>> Branch: " << (*it) << " not found!" << endl;
      continue;
    }
    nbytes += branch->GetEntry(lflag);
  }
  return nbytes;
}
// not used yet
void AnaBase::enableBranches() 
{
  chain_->SetBranchStatus("*", kFALSE); // Disable all branches
  for (auto it  = brList_.begin(); it != brList_.end(); ++it) {
    chain_->SetBranchStatus((*it).c_str(), kTRUE);
  }
}
bool AnaBase::readJob(const string& jobFile, int& nFiles)
{
  static const int BUF_SIZE = 256;

  // Open the file containing the datacards
  std::ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "Input File: " << jobFile << " could not be opened!" << endl;
    return false;
  }

  // note that you must use a pointer (reference!) to the cut map
  // in order to avoid scope related issues
  map<string, map<string, double>* > hmap;
  hmap.insert(pair<string, map<string, double>* >("vtxCutList", &vtxCutMap_));
  hmap.insert(pair<string, map<string, double>* >("electronCutList", &electronCutMap_));
  hmap.insert(pair<string, map<string, double>* >("muonCutList", &muonCutMap_));
  hmap.insert(pair<string, map<string, double>* >("tauCutList", &tauCutMap_));
  hmap.insert(pair<string, map<string, double>* >("bjetCutList", &bjetCutMap_));

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
    int vsize = tokens.size();
    assert(vsize > 1);

    string key   = tokens.at(0);
    string value = tokens.at(1);
    if (key == "dataType") {
      string vtmp(value);
      std::transform(vtmp.begin(), vtmp.end(), vtmp.begin(), ::toupper);
      vector<string> dt;
      AnaUtil::tokenize(vtmp, dt, "#");
      if (dt.size()) {
        isMC_ = (dt.at(0) == "MC") ? true : false;
        if (isMC_ && dt.size() > 1) {
          isSignal_ = (dt.at(1) == "SIGNAL") ? true : false;
        }
      }
    }
    else if (key == "readTrk") 
      readTrk_ = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "readTrigObject") 
      readTrigObject_ = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "useTrigger") 
      useTrigger_ = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "usePUWt") 
      usePUWt_ = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "usenPV") 
      usenPV_ = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "logFile")
      logFile_ = value;
    else if (key == "eventFile")
      evFile_  = value;
    else if (key == "logOption") 
      logOption_ = strtol(value.c_str(), NULL, 2);
    else if (key == "maxEvent") 
      maxEvt_ = atoi(value.c_str());
    else if (key == "histFile") 
      histFile_ = value;
    else if (key == "puHistFile") 
      puHistFile_ = value;
    else if (key == "useTrueNInt") 
      useTrueNInt_ = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "trigPathList") 
      AnaUtil::buildList(tokens, trigPathList_);
    else if (key == "inputFile") 
      AnaUtil::buildList(tokens, fileList_);
    else if (key == "eventId" && tokens.size() == 4) 
      AnaUtil::buildMap(tokens, eventIdMap_);
    else
      AnaUtil::storeCuts(tokens, hmap);
    tokens.clear();
  }
  if (!isMC_) usePUWt_ = false;

  // Close the file
  fin.close();

  // Build the chain of root files
  for (auto it  = fileList_.begin(); it != fileList_.end(); ++it) {
    string fname = *it;
    cout << ">>> INFO. Adding input file " << fname << " to TChain " << endl;
    ++nFiles;
    int nevt = setInputFile(fname.c_str());
    if (maxEvt_ > 0 && nevt >= maxEvt_) break;
  }

  if (!nFiles) {
    cerr << ">>> WARN. Input Root file list is empty! exiting ..." << endl;
    return false;
  }
  useTCHE_ = (AnaUtil::cutValue(bjetCutMap_, "useTCHE") > 0) ? true : false;
  return true;
}
void AnaBase::printJob(std::ostream& os) const
{
  os << "    datatype = " << ((isMC_) ? "mc" : "data") << endl
     << "     logFile = " << logFile_ << endl 
     << "   eventFile = " << evFile_ << endl
     << "    histFile = " << histFile_ << endl
     << "     usePUWt = " << std::boolalpha << usePUWt_ << endl
     << "     usenPV = " << std::boolalpha << usenPV_ << endl
     << "  puHistFile = " << puHistFile_ << endl
     << " useTrueNInt = " << std::boolalpha << useTrueNInt_ << endl
     << "  useTrigger = " << std::boolalpha << useTrigger_ << endl
     << "   logOption = " << logOption_ << endl
     << "    maxEvent = " << maxEvt_ << endl;

  // Trigger Path List
  if (useTrigger_) 
    AnaUtil::showList(trigPathList_, ">>> INFO. Trigger Paths used:", os);

  // InputFiles
  if (chain_) {
    TObjArray *fileElements = chain_->GetListOfFiles();
    os << ">>> INFO. nFiles: " << fileElements->GetEntries() 
       << ", Files to analyse:" 
       << endl;
    TIter next(fileElements);
    TChainElement *chEl = 0;
    while (( chEl = dynamic_cast<TChainElement*>(next()) ))
      os << chEl->GetTitle() 
         << endl;
  }
  else
    AnaUtil::showList(fileList_, ">>> INFO. inputFiles", os);

  // EventID 
  AnaUtil::showMap<string, int>(eventIdMap_, "Event List:", os);

  // Cuts
  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("vtxCutList", vtxCutMap_));
  hmap.insert(pair<string, map<string, double> >("electronCutList", electronCutMap_));
  hmap.insert(pair<string, map<string, double> >("muonCutList", muonCutMap_));
  hmap.insert(pair<string, map<string, double> >("tauCutList", tauCutMap_));
  hmap.insert(pair<string, map<string, double> >("bjetCutList", bjetCutMap_));
  AnaUtil::showCuts(hmap, os);
}
// Collect Object information
// Let's look at the event vertex  
void AnaBase::findVtxInfo(vector<Vertex>& list, Options& op, std::ostream& os) {
  if (nvertex() < 1) return;

  if (op.verbose)
    os << "=>> Vertices: " << nvertex() << endl
       << "indx     ndf     dxy       z   sumPt    chi2   ntrks  ntrkw05     sbit" 
       << endl; 
  // already in correct order, no need to sort
  int indx = 0;
  for (auto it = vertexList_->begin(); it != vertexList_->end(); ++it,indx++) {
    const Vertex& vtx = (*it);

    double dxy = sqrt(pow(vtx.x, 2) + pow(vtx.y, 2));
    int sbit = 0;
    if (vtx.ndf <= AnaUtil::cutValue(vtxCutMap_, "ndf"))    sbit |= (1 << 0);
    if (dxy >= AnaUtil::cutValue(vtxCutMap_, "dxy"))         sbit |= (1 << 1);
    if (fabs(vtx.z) >= AnaUtil::cutValue(vtxCutMap_, "z"))   sbit |= (1 << 2);

    if (op.verbose) {
      bool pp = (op.printselected && sbit) ? false : true;
      if (pp) {
        os << setprecision(2)
           << setw(4) << indx
           << setw(8) << vtx.ndf
           << setw(8) << dxy
           << setw(8) << vtx.z 
           << setw(8) << vtx.sumPt
           << setw(8) << vtx.chi2
	  //<< setw(8) << vtx.ntracks
           << setw(9) << vtx.ntracksw05;
        AnaUtil::bit_print(sbit, 8, os);
      }
    }
    if (op.usesbit && sbit) continue;
    list.push_back(vtx);
  }
}
// Let's look at the Muon collection
void AnaBase::findMuonInfo(vector<Muon>& list, double vz, Options& op, std::ostream& os) {
  if (nmuon() < 1) return;

  if (op.verbose)
    os << "=>> Muons: " << nmuon() << endl
          << "indx     Eta     Phi      Pt nChamber nMatch nStation"
          << " pixHits trkHits   gChi2      D0   D0Err      dB" 
          << " vtxDistZ   trkDz trkDzErr       dz relIso            sbit"
          << endl; 

  int indx = 0;
  for (auto it = muonList_->begin(); it != muonList_->end(); ++it,indx++) {
    const Muon& muon = (*it);

    int sbit = 0;
    if (fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap_, "eta"))                      sbit |= (1 <<  0);
    if (muon.pt < AnaUtil::cutValue(muonCutMap_, "pt"))                              sbit |= (1 <<  1);
    if (!muon.isTrackerMuon)                                                         sbit |= (1 <<  2);
    if (!muon.isGlobalMuonPromptTight)                                               sbit |= (1 <<  3);
    if (muon.nChambers < AnaUtil::cutValue(muonCutMap_, "nChambers"))                sbit |= (1 <<  4);
    if (muon.nMatches < AnaUtil::cutValue(muonCutMap_, "nMatches"))                  sbit |= (1 <<  5);
    if (muon.nMatchedStations < AnaUtil::cutValue(muonCutMap_, "nMatchedStations"))  sbit |= (1 <<  6); 
    if (muon.pixHits < AnaUtil::cutValue(muonCutMap_, "pixHits"))                    sbit |= (1 <<  7);
    if (muon.trkHits < AnaUtil::cutValue(muonCutMap_, "trkHits"))                    sbit |= (1 <<  8);
    if (muon.globalChi2 >= AnaUtil::cutValue(muonCutMap_,"globalChi2"))              sbit |= (1 <<  9);
    if (abs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap_,"trkD0"))                   sbit |= (1 << 10);
    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    //if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap_, "dz"))                sbit |= (1 << 11);
    if (fabs(muon.chargedParticleIso/muon.pt) >= AnaUtil::cutValue(muonCutMap_, "pfRelIso"))                   sbit |= (1 << 12);
    // not in baseline selection
#if 0
    if (!muon.isAllArbitrated)                                                       sbit |= (1 << 13);
    if (abs(muon.dB) >= AnaUtil::cutValue(muonCutMap_,"dB"))                         sbit |= (1 << 14);
    if (abs(muon.vtxDistZ) >= AnaUtil::cutValue(muonCutMap_, "vtxDistZ"))            sbit |= (1 << 15);
#endif
    if (op.verbose) {
      bool pp = (op.printselected && sbit) ? false : true;
      if (pp) {
        os << setprecision(2)
           << setw(4) << indx 
           << setw(8) << muon.eta
           << setw(8) << muon.phi
           << setw(8) << muon.pt
           << setw(9) << muon.nChambers
           << setw(7) << muon.nMatches
           << setw(9) << muon.nMatchedStations
           << setprecision(1)
           << setw(8) << muon.pixHits
           << setw(8) << muon.trkHits
           << setw(8) << muon.globalChi2
           << setprecision(3)
           << setw(8) << muon.trkD0
	  //<< setw(8) << muon.trkD0Error
           << setw(8) << muon.dB
           << setw(9) << muon.vtxDistZ
           << setw(8) << muon.trkDz
	  //<< setw(9) << muon.trkDzError
           << setw(9) << dz;
	  //<< setw(9) << muon.relIso;
        AnaUtil::bit_print(sbit, 15, os);
      }
    }
    // Now apply cuts
    if (op.usesbit && sbit) continue;
    list.push_back(muon);
  }
  if (list.size() > 1) 
    sort(list.begin(), list.end(), PtComparator<Muon>());
}
// Let's look at the Electron collection now
void AnaBase::findElectronInfo(vector<Electron>& list, double vz, Options& op, std::ostream& os) {
  if (nelectron() < 1) return;
  if (op.verbose)
    os << "=>> Electrons: " << nelectron() << endl
       << "indx     Eta     Phi      Pt  Energy"
       << "   scEta   dB      eleId      dz     sbit" 
       << endl; 

  int indx = 0;
  for (auto it = electronList_->begin(); it != electronList_->end(); ++it,indx++) {
    const Electron& elec = (*it);

    bool quality_EB_loose = (abs(elec.eta) <= 1.4442 
                         &&  elec.sigmaEtaEta < 0.01 
                         &&  elec.deltaEtaTrkSC < 0.007 
                         &&  elec.deltaPhiTrkSC < 0.8
                         &&  elec.hoe < 0.15);
    bool quality_EE_loose = (abs(elec.eta) >= 1.566 
                         &&  elec.sigmaEtaEta < 0.03
                         &&  elec.deltaEtaTrkSC < 0.01 
                         &&  elec.deltaPhiTrkSC < 0.7
                         &&  elec.hoe < 0.07);
    bool quality_loose = quality_EB_loose || quality_EE_loose;  

    int sbit = 0;
    if (elec.pt <= AnaUtil::cutValue(electronCutMap_, "pt"))                             sbit |= (1 << 0);
    if (abs(elec.eta) >= AnaUtil::cutValue(electronCutMap_, "eta"))                      sbit |= (1 << 1);
    if (!quality_loose)                                                                   sbit |= (1 << 2);
    if (elec.missingHits > AnaUtil::cutValue(electronCutMap_, "missingHits"))            sbit |= (1 << 3);
    //if (elec.simpleEleId95cIso <= AnaUtil::cutValue(electronCutMap_, "eleId"))         sbit |= (1 << 2);
    //if (!elec.hasGsfTrack)                                                             sbit |= (1 << 3);
    if (abs(elec.dB) >= AnaUtil::cutValue(electronCutMap_, "dB"))                        sbit |= (1 << 4);
    if ( abs(elec.scEta) >= AnaUtil::cutValue(electronCutMap_, "scEtaLow") 
      && abs(elec.scEta) <= AnaUtil::cutValue(electronCutMap_, "scEtaUp"))               sbit |= (1 << 5);

    bool isGoodVtx;
    TVector3 vele = findLeptonVtx(elec.vtxIndex, isGoodVtx);
    double dz = (isGoodVtx) ? (vele.z() - vz) : 999;
    if (abs(dz) >= AnaUtil::cutValue(electronCutMap_, "dz"))                              sbit |= (1 << 6);

    if (op.verbose) {
      bool pp = (op.printselected && sbit) ? false : true;
      if (pp) {
        os << setprecision(2)
           << setw(4) << indx 
           << setw(8) << elec.eta
           << setw(8) << elec.phi
           << setw(8) << elec.pt
           << setw(8) << elec.energy
           << setw(8) << elec.scEta
           << setw(8) << elec.dB
           << setw(8) << dz;
        AnaUtil::bit_print(sbit, 8, os);
      }
    }
    // Now apply cuts
    if (op.usesbit && sbit) continue;
    list.push_back(elec);
  }
  if (list.size() > 1) 
    sort(list.begin(), list.end(), PtComparator<Electron>());
}
// Let's look at the Tau collection
void AnaBase::findTauInfo(vector<Tau>& list, double vz, Options& op, std::ostream& os) {
  if (ntau() < 1) return;
  if (op.verbose)
    os << "=>> Taus: " << ntau() << endl
       << "indx     Eta     Phi      Pt      dz"
       << " DMF  LI  MI aMT aEL aEM eMVA    sbit" 
       << endl; 
  int indx = 0;
  for (auto it = tauList_->begin(); it != tauList_->end(); ++it,indx++) {
    const Tau& tau = (*it);

    // pre-selection
    int sbit = 0;
    if (abs(tau.eta) >= AnaUtil::cutValue(tauCutMap_, "eta"))                 sbit |= (1 << 0); 
    if (tau.pt < AnaUtil::cutValue(tauCutMap_, "pt"))                         sbit |= (1 << 1); 
    double dz = tau.zvertex - vz; 
    if (vz != -999 && abs(dz) >= AnaUtil::cutValue(tauCutMap_, "dz"))          sbit |= (1 << 2);
    if (tau.decayModeFinding <= 0.5)                                          sbit |= (1 << 3);
    if (tau.chargedIsoPtSum <= 2.0)                     sbit |= (1 << 4); 
    if (tau.againstMuonTight <= AnaUtil::cutValue(tauCutMap_, "muVeto"))      sbit |= (1 << 5);
    if (tau.againstElectronLoose <= AnaUtil::cutValue(tauCutMap_, "eleVeto")) sbit |= (1 << 6); 

    if (op.verbose) {
      bool pp = (op.printselected && sbit) ? false : true;
      if (pp) {
        os << setprecision(2)
           << setw(4) << indx 
           << setw(8) << tau.eta
           << setw(8) << tau.phi
           << setw(8) << tau.pt
           << setw(8) << dz
           << setprecision(1)
           << setw(4) << tau.decayModeFinding
	  //<< setw(4) << tau.byLooseCombinedIsolationDeltaBetaCorr
	  //<< setw(4) << tau.byMediumCombinedIsolationDeltaBetaCorr
           << setw(4) << tau.againstMuonTight
           << setw(4) << tau.againstElectronLoose
           << setw(4) << tau.againstElectronMedium
           << setw(5) << tau.againstElectronLooseMVA6;
        AnaUtil::bit_print(sbit, 8, os);
      }
    }
    if (op.usesbit && sbit) continue;
    list.push_back(tau);
  }
  if (list.size() > 1) 
    sort(list.begin(), list.end(), PtComparator<Tau>());   
}
// Let's look at the Jet collection
void AnaBase::findJetInfo(vector<Jet>& list, Options& op, std::ostream& os) {
  if (njet() < 1) return;
  if (op.verbose)
    os << "=>> Jets: " << njet() << endl
       << "indx     Eta     Phi      Pt  Energy"
       << "    TCHE    TCHP    SVHE    SVHP JPBTag JBPBTag     sbit"
       << endl; 

  int indx = 0;
  for (auto it = jetList_->begin(); it != jetList_->end(); ++it,indx++) {
    const Jet& jt = (*it);

    int sbit = 0;
    if (abs(jt.eta) >= AnaUtil::cutValue(bjetCutMap_, "eta"))                        sbit |= (1 << 0);
    if (jt.pt <= AnaUtil::cutValue(bjetCutMap_, "pt"))                               sbit |= (1 << 1);
    if (useTCHE_) {
      //if (jt.trackCountingHighEffBTag <= AnaUtil::cutValue(bjetCutMap_, "TCHECut")) sbit |= (1 << 2);
    }
    else {
      if (jt.pfCombinedInclusiveSecondaryVertexV2BJetTags <= AnaUtil::cutValue(bjetCutMap_, "CSVCut")) sbit |= (1 << 2);
    }
    if (op.verbose) {
      bool pp = (op.printselected && sbit) ? false : true;
      if (pp) {
        os << setprecision(2)
           << setw(4) << indx 
           << setw(8) << jt.eta
           << setw(8) << jt.phi
           << setw(8) << jt.pt
           << setw(8) << jt.energy;
	  //<< setprecision(1)
	  //<< setw(8) << jt.trackCountingHighEffBTag  
          //<< setw(8) << jt.trackCountingHighPurBTag
	  //<< setw(8) << jt.simpleSecondaryVertexHighEffBTag
	  //<< setw(8) << jt.simpleSecondaryVertexHighPurBTag
	  //<< setw(8) << jt.jetProbabilityBTag
	  //<< setw(8) << jt.jetBProbabilityBTag;
        AnaUtil::bit_print(sbit, 8, os);
      }
    }
    // Now apply cuts
    if (op.usesbit && sbit) continue;
    list.push_back(jt);
  }
  if (list.size() > 1) 
    sort(list.begin(), list.end(), PtComparator<Jet>());
}
void AnaBase::findTriggerObjectInfo(vector<TriggerObject>& list) {
  for (auto it = triggerObjList_->begin(); it != triggerObjList_->end(); ++it) {
    list.push_back(*it);
  }
}
TVector3 AnaBase::findLeptonVtx(int index, bool& isGoodVtx) {
  const Vertex& vtx = vertexList_->at(index);

  isGoodVtx = false;
  double dxy = sqrt(pow(vtx.x, 2) + pow(vtx.y, 2));
  if (vtx.ndf > AnaUtil::cutValue(vtxCutMap_, "ndf") && 
      dxy < AnaUtil::cutValue(vtxCutMap_, "dxy") && 
      abs(vtx.z) < AnaUtil::cutValue(vtxCutMap_, "z")) isGoodVtx = true;

  TVector3 v(vtx.x, vtx.y, vtx.z); 
  return v;
}
void AnaBase::dumpEvent(const char* optstr, std::ostream& os, bool ps) {
  // Dump original content present in the tree
  // Event
  const Event& ev = eventList_->at(0);

  os << "Event " << ev.event 
     << " Lumis " << ev.lumis 
     << " Run " << ev.run 
     << endl;

  int logOption = strtol(optstr, NULL, 2);

  // Options common for all the objects
  Options options;
  options.usesbit = false;
  options.printselected = ps;

  vector<Vertex> vList;
  options.verbose = (logOption & 0x1) ? true : false;
  findVtxInfo(vList, options, os);
  double vz = (vList.size() > 0) ? vList.at(0).z : -999;

  vector<Tau> tList;
  options.verbose = (logOption >> 1 & 0x1) ? true : false;
  findTauInfo(tList, vz, options, os);

  vector<Muon> mList;
  options.verbose = (logOption >> 2 & 0x1) ? true : false;
  findMuonInfo(mList, vz, options, os);

  vector<Electron> eList;
  options.verbose = (logOption >> 3 & 0x1) ? true : false;
  findElectronInfo(eList, vz, options, os);

  vector<Jet> bList;
  options.verbose = (logOption >> 4 & 0x1) ? true : false;
  findJetInfo(bList, options, os);

  if (isMC_ && (logOption >> 5 & 0x1)) dumpGenInfo(os); 
}
void AnaBase::findGenInfo(vector<GenParticle>& genMuonList, vector<GenParticle>& genTauList) {
  // Generator level information
  if (!ngenparticle()) return;
  for (auto jt = genParticleList_->begin(); jt != genParticleList_->end(); ++jt) {
    const GenParticle& gp = (*jt);
      
    int pdgid = abs(gp.pdgId);
    if (pdgid != 13 && pdgid != 15) continue;

    int status = gp.status;
    // looking for a Hadronically decaying Tau whose mother is Higgs
    if (pdgid == 15 && status == 2) {
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      if (index < 0 || abs(mmid) != 25) continue; // Assert Higgs

      vector<int> d = gp.daughtIndices;
      int ndau = 0; 
      for (auto it = d.begin(); it != d.end(); ++it) {
        int di = (*it);
        if (di >= ngenparticle()) continue;
        const GenParticle& dgp = genParticleList_->at(di);

        int pid = abs(dgp.pdgId);
        //if (pid == 16) continue;
        //if (pid != 11 && pid != 13 && pid != 12 && pid != 14 && pid != 22) ++ndau;
        if (pid > 100 || pid == 24) ++ndau;
      }
      if (ndau) genTauList.push_back(gp);
    }
    // Looking for Muon whose mother is W boson  
    else if (pdgid == 13 && status == 1) {
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      if (index >= 0 && abs(mmid) == 24) {
        const GenParticle& mgp = genParticleList_->at(index);
        int jmid = -1;
        int jndex = getMotherId(mgp, jmid);
        if (jndex >= 0 && abs(jmid) != 6) genMuonList.push_back(gp); // ensure W does not come from t
      }
    }
  }
  if (genMuonList.size() > 1) 
    sort(genMuonList.begin(), genMuonList.end(), PtComparator<GenParticle>());
  if (genTauList.size() > 1)
    sort(genTauList.begin(), genTauList.end(), PtComparator<GenParticle>());
}

void AnaBase::findSigDecayGenInfo(vector<GenParticle>& genHDecMuonList, vector<GenParticle>& genWDecMuonList, 
                                  vector<GenParticle>& genHDecTauList,vector<GenParticle>& genWDecTauList, 
                                  vector<GenParticle>& genHDecEleList, vector<GenParticle>& genWDecEleList) {
  // Generator level information
  if (!ngenparticle()) return;
  for (auto jt = genParticleList_->begin(); jt != genParticleList_->end(); ++jt) {
    const GenParticle& gp = (*jt);
      
    int pdgid = abs(gp.pdgId);
    if (pdgid != 11 && pdgid != 13 && pdgid != 15) continue;

    int status = gp.status;
    // looking for a Hadronically decaying Tau whose mother is Higgs or W
    if (pdgid == 15 && status == 2) {
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      if (index < 0 || (abs(mmid) != 25 && abs(mmid) != 24)) continue; // Assert Higgs or W

      vector<int> d = gp.daughtIndices;
      int ndau = 0; 
      for (auto it = d.begin(); it != d.end(); ++it) {
        int di = (*it);
        if (di >= ngenparticle()) continue;
        const GenParticle& dgp = genParticleList_->at(di);

        int pid = abs(dgp.pdgId);
        //if (pid == 16) continue;
        //if (pid != 11 && pid != 13 && pid != 12 && pid != 14 && pid != 22) ++ndau;
        if (pid > 100 || pid == 24) ++ndau;
      }
      if (ndau > 0 && abs(mmid) == 25) genHDecTauList.push_back(gp);
      else if (ndau > 0 && abs(mmid) == 24){
        const GenParticle& mgp = genParticleList_->at(index);
        int jmid = -1;
        int jndex = getMotherId(mgp, jmid);
        if (jndex >= 0 && fabs(jmid) != 6 && fabs(mmid)!=15) genWDecTauList.push_back(gp); // ensure W does not come from t
      }
    }
    // Looking for Lepton whose mother is W,H boson  
    if ((pdgid == 11 || pdgid == 13) && status == 1) {
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      if(index < 0 || (fabs(mmid) != 15 && fabs(mmid) != 24)) continue;
      const GenParticle& mgp = genParticleList_->at(index);
      int jmid = -1;
      int jndex = getMotherId(mgp, jmid);
      if(jndex <0 ) continue;
      if(fabs(mmid) ==15 && fabs(jmid) == 25){
        if(pdgid == 11) genHDecEleList.push_back(gp);
        if(pdgid == 13) genHDecMuonList.push_back(gp);
      }
      else if(jndex >= 0 && fabs(jmid) != 6 && fabs(mmid)!=15){ // ensure W does not come from t
        if(pdgid == 11) genWDecEleList.push_back(gp);
        if(pdgid == 13) genWDecMuonList.push_back(gp);
      }
    }
  }
  if (genWDecEleList.size() > 1) 
    sort(genWDecEleList.begin(), genWDecEleList.end(), PtComparator<GenParticle>());
  if (genHDecEleList.size() > 1) 
    sort(genHDecEleList.begin(), genHDecEleList.end(), PtComparator<GenParticle>());

  if (genWDecMuonList.size() > 1) 
    sort(genWDecMuonList.begin(), genWDecMuonList.end(), PtComparator<GenParticle>());
  if (genHDecMuonList.size() > 1) 
    sort(genHDecMuonList.begin(), genHDecMuonList.end(), PtComparator<GenParticle>());

  if (genWDecTauList.size() > 1) 
    sort(genWDecTauList.begin(), genWDecTauList.end(), PtComparator<GenParticle>());
  if (genHDecTauList.size() > 1) 
    sort(genHDecTauList.begin(), genHDecTauList.end(), PtComparator<GenParticle>());
}
bool AnaBase::WhHmhFilter(){
  // Generator level information
  bool WTau = false;
  bool HiggsHTau = false;
  bool HiggsMTau = false;
  if (!ngenparticle()) return false;
  for (auto jt = genParticleList_->begin(); jt != genParticleList_->end(); ++jt) {
    const GenParticle& gp = (*jt);

    int pdgid = abs(gp.pdgId);
    if (pdgid != 13 && pdgid != 15) continue;

    int status = gp.status;
    // looking for a Hadronically decaying Tau whose mother is Higgs
    if (pdgid == 15 && status >= 2) {
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      if (abs(mmid) != 25 && abs(mmid) != 24) continue; // Assert Higgs

      vector<int> d = gp.daughtIndices;
      int nHdau = 0;
      int nMdau = 0;
      for (auto it = d.begin(); it != d.end(); ++it) {
        int di = (*it);
        if (di >= ngenparticle()) continue;
	const GenParticle& dgp = genParticleList_->at(di);

        int pid = abs(dgp.pdgId);
        //if (pid == 16) continue;
        //if (pid != 11 && pid != 13 && pid != 12 && pid != 14 && pid != 22) ++ndau;
        if (pid > 100 || pid == 24) ++nHdau;
        if (pid == 13) ++nMdau;
      }
      if (abs(mmid) == 25 && nHdau > 0 )HiggsHTau = true;
      if (abs(mmid) == 24 && nHdau > 0 )WTau = true;
      if (nMdau == 1 )HiggsMTau = true;
    }
  }
  if (HiggsHTau && HiggsMTau && WTau) return true;
  else return false;
}
bool AnaBase::WmHmhFilter(){
  // Generator level information
  bool MuonOk = false;
  bool HiggsHTau = false;
  bool HiggsMTau = false;
  if (!ngenparticle()) return false;
  for (auto jt = genParticleList_->begin(); jt != genParticleList_->end(); ++jt) {
    const GenParticle& gp = (*jt);

    int pdgid = abs(gp.pdgId);
    if (pdgid != 13 && pdgid != 15) continue;

    int status = gp.status;
    // looking for a Hadronically decaying Tau whose mother is Higgs
    if (pdgid == 15 && status >= 2) {
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      if (index < 0 || abs(mmid) != 25) continue; // Assert Higgs

      vector<int> d = gp.daughtIndices;
      int nHdau = 0;
      int nMdau = 0;
      for (auto it = d.begin(); it != d.end(); ++it) {
        int di = (*it);
        if (di >= ngenparticle()) continue;
	const GenParticle& dgp = genParticleList_->at(di);

        int pid = abs(dgp.pdgId);
        //if (pid == 16) continue;
        //if (pid != 11 && pid != 13 && pid != 12 && pid != 14 && pid != 22) ++ndau;
        if (pid > 100 || pid == 24) ++nHdau;
        if (pid == 13) ++nMdau;
      }
      if (nHdau > 0 )HiggsHTau = true;
      if (nMdau == 1 )HiggsMTau = true;
    }
    // Looking for Muon whose mother is W boson
    else if (pdgid == 13 && status == 1) {
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      if (index >= 0 && abs(mmid) == 24) {
        const GenParticle& mgp = genParticleList_->at(index);
        int jmid = -1;
        int jndex = getMotherId(mgp, jmid);
        if (jndex >= 0 && abs(jmid) != 6) MuonOk = true; // ensure W does not come from t
      }
    }
  }
  if (HiggsHTau && HiggsMTau && MuonOk) return true;
  else return false;
}

void AnaBase::findZtype(vector<GenParticle>& genMuonList, vector<GenParticle>& genEleList, vector<GenParticle>& genTauList) {
  // Generator level information
  if (!ngenparticle()) return;
  for (auto jt = genParticleList_->begin(); jt != genParticleList_->end(); ++jt) {
    const GenParticle& gp = (*jt);
      
    int pdgid = abs(gp.pdgId);
    if (pdgid != 11 && pdgid != 13 && pdgid != 15) continue;

    int status = gp.status;
    // looking for a Hadronically decaying Tau whose mother is Z
    if (pdgid == 15 && status == 2) {
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      if (index < 0 || abs(mmid) != 23) continue; // Assert Z

      vector<int> d = gp.daughtIndices;
      int ndau = 0; 
      for (auto it = d.begin(); it != d.end(); ++it) {
        int di = (*it);
        if (di >= ngenparticle()) continue;
        const GenParticle& dgp = genParticleList_->at(di);

        int pid = abs(dgp.pdgId);
        //if (pid == 16) continue;
        //if (pid != 11 && pid != 13 && pid != 12 && pid != 14 && pid != 22) ++ndau;
        if (pid > 100 || pid == 24) ++ndau;
      }
      if (ndau) genTauList.push_back(gp);
    }
    // looking for a e+e- decay of Z
    if (pdgid == 11 && status == 1) {
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      if (index < 0 || abs(mmid) != 23) continue; // Assert Z

      genEleList.push_back(gp);
    }
    // looking for a mu+mu- decay of Z
    if (pdgid == 13 && status == 1) {
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      if (index < 0 || abs(mmid) != 23) continue; // Assert Z

      genMuonList.push_back(gp);
    }
  }
  if (genEleList.size() > 1) 
    sort(genEleList.begin(), genEleList.end(), PtComparator<GenParticle>());
  if (genMuonList.size() > 1) 
    sort(genMuonList.begin(), genMuonList.end(), PtComparator<GenParticle>());
  if (genTauList.size() > 1)
    sort(genTauList.begin(), genTauList.end(), PtComparator<GenParticle>());
}

int AnaBase::getMotherId(const GenParticle& gp, int& mmid) const {
  int pdgid = gp.pdgId;

  vector<int> m = gp.motherIndices;
  if (m.size() < 1) return -1;
  int indx = m.at(0);
  GenParticle& mgp = genParticleList_->at(indx);
  mmid = mgp.pdgId;
  while (mmid == pdgid) {
    m = mgp.motherIndices;
    indx = m.at(0);
    mgp = genParticleList_->at(indx);
    mmid = mgp.pdgId;
  }
  return indx;
}
void AnaBase::dumpGenInfo(std::ostream& os) const {
  if (!ngenparticle()) return;

  os << setprecision(2);
  os << "indx    status    pdgId     eta      phi      pt     energy  moIndx"
     << "      moID                   daughterID"
     << endl;
  int indx = 0;
  for (auto jt = genParticleList_->begin(); jt != genParticleList_->end(); ++jt,indx++) {
    const GenParticle& gp = (*jt);

    std::ostringstream mID;
    vector<int> m = gp.motherIndices;
    for (auto it = m.begin(); it != m.end(); ++it) {
      int mi = (*it);
      if (mi >= ngenparticle()) continue;
      const GenParticle& mgp = genParticleList_->at(mi);
      mID << " " << mgp.pdgId; 
    }
    string ms = mID.str();
    if (!ms.length()) ms = " -";
    
    std::ostringstream dID;
    vector<int> d = gp.daughtIndices;
    for (auto it = d.begin(); it != d.end(); ++it) {
      int di = (*it);
      if (di >= ngenparticle()) continue;
      const GenParticle& dgp = genParticleList_->at(di);
      double energy = dgp.energy;
      int pdgid = dgp.pdgId;
      if (abs(pdgid) == 21 && energy <= 10) continue;
      dID << " " << dgp.pdgId; 
    }
    string ds = dID.str();
    if (!ds.length()) ds = " -";

    os << setw(4)  << indx
       << setw(8)  << gp.status
       << setw(10) << gp.pdgId
       << setw(10) << gp.eta
       << setw(9)  << gp.phi
       << setw(9)  << gp.pt
       << setw(9)  << gp.energy
       << setw(8)  << gp.motherIndex 
       << setw(10) << ms 
       << ds
       << endl;
  }
}


bool AnaBase::readPileUpHist() {
  size_t found = puHistFile_.find(".root");
  if (found == string::npos) {
    cerr << ">>> Warning: <<" << puHistFile_ << ">> does not have .root extension!!" << endl;
    return false;
  }

  const char* fname = gSystem->ExpandPathName(puHistFile_.c_str());
  if (gSystem->AccessPathName(fname)) {
    cerr << ">>> Warning: File <<" << puHistFile_ << ">> was not found!!" << endl;
    return false;
  }

  TFile file(fname);
  TH1D *h = dynamic_cast<TH1D*>(file.Get("pileup"));

  int nx = h->GetXaxis()->GetNbins();
  for (int i = 0; i < nx; ++i) {

    double wt = h->GetBinContent(i);
    puWtList_.push_back(wt);
  }
  return true;
}

double AnaBase::  wtPileUp(int& nPU) const {
  if (!usePUWt_) return 1.0;;
  nPU = 0;

size_t found = puHistFile_.find(".root");
  if (found == string::npos) {
    cerr << ">>> Warning: <<" << puHistFile_ << ">> does not have .root extension!!" << endl;
    return false;
  }

  const char* fname = gSystem->ExpandPathName(puHistFile_.c_str());
  if (gSystem->AccessPathName(fname)) {
    cerr << ">>> Warning: File <<" << puHistFile_ << ">> was not found!!" << endl;
    return false;
  }

  TFile file(fname);
  TH1D *h = dynamic_cast<TH1D*>(file.Get("pileup"));

  const Event& evt = eventList_->at(0);

  auto list = (useTrueNInt_) ? evt.trueNInt : evt.nPU;
  if (!list.size()) return 1.0;


/*
  int nbins = h->GetNbinsX();
  //Important, to use the nPV, clean this up 
  //if (!usenPV_) nPU = list.at(0)
   nPU = list.at(0);
  if (nPU < 0) nPU = 0;
  if (nPU >= nbins) nPU = nbins - 1;
  if (nPU > h->GetBinCenter(nbins)) 
  
  nPU = h->GetBinCenter(nbins);
  
  double bin = h->GetXaxis()->FindBin(nPU);
  double wt  = h->GetBinContent(bin);               
  //std::cout << "nPU: " << nPU << ", wt: " << wt << ", PUHistBin No: " << bin  <<  std::endl;
  return wt;
  
*/

int nbins = puWtList_.size();

nPU = list.at(0);
if (nPU < 0 ) nPU = 0 ;
if (nPU >= nbins) nPU = nbins - 1 ;
return puWtList_.at(nPU);
 



}

/*
double AnaBase::wtPileUp(int& nPU) const {
  if (!usePUWt_) return 1.0;;
  nPU = 0;

  const Event& evt = eventList_->at(0);

  auto list = (useTrueNInt_) ? evt.trueNInt : evt.nPU;
  if (!list.size()) return 1.0;

  int nbins = puWtList_.size();
  //Important, to use the nPV, clean this up 
  //if (!usenPV_) nPU = list.at(0);

  nPU = list.at(0);
  if (nPU < 0) nPU = 0;
//  if (nPU >= nbins) nPU = nbins - 1;
  cout <<"wt in the nPU "<<puWtList_.at(nPU)<<endl;
  return puWtList_.at(nPU);

}
*/

bool AnaBase::isTriggered(bool check_prescale, bool verbose) const {
  bool flag = false;
  for (size_t i = 0; i < hltpaths_->size(); ++i) {
    string str = (*hltpaths_).at(i);
    bool found = false;
    for (auto it = trigPathList_.begin(); it != trigPathList_.end(); ++it) {
      if (str.find(*it) != string::npos) {
	found = true;
        break;
      }
    }
    if (!found) continue;

    int prescl = (check_prescale) ? (*hltprescales_).at(i) : 1;
    int result = (*hltresults_).at(i);
    if (verbose) cout << ">>> HLT Path = " << str 
                      << ", fired=" << result
                      << ", prescale=" << prescl
                      << endl;
    if (result == 1 && prescl == 1) {
      flag = true;
      break;
    }
  }
  return flag;
}
void AnaBase::dumpTriggerPaths(std::ostream& os, bool check_prescale) const 
{
  os << "=> Trigger paths" << endl;
  os << setw(96) << "Path"
     << setw(8) << "prescl"
     << setw(8) << "result"
     << endl;
  for (uint i = 0; i < hltpaths_->size(); ++i) {
    string path_name = (*hltpaths_).at(i);
    int prescale     = (*hltprescales_).at(i);  
    int result       = (*hltresults_).at(i);  
    if ((check_prescale && prescale != 1) || result != 1) continue;
    os << setw(96) << path_name 
       << setw(8) << prescale 
       << setw(8) << result
       << endl;
  }
}
void AnaBase::dumpTriggerObjectInfo(const vector<TriggerObject>& list, std::ostream& os) const
{
  os << setprecision(2);
  os << "=> TriggerObjects: " << list.size() << endl;
  os << "Indx     Eta     Phi      Pt  Energy            =Trigger path list=" << endl;
  int indx = 0;
  for (auto it  = list.begin(); it != list.end(); ++it,++indx) {
    const TriggerObject& tobj = (*it);
    os << setw(4) << indx 
       << setw(8) << tobj.eta
       << setw(8) << tobj.phi
       << setw(8) << tobj.pt
       << setw(8) << tobj.energy
       << endl;
    map<string, uint> path_list = tobj.pathList;
    for (auto jt = path_list.begin(); jt != path_list.end(); ++jt) {
      os << "\t\t\t\t\t" << jt->first << " " << jt->second << endl;
    }
  }
}
bool AnaBase::matchTriggerPath(const vector<string>& pathList, const string& path) const {
  bool result = false;
  for (auto it = pathList.begin(); it != pathList.end(); ++it) {
    string lname = (*it);
    if (path.find(lname) != string::npos) {
      result = true;
      break;
    }
  }
  return result;
}
double AnaBase::matchTriggerObject(const vector<TriggerObject>& trigObjList, 
                                   const TLorentzVector& obj, 
                                   const string& trigPath, 
                                   int trig_skip, 
                                   double maxPtDiff,
                                   int& trig_indx) const
{
  double dRmin = 999; 
  trig_indx = -1;
  double obj_pt = obj.Pt();
  int indx = 0;

  for (auto it = trigObjList.begin(); it != trigObjList.end(); ++it,++indx) {
    if (indx == trig_skip) continue;
    const TriggerObject& trigObj = (*it);
    const map<string, uint>& path_list = trigObj.pathList;
    bool matched = false;
    for (auto istr = path_list.begin(); istr != path_list.end(); ++istr) {
      string path = istr->first;
      int flag = istr->second;
      if (path.find(trigPath) != string::npos && flag >= 1) {
        matched = true;
        break;
      }
    }
    if (!matched) continue;

    // check deltaPt
    if (fabs(trigObj.pt - obj_pt) > maxPtDiff) continue;

    TLorentzVector trigTL;
    trigTL.SetPtEtaPhiE(trigObj.pt, trigObj.eta, trigObj.phi, trigObj.energy);
    double dR = AnaUtil::deltaR(obj, trigTL);
    if (dR < dRmin) {
      dRmin = dR;
      trig_indx = indx;
    }
  }
  return dRmin;
}

std::tuple<double, int> AnaBase::matchTriggerObject(const vector<TriggerObject>& trigObjList, 
                                   const TLorentzVector& obj, 
                                   const vector<string>& trigLegList, 
                                   int trig_skip, 
                                   double maxPtDiff,
                                   int& trig_indx, 
                                   uint& flag) const
{
 double dRmin = 999; 
  trig_indx = -1;
  double obj_pt = obj.Pt();
  int indx = 0;
  for (vector<TriggerObject>::const_iterator it  = trigObjList.begin(); 
                                             it != trigObjList.end(); ++it,indx++) {
    if (indx == trig_skip) continue;
    const TriggerObject& trigObj = (*it);
    const map<string, uint>& path_list = trigObj.pathList;
    bool matched = false;
    flag = 0;
    for (map<string, uint>::const_iterator istr  = path_list.begin(); 
                                           istr != path_list.end(); ++istr) {
      string pname = istr->first;
      flag = istr->second;
      if (matchTriggerPath(trigLegList, pname)) {
        matched = true;
        break;
      }
    }
    if (!matched) continue;

    // check deltaPt
    if (abs(trigObj.pt - obj_pt) > maxPtDiff) continue;

    TLorentzVector trigTL;
    trigTL.SetPtEtaPhiE(trigObj.pt, trigObj.eta, trigObj.phi, trigObj.energy);
    double dR = AnaUtil::deltaR(obj, trigTL);
    if (dR <dRmin) {
      dRmin = dR;
      trig_indx = indx;
    }
  }
  return std::make_tuple(dRmin,trig_indx);
}





int AnaBase::vetoMuon(double zvTau, double vetoPtCut, double dzTauCut) {
  int nm = 0;
  for (auto it = muonList_->begin(); it != muonList_->end(); ++it) {
    const Muon& muon = (*it);

    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double delz = (isGoodVtx) ? abs(zvTau - vmu.z()) : 999; 
    
    if (!(muon.isTrackerMuon) ||
	!(muon.isGlobalMuonPromptTight) ||
        std::abs(muon.eta) >= AnaUtil::cutValue(muonCutMap_, "eta") ||
        muon.pt <= vetoPtCut ||  
        delz >= dzTauCut) continue;
    ++nm;
  }
  return nm;
}
int AnaBase::vetoMuon_wPt(double vetoPtCut) {
  int nm = 0;
  for (auto it = muonList_->begin(); it != muonList_->end(); ++it) {
    const Muon& muon = (*it);

    if (!(muon.isTrackerMuon) ||
	!(muon.isGlobalMuonPromptTight) ||
        std::abs(muon.eta) >= AnaUtil::cutValue(muonCutMap_, "eta") ||
        muon.pt <= vetoPtCut) continue;
    ++nm;
  }
  return nm;
}
int AnaBase::vetoElectron(double zvTau, double vetoPtCut, double dzTauCut) {
  int nel = 0;
  for (auto it = electronList_->begin(); it != electronList_->end(); ++it) {
    const Electron& ele = (*it);

    bool isGoodVtx;
    TVector3 vele = findLeptonVtx(ele.vtxIndex, isGoodVtx);
    double delz = (isGoodVtx) ? abs(zvTau - vele.z()) : 999; 
    
    if (fabs(ele.eta) >= AnaUtil::cutValue(electronCutMap_, "eta")) continue;
    if ( (fabs(ele.eta) >= AnaUtil::cutValue(electronCutMap_, "etaLow") && 
          fabs(ele.eta) <= AnaUtil::cutValue(electronCutMap_, "etaUp")) ||
          ele.pt <= vetoPtCut || 
	 //          ele.simpleEleId95cIso != AnaUtil::cutValue(electronCutMap_, "eleId") ||
          delz >= dzTauCut) continue;
         
    if (!eleId(ele, 25)) continue;
    ++nel;
  }
  return nel;
}
int AnaBase::vetoElectron_wPt(double vetoPtCut) {
  int nel = 0;
  for (auto it = electronList_->begin(); it != electronList_->end(); ++it) {
    const Electron& ele = (*it);

    if (fabs(ele.eta) >= AnaUtil::cutValue(electronCutMap_, "eta")) continue;
    if ( (fabs(ele.eta) >= AnaUtil::cutValue(electronCutMap_, "etaLow") && 
          fabs(ele.eta) <= AnaUtil::cutValue(electronCutMap_, "etaUp")) ||
          ele.pt <= vetoPtCut 
       ) continue;
    if (!eleId(ele, 25)) continue;
    ++nel;
  }
  return nel;
}
bool AnaBase::TightEleId(const Electron&  ele){
  bool elIdtight = false ;
  bool quality_EB_tight = fabs(ele.eta) <= 1.479
      && ele.deltaEtaTrkSC < 0.009
      && ele.deltaPhiTrkSC < 0.03
      && ele.sigmaIEtaIEta < 0.01
      && ele.hoe < 0.12
      && fabs(1/ele.caloEnergy - ele.eop/ele.caloEnergy) < 0.05
      && abs(ele.vtxDistZ) < 0.1
      && fabs(ele.trkD0) < 0.02
      && !ele.hasMatchedConv
      && ele.missingHits <= 0
      && ele.pfRelIso < 0.18;

  bool quality_EE_tight = 1.479 < fabs(ele.eta) && fabs(ele.eta) < 2.5
      && ele.deltaEtaTrkSC < 0.01
      && ele.deltaPhiTrkSC < 0.02
      && ele.sigmaIEtaIEta < 0.031
      && ele.hoe < 0.12
      && fabs(1/ele.caloEnergy - ele.eop/ele.caloEnergy) < 0.05
      && abs(ele.vtxDistZ) < 0.1
      && fabs(ele.trkD0) < 0.02
      && !ele.hasMatchedConv
      && ele.missingHits <= 0
      && ele.pfRelIso < 0.16;

  elIdtight = (quality_EB_tight || quality_EE_tight);

  return elIdtight;
}

bool AnaBase::eleId(const vhtm::Electron& ele, int bx=25){
  bool elIdloose=false;
  if(bx==50) {
    bool quality_EB_loose = fabs(ele.eta) <= 1.479
      && ele.deltaEtaTrkSC < 0.012
      && ele.deltaPhiTrkSC < 0.15
      && ele.sigmaIEtaIEta < 0.012
      && ele.hoe < 0.12
      && fabs(1/ele.caloEnergy - ele.eop/ele.caloEnergy) < 0.05
      && !ele.hasMatchedConv // corrected for MiniAOD!!!
      && ele.missingHits <= 1;
    
    bool quality_EE_loose = 1.479 < fabs(ele.eta) && fabs(ele.eta) < 2.5
      && ele.deltaEtaTrkSC < 0.021
      && ele.deltaPhiTrkSC < 0.10
      && ele.sigmaIEtaIEta < 0.033
      && ele.hoe < 0.12
      && fabs(1/ele.caloEnergy - ele.eop/ele.caloEnergy) < 0.05
      && !ele.hasMatchedConv
      && ele.missingHits <= 1;
    elIdloose = quality_EB_loose || quality_EE_loose;
  }
  
  else if ( bx==25) {
    bool quality_EB_loose = fabs(ele.eta) <= 1.479
      && ele.deltaEtaTrkSC < 0.012
      && ele.deltaPhiTrkSC < 0.15
      && ele.sigmaIEtaIEta < 0.01
      && ele.hoe < 0.12
      && fabs(1/ele.caloEnergy - ele.eop/ele.caloEnergy) < 0.05
      && !ele.hasMatchedConv
      && ele.missingHits <= 1;
    
    bool quality_EE_loose = 1.479 < fabs(ele.eta) && fabs(ele.eta) < 2.5
      && ele.deltaEtaTrkSC < 0.014
      && ele.deltaPhiTrkSC < 0.10
      && ele.sigmaIEtaIEta < 0.033
      && ele.hoe < 0.12
      && fabs(1/ele.caloEnergy - ele.eop/ele.caloEnergy) < 0.05
      && !ele.hasMatchedConv
      && ele.missingHits <= 1;
    
    elIdloose = quality_EB_loose || quality_EE_loose;
  }
  else {
    std::cout << ">>> Invalid bunch crossing value!!" << std::endl;
  }
  return elIdloose;
}

int AnaBase::vetoElectronSpl(double vetoPtCut) {
  int nel = 0;
#if 0
  for (auto it = electronList_->begin(); it != electronList_->end(); ++it) {
    const Electron& ele = (*it);

    bool idMVAloose = false;
    if(ele.pt < 20){
      idMVAloose = (abs(ele.eta) < 0.8 && ele.mvaPOGNonTrig > 0.925) ||
	(abs(ele.eta) > 0.8 && abs(ele.eta) < 1.479 && ele.mvaPOGNonTrig > 0.915) ||
	(abs(ele.eta) > 1.479 && ele.mvaPOGNonTrig > 0.965);
    }
    else 
      {
	idMVAloose = (abs(ele.eta) < 0.8 && ele.mvaPOGNonTrig > 0.905) ||
	  (abs(ele.eta) > 0.8 && abs(ele.eta) < 1.479 && ele.mvaPOGNonTrig > 0.955) ||
	  (abs(ele.eta) > 1.479 && ele.mvaPOGNonTrig > 0.975);
      }

    if (!idMVAloose ||
        ele.pt <= vetoPtCut ||
        (abs(ele.eta) >= 1.4442 && abs(ele.eta) <= 1.566) ||
        abs(ele.eta) >= 2.5) continue;
    ++nel;
  }
#endif
  return nel;
}
