#include <iostream>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "MVASkim.h"

using std::string;
using std::cout;
using std::endl;

MVASkim::MVASkim(const string& filename) {
  _mvaFile = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  _tree = new TTree("RTree", "RTree");
  _tree->Branch("lepEta",      &_varList.lepEta,      "lepEta/F");
  _tree->Branch("lepPt",       &_varList.lepPt,       "lepPt/F");
  _tree->Branch("tauOSEta",    &_varList.tauOSEta,    "tauOSEta/F");
  _tree->Branch("tauOSPt",     &_varList.tauOSPt,     "tauOSPt/F");
  _tree->Branch("tauSSEta",    &_varList.tauSSEta,    "tauSSEta/F");
  _tree->Branch("tauSSPt",     &_varList.tauSSPt,     "tauSSPt/F");
  _tree->Branch("diTaudEta",   &_varList.diTaudEta,   "diTaudEta/F");
  _tree->Branch("diTaudPhi",   &_varList.diTaudPhi,   "diTaudPhi/F");
  _tree->Branch("diTauPt",    &_varList.diTauPt,    "diTauPt/F");
  _tree->Branch("diTauMass",    &_varList.diTauMass,    "diTauMass/F");
  _tree->Branch("dphilepTauOS", &_varList.dphilepTauOS, "dphilepTauOS/F");
  _tree->Branch("dphilepTauSS", &_varList.dphilepTauSS, "dphilepTauSS/F");
  _tree->Branch("met",        &_varList.met,        "met/F");
  _tree->Branch("alpha",      &_varList.alpha,      "alpha/F");
  _tree->Branch("alphatt",      &_varList.alphatt,      "alphatt/F");
  _tree->Branch("alphalditau",      &_varList.alphalditau,      "alphalditau/F");
  _tree->Branch("alphalprod",      &_varList.alphalprod,      "alphalprod/F");
  _tree->Branch("acop",      &_varList.alpha,      "acop/F");
  _tree->Branch("dphilepDiTau", &_varList.dphilepDiTau, "dphilepDiTau/F");
  _tree->Branch("DeltaRDiTau",   &_varList.DeltaRDiTau,   "DeltaRDiTau/F");
  _tree->Branch("dzTau12Vtx",   &_varList.dzTau12Vtx,   "dzTau12Vtx/F");
  _mvaFile->ls();
}
MVASkim::~MVASkim() {
  delete _mvaFile;
}
void MVASkim::fill(const TreeVariables& varList) {
  memcpy(&_varList, &varList, sizeof(varList));
  _mvaFile->cd();
  _tree->Fill();  
}
void MVASkim::close() {
  _mvaFile->cd();
  _tree->Print();
  _tree->Write();
  _mvaFile->Write();
  _mvaFile->Close();
}
