#include <iostream>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "MVASkim_BkgkNN.h"

using std::string;
using std::cout;
using std::endl;

MVASkim_BkgkNN::MVASkim_BkgkNN(const string& filename) {
  _mvaFile = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  _tree = new TTree("RTree", "RTree");
  _tree->Branch("muEta",      &_varList.muEta,      "muEta/F");
  _tree->Branch("muPt",       &_varList.muPt,       "muPt/F");
  _tree->Branch("nJet",       &_varList.nJet,       "nJet/F");
  _mvaFile->ls();
}
MVASkim_BkgkNN::~MVASkim_BkgkNN() {
  delete _mvaFile;
}
void MVASkim_BkgkNN::fill(const BkgVariables& varList) {
  memcpy(&_varList, &varList, sizeof(varList));
  _mvaFile->cd();
  _tree->Fill();  
}
void MVASkim_BkgkNN::close() {
  _mvaFile->cd();
  _tree->Print();
  _tree->Write();
  _mvaFile->Write();
  _mvaFile->Close();
}
