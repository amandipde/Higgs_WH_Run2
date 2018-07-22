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
#include "TH3.h"
#include "TProfile.h"

#include "FakeRate.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"

#include "MVASkim_SigkNN.h"
#include "MVASkim_BkgkNN.h"

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
FakeRate::FakeRate()
  : AnaBase(),
    _createMVATree(false),
    _readMVA(false),
    _readMVAFK(false),
    _mvaInputFileSigkNN(""),
    _mvaInputFileBkgkNN(""),
    _MVAdisFile(""),
    _MVAFKdisFile(""),
    _skimSig(0),
    _skimBkg(0)

{}

// ----------
// Destructor
// ----------
FakeRate::~FakeRate()
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool FakeRate::beginJob() 
{ 
  AnaBase::beginJob();


  histf()->cd();
  bookHistograms();

  // Optionally write selected events to another tree                                                                                                $
  if (_createMVATree) _skimSig = new MVASkim_SigkNN(_mvaInputFileSigkNN);
  if (_createMVATree) _skimBkg = new MVASkim_BkgkNN(_mvaInputFileBkgkNN);

  if (_readMVAFK) {
    reader1 = new TMVA::Reader("Silent");
    reader1->AddVariable("LeadPt", &leadPt);
    reader1->AddVariable("SubPt", &subPt);
    reader1->AddVariable("Met", &met);
    reader1->AddVariable("DeltaRDiTau", &deltaRDiTau);
    reader1->AddVariable("PtRatio", &ptRatio);

    //reader->BookMVA("MLP::MLPBNN", "TMVAClassification_MLPBNN.weights.xml");
    //reader->BookMVA("MLPBNN", _MVAdisFile);
    reader1->BookMVA("BDT8", _MVAFKdisFile);
  }

  if (_readMVA) {
    reader = new TMVA::Reader("Silent");
    reader->AddVariable("muEta", &muEta);
    reader->AddVariable("muPt", &muPt);
    reader->AddVariable("tau1Eta", &tau1Eta);
    reader->AddVariable("tau1Pt", &tau1Pt);
    //reader->AddVariable("tau2Eta", &tau2Eta);
    reader->AddVariable("tau2Pt", &tau2Pt);
    //reader->AddVariable("diTaudR", &diTaudR);
    //reader->AddVariable("dphiMuTau1", &dphiMuTau1);
    //reader->AddVariable("dphiMuDiTau", &dphiMuDiTau);
    reader->AddVariable("met", &met);
    reader->AddVariable("diTauPt/(tau1Pt+tau2Pt)", &ptRatio);

    //reader->BookMVA("MLP::MLPBNN", "TMVAClassification_MLPBNN.weights.xml");
    //reader->BookMVA("MLPBNN", _MVAdisFile);
    reader->BookMVA("MLP", _MVAdisFile);
  }

  return true;
}
// ---------------
// Book histograms
// ---------------
void FakeRate::bookHistograms() 
{



  new TH1F("npu", "npu distrubution//MC only", 80, 0, 80);
  new TH1F("trueNInt", "no of true int distrubution//MC only", 80, 0, 80);
  new TH1F("npv", "npv distrubution", 80, 0, 80);
  new TH1F("nPVCorrected", "npv distrubution", 80, 0, 80);
  new TH1F("npv_PU", "npv distrubution with PU correction applied//MC only", 80, 0, 80);
  new TH1F("puwt", "PU weight factor", 80, 0, 80);
 
  new TH1D("getnEvents", "No of events to analyze ", 500001, -250000.5, 250000.5);
  new TH1D("evcount", "Selected event counter", 12, 0, 12);


  new TH1D("eventWeight", "Event Weight", 500001, -250000.5, 250000.5);
  new TH1D("eventWeight_initial", "Event Weight", 1, 0.5, 1.5);
  new TH1D("eventWeight_fakeCR_mtt", "Event Weight", 1, 0.5, 1.5);
  new TH1D("eventWeight_CR_mtt", "Event Weight", 1, 0.5, 1.5);
  new TH1D("eventWeight_BL_mtt", "Event Weight", 1, 0.5, 1.5);
  new TH1D("eventWeight_AA_mtt", "Event Weight", 1, 0.5, 1.5);
  new TH1D("eventWeight_BL_mmt", "Event Weight", 1, 0.5, 1.5);
  new TH1D("eventWeight_leadFakeCR_mmt", "Event Weight", 1, 0.5, 1.5);
  new TH1D("eventWeight_slFakeCR_mmt", "Event Weight", 1, 0.5, 1.5);
  new TH1D("eventWeight_bothFakeCR_mmt", "Event Weight", 1, 0.5, 1.5);

  new TH1F("fakemu_deno", "Denominator of the JetToMu Fake Function", 140, 0, 140);
  new TH1F("fakemu_nume", "Numerator of the JetToMu Fake Function", 140, 0, 140);

  new TH1F("faketaupt_deno_DYAll", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_DYC", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_DYI", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_DYF", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_DYAll", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_DYC", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_DYI", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_DYF", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);

  new TH1F("faketaupt_deno_WJetsAll_SSSS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsC_SSSS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsI_SSSS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsF_SSSS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsAll_SSSS", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsC_SSSS", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsI_SSSS", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsF_SSSS", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);

  new TH1F("faketaupt_deno_WJetsAll_OSSS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsC_OSSS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsI_OSSS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsF_OSSS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsAll_OSSS", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsC_OSSS", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsI_OSSS", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsF_OSSS", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);

  new TH1F("faketaupt_deno_WJetsAll_SSOS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsC_SSOS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsI_SSOS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsF_SSOS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsAll_SSOS", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsC_SSOS", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsI_SSOS", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsF_SSOS", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);

  new TH1F("faketaupt_deno_WJetsAll_NoW", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsC_NoW", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsI_NoW", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_WJetsF_NoW", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsAll_NoW", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsC_NoW", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsI_NoW", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_WJetsF_NoW", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);

  new TH1F("faketaupt_deno_TTJetsAll_2j", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_TTJetsC_2j", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_TTJetsI_2j", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_TTJetsF_2j", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsAll_2j", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsC_2j", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsI_2j", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsF_2j", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);

  new TH1F("faketaupt_deno_TTJetsAll_2j_ss", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_TTJetsC_2j_ss", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_TTJetsI_2j_ss", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_TTJetsF_2j_ss", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsAll_2j_ss", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsC_2j_ss", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsI_2j_ss", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsF_2j_ss", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);

  new TH1F("faketaupt_deno_TTJetsAll_1j", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_TTJetsC_1j", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_TTJetsI_1j", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_TTJetsF_1j", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsAll_1j", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsC_1j", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsI_1j", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_TTJetsF_1j", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);

  new TH1F("faketaupt_deno_QCDAll_SS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_QCDC_SS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_QCDI_SS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_deno_QCDF_SS", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_QCDAll_SS", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_QCDC_SS", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_QCDI_SS", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_nume_QCDF_SS", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);

  new TH1F("faketaupt_Isolated_WJetsAll", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_AntiIsolated_WJetsAll", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_Isolated_DYAll", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_AntiIsolated_DYAll", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);

  new TH1F("faketaupt_antianti_WJetsOSS", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketaupt_antianti_WJetsSSS", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);

  new TH1F("ZMass_mu", "Z mass of the MuMu pair after selecting the Denominator of the JetToMu Fake function", 200, 0, 200);


  new TH1F("met_DYjets", "met distribution of the event in DYjets Region", 200, 0, 200);
  new TH1F("pt_mu1_DY_CR", "Mu1 Pt in the DY Region", 140, 0, 140);
  new TH1F("pt_mu2_DY_CR", "Mu2 Pt in the DY Region", 140, 0, 140);
  new TH1F("mu1_eta", "Mu1 eta in the DY Region", 50, -2.5, 2.5);
  new TH1F("mu1_phi", "Mu1 phi in the DY Region", 140, -3.5, 3.5);
  new TH1F("Delta_R_M1_M2", "deltaR angle between mu1 and mu2 after cuts of DY Region", 100, 0, 10);
  new TH1F("npv_DY", "npv distribution in the DY Region", 80, 0, 80);
  new TH1F("Relative Isolation of muons", "Relative isolation of muons in DY region", 100, 0, 10);


  new TH1F("MT_Mu_Met_Wjet_before_tau_sel","transverse mass plot for Muon and MET in Wjet Region before tau selections", 140, 0, 140);
  new TH1F("MT_Mu_Met_Wjet_after_tau_sel","transverse mass plot for Muon and MET in Wjet Region after tau selections", 140, 0, 140);
  new TH1F("met", "met distribution of the event after signal selection", 200, 0, 200);
  new TH1F("met_Wjets", "met distribution of the event in Wjets Region", 200, 0, 200);
  new TH1F("npv_Wjets", "npv distribution in the Wjets  Region", 80, 0, 80);
  
  new TH1F("tau2fake_pt", "pt distribution of antiisolated tau2::same sign to the Muon", 140, 0, 140);
  new TH1F("tau2fake_pt_mva1", "pt distribution of antiisolated tau2::same sign to the Muon", 140, 0, 140);
  new TH1F("tau2fake_pt_mva2", "pt distribution of antiisolated tau2::same sign to the Muon", 140, 0, 140);
  new TH1F("tau2fake_pt_mva3", "pt distribution of antiisolated tau2::same sign to the Muon", 140, 0, 140);
  new TH1F("tau2fake_pt_mva4", "pt distribution of antiisolated tau2::same sign to the Muon", 140, 0, 140);
  new TH1F("tau2fake_pt_mva5", "pt distribution of antiisolated tau2::same sign to the Muon", 140, 0, 140);
  new TH1F("tau2fake_pt_mva6", "pt distribution of antiisolated tau2::same sign to the Muon", 140, 0, 140);
  new TH1F("tau2fake_pt_mva7", "pt distribution of antiisolated tau2::same sign to the Muon", 140, 0, 140);

  new TH1F("tautauMass_fakeCR", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2F("tautauMass_fakeCR_p", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320, 140, 0, 140);
  new TH2F("tau1Pt_fakeCR_p", "Tau OS Pt in the SS-AntiIsolated control region", 140, 0, 140, 140, 0, 140);
  new TH2F("tau1Eta_fakeCR_p", "Tau OS Eta in the SS-AntiIsolated control region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("tau1Phi_fakeCR_p", "Tau OS Phi in the SS-AntiIsolated control region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("tau2Pt_fakeCR_p", "Tau SS Pt in the SS-AntiIsolated control region", 140, 0, 140, 140, 0, 140);
  new TH2F("tau2Eta_fakeCR_p", "Tau SS Eta in the SS-AntiIsolated control region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("tau2Phi_fakeCR_p", "Tau SS Phi in the SS-AntiIsolated control region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("muPt_fakeCR_p", "Muon Pt in the SS-AntiIsolated control region", 140, 0, 140, 140, 0, 140);
  new TH2F("muEta_fakeCR_p", "Muon Eta in the SS-AntiIsolated control region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("muPhi_fakeCR_p", "Muon Phi in the SS-AntiIsolated control region", 140, -3.5, 3.5, 140, 0, 140);

  new TH1F("tautauMass_BL", "tautau Mass in the Signal Region", 320, 0, 320);
  new TH1F("tau1Pt_BL", "Tau OS Pt in the Signal Region", 140, 0, 140);
  new TH1F("tau1Eta_BL", "Tau OS Eta in the Signal Region", 50, -2.5, 2.5);
  new TH1F("tau1Phi_BL", "Tau OS Phi in the Signal Region", 140, -3.5, 3.5);
  new TH1F("tau2Pt_BL", "Tau SS Pt in the Signal Region", 140, 0, 140);
  new TH1F("tau2Eta_BL", "Tau SS Eta in the Signal Region", 50, -2.5, 2.5);
  new TH1F("tau2Phi_BL", "Tau SS Phi in the Signal Region", 140, -3.5, 3.5);
  new TH1F("muPt_BL", "Muon Pt in the Signal Region", 140, 0, 140);
  new TH1F("muEta_BL", "Muon Eta in the Signal Region", 50, -2.5, 2.5);
  new TH1F("muPhi_BL", "Muon Phi in the Signal Region", 140, -3.5, 3.5);
  new TH3F("tautauMass_fakeCR_tau2fake_pe", "tautau Mass in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 50, -2.5, 2.5, 320, 0, 320);

  new TH1F("mvaoutput_BL", "mvaoutput in the signal region", 100, -2, 2);
  new TH3F("mvaoutput_fakeCR_tau2fake", "mvaoutput in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 50, -2.5, 2.5, 100, -2, 2);

  new TH1F("tau1Pt_CR", "Tau OS Pt in the OS-anti Control Region", 140, 0, 140);
  new TH1F("tau1Eta_CR", "Tau OS Eta in the OS-anti Control Region", 50, -2.5, 2.5);
  new TH1F("tau1Phi_CR", "Tau OS Phi in the OS-anti Control Region", 140, -3.5, 3.5);
  new TH1F("tau2Pt_CR", "Tau SS Pt in the OS-anti Control Region", 140, 0, 140);
  new TH1F("tau2Eta_CR", "Tau SS Eta in the OS-anti Control Region", 50, -2.5, 2.5);
  new TH1F("tau2Phi_CR", "Tau SS Phi in the OS-anti Control Region", 140, -3.5, 3.5);
  new TH1F("muPt_CR", "Muon Pt in the OS-anti Control Region", 140, 0, 140);
  new TH1F("muEta_CR", "Muon Eta in the OS-anti Control Region", 50, -2.5, 2.5);
  new TH1F("muPhi_CR", "Muon Phi in the OS-anti Control Region", 140, -3.5, 3.5);
  new TH1F("tautauMass_CR", "tautau Mass in the OS-AntiIsolated control region", 320, 0, 320);
  new TH1F("mvaoutput_CR", "mvaoutput in the signal region", 100, -2, 2);
  new TH1F("mvaOutputFK_CR", "mvaoutput in the signal region after BL", 100, -1, 1);

  new TH2F("tau1Pt_CR_e", "Tau OS Pt in the OS-anti Control Region", 140, 0, 140, 50, -2.5, 2.5);
  new TH2F("tau1Eta_CR_e", "Tau OS Eta in the OS-anti Control Region", 50, -2.5, 2.5, 50, -2.5, 2.5);
  new TH2F("tau1Phi_CR_e", "Tau OS Phi in the OS-anti Control Region", 140, -3.5, 3.5, 50, -2.5, 2.5);
  new TH2F("tau2Pt_CR_e", "Tau SS Pt in the OS-anti Control Region", 140, 0, 140, 50, -2.5, 2.5);
  new TH2F("tau2Eta_CR_e", "Tau SS Eta in the OS-anti Control Region", 50, -2.5, 2.5, 50, -2.5, 2.5);
  new TH2F("tau2Phi_CR_e", "Tau SS Phi in the OS-anti Control Region", 140, -3.5, 3.5, 50, -2.5, 2.5);
  new TH2F("muPt_CR_e", "Muon Pt in the OS-anti Control Region", 140, 0, 140, 50, -2.5, 2.5);
  new TH2F("muEta_CR_e", "Muon Eta in the OS-anti Control Region", 50, -2.5, 2.5, 50, -2.5, 2.5);
  new TH2F("muPhi_CR_e", "Muon Phi in the OS-anti Control Region", 140, -3.5, 3.5, 50, -2.5, 2.5);
  new TH2F("tautauMass_CR_e", "tautau Mass in the OS-AntiIsolated control region", 320, 0, 320, 50, -2.5, 2.5);

  new TH3F("mvaOutputFK_AA_tau2fake_pe", "mvaoutput in the anti anti region", 140, 0, 140, 50, -2.5, 2.5, 100, -1, 1);
  new TH3F("tautauMass_AA_tau2fake_pe", "tautau Mass in the anti anti control region", 140, 0, 140, 50, -2.5, 2.5, 320, 0, 320);
  new TH3F("mvaoutput_AA_tau2fake_pe", "mvaoutput in the anti anti control region", 140, 0, 140, 50, -2.5, 2.5, 100, -2, 2);
  new TH3F("tau1Pt_AA_tau2fake_pe", "tau1Pt in the anti anti control region", 140, 0, 140, 50, -2.5, 2.5, 140, 0, 140);
  new TH3F("tau2Pt_AA_tau2fake_pe", "tau2Pt in the anti anti control region", 140, 0, 140, 50, -2.5, 2.5, 140, 0, 140);
  new TH3F("muPt_AA_tau2fake_pe", "muPt in the anti anti control region", 140, 0, 140, 50, -2.5, 2.5, 140, 0, 140);

  new TH2F("mvaOutputFK_AA_tau2fake_p", "mvaoutput in the anti anti region", 100, -1, 1, 140, 0, 140);
  new TH2F("tautauMass_AA_tau2fake_p", "tautau Mass in the anti anti control region", 320, 0, 320, 140, 0, 140);
  new TH2F("tau1Pt_AA_tau2fake_p", "Tau OS Pt in the anti-anti region", 140, 0, 140, 140, 0, 140);
  new TH2F("tau1Eta_AA_tau2fake_p", "Tau OS Eta in the anti-anti region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("tau1Phi_AA_tau2fake_p", "Tau OS Phi in the anti-anti region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("tau2Pt_AA_tau2fake_p", "Tau SS Pt in the anti-anti region", 140, 0, 140, 140, 0, 140);
  new TH2F("tau2Eta_AA_tau2fake_p", "Tau SS Eta in the anti-anti region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("tau2Phi_AA_tau2fake_p", "Tau SS Phi in the anti-anti region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("muPt_AA_tau2fake_p", "Muon Pt in the anti-anti region", 140, 0, 140, 140, 0, 140);
  new TH2F("muEta_AA_tau2fake_p", "Muon Eta in the anti-anti region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("muPhi_AA_tau2fake_p", "Muon Phi in the anti-anti region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("mvaoutput_AA_tau2fake_p", "mvaoutput in the anti anti control region", 100, -2, 2, 140, 0, 140);

  new TH1F("tautauMass_AA", "tautau Mass in the anti anti control region", 320, 0, 320);
  new TH1F("tau1Pt_AA", "tau1Pt in the anti anti control region", 140, 0, 140);
  new TH1F("tau2Pt_AA", "tau2Pt in the anti anti control region", 140, 0, 140);
  new TH1F("muPt_AA", "muPt in the anti anti control region", 140, 0, 140);

  new TH1F("mvaOutputFK_BL", "mvaoutput in the signal region after BL", 100, -1, 1);
  new TH3F("mvaOutputFK_fakeCR_tau2fake", "mvaoutput in the fake CR region after BL", 140, 0, 140, 50, -2.5, 2.5, 100, -1, 1);

  new TH1F("mvaoutput_mva1", "mvaoutput in the signal region after mva1", 100, -2, 2);
  new TH1F("mvaoutput_mva2", "mvaoutput in the signal region after mva2", 100, -2, 2);
  new TH1F("mvaoutput_mva3", "mvaoutput in the signal region after mva3", 100, -2, 2);
  new TH1F("mvaoutput_mva4", "mvaoutput in the signal region after mva4", 100, -2, 2);
  new TH1F("mvaoutput_mva5", "mvaoutput in the signal region after mva5", 100, -2, 2);
  new TH1F("mvaoutput_mva6", "mvaoutput in the signal region after mva6", 100, -2, 2);
  new TH1F("mvaoutput_mva7", "mvaoutput in the signal region after mva7", 100, -2, 2);

  new TH1F("tautauMass_mva1", "tautau Mass in the signal region after BL and MVA1", 320, 0, 320);
  new TH1F("tautauMass_mva2", "tautau Mass in the signal region after BL and MVA2", 320, 0, 320);
  new TH1F("tautauMass_mva3", "tautau Mass in the signal region after BL and MVA3", 320, 0, 320);
  new TH1F("tautauMass_mva4", "tautau Mass in the signal region after BL and MVA4", 320, 0, 320);
  new TH1F("tautauMass_mva5", "tautau Mass in the signal region after BL and MVA5", 320, 0, 320);
  new TH1F("tautauMass_mva6", "tautau Mass in the signal region after BL and MVA6", 320, 0, 320);
  new TH1F("tautauMass_mva7", "tautau Mass in the signal region after BL and MVA7", 320, 0, 320);

  new TH1F("tautauMass_fakeCR_mva1", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH1F("tautauMass_fakeCR_mva2", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH1F("tautauMass_fakeCR_mva3", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH1F("tautauMass_fakeCR_mva4", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH1F("tautauMass_fakeCR_mva5", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH1F("tautauMass_fakeCR_mva6", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH1F("tautauMass_fakeCR_mva7", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);

  new TH2F("tautauMass_fakeCR_tau2fake_mva1", "tautau Mass in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 320, 0, 320);
  new TH2F("tautauMass_fakeCR_tau2fake_mva2", "tautau Mass in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 320, 0, 320);
  new TH2F("tautauMass_fakeCR_tau2fake_mva3", "tautau Mass in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 320, 0, 320);
  new TH2F("tautauMass_fakeCR_tau2fake_mva4", "tautau Mass in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 320, 0, 320);
  new TH2F("tautauMass_fakeCR_tau2fake_mva5", "tautau Mass in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 320, 0, 320);
  new TH2F("tautauMass_fakeCR_tau2fake_mva6", "tautau Mass in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 320, 0, 320);
  new TH2F("tautauMass_fakeCR_tau2fake_mva7", "tautau Mass in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 320, 0, 320);

  new TH2F("mvaoutput_fakeCR_tau2fake_mva1", "mvaoutput of IrrMVA in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 100, -2, 2);
  new TH2F("mvaoutput_fakeCR_tau2fake_mva2", "mvaoutput of IrrMVA in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 100, -2, 2);
  new TH2F("mvaoutput_fakeCR_tau2fake_mva3", "mvaoutput of IrrMVA in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 100, -2, 2);
  new TH2F("mvaoutput_fakeCR_tau2fake_mva4", "mvaoutput of IrrMVA in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 100, -2, 2);
  new TH2F("mvaoutput_fakeCR_tau2fake_mva5", "mvaoutput of IrrMVA in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 100, -2, 2);
  new TH2F("mvaoutput_fakeCR_tau2fake_mva6", "mvaoutput of IrrMVA in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 100, -2, 2);
  new TH2F("mvaoutput_fakeCR_tau2fake_mva7", "mvaoutput of IrrMVA in the SS-AntiIsolated control region when Tau2 is fake", 140, 0, 140, 100, -2, 2);

  new TH1D("EventCount", "no of events after all cut, in signal region selection", 2, 0.5, 2.5);
  new TH1D("EventCount_mva1", "no of events after all cut, in signal region selection", 2, 0.5, 2.5);
  new TH1D("EventCount_mva2", "no of events after all cut, in signal region selection", 2, 0.5, 2.5);
  new TH1D("EventCount_mva3", "no of events after all cut, in signal region selection", 2, 0.5, 2.5);
  new TH1D("EventCount_mva4", "no of events after all cut, in signal region selection", 2, 0.5, 2.5);
  new TH1D("EventCount_mva5", "no of events after all cut, in signal region selection", 2, 0.5, 2.5);
  new TH1D("EventCount_mva6", "no of events after all cut, in signal region selection", 2, 0.5, 2.5);
  new TH1D("EventCount_mva7", "no of events after all cut, in signal region selection", 2, 0.5, 2.5);
  
  new TH1F("test", "test hist", 9, -0.5, 8.5);

  new TH1F("isoTau", "tauPt-genPt//genPt, for Isolated Taus", 100, -1, 1);
  new TH1F("antiTau", "tauPt-genPt//genPt, for AntiIsolated Taus", 100, -1, 1);
  new TH1F("antiJet", "jetPt-genPt//genPt, for AntiIsolated Taus", 100, -1, 1);
  new TH1F("tvsj", "tau-jet Pt in Gen Level", 100, -50, 50);

  // Histograms for fake function for MuMuTau
  new TH1F("JetToMuFake_WJets_mmt_deno_dummy", "Mu Pt distribution for denominator", 140, 0, 140);
  new TH1F("JetToMuFake_WJets_mmt_deno", "Mu Pt distribution for denominator", 140, 0, 140);
  new TH1F("JetToMuFake_WJets_mmt_nume_lead", "Mu Pt distribution for lead condition", 140, 0, 140);
  new TH1F("JetToMuFake_WJets_mmt_nume_sublead", "Mu Pt distribution for sublead condition", 140, 0, 140);

  new TH1F("JetToMuFake_QCD_mmt_deno", "Mu Pt distribution for denominator", 140, 0, 140);
  new TH1F("JetToMuFake_QCD_mmt_nume_lead", "Mu Pt distribution for lead condition", 140, 0, 140);
  new TH1F("JetToMuFake_QCD_mmt_nume_sublead", "Mu Pt distribution for sublead condition", 140, 0, 140);

  new TH1F("JetToMuFake_ZJets_mmt_deno", "Mu Pt distribution for denominator", 140, 0, 140);
  new TH1F("JetToMuFake_ZJets_mmt_nume_lead", "Mu Pt distribution for lead condition", 140, 0, 140);
  new TH1F("JetToMuFake_ZJets_mmt_nume_sublead", "Mu Pt distribution for sublead condition", 140, 0, 140);

  new TH1F("TagMet_MT_emt_WJets", "MT for WJets selection", 140, 0, 140);
  new TH1F("TagpfMet_MT_mmt_WJets", "MT for WJets selection", 140, 0, 140);
  new TH1F("TagmvaMet_MT_mmt_WJets", "MT for WJets selection", 140, 0, 140);
  new TH1F("Tag_ZMass_mmt_ZJets", "Z Mass for ZJets selection", 140, 0, 140);
  new TH1F("testPt", "pt of the ZJets probeMuon denominator", 1000, 0, 10000);

  //Apply Fake Rate in MMT
  new TH1F("mtMass_mmt_signal", "MuTau Mass (Higgs) in signal region", 320, 0, 320);
  new TH1F("mu1Pt_mmt_signal", "Lead Muon Pt in signal region", 140, 0, 140);
  new TH1F("mu1Eta_mmt_signal", "Lead Muon Eta in signal region", 50, -2.5, 2.5);
  new TH1F("mu1Phi_mmt_signal", "Lead Muon Phi in signal region", 140, -3.5, 3.5);
  new TH1F("mu2Pt_mmt_signal", "subLead Muon Pt in signal region", 140, 0, 140);
  new TH1F("mu2Eta_mmt_signal", "subLead Muon Eta in signal region", 50, -2.5, 2.5);
  new TH1F("mu2Phi_mmt_signal", "subLead Muon Phi in signal region", 140, -3.5, 3.5);
  new TH1F("tauPt_mmt_signal", "Tau Pt in signal region", 140, 0, 140);
  new TH1F("tauEta_mmt_signal", "Tau Eta in signal region", 50, -2.5, 2.5);
  new TH1F("tauPhi_mmt_signal", "Tau Phi in signal region", 140, -3.5, 3.5);

  new TH1F("mu1Pt_leadFakeCR", "pt of lead muon in lead fake region", 140, 0, 140);
  new TH2F("mtMass_mmt_leadFakeCR_p", "Pt of Lead muon in Lead Fake region", 320, 0, 320, 140, 0, 140);
  new TH2F("mu1Pt_mmt_leadFakeCR_p", "Pt of Lead muon in Lead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("mu1Eta_mmt_leadFakeCR_p", "Eta of Lead muon in Lead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("mu1Phi_mmt_leadFakeCR_p", "Phi of Lead muon in Lead Fake region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("mu2Pt_mmt_leadFakeCR_p", "Pt of subLead muon in Lead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("mu2Eta_mmt_leadFakeCR_p", "Eta of subLead muon in Lead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("mu2Phi_mmt_leadFakeCR_p", "Phi of subLead muon in Lead Fake region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("tauPt_mmt_leadFakeCR_p", "Pt of Tau in Lead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("tauEta_mmt_leadFakeCR_p", "Eta of Tau in Lead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("tauPhi_mmt_leadFakeCR_p", "Phi of Tau in Lead Fake region", 140, -3.5, 3.5, 140, 0, 140);

  new TH1F("mu2Pt_slFakeCR", "pt of sub lead muon in sub lead fake region", 140, 0, 140);
  new TH2F("mtMass_mmt_slFakeCR_p", "Pt of Lead muon in subLead Fake region", 320, 0, 320, 140, 0, 140);
  new TH2F("mu1Pt_mmt_slFakeCR_p", "Pt of Lead muon in subLead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("mu1Eta_mmt_slFakeCR_p", "Eta of Lead muon in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("mu1Phi_mmt_slFakeCR_p", "Phi of Lead muon in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("mu2Pt_mmt_slFakeCR_p", "Pt of subLead muon in subLead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("mu2Eta_mmt_slFakeCR_p", "Eta of subLead muon in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("mu2Phi_mmt_slFakeCR_p", "Phi of subLead muon in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("tauPt_mmt_slFakeCR_p", "Pt of Tau in subLead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("tauEta_mmt_slFakeCR_p", "Eta of Tau in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("tauPhi_mmt_slFakeCR_p", "Phi of Tau in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140);

  new TH1F("mu1Pt_bothFakeCR", "pt of lead muon in both fake region", 140, 0, 140);
  new TH1F("mu2Pt_bothFakeCR", "pt of sublead muon in both fake region", 140, 0, 140);
  new TH3F("mtMass_mmt_bothFakeCR_p", "Pt of Lead muon in subLead Fake region", 320, 0, 320, 140, 0, 140, 140, 0, 140);
  new TH3F("mu1Pt_mmt_bothFakeCR_p", "Pt of Lead muon in subLead Fake region", 140, 0, 140, 140, 0, 140, 140, 0, 140);
  new TH3F("mu1Eta_mmt_bothFakeCR_p", "Eta of Lead muon in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140, 140, 0, 140);
  new TH3F("mu1Phi_mmt_bothFakeCR_p", "Phi of Lead muon in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140, 140, 0, 140);
  new TH3F("mu2Pt_mmt_bothFakeCR_p", "Pt of subLead muon in subLead Fake region", 140, 0, 140, 140, 0, 140, 140, 0, 140);
  new TH3F("mu2Eta_mmt_bothFakeCR_p", "Eta of subLead muon in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140, 140, 0, 140);
  new TH3F("mu2Phi_mmt_bothFakeCR_p", "Phi of subLead muon in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140, 140, 0, 140);
  new TH3F("tauPt_mmt_bothFakeCR_p", "Pt of Tau in subLead Fake region", 140, 0, 140, 140, 0, 140, 140, 0, 140);
  new TH3F("tauEta_mmt_bothFakeCR_p", "Eta of Tau in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140, 140, 0, 140);
  new TH3F("tauPhi_mmt_bothFakeCR_p", "Phi of Tau in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140, 140, 0, 140);


  new TH1F("mtMass_mmt_F3", "MuTau Mass (Higgs) in F3 region", 320, 0, 320);
  new TH1F("mu1Pt_mmt_F3", "Lead Muon Pt in F3 region", 140, 0, 140);
  new TH1F("mu1Eta_mmt_F3", "Lead Muon Eta in F3 region", 50, -2.5, 2.5);
  new TH1F("mu1Phi_mmt_F3", "Lead Muon Phi in F3 region", 140, -3.5, 3.5);
  new TH1F("mu2Pt_mmt_F3", "subLead Muon Pt in F3 region", 140, 0, 140);
  new TH1F("mu2Eta_mmt_F3", "subLead Muon Eta in F3 region", 50, -2.5, 2.5);
  new TH1F("mu2Phi_mmt_F3", "subLead Muon Phi in F3 region", 140, -3.5, 3.5);
  new TH1F("tauPt_mmt_F3", "Tau Pt in F3 region", 140, 0, 140);
  new TH1F("tauEta_mmt_F3", "Tau Eta in F3 region", 50, -2.5, 2.5);
  new TH1F("tauPhi_mmt_F3", "Tau Phi in F3 region", 140, -3.5, 3.5);

  new TH2F("mtMass_mmt_leadFakeCR_F3_p", "Pt of Lead muon in Lead Fake region", 320, 0, 320, 140, 0, 140);
  new TH2F("mu1Pt_mmt_leadFakeCR_F3_p", "Pt of Lead muon in Lead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("mu1Eta_mmt_leadFakeCR_F3_p", "Eta of Lead muon in Lead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("mu1Phi_mmt_leadFakeCR_F3_p", "Phi of Lead muon in Lead Fake region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("mu2Pt_mmt_leadFakeCR_F3_p", "Pt of subLead muon in Lead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("mu2Eta_mmt_leadFakeCR_F3_p", "Eta of subLead muon in Lead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("mu2Phi_mmt_leadFakeCR_F3_p", "Phi of subLead muon in Lead Fake region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("tauPt_mmt_leadFakeCR_F3_p", "Pt of Tau in Lead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("tauEta_mmt_leadFakeCR_F3_p", "Eta of Tau in Lead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("tauPhi_mmt_leadFakeCR_F3_p", "Phi of Tau in Lead Fake region", 140, -3.5, 3.5, 140, 0, 140);

  new TH2F("mtMass_mmt_slFakeCR_F3_p", "Pt of Lead muon in subLead Fake region", 320, 0, 320, 140, 0, 140);
  new TH2F("mu1Pt_mmt_slFakeCR_F3_p", "Pt of Lead muon in subLead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("mu1Eta_mmt_slFakeCR_F3_p", "Eta of Lead muon in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("mu1Phi_mmt_slFakeCR_F3_p", "Phi of Lead muon in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("mu2Pt_mmt_slFakeCR_F3_p", "Pt of subLead muon in subLead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("mu2Eta_mmt_slFakeCR_F3_p", "Eta of subLead muon in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("mu2Phi_mmt_slFakeCR_F3_p", "Phi of subLead muon in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140);
  new TH2F("tauPt_mmt_slFakeCR_F3_p", "Pt of Tau in subLead Fake region", 140, 0, 140, 140, 0, 140);
  new TH2F("tauEta_mmt_slFakeCR_F3_p", "Eta of Tau in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140);
  new TH2F("tauPhi_mmt_slFakeCR_F3_p", "Phi of Tau in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140);

  new TH3F("mtMass_mmt_bothFakeCR_F3_p", "Pt of Lead muon in subLead Fake region", 320, 0, 320, 140, 0, 140, 140, 0, 140);
  new TH3F("mu1Pt_mmt_bothFakeCR_F3_p", "Pt of Lead muon in subLead Fake region", 140, 0, 140, 140, 0, 140, 140, 0, 140);
  new TH3F("mu1Eta_mmt_bothFakeCR_F3_p", "Eta of Lead muon in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140, 140, 0, 140);
  new TH3F("mu1Phi_mmt_bothFakeCR_F3_p", "Phi of Lead muon in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140, 140, 0, 140);
  new TH3F("mu2Pt_mmt_bothFakeCR_F3_p", "Pt of subLead muon in subLead Fake region", 140, 0, 140, 140, 0, 140, 140, 0, 140);
  new TH3F("mu2Eta_mmt_bothFakeCR_F3_p", "Eta of subLead muon in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140, 140, 0, 140);
  new TH3F("mu2Phi_mmt_bothFakeCR_F3_p", "Phi of subLead muon in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140, 140, 0, 140);
  new TH3F("tauPt_mmt_bothFakeCR_F3_p", "Pt of Tau in subLead Fake region", 140, 0, 140, 140, 0, 140, 140, 0, 140);
  new TH3F("tauEta_mmt_bothFakeCR_F3_p", "Eta of Tau in subLead Fake region", 50, -2.5, 2.5, 140, 0, 140, 140, 0, 140);
  new TH3F("tauPhi_mmt_bothFakeCR_F3_p", "Phi of Tau in subLead Fake region", 140, -3.5, 3.5, 140, 0, 140, 140, 0, 140);

  //new TH1F("mtMass_mmt_F3", "MuTau Mass (Higgs) in F3 region", 320, 0, 320);
  new TH1F("mu1Pt_leadFakeCR_F3", "pt of lead muon in lead fake region F3", 140, 0, 140);
  new TH1F("mu2Pt_slFakeCR_F3", "pt of sub lead muon in sub lead fake region F3", 140, 0, 140);
  new TH1F("mu1Pt_bothFakeCR_F3", "pt of lead muon in both fake region F3", 140, 0, 140);
  new TH1F("mu2Pt_bothFakeCR_F3", "pt of sublead muon in both fake region F3", 140, 0, 140);

  new TH1F("JetToMuFake_WJets_emt_deno", "Mu Pt distribution for denominator", 140, 0, 140);
  new TH1F("JetToMuFake_WJets_emt_nume", "Mu Pt distribution for numerator", 140, 0, 140);
  new TH1F("JetToMuFake_ZJets_emt_deno", "Mu Pt distribution for denominator", 140, 0, 140);
  new TH1F("JetToMuFake_ZJets_emt_nume", "Mu Pt distribution for lead condition", 140, 0, 140);

  new TH2F("mu1Pt_mu2Pt_bothFakeCR_F3", "Mu1 Pt vs Mu2 Pt in bothfakeCR F3", 140, 0, 140, 140, 0, 140);
  new TH2F("mu1Pt_mu2Pt_bothFakeCR", "Mu1 Pt vs Mu2 Pt in bothfake CR", 140, 0, 140, 140, 0, 140);
}
// -------------------
// The main event loop
// -------------------
void FakeRate::clearLists() {
  vtxList.clear();
  muoList.clear();
  eleList.clear();
  tauList.clear();
  bjetList.clear();

  selectTauList.clear();
  selectTightMuonList.clear();
  selectLooseMuonList.clear();
  selectTightElectronList.clear();
}
void FakeRate::eventLoop() 
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
  // -------------------
int noofEvents= nEvents();   

 string lastFile;
  for (int ev = 0; ev < nEvents(); ++ev) {
    clearEvent();
    clearLists();
 
    int lflag = chain()->LoadTree(ev); 
    int nbytes = getEntry(lflag);    // returns total bytes read

    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName())); 

    const Event& evt = eventColl()->at(0);

    histf()->cd();

    //Defined here since I may also need nPV for PU Reweight
    findVtxInfo(vtxList, op, fLog());
    int nvtx = vtxList.size();
   

    // PileUP weight
    puevWt_ = 1; // for data
    if (isMC()) {
      //If need	nPV for	Reweight, pass it here instead of npu=0    
      int npu = 0;
      puevWt_ = wtPileUp(npu);
      
      histf()->cd();
      //cout<<"pile up weight "<<puevWt_<<endl;
      //cout<<"npu in the Fakerate "<<npu<<endl;     
      AnaUtil::fillHist1D("npu", npu, 1.0);
      AnaUtil::fillHist1D("puwt", puevWt_, 1.0);


      const GenEvent& genevt = genEventColl()->at(0);
      eventWt_ = genevt.weightevt;
      AnaUtil::fillHist1D("eventWeight", eventWt_, puevWt_);
      AnaUtil::fillHist1D("eventWeight_initial", 1.0, eventWt_*puevWt_);
      //keep it closed, not needed LO 
      //puevWt_ *= eventWt_;

      vector<int> list = evt.trueNInt;
      AnaUtil::fillHist1D("trueNInt", (list.size() ? list.at(0) : 0), 1.0);


    }

    AnaUtil::fillHist1D("getnEvents", noofEvents, puevWt_ );     
 
    AnaUtil::fillHist1D("evcount", 1, puevWt_);

    AnaUtil::fillHist1D("npv_PU", nvtx, puevWt_);
    AnaUtil::fillHist1D("npv", nvtx, 1);
       
    AnaUtil::fillHist1D("EventCount", 1.0, puevWt_);    
    AnaUtil::fillHist1D("EventCount_mva1", 1.0, puevWt_);    
    AnaUtil::fillHist1D("EventCount_mva2", 1.0, puevWt_);    
    AnaUtil::fillHist1D("EventCount_mva3", 1.0, puevWt_);    
    AnaUtil::fillHist1D("EventCount_mva4", 1.0, puevWt_);    
    AnaUtil::fillHist1D("EventCount_mva5", 1.0, puevWt_);    
    AnaUtil::fillHist1D("EventCount_mva6", 1.0, puevWt_);    

    AnaUtil::fillHist1D("test", 0, puevWt_);    
    
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

    // Trigger selection
    if (useTrigger() && !isTriggered()) continue;

     AnaUtil::fillHist1D("evcount", 2, puevWt_);

    AnaUtil::fillHist1D("nPVCorrected", nvtx, puevWt_);

    op.verbose = (logOption() >> 2 & 0x1); 
    findVtxInfo(vtxList, op, fLog());

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
    if (vtxList.size() < 1) continue;
    vz = ( nvtx > 0 ) ? vtxList.at(0).z : -999;

    AnaUtil::fillHist1D("test", 1.0, puevWt_);    

    op.verbose = (logOption() >> 3 & 0x1);
    //findMuonInfo(muoList, vz, op, fLog());


     AnaUtil::fillHist1D("evcount", 3, puevWt_);

    selectEvent();
  }  
  // Analysis is over
  endJob();
}
void FakeRate::selectEvent() {
  //Object Selection: will be used in later modules
  selectTau();
  selectTightMuon();
  selectTightElectron();
  selectLooseMuon();

  //investigation purposes
  //investigateFake();

  //Measure Fake Rate for Fully Hadronic Channel
  JetToFakeDYJets();                       // Here We Calculate Jet-> Mu, Tau fake rate in DYToMuMu + Jets Control region
  WjetsCR();                               // Here we want to calculate Mt and other variables in Wjets CR 
  JetToTauFakeWJets();                     // Here We Calculate Jet-> Tau fake function in W+Jets Control Region 
  JetToTauFakeTTJets_2j();                     // Here We Calculate Jet-> Tau fake function in W+Jets Control Region 
  JetToTauFakeTTJets_2j_ss();                     // Here We Calculate Jet-> Tau fake function in W+Jets Control Region 
  JetToTauFakeTTJets_1j();                     // Here We Calculate Jet-> Tau fake function in W+Jets Control Region 
  JetToTauFakeQCD();                     // Here We Calculate Jet-> Tau fake function in W+Jets Control Region 
  ApplyFR();

  //Measure Fake Rate and Apply Fake Rate Modules for Semi-Leptonic Channel
  MuMuTauApplyFR();
  MuMuTauFakeQCD();
  MuMuTauFakeZJets();
  //EleMuTauFakeZJets();
  //EleMuTauFakeWJets();

  //Keep it in last, for kNN MVA skim filling, since the problem of switching between skim hist and normal hist
  MuMuTauFakeWJets();
}

void FakeRate::selectTau()                
{
  for (auto it = tauColl()->begin(); it != tauColl()->end(); ++it) {
    const Tau& tau = (*it);
    if ((fabs(tau.eta) >= AnaUtil::cutValue(_tau1CutMap, "eta")) ||
        (tau.decayModeFinding <= 0.5) ||
        (tau.againstMuonTight3 <= 0.5) ||
        (tau.againstElectronTightMVA6 <= 0.5) ||
        (tau.byVLooseIsolationMVArun2v1DBoldDMwLT <= 0.5)
       ) continue; 

    selectTauList.push_back(tau);
  }
  if (selectTauList.size() > 1)
    sort(selectTauList.begin(), selectTauList.end(), PtComparator<Tau>());
}



//Without Pt and Isolation
void FakeRate::selectTightMuon()
{
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    if (fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
#if 0
        (!muon.isTrackerMuon) ||
        (!muon.isGlobalMuonPromptTight) ||
        (muon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) ||
        (muon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches")) ||
        (muon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) ||
        (muon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) ||
        (muon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) ||
        (muon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) ||
        (abs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0"))
#endif
        (!muon.isTightMuon)
       ) continue;

    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || fabs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;

    selectTightMuonList.push_back(muon);
  }
  if (selectTightMuonList.size() > 1)
    sort(selectTightMuonList.begin(), selectTightMuonList.end(), PtComparator<Muon>());
}



void FakeRate::selectLooseMuon()
{
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    if (fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
        (!muon.isTrackerMuon) ||
        (!muon.isGlobalMuonPromptTight) ||
        (muon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) 
       ) continue;

    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;

    selectLooseMuonList.push_back(muon);
  }
  if (selectLooseMuonList.size() > 1)
    sort(selectLooseMuonList.begin(), selectLooseMuonList.end(), PtComparator<Muon>());
}


//With Pt and Isolation
void FakeRate::selectTightElectron()
{
  for (auto it = electronColl()->begin(); it != electronColl()->end(); ++it) {
    const Electron& ele = (*it);
    if (!TightEleId(ele)) continue;
    if (ele.pt <= 20.0) continue;
    selectTightElectronList.push_back(ele);
  }
  if (selectTightElectronList.size() > 1)
    sort(selectTightElectronList.begin(), selectTightElectronList.end(), PtComparator<Electron>());
}



void FakeRate::investigateFake()                
{
#if 0
  vector<Tau> list;
  for (auto it = tauColl()->begin(); it != tauColl()->end(); ++it) {
    const Tau& tau = (*it);
    if ((abs(tau.eta) >= AnaUtil::cutValue(_tau1CutMap, "eta")) ||
        (tau.decayModeFinding <= 0.5) ||
        (tau.againstMuonTight <= 0.5) ||
        (tau.againstElectronTight <= 0.5) ||
        (tau.pt < 30)
       ) continue; 

    list.push_back(tau);
  }
  if (list.size() > 1)
    sort(list.begin(), list.end(), PtComparator<Tau>());

  cout << "list size = " << list.size() 
       << ", genparticle = " << ngenparticle() 
       << endl;
  if (!list.size() || ngenparticle()) return;

  for (auto it = list.begin(); it != list.end(); ++it) {
    const Tau& tau = (*it);
    TLorentzVector ft;
    ft.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
    vector<GenParticle> candTau, candJet;
    for (auto jt = genParticleColl()->begin(); jt != genParticleColl()->end(); ++jt) {
      const GenParticle& gp = (*jt);
      if (abs(gp.pdgId) == 15) continue;
      TLorentzVector g;
      g.SetPtEtaPhiE (gp.pt, gp.eta, gp.phi, gp.energy);
      if (AnaUtil::deltaR(ft, g) < 0.3)
        candTau.push_back(gp);
      double dR = sqrt(pow(abs(gp.eta - tau.jetEta),2) + pow(AnaUtil::deltaPhi(gp.phi, tau.jetPhi) ,2));
      if (dR < 0.3)
        candJet.push_back(gp);
    }
    if (candTau.size() > 1)
      sort(candTau.begin(), candTau.end(), PtComparator<GenParticle>());

    if (candJet.size() > 1)
      sort(candJet.begin(), candJet.end(), PtComparator<GenParticle>());

    // Isolated Fake Tau
    if (candTau.size() > 0) {
      if (tau.byTightCombinedIsolationDeltaBetaCorr >= 0.5) {
        double isoTau = (tau.pt - candTau.at(0).pt)/candTau.at(0).pt;
        AnaUtil::fillHist1D("isoTau", isoTau, puevWt_);
      }
      else {
        double antiisoTau = (tau.pt - candTau.at(0).pt)/candTau.at(0).pt;
        AnaUtil::fillHist1D("antiTau", antiisoTau, puevWt_);
      }
    }

    // AntiIsolated Fake Jet
    if (candJet.size() > 0 && tau.byTightCombinedIsolationDeltaBetaCorr < 0.5) {
      double antiisoJet = (tau.jetPt - candJet.at(0).pt)/candJet.at(0).pt;
      AnaUtil::fillHist1D("antiJet", antiisoJet, puevWt_);
    }
    if (candTau.size() > 0 && candJet.size() > 0)
      AnaUtil::fillHist1D("tvsj", candTau.at(0).pt - candJet.at(0).pt, puevWt_);
  }
#endif
}
void FakeRate::JetToFakeDYJets()
{
  double vz = vtxList.at(0).z;
  vector<Muon> fMuoList, frMuoList;

  AnaUtil::fillHist1D("evcount", 4, puevWt_);
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    if (fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||

#if 0
        (!muon.isTrackerMuon) || 
        (!muon.isGlobalMuonPromptTight) || 
        (muon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) || 
        (muon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches")) ||
        (muon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) || 
        (muon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) || 
        (muon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) || 
        (muon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) || 
        (abs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0"))
#endif
        (!muon.isTightMuon)
       ) continue;

    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;

    if (fMuoList.size() >= 2)

      frMuoList.push_back(muon); // This is to Catch the 3rd Muon // Jet->Muon Fake

    AnaUtil::fillHist1D("Relative Isolation of muons",fabs(muon.pfRelIsoDB03),puevWt_);

    if (fabs(muon.pfRelIsoDB03) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso") ||
	muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt")) continue;     

   fMuoList.push_back(muon);
  }
  // atleast 2 good muon
  if (fMuoList.size() < 2) return;
  const Muon& muo1 = fMuoList.at(0);
  const Muon& muo2 = fMuoList.at(1);
  
  const MET& mt = metColl()->at(0);

AnaUtil::fillHist1D("evcount", 5, puevWt_);

  //---------------------------------------------
  //
  //              Electron Veto
  //
  //---------------------------------------------
  if (vetoElectron_wPt(15) != 0) return;

AnaUtil::fillHist1D("evcount", 6, puevWt_);


  //---------------------------------------------
  //
  //              b-tagged Jets Veto
  //
  //---------------------------------------------
  if (bjetList.size() > 0) return; 
  
  AnaUtil::fillHist1D("evcount", 7, puevWt_);

if (muo1.pt <= 24 || muo2.pt <= 24) return;               
  
AnaUtil::fillHist1D("evcount", 8, puevWt_);

  TLorentzVector M1, M2,T1,sum;
  M1.SetPtEtaPhiE (muo1.pt, muo1.eta, muo1.phi, muo1.energy);
  M2.SetPtEtaPhiE (muo2.pt, muo2.eta, muo2.phi, muo2.energy);
  sum = M1 + M2;
  double drM1M2 = AnaUtil::deltaR(M1, M2);
  if (AnaUtil::deltaR(M1, M2) < 0.4) return;


AnaUtil::fillHist1D("evcount", 9, puevWt_);


  if ((muo1.charge * muo2.charge) > 0) return;

AnaUtil::fillHist1D("evcount", 10, puevWt_);

  if (sum.M() < 60 ) return;


AnaUtil::fillHist1D("evcount", 11, puevWt_);


  AnaUtil::fillHist1D("met_DYjets", mt.met, puevWt_);
  AnaUtil::fillHist1D("ZMass_mu", sum.M(), puevWt_);  
  AnaUtil::fillHist1D("pt_mu1_DY_CR", muo1.pt, puevWt_);
  AnaUtil::fillHist1D("pt_mu2_DY_CR", muo2.pt, puevWt_);
  AnaUtil::fillHist1D("mu1_eta", muo1.eta, puevWt_);
  AnaUtil::fillHist1D("mu1_phi", muo1.phi, puevWt_);
  AnaUtil::fillHist1D("Delta_R_M1_M2", drM1M2, puevWt_);
  int nvtx = vtxList.size();
  AnaUtil::fillHist1D("npv_DY", nvtx, puevWt_);




















if (sum.M() < 70 || sum.M() > 110) return;
//const MET& mt = metColl()->at(0);

  ///////////////////////////////////////////////////// Muon Fake Function for DYJets Control Region ////////////////////////////////////////////////////////////

  if (frMuoList.size() >= 1) {
    for (unsigned int indx = 0; indx < frMuoList.size(); ++indx) {
      const Muon& fake = frMuoList[indx];
      double mass1 = sqrt(2*fake.pt*mt.met*(1-cos(AnaUtil::deltaPhi(fake.phi, mt.metphi))));
      if (mass1 >= 20) continue;                                                          // To Remove WZ Contamination

      AnaUtil::fillHist1D("fakemu_deno", fake.pt, puevWt_);
      if (fabs(fake.pfRelIsoDB03) < 0.15)
        AnaUtil::fillHist1D("fakemu_nume", fake.pt, puevWt_);
    }
  }
  ////////////////////////////////////////////////////// Jet to Tau Fake Rate for DYJets Control Region///////////////////////////////////////////////////////
  if (!selectTauList.size()) return;      

  int closureCountDY = 0;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& fake = selectTauList[indx];
    T1.SetPtEtaPhiE(fake.pt, fake.eta, fake.phi, fake.energy);
    double dr1 = AnaUtil::deltaR(M1, T1);
    double dr2 = AnaUtil::deltaR(M2, T1);
    if ((dr1 < AnaUtil::cutValue(_tau1CutMap, "drMuTau")) || 
        (dr2 < AnaUtil::cutValue(_tau1CutMap, "drMuTau"))) continue;
  
    //double mass1 = sqrt(2*fake.pt*mt->met*(1-cos(AnaUtil::deltaPhi(fake.phi, mt->metphi))));
    //if (mass1 >= 20) continue;                                                                          // To Remove WZ Contamination
    if (fabs(fake.eta) < 0.8) { 
      AnaUtil::fillHist1D("faketaupt_deno_DYC", fake.pt, puevWt_);
      if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
        AnaUtil::fillHist1D("faketaupt_nume_DYC", fake.pt, puevWt_);
    }
    if (fabs(fake.eta) >= 0.8 && fabs(fake.eta) < 1.6) { 
      AnaUtil::fillHist1D("faketaupt_deno_DYI", fake.pt, puevWt_);
      if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
        AnaUtil::fillHist1D("faketaupt_nume_DYI", fake.pt, puevWt_);
    }
    if (fabs(fake.eta) >= 1.6) { 
      AnaUtil::fillHist1D("faketaupt_deno_DYF", fake.pt, puevWt_);
      if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
        AnaUtil::fillHist1D("faketaupt_nume_DYF", fake.pt, puevWt_);
    }
    if (fabs(fake.eta) <= 2.3) {
      closureCountDY += 1.0; 
      AnaUtil::fillHist1D("faketaupt_deno_DYAll", fake.pt, puevWt_);
      if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
        AnaUtil::fillHist1D("faketaupt_nume_DYAll", fake.pt, puevWt_);
      if (closureCountDY == 1.0){
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5) AnaUtil::fillHist1D("faketaupt_Isolated_DYAll", fake.pt, puevWt_);
        else AnaUtil::fillHist1D("faketaupt_AntiIsolated_DYAll", fake.pt, puevWt_);
      }
    }
  }
}

void FakeRate::WjetsCR()
{
  double vz = vtxList.at(0).z;
  vector<Muon> fMuoList;

for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);
    if (fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") || (fabs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0")) || (!muon.isTightMuon)) 
    continue;
 bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;

    if (fabs(muon.pfRelIsoDB03) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso") ||
        muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt")) continue;

    fMuoList.push_back(muon);
  }
if (fMuoList.size() != 1) return; //Only one good Muon
  const Muon& muo1 = fMuoList.at(0);
  TLorentzVector M1;
  M1.SetPtEtaPhiE (muo1.pt, muo1.eta, muo1.phi, muo1.energy);
  const MET& mt = metColl()->at(0);
  double mass1 = sqrt(2*muo1.pt*mt.met*(1-cos(AnaUtil::deltaPhi(muo1.phi, mt.metphi))));
 
  AnaUtil::fillHist1D("MT_Mu_Met_Wjet_before_tau_sel", mass1, puevWt_);

if (mass1 <= 40) return;                                                          
 if (vetoElectron_wPt(15) != 0) return;
 if (bjetList.size() > 0) return;
 if (vetoMuon_wPt(10) != 1) return;

  TLorentzVector T1;

  if (selectTauList.size() < 2) return;      // No Tau1 means no Tau2, so we can safely reject the event
  vector<Tau> tau1list;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& fake = selectTauList[indx]; 
    T1.SetPtEtaPhiE(fake.pt, fake.eta, fake.phi, fake.energy);
    double dr1 = AnaUtil::deltaR(M1, T1);
    if (dr1 < AnaUtil::cutValue(_tau1CutMap, "drMuTau")) continue; 
    tau1list.push_back(fake);
  } 

  if (tau1list.size() < 2) return;
  const Tau& ta1 = tau1list.at(0);
  const Tau& ta2 = tau1list.at(1);
if ((ta1.charge == ta2.charge) && (muo1.charge == ta1.charge)) { 

AnaUtil::fillHist1D("MT_Mu_Met_Wjet_after_tau_sel", mass1, puevWt_);
int nvtx = vtxList.size();
AnaUtil::fillHist1D("npv_Wjets", nvtx, puevWt_);
  }
}







void FakeRate::JetToTauFakeWJets()
{
  double vz = vtxList.at(0).z;

  // atleast 1 good muon
  //if (muoList.size() < 1) return;
  vector<Muon> fMuoList;
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    if (fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
#if 0
        (!muon.isTrackerMuon) ||
        (!muon.isGlobalMuonPromptTight) ||
        (muon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) ||
        (muon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches")) ||
        (muon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) ||
        (muon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) ||
        (muon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) ||
	(muon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) ||

#endif
        (fabs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0")) ||
        (!muon.isTightMuon)
       ) continue;

    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;

    if (fabs(muon.pfRelIsoDB03) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso") ||
        muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt")) continue;

    fMuoList.push_back(muon);
  }
  if (fMuoList.size() != 1) return; //Only one good Muon
  const Muon& muo1 = fMuoList.at(0);
  TLorentzVector M1;
  M1.SetPtEtaPhiE (muo1.pt, muo1.eta, muo1.phi, muo1.energy);
  const MET& mt = metColl()->at(0);
  double mass1 = sqrt(2*muo1.pt*mt.met*(1-cos(AnaUtil::deltaPhi(muo1.phi, mt.metphi))));

  if (mass1 <= 40) return;                                                          
  AnaUtil::fillHist1D("met_Wjets", mt.met, puevWt_);
  //---------------------------------------------
  //
  //              Electron Veto
  //
  //---------------------------------------------
  if (vetoElectron_wPt(15) != 0) return;

  //---------------------------------------------
  //
  //              b-tagged Jets Veto
  //
  //---------------------------------------------
  if (bjetList.size() > 0) return; 

/*  ------ For Jet---
  //Require 1 Jet to mimic the presence of a Tau object
  vector<Jet> JetList;
  float nJet = 0;
  TLorentzVector J1,J2;
  bool JetFlag = false;
  for (auto it = jetColl()->begin(); it != jetColl()->end(); ++it) {
    const Jet& jt = (*it);
    if (fabs(jt.eta) >= 2.4) continue;
    nJet += 1;
    if (JetList.size() >= 2)
      JetList.push_back(jt);
}
if (JetList.size() < 2) return;
  const Jet& jt1 = JetList.at(0);
  const Jet& jt2 = JetList.at(1);

    J1.SetPtEtaPhiE(jt1.pt, jt1.eta, jt1.phi, jt1.energy);
    J2.SetPtEtaPhiE(jt2.pt, jt2.eta, jt2.phi, jt2.energy);
    if (AnaUtil::deltaR(M1, J1) < 0.5 || AnaUtil::deltaR(M1, J2) < 0.5) continue;
    JetFlag = true;
    //break; //just need 1 non-overlapping jet
  if (!JetFlag) return;

 --------------
*/



  if (vetoMuon_wPt(10) != 1) return;
  /*
  int nm = 0;

  if (nmuon() < 1) return;
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    if (muon.pt <= 15 ||
        fabs(muon.eta) >= 2.1 ||
       	!muon.isTrackerMuon ||
        !muon.isGlobalMuonPromptTight) continue;
    ++nm;
  }
  if (nm != 1) return;
  */

  TLorentzVector T1;

  if (selectTauList.size() < 2) return;      // No Tau1 means no Tau2, so we can safely reject the event
  vector<Tau> tau1list;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& fake = selectTauList[indx];
    T1.SetPtEtaPhiE(fake.pt, fake.eta, fake.phi, fake.energy);
    double dr1 = AnaUtil::deltaR(M1, T1);
    if (dr1 < AnaUtil::cutValue(_tau1CutMap, "drMuTau")) continue; 
    tau1list.push_back(fake);
  } 

  if (tau1list.size() < 2) return;
  const Tau& ta1 = tau1list.at(0);
  const Tau& ta2 = tau1list.at(1);

  //W and TauTau charge correlation
  if ((ta1.charge == ta2.charge) && (muo1.charge == ta1.charge)) { 
    for (unsigned int indx = 0; indx < 2; ++indx) {
      const Tau& fake = tau1list.at(indx);
      if (fabs(fake.eta) < 0.8) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsC_SSSS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsC_SSSS", fake.pt, puevWt_);
      }

      if (fabs(fake.eta) >= 0.8 && fabs(fake.eta) < 1.6) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsI_SSSS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsI_SSSS", fake.pt, puevWt_);
      }

      if (fabs(fake.eta) >= 1.6 && fabs(fake.eta) < 2.3) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsF_SSSS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsF_SSSS", fake.pt, puevWt_);
      }
      if (fabs(fake.eta) < 2.3) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsAll_SSSS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5){
          AnaUtil::fillHist1D("faketaupt_nume_WJetsAll_SSSS", fake.pt, puevWt_);
          AnaUtil::fillHist1D("faketaupt_Isolated_WJetsAll", fake.pt, puevWt_);
        }
        else AnaUtil::fillHist1D("faketaupt_AntiIsolated_WJetsAll", fake.pt, puevWt_);
      }
    }
  }

  //W and TauTau charge correlation
  if (ta1.charge != ta2.charge) { 
    for (unsigned int indx = 0; indx < 2; ++indx) {
      const Tau& fake = tau1list.at(indx);
      if ((muo1.charge + fake.charge) == 0) continue; //Ensure W to SS 
      if (fabs(fake.eta) < 0.8) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsC_OSSS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsC_OSSS", fake.pt, puevWt_);
      }

      if (fabs(fake.eta) >= 0.8 && fabs(fake.eta) < 1.6) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsI_OSSS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsI_OSSS", fake.pt, puevWt_);
      }

      if (fabs(fake.eta) >= 1.6 && fabs(fake.eta) < 2.3) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsF_OSSS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsF_OSSS", fake.pt, puevWt_);
      }
      if (fabs(fake.eta) < 2.3) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsAll_OSSS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsAll_OSSS", fake.pt, puevWt_);
      }
    }
  }

  //W and TauTau charge correlation
  if ((ta1.charge == ta2.charge) && (muo1.charge != ta1.charge)) { 
    for (unsigned int indx = 0; indx < 2; ++indx) {
      const Tau& fake = tau1list.at(indx);
      if (fabs(fake.eta) < 0.8) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsC_SSOS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsC_SSOS", fake.pt, puevWt_);
      }

      if (fabs(fake.eta) >= 0.8 && fabs(fake.eta) < 1.6) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsI_SSOS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsI_SSOS", fake.pt, puevWt_);
      }

      if (fabs(fake.eta) >= 1.6 && fabs(fake.eta) < 2.3) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsF_SSOS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsF_SSOS", fake.pt, puevWt_);
      }
      if (fabs(fake.eta) < 2.3) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsAll_SSOS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsAll_SSOS", fake.pt, puevWt_);
      }
    }
  }

  //W and TauTau charge correlation
  if (ta1.charge == ta2.charge) { 
    for (unsigned int indx = 0; indx < 2; ++indx) {
      const Tau& fake = tau1list.at(indx);
      if (fabs(fake.eta) < 0.8) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsC_NoW", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsC_NoW", fake.pt, puevWt_);
      }

      if (fabs(fake.eta) >= 0.8 && fabs(fake.eta) < 1.6) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsI_NoW", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsI_NoW", fake.pt, puevWt_);
      }

      if (fabs(fake.eta) >= 1.6 && fabs(fake.eta) < 2.3) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsF_NoW", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsF_NoW", fake.pt, puevWt_);
      }
      if (fabs(fake.eta) < 2.3) {
        AnaUtil::fillHist1D("faketaupt_deno_WJetsAll_NoW", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_WJetsAll_NoW", fake.pt, puevWt_);
      }
    }
  }

  //W and TauTau charge correlation
  if ((ta1.charge != ta2.charge) && ta1.byTightIsolationMVArun2v1DBoldDMwLT < 0.5 && ta2.byTightIsolationMVArun2v1DBoldDMwLT < 0.5) { 
       if (muo1.charge == ta1.charge)
         AnaUtil::fillHist1D("faketaupt_antianti_WJetsOSS", ta1.pt, puevWt_);
       else if (muo1.charge == ta2.charge)
         AnaUtil::fillHist1D("faketaupt_antianti_WJetsOSS", ta2.pt, puevWt_);
  }
  if ((ta1.charge == ta2.charge) && (muo1.charge == ta1.charge) && 
  //if ((ta1.charge == ta2.charge) && 
      ta1.byTightIsolationMVArun2v1DBoldDMwLT < 0.5 && 
      ta2.byTightIsolationMVArun2v1DBoldDMwLT < 0.5) { 
        AnaUtil::fillHist1D("faketaupt_antianti_WJetsSSS", ta1.pt, puevWt_);
        AnaUtil::fillHist1D("faketaupt_antianti_WJetsSSS", ta2.pt, puevWt_);
  }
}
void FakeRate::JetToTauFakeTTJets_2j()
{
  double vz = vtxList.at(0).z;

  // atleast 1 good muon
  //if (muoList.size() < 1) return;
  vector<Muon> fMuoList;
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    if (fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
        (fabs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0")) ||
        (!muon.isTightMuon)
       ) continue;

    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;

    if (fabs(muon.pfRelIsoDB03) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso") ||
        muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt")) continue;

    fMuoList.push_back(muon);
  }
  if (fMuoList.size() != 1) return; //only one good Muon
  const Muon& muo1 = fMuoList.at(0);
  TLorentzVector M1;
  M1.SetPtEtaPhiE (muo1.pt, muo1.eta, muo1.phi, muo1.energy);
  const MET& mt = metColl()->at(0);
  double mass1 = sqrt(2*muo1.pt*mt.met*(1-cos(AnaUtil::deltaPhi(muo1.phi, mt.metphi))));
  if (mass1 <= 40.0) return;                                                          
  if (mt.met <= 50.0) return;

  //---------------------------------------------
  //
  //              Electron Veto
  //
  //---------------------------------------------
  if (vetoElectron_wPt(15) != 0) return;

  //---------------------------------------------
  //
  //              b-tagged Jets Veto
  //
  //---------------------------------------------
  if (bjetList.size() < 2) return; 

  if (vetoMuon_wPt(10) != 1) return;

  TLorentzVector T1;

  if (selectTauList.size() < 2) return;      // No Tau1 means no Tau2, so we can safely reject the event

  vector<Tau> taulist;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& sig = selectTauList[indx];
    T1.SetPtEtaPhiE(sig.pt, sig.eta, sig.phi, sig.energy);
    double dr1 = AnaUtil::deltaR(M1, T1);
    if (dr1 < AnaUtil::cutValue(_tau1CutMap, "drMuTau")) continue;
    if (sig.byTightIsolationMVArun2v1DBoldDMwLT < 0.5) continue;
    if ((muo1.charge + sig.charge) != 0) continue;
    taulist.push_back(sig);
  } 
  if (taulist.size() > 0) return;

  vector<Tau> tau1list;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& fake = selectTauList[indx];
    T1.SetPtEtaPhiE(fake.pt, fake.eta, fake.phi, fake.energy);
    double dr1 = AnaUtil::deltaR(M1, T1);
    if (dr1 < AnaUtil::cutValue(_tau1CutMap, "drMuTau")) continue; 
    if ((muo1.charge + fake.charge) == 0) continue;
    tau1list.push_back(fake);
  } 

  if (tau1list.size() < 1) return;

  const Tau& fake = tau1list.at(0);

  if (fabs(fake.eta) < 0.8) {
    AnaUtil::fillHist1D("faketaupt_deno_TTJetsC_2j", fake.pt, puevWt_);
    if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
      AnaUtil::fillHist1D("faketaupt_nume_TTJetsC_2j", fake.pt, puevWt_);
  }
  if (fabs(fake.eta) >= 0.8 && fabs(fake.eta) < 1.6) {
    AnaUtil::fillHist1D("faketaupt_deno_TTJetsI_2j", fake.pt, puevWt_);
    if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
      AnaUtil::fillHist1D("faketaupt_nume_TTJetsI_2j", fake.pt, puevWt_);
  }

  if (fabs(fake.eta) >= 1.6 && fabs(fake.eta) < 2.3) {
    AnaUtil::fillHist1D("faketaupt_deno_TTJetsF_2j", fake.pt, puevWt_);
    if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
      AnaUtil::fillHist1D("faketaupt_nume_TTJetsF_2j", fake.pt, puevWt_);
  }
  if (fabs(fake.eta) < 2.3) {
    AnaUtil::fillHist1D("faketaupt_deno_TTJetsAll_2j", fake.pt, puevWt_);
    if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
      AnaUtil::fillHist1D("faketaupt_nume_TTJetsAll_2j", fake.pt, puevWt_);
  }
}
void FakeRate::JetToTauFakeTTJets_2j_ss()
{
  double vz = vtxList.at(0).z;

  // atleast 1 good muon
  //if (muoList.size() < 1) return;
  vector<Muon> fMuoList;
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    if (fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
        (fabs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0")) ||
        (!muon.isTightMuon)
       ) continue;

    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;

    if (fabs(muon.pfRelIsoDB03) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso") ||
        muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt")) continue;

    fMuoList.push_back(muon);
  }
  if (fMuoList.size() != 1) return; //only one good Muon
  const Muon& muo1 = fMuoList.at(0);
  TLorentzVector M1;
  M1.SetPtEtaPhiE (muo1.pt, muo1.eta, muo1.phi, muo1.energy);
  const MET& mt = metColl()->at(0);
  double mass1 = sqrt(2*muo1.pt*mt.met*(1-cos(AnaUtil::deltaPhi(muo1.phi, mt.metphi))));
  if (mass1 <= 40.0) return;                                                          
  if (mt.met <= 50.0) return;

  //---------------------------------------------
  //
  //              Electron Veto
  //
  //---------------------------------------------
  if (vetoElectron_wPt(15) != 0) return;

  //---------------------------------------------
  //
  //              b-tagged Jets Veto
  //
  //---------------------------------------------
  if (bjetList.size() < 2) return; 

  if (vetoMuon_wPt(10) != 1) return;

  TLorentzVector T1;

  if (selectTauList.size() < 2) return;

  vector<Tau> tau1list;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& fake = selectTauList[indx];
    T1.SetPtEtaPhiE(fake.pt, fake.eta, fake.phi, fake.energy);
    double dr1 = AnaUtil::deltaR(M1, T1);
    if (dr1 < AnaUtil::cutValue(_tau1CutMap, "drMuTau")) continue; 
    if ((muo1.charge + fake.charge) == 0) continue;
    tau1list.push_back(fake);
  } 

  if (tau1list.size() < 2) return;


  for (unsigned int indx = 0; indx < tau1list.size(); ++indx) {
    const Tau& fake = tau1list.at(indx);
    if (fabs(fake.eta) < 0.8) {
      AnaUtil::fillHist1D("faketaupt_deno_TTJetsC_2j_ss", fake.pt, puevWt_);
      if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
        AnaUtil::fillHist1D("faketaupt_nume_TTJetsC_2j_ss", fake.pt, puevWt_);
    }
    if (fabs(fake.eta) >= 0.8 && fabs(fake.eta) < 1.6) {
      AnaUtil::fillHist1D("faketaupt_deno_TTJetsI_2j_ss", fake.pt, puevWt_);
      if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
        AnaUtil::fillHist1D("faketaupt_nume_TTJetsI_2j_ss", fake.pt, puevWt_);
    }

    if (fabs(fake.eta) >= 1.6 && fabs(fake.eta) < 2.3) {
      AnaUtil::fillHist1D("faketaupt_deno_TTJetsF_2j_ss", fake.pt, puevWt_);
      if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
        AnaUtil::fillHist1D("faketaupt_nume_TTJetsF_2j_ss", fake.pt, puevWt_);
    }
    if (fabs(fake.eta) < 2.3) {
      AnaUtil::fillHist1D("faketaupt_deno_TTJetsAll_2j_ss", fake.pt, puevWt_);
      if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
        AnaUtil::fillHist1D("faketaupt_nume_TTJetsAll_2j_ss", fake.pt, puevWt_);
    }
  }
}
void FakeRate::JetToTauFakeTTJets_1j()
{
  double vz = vtxList.at(0).z;

  vector<Muon> fMuoList;
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    if (fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
        (fabs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0")) ||
        (!muon.isTightMuon)
       ) continue;

    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;

    if (fabs(muon.pfRelIsoDB03) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso") ||
        muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt")) continue;

    fMuoList.push_back(muon);
  }
  if (fMuoList.size() != 2) return; //only two good Muon
  const Muon& muo1 = fMuoList.at(0);
  const Muon& muo2 = fMuoList.at(1);
  if ((muo1.charge + muo2.charge) != 0) return;
  TLorentzVector M1, M2, zMass;
  M1.SetPtEtaPhiE (muo1.pt, muo1.eta, muo1.phi, muo1.energy);
  M2.SetPtEtaPhiE (muo2.pt, muo2.eta, muo2.phi, muo2.energy);
  zMass = M1 + M2;
  const MET& mt = metColl()->at(0);
  if (fabs(zMass.M() - 91) < 15) return; // Z Mass veto
  if (mt.met <= 50.0) return;

  //---------------------------------------------
  //
  //              Electron Veto
  //
  //---------------------------------------------
  if (vetoElectron_wPt(15) != 0) return;

  //---------------------------------------------
  //
  //              b-tagged Jets Veto
  //
  //---------------------------------------------
  if (bjetList.size() < 2) return; 

  if (vetoMuon_wPt(10) != 2) return;

  TLorentzVector T1;

  if (selectTauList.size() < 1) return;

  vector<Tau> tau1list;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& fake = selectTauList[indx];
    T1.SetPtEtaPhiE(fake.pt, fake.eta, fake.phi, fake.energy);
    double dr1 = AnaUtil::deltaR(M1, T1);
    double dr2 = AnaUtil::deltaR(M2, T1);
    if (dr1 < AnaUtil::cutValue(_tau1CutMap, "drMuTau")) continue; 
    if (dr2 < AnaUtil::cutValue(_tau1CutMap, "drMuTau")) continue; 
    tau1list.push_back(fake);
  } 

  if (tau1list.size() < 1) return;

  const Tau& fake = tau1list.at(0);

  if (fabs(fake.eta) < 0.8) {
    AnaUtil::fillHist1D("faketaupt_deno_TTJetsC_1j", fake.pt, puevWt_);
    if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
      AnaUtil::fillHist1D("faketaupt_nume_TTJetsC_1j", fake.pt, puevWt_);
  }
  if (fabs(fake.eta) >= 0.8 && fabs(fake.eta) < 1.6) {
    AnaUtil::fillHist1D("faketaupt_deno_TTJetsI_1j", fake.pt, puevWt_);
    if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
      AnaUtil::fillHist1D("faketaupt_nume_TTJetsI_1j", fake.pt, puevWt_);
  }

  if (fabs(fake.eta) >= 1.6 && fabs(fake.eta) < 2.3) {
    AnaUtil::fillHist1D("faketaupt_deno_TTJetsF_1j", fake.pt, puevWt_);
    if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
      AnaUtil::fillHist1D("faketaupt_nume_TTJetsF_1j", fake.pt, puevWt_);
  }
  if (fabs(fake.eta) < 2.3) {
    AnaUtil::fillHist1D("faketaupt_deno_TTJetsAll_1j", fake.pt, puevWt_);
    if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
      AnaUtil::fillHist1D("faketaupt_nume_TTJetsAll_1j", fake.pt, puevWt_);
  }
}

void FakeRate::JetToTauFakeQCD()
{
  double vz = vtxList.at(0).z;

  vector<Muon> fMuoList;
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    if (fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
        (fabs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0")) ||
        (!muon.isTightMuon)
       ) continue;

    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;

    //if (fabs(muon.pfRelIsoDB03) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso") ||
    if (fabs(muon.pfRelIsoDB03) < 0.30 ||
        muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt")) continue;

    fMuoList.push_back(muon);
  }
  if (fMuoList.size() != 1) return; //only two good Muon
  const Muon& muo1 = fMuoList.at(0);
  TLorentzVector M1;
  M1.SetPtEtaPhiE (muo1.pt, muo1.eta, muo1.phi, muo1.energy);

  const MET& mt = metColl()->at(0);
  if (mt.met > 20.0) return;

  //---------------------------------------------
  //
  //              Electron Veto
  //
  //---------------------------------------------
  if (vetoElectron_wPt(15) != 0) return;

  //---------------------------------------------
  //
  //              b-tagged Jets Veto
  //
  //---------------------------------------------
  if (bjetList.size() > 0) return; 

  if (vetoMuon_wPt(10) != 1) return;

  TLorentzVector T1;

  vector<Tau> tau1list;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& fake = selectTauList[indx];
    T1.SetPtEtaPhiE(fake.pt, fake.eta, fake.phi, fake.energy);
    double dr1 = AnaUtil::deltaR(M1, T1);
    if (dr1 < AnaUtil::cutValue(_tau1CutMap, "drMuTau")) continue; 
    tau1list.push_back(fake);
  } 

  if (tau1list.size() < 2) return;
  const Tau& ta1 = tau1list.at(0);
  const Tau& ta2 = tau1list.at(1);

  //W and TauTau charge correlation
  if (ta1.charge == ta2.charge) { 
    for (unsigned int indx = 0; indx < 2; ++indx) {
      const Tau& fake = tau1list.at(indx);
      if (fabs(fake.eta) < 0.8) {
        AnaUtil::fillHist1D("faketaupt_deno_QCDC_SS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_QCDC_SS", fake.pt, puevWt_);
      }

      if (fabs(fake.eta) >= 0.8 && fabs(fake.eta) < 1.6) {
        AnaUtil::fillHist1D("faketaupt_deno_QCDI_SS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_QCDI_SS", fake.pt, puevWt_);
      }

      if (fabs(fake.eta) >= 1.6 && fabs(fake.eta) < 2.3) {
        AnaUtil::fillHist1D("faketaupt_deno_QCDF_SS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_QCDF_SS", fake.pt, puevWt_);
      }
      if (fabs(fake.eta) < 2.3) {
        AnaUtil::fillHist1D("faketaupt_deno_QCDAll_SS", fake.pt, puevWt_);
        if (fake.byTightIsolationMVArun2v1DBoldDMwLT > 0.5)
          AnaUtil::fillHist1D("faketaupt_nume_QCDAll_SS", fake.pt, puevWt_);
      }
    }
  }
}

void FakeRate::ApplyFR()
{
  if (nmuon() < 1) return;
  double vz = vtxList.at(0).z;

/*
  // atleast 1 good muon
  if (muoList.size() < 1) return;
*/
  //To x-check whether the discrepancy was really coming from muon
  vector<Muon> fMuoList;
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    if (muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt")
      || fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta")) continue;
#if 0
    if (!muon.isTrackerMuon || !muon.isGlobalMuonPromptTight) continue;
    if (muon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) continue;
    if (muon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches") ||
        muon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) continue;
    if (muon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) continue;
    if (muon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) continue;
    if (muon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) continue;
    if (abs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0")) continue;
    //if (abs(muon.vtxDistZ) >= AnaUtil::cutValue(muonCutMap(), "vtxDistZ")) continue;
    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;
#endif
    if (!muon.isTightMuon) continue;
    //if (!muon.isMediumMuon) continue;
    //if (fabs(muon.chargedHadronIso/muon.pt) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso")) continue;
    if (fabs(muon.pfRelIsoDB03) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso")) continue;
    fMuoList.push_back(muon);
  }

  // atleast 1 good muon
  if (fMuoList.size() < 1) return;

  AnaUtil::fillHist1D("test", 2.0, puevWt_);    

  const Muon& muo = fMuoList.at(0);
  TLorentzVector M1,T;
  M1.SetPtEtaPhiE (muo.pt, muo.eta, muo.phi, muo.energy);

  if (selectTauList.size() < 2) return;

  vector<Tau> tau1List, tau2List, antitau1List, antitau2List, tauOSList, tauSSList;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& tau = selectTauList[indx];
    T.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
    double dr1 = AnaUtil::deltaR(M1, T);
    double vtxdz = fabs(tau.zvertex - vz);
    bool decision1 = dr1 >= AnaUtil::cutValue(_tau1CutMap, "drMuTau") &&
                     vtxdz < AnaUtil::cutValue(_tau1CutMap, "dz") &&
                     tau.pt > AnaUtil::cutValue(_tau2CutMap, "pt") &&
                     (tau.charge + muo.charge) == 0;

    if (decision1) tauOSList.push_back(tau);
    //if (decision1 &&  tau.byTightIsolationMVArun2v1DBoldDMwLT > 0.5) tau1List.push_back(tau);
  }

  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& tau = selectTauList[indx];
    T.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
    double dr1 = AnaUtil::deltaR(M1, T);
    double vtxdz = fabs(tau.zvertex - vz);

    bool decision2 = dr1 >= AnaUtil::cutValue(_tau2CutMap, "drMuTau") &&
                   tau.pt > AnaUtil::cutValue(_tau2CutMap, "pt") &&
                   vtxdz < AnaUtil::cutValue(_tau2CutMap, "dz") &&
      //tau.againstElectronMediumMVA6 >= 0.5 &&
                   (tau.charge + muo.charge) != 0;

    if (decision2) tauSSList.push_back(tau);
    //if (decision2 && tau.byTightIsolationMVArun2v1DBoldDMwLT > 0.5) tau2List.push_back(tau);
  }


  if (tauOSList.size() < 1 || tauSSList.size() < 1) return;    //No GOOD Tau1 and antiTau1, so cancel the event
  if (tauOSList[0].pt <= AnaUtil::cutValue(_tau1CutMap, "pt") && tauSSList[0].pt <= AnaUtil::cutValue(_tau1CutMap, "pt")) return;


  for (unsigned int indx = 0; indx < tauOSList.size(); ++indx) {
    const Tau& tau = tauOSList[indx];

    if (tau.byMediumIsolationMVArun2v1DBoldDMwLT > 0.5) tau1List.push_back(tau);
    else antitau1List.push_back(tau);
  }
  for (unsigned int indx = 0; indx < tauSSList.size(); ++indx) {
    const Tau& tau = tauSSList[indx];

    if (tau.byTightIsolationMVArun2v1DBoldDMwLT > 0.5) tau2List.push_back(tau);
    else antitau2List.push_back(tau);
  }

  int tau1_indx, tau2_indx;
  double tauazvtx, taubzvtx;
  TLorentzVector T1, T2;

  double tau1charge;
  bool foundtaupair = false;
  bool signal = false;
  if (tau1List.size() > 0 && tau2List.size() > 0) {
    for (unsigned int indx = 0; indx < tau1List.size(); ++indx) {
      const Tau& tau1 = tau1List[indx];
      T1.SetPtEtaPhiE(tau1.pt, tau1.eta, tau1.phi, tau1.energy);
      for (unsigned int jndx = 0; jndx < tau2List.size(); ++jndx) {
        const Tau& tau2 = tau2List[jndx];
        T2.SetPtEtaPhiE(tau2.pt, tau2.eta, tau2.phi, tau2.energy);
        if (AnaUtil::deltaR(T1, T2) >= 0.5 && (tau1.pt > AnaUtil::cutValue(_tau1CutMap, "pt") || tau2.pt > AnaUtil::cutValue(_tau1CutMap, "pt"))) {
          tau2_indx = jndx;
          taubzvtx = tau2.zvertex;
          foundtaupair = true;
          break;
        }
      }
      if (foundtaupair) { 
        tau1_indx = indx;
        tauazvtx = tau1.zvertex;
        signal = true;
        tau1charge = tau1.charge;
        break;
      }
    }
  }

  bool OSantiSSanti = false;
  if (!foundtaupair && tau1List.size() < 1 && tau2List.size() < 1 && antitau1List.size() > 0 && antitau2List.size() > 0) {
    for (unsigned int indx = 0; indx < antitau1List.size(); ++indx) {
      const Tau& tau1 = antitau1List[indx];
      T1.SetPtEtaPhiE(tau1.pt, tau1.eta, tau1.phi, tau1.energy);
      for (unsigned int jndx = 0; jndx < antitau2List.size(); ++jndx) {
        const Tau& tau2 = antitau2List[jndx];
        T2.SetPtEtaPhiE(tau2.pt, tau2.eta, tau2.phi, tau2.energy);
        if (AnaUtil::deltaR(T1, T2) >= 0.5 && (tau1.pt > AnaUtil::cutValue(_tau1CutMap, "pt") || tau2.pt > AnaUtil::cutValue(_tau1CutMap, "pt"))) {
          tau2_indx = jndx;
          taubzvtx = tau2.zvertex;
          OSantiSSanti = true;
          break;
        }
      }
      if (OSantiSSanti) {
        tau1_indx = indx;
        tauazvtx = tau1.zvertex;
        break;
      }
    }
  }

  bool OSantiCR = false;
  if (!foundtaupair && !OSantiSSanti && tau2List.size() > 0 && antitau1List.size() > 0) {
    for (unsigned int indx = 0; indx < antitau1List.size(); ++indx) {
      const Tau& tau1 = antitau1List[indx];
      T1.SetPtEtaPhiE(tau1.pt, tau1.eta, tau1.phi, tau1.energy);
      for (unsigned int jndx = 0; jndx < tau2List.size(); ++jndx) {
        const Tau& tau2 = tau2List[jndx];
        T2.SetPtEtaPhiE(tau2.pt, tau2.eta, tau2.phi, tau2.energy);
        if (AnaUtil::deltaR(T1, T2) >= 0.5 && (tau1.pt > AnaUtil::cutValue(_tau1CutMap, "pt") || tau2.pt > AnaUtil::cutValue(_tau1CutMap, "pt"))) {
          tau2_indx = jndx;
          taubzvtx = tau2.zvertex;
          OSantiCR = true;
          break;
        }
      }
      if (OSantiCR) {
        tau1_indx = indx;
        tauazvtx = tau1.zvertex;
        break;
      }
    }
  }

  bool JetToTau2Fake = false;
  //This is the case of JetToTau2 Fake//SS
  if (!foundtaupair && !OSantiSSanti && !OSantiCR && tau1List.size() > 0 && antitau2List.size() > 0) {
    for (unsigned int indx = 0; indx < tau1List.size(); ++indx) {
      const Tau& tau1 = tau1List[indx];
      T1.SetPtEtaPhiE(tau1.pt, tau1.eta, tau1.phi, tau1.energy);
      for (unsigned int jndx = 0; jndx < antitau2List.size(); ++jndx) {
        const Tau& tau2 = antitau2List[jndx];
        T2.SetPtEtaPhiE(tau2.pt, tau2.eta, tau2.phi, tau2.energy);
        if (AnaUtil::deltaR(T1, T2) >= 0.5 && (tau1.pt > AnaUtil::cutValue(_tau1CutMap, "pt") || tau2.pt > AnaUtil::cutValue(_tau1CutMap, "pt"))) {
          tau2_indx = jndx;
          taubzvtx = tau2.zvertex;
          JetToTau2Fake = true;
          break;
        }
      }
      if (JetToTau2Fake) {
        tauazvtx = tau1.zvertex;
        tau1_indx = indx;
        break;
      }
    }    
  }
  

  //if (!signal && !JetToTau2Fake && !antiTau2 && !OSantiCR && !OSantiSSanti) return;
  if (fabs(tauazvtx - taubzvtx) >= 0.14) return;

  bool isGoodVtx;
  TVector3 vtxmu = findLeptonVtx(muo.vtxIndex, isGoodVtx);
  double mutaudz = fabs(tauazvtx - vtxmu.z());
  if (!isGoodVtx || mutaudz >= 0.14) return;

  //---------------------------------------------
  //
  //              Muon Veto
  //
  //---------------------------------------------
  if (vetoMuon(tauazvtx, 10, 0.14) != 1) return;

  if (signal)
    AnaUtil::fillHist1D("test", 3.0, puevWt_);    


  //---------------------------------------------
  //
  //              Electron Veto
  //
  //---------------------------------------------
  if (vetoElectron(tauazvtx, 15, 0.14) != 0) return;

  if (signal)
    AnaUtil::fillHist1D("test", 4.0, puevWt_);    

  //No b-Jet Veto
  if (bjetList.size() > 0) return;

  if (signal)
    AnaUtil::fillHist1D("test", 5.0, puevWt_);    

  const MET& mt = metColl()->at(0);
  //if (mt->met < 20) return;

  if (signal)
    AnaUtil::fillHist1D("test", 6.0, puevWt_);    

  double mass1 = sqrt(2*muo.pt*mt.met*(1-cos(AnaUtil::deltaPhi(muo.phi, mt.metphi))));
  //if (mass1 <= 30) return;
  if (signal)
    AnaUtil::fillHist1D("met", mt.met, puevWt_); //after choosing signal

  if (signal)
    AnaUtil::fillHist1D("test", 7.0, puevWt_);    

  TLorentzVector zMuTau, diTau;
  diTau = T1 + T2;

  double mvaOutputFK = -999;      //Set at high negative value, as MvaOutput > -0.106
  if (_readMVAFK) {
    leadPt     = T1.Pt();
    subPt     = T2.Pt();
    met        = mt.met;
    deltaRDiTau   = AnaUtil::deltaR(T1, T2);
    ptRatio = diTau.Pt()/(T1.Pt() + T2.Pt());

    mvaOutputFK = reader1->EvaluateMVA("BDT8");
    //if (mvaOutputFK <= -0.106) return;
  }

  double mvaOutput = -999;	//Set at high negative value, as MvaOutput > -0.106
  if (_readMVA) {
    muEta      = M1.Eta();
    muPt       = M1.Pt();
    tau1Eta    = T1.Eta();
    tau1Pt     = T1.Pt();
    //tau2Eta    = T2.Eta();
    tau2Pt     = T2.Pt();
    //diTaudR   = AnaUtil::deltaR(T1, T2);
    //dphiMuTau1 = AnaUtil::deltaPhi(M1, T1);
    //dphiMuDiTau = AnaUtil::deltaPhi(M1, diTau);
    met        = mt.met;
    ptRatio = diTau.Pt()/(T1.Pt() + T2.Pt());

    //mvaOutput = reader->EvaluateMVA("MLPBNN");
    mvaOutput = reader->EvaluateMVA("MLP");

    //for the time being returning at this value
    //if (mvaOutput <= 0.65) return;     //This value has to be optimized
  }

  if (JetToTau2Fake) {
    zMuTau = T1 + M1;
    if (mt.met < 20) return;
    if (mass1 <= 30) return;
    if (zMuTau.M() >= 80 || diTau.Pt() >= 50) {
      if (isMC()) AnaUtil::fillHist1D("eventWeight_fakeCR_mtt", 1.0, eventWt_*puevWt_);

      AnaUtil::fillHist1D("tau2fake_pt", T2.Pt(), puevWt_);
      AnaUtil::fillHist1D("tautauMass_fakeCR", diTau.M(),  puevWt_);
      AnaUtil::fillHist2D("tautauMass_fakeCR_p", diTau.M(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tau1Pt_fakeCR_p", T1.Pt(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tau2Pt_fakeCR_p", T2.Pt(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tau1Eta_fakeCR_p", T1.Eta(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tau2Eta_fakeCR_p", T2.Eta(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tau1Phi_fakeCR_p", T1.Phi(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tau2Phi_fakeCR_p", T2.Phi(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("muPt_fakeCR_p", M1.Pt(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("muEta_fakeCR_p", M1.Eta(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("muPhi_fakeCR_p", M1.Phi(), T2.Pt(), puevWt_);
      AnaUtil::fillHist3D("tautauMass_fakeCR_tau2fake_pe", T2.Pt(), T2.Eta(), diTau.M(), puevWt_);
      //AnaUtil::fillHist3D("mvaoutput_fakeCR_tau2fake_pe", T2.Pt(), T2.Eta(), mvaOutput, puevWt_);
      //AnaUtil::fillHist3D("mvaOutputFK_fakeCR_tau2fake_pe", T2.Pt(), T2.Eta(), mvaOutputFK, puevWt_);
      if (mvaOutputFK > -0.40) {
        AnaUtil::fillHist1D("tau2fake_pt_mva1", T2.Pt(), puevWt_);
        AnaUtil::fillHist1D("tautauMass_fakeCR_mva1", diTau.M(), puevWt_);
        AnaUtil::fillHist2D("tautauMass_fakeCR_tau2fake_mva1", T2.Pt(), diTau.M(), puevWt_);
        AnaUtil::fillHist2D("mvaoutput_fakeCR_tau2fake_mva1", T2.Pt(), mvaOutput, puevWt_);
        AnaUtil::fillHist1D("tau2fake_eta", T2.Eta(), puevWt_);
        AnaUtil::fillHist1D("tau2fake_phi", T2.Phi(), puevWt_);
        AnaUtil::fillHist1D("tau2fake_ene", T2.Energy(), puevWt_);
      }
      if (mvaOutputFK > -0.35) {
        AnaUtil::fillHist1D("tau2fake_pt_mva2", T2.Pt(), puevWt_);
        AnaUtil::fillHist1D("tautauMass_fakeCR_mva2", diTau.M(), puevWt_);
        AnaUtil::fillHist2D("tautauMass_fakeCR_tau2fake_mva2", T2.Pt(), diTau.M(), puevWt_);
        AnaUtil::fillHist2D("mvaoutput_fakeCR_tau2fake_mva2", T2.Pt(), mvaOutput, puevWt_);
      }
      if (mvaOutputFK > -0.30) {
        AnaUtil::fillHist1D("tau2fake_pt_mva3", T2.Pt(), puevWt_);
        AnaUtil::fillHist1D("tautauMass_fakeCR_mva3", diTau.M(), puevWt_);
        AnaUtil::fillHist2D("tautauMass_fakeCR_tau2fake_mva3", T2.Pt(), diTau.M(), puevWt_);
        AnaUtil::fillHist2D("mvaoutput_fakeCR_tau2fake_mva3", T2.Pt(), mvaOutput, puevWt_);
      }
      if (mvaOutputFK > -0.25) {
        AnaUtil::fillHist1D("tau2fake_pt_mva4", T2.Pt(), puevWt_); 
        AnaUtil::fillHist1D("tautauMass_fakeCR_mva4", diTau.M(), puevWt_);
        AnaUtil::fillHist2D("tautauMass_fakeCR_tau2fake_mva4", T2.Pt(), diTau.M(), puevWt_);
        AnaUtil::fillHist2D("mvaoutput_fakeCR_tau2fake_mva4", T2.Pt(), mvaOutput, puevWt_);
      }
      if (mvaOutputFK > -0.20) {
        AnaUtil::fillHist1D("tau2fake_pt_mva5", T2.Pt(), puevWt_);
        AnaUtil::fillHist1D("tautauMass_fakeCR_mva5", diTau.M(), puevWt_);
        AnaUtil::fillHist2D("tautauMass_fakeCR_tau2fake_mva5", T2.Pt(), diTau.M(), puevWt_);
        AnaUtil::fillHist2D("mvaoutput_fakeCR_tau2fake_mva5", T2.Pt(), mvaOutput, puevWt_);
      }
      if (mvaOutputFK > -0.15) {
        AnaUtil::fillHist1D("tau2fake_pt_mva6", T2.Pt(), puevWt_);
        AnaUtil::fillHist1D("tautauMass_fakeCR_mva6", diTau.M(), puevWt_);
        AnaUtil::fillHist2D("tautauMass_fakeCR_tau2fake_mva6", T2.Pt(), diTau.M(), puevWt_);
        AnaUtil::fillHist2D("mvaoutput_fakeCR_tau2fake_mva6", T2.Pt(), mvaOutput, puevWt_);
      }
      if (mvaOutputFK > -0.10) {
        AnaUtil::fillHist1D("tau2fake_pt_mva7", T2.Pt(), puevWt_);
        AnaUtil::fillHist1D("tautauMass_fakeCR_mva7", diTau.M(), puevWt_);
        AnaUtil::fillHist2D("tautauMass_fakeCR_tau2fake_mva7", T2.Pt(), diTau.M(), puevWt_);
        AnaUtil::fillHist2D("mvaoutput_fakeCR_tau2fake_mva7", T2.Pt(), mvaOutput, puevWt_);
      }
    }
  }

  if (signal) {
    zMuTau =T1 + M1;
    if (mt.met < 20) return;
    if (mass1 <= 30) return;
    if (zMuTau.M() >= 80 || diTau.Pt() >= 50) {
      if (isMC()) AnaUtil::fillHist1D("eventWeight_BL_mtt", 1.0, eventWt_*puevWt_);
      AnaUtil::fillHist1D("test", 8.0, puevWt_);    
      AnaUtil::fillHist1D("EventCount", 2.0, puevWt_);
      AnaUtil::fillHist1D("tautauMass_BL", diTau.M(), puevWt_);    
      AnaUtil::fillHist1D("tau1Pt_BL", T1.Pt(), puevWt_);    
      AnaUtil::fillHist1D("tau2Pt_BL", T2.Pt(), puevWt_);    
      AnaUtil::fillHist1D("tau1Eta_BL", T1.Eta(), puevWt_);    
      AnaUtil::fillHist1D("tau2Eta_BL", T2.Eta(), puevWt_);    
      AnaUtil::fillHist1D("tau1Phi_BL", T1.Phi(), puevWt_);    
      AnaUtil::fillHist1D("tau2Phi_BL", T2.Phi(), puevWt_);    
      AnaUtil::fillHist1D("muPt_BL", M1.Pt(), puevWt_);    
      AnaUtil::fillHist1D("muEta_BL", M1.Eta(), puevWt_);    
      AnaUtil::fillHist1D("muPhi_BL", M1.Phi(), puevWt_);    
      AnaUtil::fillHist1D("mvaoutput_BL", mvaOutput, puevWt_);    
      AnaUtil::fillHist1D("mvaOutputFK_BL", mvaOutputFK, puevWt_);
      if (mvaOutputFK > -0.40) {
        AnaUtil::fillHist1D("EventCount_mva1", 2.0, puevWt_);
        AnaUtil::fillHist1D("tautauMass_mva1", diTau.M(), puevWt_);    
        AnaUtil::fillHist1D("mvaoutput_mva1", mvaOutput, puevWt_);
      }
      if (mvaOutputFK > -0.35) {
        AnaUtil::fillHist1D("EventCount_mva2", 2.0, puevWt_);
        AnaUtil::fillHist1D("tautauMass_mva2", diTau.M(), puevWt_);    
        AnaUtil::fillHist1D("mvaoutput_mva2", mvaOutput, puevWt_);
      }
      if (mvaOutputFK > -0.30) {
        AnaUtil::fillHist1D("EventCount_mva3", 2.0, puevWt_);
        AnaUtil::fillHist1D("tautauMass_mva3", diTau.M(), puevWt_);    
        AnaUtil::fillHist1D("mvaoutput_mva3", mvaOutput, puevWt_);
      }
      if (mvaOutputFK > -0.25) {
        AnaUtil::fillHist1D("EventCount_mva4", 2.0, puevWt_);
        AnaUtil::fillHist1D("mvaoutput_mva4", mvaOutput, puevWt_);
        AnaUtil::fillHist1D("tautauMass_mva4", diTau.M(), puevWt_);            
      }
      if (mvaOutputFK > -0.20) {
        AnaUtil::fillHist1D("EventCount_mva5", 2.0, puevWt_);
        AnaUtil::fillHist1D("tautauMass_mva5", diTau.M(), puevWt_);    
        AnaUtil::fillHist1D("mvaoutput_mva5", mvaOutput, puevWt_);
      }
      if (mvaOutputFK > -0.15) {
        AnaUtil::fillHist1D("EventCount_mva6", 2.0, puevWt_);
        AnaUtil::fillHist1D("tautauMass_mva6", diTau.M(), puevWt_);    
        AnaUtil::fillHist1D("mvaoutput_mva6", mvaOutput, puevWt_);
      }
      if (mvaOutputFK > -0.10) {
        AnaUtil::fillHist1D("EventCount_mva7", 2.0, puevWt_);
        AnaUtil::fillHist1D("tautauMass_mva7", diTau.M(), puevWt_);    
        AnaUtil::fillHist1D("mvaoutput_mva7", mvaOutput, puevWt_);
      }
    }
  }

  if (OSantiCR) {
    zMuTau =T1 + M1;
    //if (zMuTau.M() >= 80 || diTau.Pt() >= 50) {
      if (isMC()) AnaUtil::fillHist1D("eventWeight_CR_mtt", 1.0, eventWt_*puevWt_);
      AnaUtil::fillHist1D("tau1Pt_CR", T1.Pt(), puevWt_);
      AnaUtil::fillHist1D("tau1Eta_CR", T1.Eta(), puevWt_);
      AnaUtil::fillHist1D("tau1Phi_CR", T1.Phi(), puevWt_);
      AnaUtil::fillHist1D("tau2Pt_CR", T2.Pt(), puevWt_);
      AnaUtil::fillHist1D("tau2Eta_CR", T2.Eta(), puevWt_);
      AnaUtil::fillHist1D("tau2Phi_CR", T2.Phi(), puevWt_);
      AnaUtil::fillHist1D("muPt_CR", M1.Pt(), puevWt_);
      AnaUtil::fillHist1D("muEta_CR", M1.Eta(), puevWt_);
      AnaUtil::fillHist1D("muPhi_CR", M1.Phi(), puevWt_);
      AnaUtil::fillHist1D("tautauMass_CR", diTau.M(), puevWt_);
      AnaUtil::fillHist1D("mvaoutput_CR", mvaOutput, puevWt_);    
      AnaUtil::fillHist1D("mvaOutputFK_CR", mvaOutputFK, puevWt_);

      AnaUtil::fillHist2D("tau1Pt_CR_e", T1.Pt(), T2.Eta(), puevWt_);
      AnaUtil::fillHist2D("tau1Eta_CR_e", T1.Eta(), T2.Eta(), puevWt_);
      AnaUtil::fillHist2D("tau1Phi_CR_e", T1.Phi(), T2.Eta(), puevWt_);
      AnaUtil::fillHist2D("tau2Pt_CR_e", T2.Pt(), T2.Eta(), puevWt_);
      AnaUtil::fillHist2D("tau2Eta_CR_e", T2.Eta(), T2.Eta(), puevWt_);
      AnaUtil::fillHist2D("tau2Phi_CR_e", T2.Phi(), T2.Eta(), puevWt_);
      AnaUtil::fillHist2D("muPt_CR_e", M1.Pt(), T2.Eta(), puevWt_);
      AnaUtil::fillHist2D("muEta_CR_e", M1.Eta(), T2.Eta(), puevWt_);
      AnaUtil::fillHist2D("muPhi_CR_e", M1.Phi(), T2.Eta(), puevWt_);
      AnaUtil::fillHist2D("tautauMass_CR_e", diTau.M(), T2.Eta(), puevWt_);
      //}
  }
  if (OSantiSSanti) {
    zMuTau =T1 + M1;
    //if (zMuTau.M() >= 80 || diTau.Pt() >= 50) {
      if (isMC()) AnaUtil::fillHist1D("eventWeight_AA_mtt", 1.0, eventWt_*puevWt_);
      AnaUtil::fillHist3D("tau1Pt_AA_tau2fake_pe", T2.Pt(), T2.Eta(), T1.Pt(), puevWt_);
      AnaUtil::fillHist3D("tau2Pt_AA_tau2fake_pe", T2.Pt(), T2.Eta(), T2.Pt(), puevWt_);
      AnaUtil::fillHist3D("muPt_AA_tau2fake_pe", T2.Pt(), T2.Eta(), M1.Pt(), puevWt_);
      AnaUtil::fillHist3D("tautauMass_AA_tau2fake_pe", T2.Pt(), T2.Eta(), diTau.M(), puevWt_);
      AnaUtil::fillHist3D("mvaoutput_AA_tau2fake_pe", T2.Pt(), T2.Eta(), mvaOutput, puevWt_);
      AnaUtil::fillHist3D("mvaOutputFK_AA_tau2fake_pe", T2.Pt(), T2.Eta(), mvaOutputFK, puevWt_);

      AnaUtil::fillHist2D("tau1Pt_AA_tau2fake_p", T1.Pt(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tau1Eta_AA_tau2fake_p", T1.Eta(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tau1Phi_AA_tau2fake_p", T1.Phi(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tau2Pt_AA_tau2fake_p", T2.Pt(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tau2Eta_AA_tau2fake_p", T2.Eta(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tau2Phi_AA_tau2fake_p", T2.Phi(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("muPt_AA_tau2fake_p", M1.Pt(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("muEta_AA_tau2fake_p", M1.Eta(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("muPhi_AA_tau2fake_p", M1.Phi(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("tautauMass_AA_tau2fake_p", diTau.M(), T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("mvaoutput_AA_tau2fake_p", mvaOutput, T2.Pt(), puevWt_);
      AnaUtil::fillHist2D("mvaOutputFK_AA_tau2fake_p", mvaOutputFK, T2.Pt(), puevWt_);

      AnaUtil::fillHist1D("tau1Pt_AA", T1.Pt(), puevWt_);
      AnaUtil::fillHist1D("tau2Pt_AA", T2.Pt(), puevWt_);
      AnaUtil::fillHist1D("muPt_AA", M1.Pt(), puevWt_);
      AnaUtil::fillHist1D("tautauMass_AA", diTau.M(), puevWt_);
      //}
  }
}

/////////////////////////////////////////////////////////////////////
///////////////Here comes the fake estimation for MuMuTau channel////
/////////////////////////////////////////////////////////////////////
void FakeRate::MuMuTauFakeWJets()
{
  if (selectTightMuonList.size() < 1) return;
  TLorentzVector M1;
  float tagIndx = -1;
  for (unsigned int indx = 0; indx < selectTightMuonList.size(); ++indx) {
    const Muon& tag = selectTightMuonList[indx];
    M1.SetPtEtaPhiE(tag.pt, tag.eta, tag.phi, tag.energy);
    if (fabs(tag.pfRelIsoDB03) >= 0.15 ||
        tag.pt <= AnaUtil::cutValue(muonCutMap(), "ptTag")) continue;
    tagIndx = indx;
    break;
  }
  if (tagIndx < 0) return; //No Tag Muon Found
  const Muon& tagMuon = selectTightMuonList[tagIndx];

  TLorentzVector M2;
  float probeIndx = -1;
  for (unsigned int indx = 0; indx < selectLooseMuonList.size(); ++indx) {
    const Muon& probe = selectLooseMuonList[indx];
    M2.SetPtEtaPhiE(probe.pt, probe.eta, probe.phi, probe.energy);
    if (AnaUtil::deltaR(M1, M2) < 0.5) continue;
    if ((tagMuon.charge + probe.charge) == 0) continue;  //To avoid Z/Gamma jets//Do Not need for WJets MC
    if (probe.pt <= AnaUtil::cutValue(muonCutMap(), "ptProbe")) continue;
    probeIndx = indx;
    break;
  }
  if (probeIndx < 0) return;
  const Muon& probeMuon = selectLooseMuonList[probeIndx];

  const MET& pfmt = metColl()->at(0);
  //const MET& mvamt = mvametColl()->at(0);
  float TagMet_MT_pf  = sqrt(2*tagMuon.pt*pfmt.met*(1-cos(AnaUtil::deltaPhi(tagMuon.phi, pfmt.metphi))));
  //float TagMet_MT_mva = sqrt(2*tagMuon.pt*mvamt.met*(1-cos(AnaUtil::deltaPhi(tagMuon.phi, mvamt.metphi))));

  AnaUtil::fillHist1D("TagpfMet_MT_mmt_WJets", TagMet_MT_pf, puevWt_);
  //AnaUtil::fillHist1D("TagmvaMet_MT_mmt_WJets", TagMet_MT_mva, puevWt_);
  if (TagMet_MT_pf < 35.0) return;

  float ProbeMet_MT = sqrt(2*probeMuon.pt*pfmt.met*(1-cos(AnaUtil::deltaPhi(probeMuon.phi, pfmt.metphi))));
  if (ProbeMet_MT > 35.0) return;

  if (vetoMuon_wPt(10.0) != 2) return; //No Extra Muons
  if (vetoElectron_wPt(10.0) != 0) return; // No Extra Electrons
  if (bjetList.size() > 0) return; //No bJets
  if (selectTauList.size() > 0 ){
    if ((selectTauList[0].pt > 20.0) && (selectTauList[0].chargedIsoPtSum < 1.0)) return; //No Selected Taus
  }
  AnaUtil::fillHist1D("JetToMuFake_WJets_mmt_deno_dummy", probeMuon.pt, puevWt_);  
  //Require 1 Jet to mimic the presence of a Tau object
  float nJet = 0;
  TLorentzVector J1;
  bool JetFlag = false;
  for (auto it = jetColl()->begin(); it != jetColl()->end(); ++it) {
    const Jet& jt = (*it);
    if (fabs(jt.eta) >= 2.4 || jt.pt <= 20.0) continue;
    nJet += 1;
    J1.SetPtEtaPhiE(jt.pt, jt.eta, jt.phi, jt.energy);
    if (AnaUtil::deltaR(M1, J1) < 0.5 || AnaUtil::deltaR(M2, J1) < 0.5) continue;
    JetFlag = true;
    //break; //just need 1 non-overlapping jet
  }
  if (!JetFlag) return;

  //Fill Here the isolation denominator and numerator
  AnaUtil::fillHist1D("JetToMuFake_WJets_mmt_deno", probeMuon.pt, puevWt_);  

  //Check if Probe Muon Passes the Tight ID
  bool probeTID = true;
  if (fabs(probeMuon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
#if 0
     (!probeMuon.isTrackerMuon) ||
     (!probeMuon.isGlobalMuonPromptTight) ||
     (probeMuon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) ||
     (probeMuon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches")) ||
     (probeMuon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) ||
     (probeMuon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) ||
     (probeMuon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) ||
     (probeMuon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) ||
     (fabs(probeMuon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0"))
#endif
     (!probeMuon.isTightMuon)
  ) probeTID = false;

  if (fabs(probeMuon.eta) < 1.479){
    if (probeMuon.pfRelIsoDB03 < 0.15 && probeTID) AnaUtil::fillHist1D("JetToMuFake_WJets_mmt_nume_lead", probeMuon.pt, puevWt_);
    if (probeMuon.pfRelIsoDB03 < 0.20 && probeTID) AnaUtil::fillHist1D("JetToMuFake_WJets_mmt_nume_sublead", probeMuon.pt, puevWt_);
  }
  else{
    if (probeMuon.pfRelIsoDB03 < 0.10 && probeTID) AnaUtil::fillHist1D("JetToMuFake_WJets_mmt_nume_lead", probeMuon.pt, puevWt_);
    if (probeMuon.pfRelIsoDB03 < 0.15 && probeTID) AnaUtil::fillHist1D("JetToMuFake_WJets_mmt_nume_sublead", probeMuon.pt, puevWt_);
  }

  //Fill histogram for kNN training, this should be called only once
  //skimForKNN(probeMuon, probeTID, nJet);
}

void FakeRate::MuMuTauFakeQCD()
{
  if (selectTightMuonList.size() < 1) return;
  TLorentzVector M1;
  float tagIndx = -1;
  for (unsigned int indx = 0; indx < selectTightMuonList.size(); ++indx) {
    const Muon& tag = selectTightMuonList[indx];
    M1.SetPtEtaPhiE(tag.pt, tag.eta, tag.phi, tag.energy);
    //Anti-Isolation for QCD
    if (fabs(tag.pfRelIsoDB03) < 0.30 ||
        tag.pt <= AnaUtil::cutValue(muonCutMap(), "ptTag")) continue;
    tagIndx = indx;
    break;
  }
  if (tagIndx < 0) return; //No Tag Muon Found
  const Muon& tagMuon = selectTightMuonList[tagIndx];

  TLorentzVector M2;
  float probeIndx = -1;
  for (unsigned int indx = 0; indx < selectLooseMuonList.size(); ++indx) {
    const Muon& probe = selectLooseMuonList[indx];
    M2.SetPtEtaPhiE(probe.pt, probe.eta, probe.phi, probe.energy);
    if (AnaUtil::deltaR(M1, M2) < 0.5) continue;
    if ((tagMuon.charge + probe.charge) == 0) continue;  //To avoid Z/Gamma jets
    if (probe.pt <= AnaUtil::cutValue(muonCutMap(), "ptProbe")) continue;
    probeIndx = indx;
    break;
  }
  if (probeIndx < 0) return;
  const Muon& probeMuon = selectLooseMuonList[probeIndx];

  const MET& mt = metColl()->at(0);
  if (mt.met > 20.0) return;

  if (vetoMuon_wPt(5.0) != 2) return; //No Extra Muons
  if (vetoElectron_wPt(10.0) != 0) return; // No Extra Electrons
  if (bjetList.size() > 0) return; //No bJets
  if (selectTauList.size() > 0 ){
    if ((selectTauList[0].pt > 20.0) && (selectTauList[0].chargedIsoPtSum < 2.0)) return; //No Selected Taus
  }
  //Require 1 Jet to mimic the presence of a Tau object
  TLorentzVector J1;
  bool JetFlag = false;
  for (auto it = jetColl()->begin(); it != jetColl()->end(); ++it) {
    const Jet& jt = (*it);
    if (fabs(jt.eta) >= 2.4 || jt.pt <= 20.0) continue;
    J1.SetPtEtaPhiE(jt.pt, jt.eta, jt.phi, jt.energy);
    if (AnaUtil::deltaR(M1, J1) < 0.5 || AnaUtil::deltaR(M2, J1) < 0.5) continue;
    JetFlag = true;
    break; //just need 1 non-overlapping jet
  }
  if (!JetFlag) return;

  //Fill Here the isolation denominator and numerator
  AnaUtil::fillHist1D("JetToMuFake_QCD_mmt_deno", probeMuon.pt, puevWt_);  

  //Check if Probe Muon Passes the Tight ID
  if (fabs(probeMuon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
#if 0
     (!probeMuon.isTrackerMuon) ||
     (!probeMuon.isGlobalMuonPromptTight) ||
     (probeMuon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) ||
     (probeMuon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches")) ||
     (probeMuon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) ||
     (probeMuon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) ||
     (probeMuon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) ||
     (probeMuon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) ||
     (fabs(probeMuon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0"))
#endif
     (!probeMuon.isMediumMuon)
  ) return;

  if (fabs(probeMuon.eta) < 1.479){
    if (probeMuon.pfRelIsoDB03 < 0.15) AnaUtil::fillHist1D("JetToMuFake_QCD_mmt_nume_lead", probeMuon.pt, puevWt_);
    if (probeMuon.pfRelIsoDB03 < 0.20) AnaUtil::fillHist1D("JetToMuFake_QCD_mmt_nume_sublead", probeMuon.pt, puevWt_);
  }
  else{
    if (probeMuon.pfRelIsoDB03 < 0.10) AnaUtil::fillHist1D("JetToMuFake_QCD_mmt_nume_lead", probeMuon.pt, puevWt_);
    if (probeMuon.pfRelIsoDB03 < 0.15) AnaUtil::fillHist1D("JetToMuFake_QCD_mmt_nume_sublead", probeMuon.pt, puevWt_);
  }
}

void FakeRate::MuMuTauFakeZJets()
{
  if (selectTightMuonList.size() < 1) return;
  TLorentzVector M1;
  float tagIndx = -1;
  for (unsigned int indx = 0; indx < selectTightMuonList.size(); ++indx) {
    const Muon& tag = selectTightMuonList[indx];
    M1.SetPtEtaPhiE(tag.pt, tag.eta, tag.phi, tag.energy);
    if (fabs(tag.pfRelIsoDB03) >= 0.15 ||
        tag.pt <= AnaUtil::cutValue(muonCutMap(), "ptTag")) continue;
    tagIndx = indx;
    break;
  }
  if (tagIndx < 0) return; //No Tag Muon Found
  const Muon& tagMuon = selectTightMuonList[tagIndx];

  TLorentzVector M2;
  float tag2Indx = -1;
  for (unsigned int indx = 0; indx < selectTightMuonList.size(); ++indx) {
    const Muon& tag2 = selectTightMuonList[indx];
    M2.SetPtEtaPhiE(tag2.pt, tag2.eta, tag2.phi, tag2.energy);
    if (AnaUtil::deltaR(M1, M2) < 0.5) continue;
    if ((tagMuon.charge + tag2.charge) != 0) continue;  //To select Z-jets
    if (tag2.pt <= AnaUtil::cutValue(muonCutMap(), "ptTag")) continue;
    //if (fabs(tag2.chargedHadronIso/tag2.pt) >= 0.15 || tag2.pt <= AnaUtil::cutValue(muonCutMap(), "ptTag")) continue;
    TLorentzVector Mass;
    Mass = M1 + M2;
    AnaUtil::fillHist1D("Tag_ZMass_mmt_ZJets", Mass.M(), puevWt_);
    if (Mass.M() > 110.0 || Mass.M() < 70.0) continue;  // Z window
    tag2Indx = indx;
    break;
  }
  if (tag2Indx < 0) return;

  TLorentzVector M3;
  float probeIndx = -1;
  for (unsigned int indx = 0; indx < selectLooseMuonList.size(); ++indx) {
    const Muon& probe = selectLooseMuonList[indx];
    M3.SetPtEtaPhiE(probe.pt, probe.eta, probe.phi, probe.energy);
    if (AnaUtil::deltaR(M1, M3) < 0.5) continue;
    if (AnaUtil::deltaR(M2, M3) < 0.5) continue;
    if (probe.pt <= AnaUtil::cutValue(muonCutMap(), "ptProbe")) continue;
    probeIndx = indx;
    break;
  }
  if (probeIndx < 0) return;
  const Muon& probeMuon = selectLooseMuonList[probeIndx];

  const MET& mt = metColl()->at(0);
  if (mt.met >=20.0) return;

  if (vetoMuon_wPt(10.0) != 3) return; //No Extra Muons
  if (vetoElectron_wPt(10.0) != 0) return; // No Extra Electrons
  if (bjetList.size() > 0) return; //No bJets

  AnaUtil::fillHist1D("testPt", probeMuon.pt, puevWt_);
  if (selectTauList.size() > 0 ){
    if ((selectTauList[0].pt > 20.0) && (selectTauList[0].chargedIsoPtSum < 1.0)) return; //No Selected Taus
  }

  float ProbeMet_MT = sqrt(2*probeMuon.pt*mt.met*(1-cos(AnaUtil::deltaPhi(probeMuon.phi, mt.metphi))));
  //if (ProbeMet_MT > 20.0) return;  // To supress WZ contamination

  //Fill Here the isolation denominator and numerator
  AnaUtil::fillHist1D("JetToMuFake_ZJets_mmt_deno", probeMuon.pt, puevWt_);  


  //Check if Probe Muon Passes the Tight ID
  if (fabs(probeMuon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
#if 0
     (!probeMuon.isTrackerMuon) ||
     (!probeMuon.isGlobalMuonPromptTight) ||
     (probeMuon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) ||
     (probeMuon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches")) ||
     (probeMuon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) ||
     (probeMuon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) ||
     (probeMuon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) ||
     (probeMuon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) ||
     (fabs(probeMuon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0"))
#endif
     (!probeMuon.isMediumMuon)
  ) return;

  if (fabs(probeMuon.eta) < 1.479){
    if (probeMuon.pfRelIsoDB03 < 0.15) AnaUtil::fillHist1D("JetToMuFake_ZJets_mmt_nume_lead", probeMuon.pt, puevWt_);
    if (probeMuon.pfRelIsoDB03 < 0.20) AnaUtil::fillHist1D("JetToMuFake_ZJets_mmt_nume_sublead", probeMuon.pt, puevWt_);
  }
  else{
    if (probeMuon.pfRelIsoDB03 < 0.10) AnaUtil::fillHist1D("JetToMuFake_ZJets_mmt_nume_lead", probeMuon.pt, puevWt_);
    if (probeMuon.pfRelIsoDB03 < 0.15) AnaUtil::fillHist1D("JetToMuFake_ZJets_mmt_nume_sublead", probeMuon.pt, puevWt_);
  }
}

void FakeRate::MuMuTauApplyFR()
{
  double vz = vtxList.at(0).z;

  vector<Muon> fMuoList;
  for (unsigned int indx = 0; indx < selectLooseMuonList.size(); ++indx) {
    const Muon& muon = selectLooseMuonList[indx];

    if (muon.pt <= 10) continue; 
    if (!muon.isMediumMuon) continue; 
    //if (fabs(muon.chargedHadronIso/muon.pt) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso")) continue;

    fMuoList.push_back(muon);
  }
  // atleast 2 good muons
  if (fMuoList.size() < 2) return;
  const Muon& muo1 = fMuoList.at(0);
  if (muo1.pt <= 22.0) return;
  const Muon& muo2 = fMuoList.at(1);

  //SS Muons
  if ((muo1.charge + muo2.charge) == 0) return;

  TLorentzVector M1, M2, T1;
  M1.SetPtEtaPhiE(muo1.pt, muo1.eta, muo1.phi, muo1.energy);
  M2.SetPtEtaPhiE(muo2.pt, muo2.eta, muo2.phi, muo2.energy);

  vector<Tau> fTauList;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& tau = selectTauList[indx];
    if (tau.pt <= 0.0) continue;

    T1.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
    float dr1 = AnaUtil::deltaR(M1, T1);
    float dr2 = AnaUtil::deltaR(M2, T1);
    if ((dr1 < AnaUtil::cutValue(_tau1CutMap, "drMuTau")) || (dr2 < AnaUtil::cutValue(_tau1CutMap, "drMuTau"))) continue;
    if (tau.pt <= AnaUtil::cutValue(_tau1CutMap, "pt") || 
        fabs(tau.eta) >= AnaUtil::cutValue(_tau1CutMap, "eta")) continue;  

    //if (tau.chargedIsoPtSum >= 2.0) continue;                     
    if (fabs(tau.zvertex - vz) >= AnaUtil::cutValue(_tau1CutMap, "dz")) continue;   
    if ((muo1.charge + tau.charge) != 0) continue;
    fTauList.push_back(tau);
  }
  
  if (fTauList.size() < 1) return; 
  const Tau& tau = fTauList.at(0);
  T1.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);

  // Here comes MET
  //const MET& mvamet = metColl()->at(0);
  const MET& pfmet = metColl()->at(0);

  const Tau& taua = fTauList.at(0);
  bool tauIso = (taua.byMediumIsolationMVArun2v1DBoldDMwLT > 0.5);

  bool isGoodVtx;
  TVector3 vtxmu = findLeptonVtx(muo1.vtxIndex, isGoodVtx);
  double mutaudz = abs(taua.zvertex - vtxmu.z());
  if (!isGoodVtx || mutaudz >= 0.14) return;
  if (vetoMuon_wPt(10.0) != 2) return;
  if (vetoElectron_wPt(10.0) != 0) return;
  if (bjetList.size() > 0) return; 

  float sumPt = M1.Pt() + M2.Pt() + T1.Pt();
  //if (sumPt <= 80) return;
  
  bool Signal = (
                 (fabs(muo1.eta) <  1.479 && (muo1.pfRelIsoDB03) < 0.15) || 
                 (fabs(muo1.eta) >= 1.479 && (muo1.pfRelIsoDB03) < 0.10)
                ) &&
                (
                 (fabs(muo2.eta) <  1.479 && (muo2.pfRelIsoDB03) < 0.20) || 
                 (fabs(muo2.eta) >= 1.479 && (muo2.pfRelIsoDB03) < 0.15)
                );

  bool leadFakeCR = (
                 (fabs(muo1.eta) <  1.479 && (muo1.pfRelIsoDB03) >= 0.15) || 
                 (fabs(muo1.eta) >= 1.479 && (muo1.pfRelIsoDB03) >= 0.10)
                ) &&
                (
                 (fabs(muo2.eta) <  1.479 && (muo2.pfRelIsoDB03) < 0.20) || 
                 (fabs(muo2.eta) >= 1.479 && (muo2.pfRelIsoDB03) < 0.15)
                );

  bool slFakeCR = (
                 (fabs(muo1.eta) <  1.479 && (muo1.pfRelIsoDB03) < 0.15) || 
                 (fabs(muo1.eta) >= 1.479 && (muo1.pfRelIsoDB03) < 0.10)
                ) &&
                (
                 (fabs(muo2.eta) <  1.479 && (muo2.pfRelIsoDB03) >= 0.20) || 
                 (fabs(muo2.eta) >= 1.479 && (muo2.pfRelIsoDB03) >= 0.15)
                );


  bool bothFakeCR = (
                 (fabs(muo1.eta) <  1.479 && (muo1.pfRelIsoDB03) >= 0.15) || 
                 (fabs(muo1.eta) >= 1.479 && (muo1.pfRelIsoDB03) >= 0.10)
                ) &&
                (
                 (fabs(muo2.eta) <  1.479 && (muo2.pfRelIsoDB03) >= 0.20) || 
                 (fabs(muo2.eta) >= 1.479 && (muo2.pfRelIsoDB03) >= 0.15)
                );

  TLorentzVector Higgs;
  Higgs = M2+T1;
  if (Signal && tauIso){
    if (isMC()) AnaUtil::fillHist1D("eventWeight_BL_mmt", 1.0, eventWt_*puevWt_);
    AnaUtil::fillHist1D("mtMass_mmt_signal", Higgs.M(), puevWt_);
    AnaUtil::fillHist1D("mu1Pt_mmt_signal", M1.Pt(), puevWt_);
    AnaUtil::fillHist1D("mu1Eta_mmt_signal", M1.Eta(), puevWt_);
    AnaUtil::fillHist1D("mu1Phi_mmt_signal", M1.Phi(), puevWt_);
    AnaUtil::fillHist1D("mu2Pt_mmt_signal", M2.Pt(), puevWt_);
    AnaUtil::fillHist1D("mu2Eta_mmt_signal", M2.Eta(), puevWt_);
    AnaUtil::fillHist1D("mu2Phi_mmt_signal", M2.Phi(), puevWt_);
    AnaUtil::fillHist1D("tauPt_mmt_signal", T1.Pt(), puevWt_);
    AnaUtil::fillHist1D("tauEta_mmt_signal", T1.Eta(), puevWt_);
    AnaUtil::fillHist1D("tauPhi_mmt_signal", T1.Phi(), puevWt_);
  }
  if (leadFakeCR && tauIso){
    if (isMC()) AnaUtil::fillHist1D("eventWeight_leadFakeCR_mmt", 1.0, eventWt_*puevWt_);
    AnaUtil::fillHist1D("mu1Pt_leadFakeCR", M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mtMass_mmt_leadFakeCR_p", Higgs.M(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Pt_mmt_leadFakeCR_p", M1.Pt(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Eta_mmt_leadFakeCR_p", M1.Eta(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Phi_mmt_leadFakeCR_p", M1.Phi(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Pt_mmt_leadFakeCR_p", M2.Pt(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Eta_mmt_leadFakeCR_p", M2.Eta(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Phi_mmt_leadFakeCR_p", M2.Phi(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauPt_mmt_leadFakeCR_p", T1.Pt(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauEta_mmt_leadFakeCR_p", T1.Eta(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauPhi_mmt_leadFakeCR_p", T1.Phi(), M1.Pt(), puevWt_);
  }
  if (slFakeCR && tauIso){
    if (isMC()) AnaUtil::fillHist1D("eventWeight_slFakeCR_mmt", 1.0, eventWt_*puevWt_);
    AnaUtil::fillHist1D("mu2Pt_slFakeCR", M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mtMass_mmt_slFakeCR_p", Higgs.M(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Pt_mmt_slFakeCR_p", M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Eta_mmt_slFakeCR_p", M1.Eta(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Phi_mmt_slFakeCR_p", M1.Phi(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Pt_mmt_slFakeCR_p", M2.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Eta_mmt_slFakeCR_p", M2.Eta(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Phi_mmt_slFakeCR_p", M2.Phi(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauPt_mmt_slFakeCR_p", T1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauEta_mmt_slFakeCR_p", T1.Eta(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauPhi_mmt_slFakeCR_p", T1.Phi(), M2.Pt(), puevWt_);
  }
  if (bothFakeCR && tauIso){
    if (isMC()) AnaUtil::fillHist1D("eventWeight_bothFakeCR_mmt", 1.0, eventWt_*puevWt_);
    AnaUtil::fillHist1D("mu1Pt_bothFakeCR", M1.Pt(), puevWt_);
    AnaUtil::fillHist1D("mu2Pt_bothFakeCR", M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Pt_mu2Pt_bothFakeCR", M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mtMass_mmt_bothFakeCR_p", Higgs.M(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu1Pt_mmt_bothFakeCR_p", M1.Pt(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu1Eta_mmt_bothFakeCR_p", M1.Eta(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu1Phi_mmt_bothFakeCR_p", M1.Phi(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu2Pt_mmt_bothFakeCR_p", M2.Pt(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu2Eta_mmt_bothFakeCR_p", M2.Eta(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu2Phi_mmt_bothFakeCR_p", M2.Phi(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("tauPt_mmt_bothFakeCR_p", T1.Pt(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("tauEta_mmt_bothFakeCR_p", T1.Eta(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("tauPhi_mmt_bothFakeCR_p", T1.Phi(), M1.Pt(), M2.Pt(), puevWt_);

  }

  if (Signal && !tauIso){
    AnaUtil::fillHist1D("mtMass_mmt_F3", Higgs.M(), puevWt_);
    AnaUtil::fillHist1D("mu1Pt_mmt_F3", M1.Pt(), puevWt_);
    AnaUtil::fillHist1D("mu1Eta_mmt_F3", M1.Eta(), puevWt_);
    AnaUtil::fillHist1D("mu1Phi_mmt_F3", M1.Phi(), puevWt_);
    AnaUtil::fillHist1D("mu2Pt_mmt_F3", M2.Pt(), puevWt_);
    AnaUtil::fillHist1D("mu2Eta_mmt_F3", M2.Eta(), puevWt_);
    AnaUtil::fillHist1D("mu2Phi_mmt_F3", M2.Phi(), puevWt_);
    AnaUtil::fillHist1D("tauPt_mmt_F3", T1.Pt(), puevWt_);
    AnaUtil::fillHist1D("tauEta_mmt_F3", T1.Eta(), puevWt_);
    AnaUtil::fillHist1D("tauPhi_mmt_F3", T1.Phi(), puevWt_);
  }
  if (leadFakeCR && !tauIso){
    AnaUtil::fillHist1D("mu1Pt_leadFakeCR_F3", M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mtMass_mmt_leadFakeCR_F3_p", Higgs.M(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Pt_mmt_leadFakeCR_F3_p", M1.Pt(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Eta_mmt_leadFakeCR_F3_p", M1.Eta(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Phi_mmt_leadFakeCR_F3_p", M1.Phi(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Pt_mmt_leadFakeCR_F3_p", M2.Pt(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Eta_mmt_leadFakeCR_F3_p", M2.Eta(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Phi_mmt_leadFakeCR_F3_p", M2.Phi(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauPt_mmt_leadFakeCR_F3_p", T1.Pt(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauEta_mmt_leadFakeCR_F3_p", T1.Eta(), M1.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauPhi_mmt_leadFakeCR_F3_p", T1.Phi(), M1.Pt(), puevWt_);

  }
  if (slFakeCR && !tauIso){
    AnaUtil::fillHist1D("mu2Pt_slFakeCR_F3", M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mtMass_mmt_slFakeCR_F3_p", Higgs.M(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Pt_mmt_slFakeCR_F3_p", M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Eta_mmt_slFakeCR_F3_p", M1.Eta(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Phi_mmt_slFakeCR_F3_p", M1.Phi(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Pt_mmt_slFakeCR_F3_p", M2.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Eta_mmt_slFakeCR_F3_p", M2.Eta(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu2Phi_mmt_slFakeCR_F3_p", M2.Phi(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauPt_mmt_slFakeCR_F3_p", T1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauEta_mmt_slFakeCR_F3_p", T1.Eta(), M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauPhi_mmt_slFakeCR_F3_p", T1.Phi(), M2.Pt(), puevWt_);

  }
  if (bothFakeCR && !tauIso){
    AnaUtil::fillHist1D("mu1Pt_bothFakeCR_F3", M1.Pt(), puevWt_);
    AnaUtil::fillHist1D("mu2Pt_bothFakeCR_F3", M2.Pt(), puevWt_);
    AnaUtil::fillHist2D("mu1Pt_mu2Pt_bothFakeCR_F3", M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mtMass_mmt_bothFakeCR_F3_p", Higgs.M(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu1Pt_mmt_bothFakeCR_F3_p", M1.Pt(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu1Eta_mmt_bothFakeCR_F3_p", M1.Eta(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu1Phi_mmt_bothFakeCR_F3_p", M1.Phi(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu2Pt_mmt_bothFakeCR_F3_p", M2.Pt(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu2Eta_mmt_bothFakeCR_F3_p", M2.Eta(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("mu2Phi_mmt_bothFakeCR_F3_p", M2.Phi(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("tauPt_mmt_bothFakeCR_F3_p", T1.Pt(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("tauEta_mmt_bothFakeCR_F3_p", T1.Eta(), M1.Pt(), M2.Pt(), puevWt_);
    AnaUtil::fillHist3D("tauPhi_mmt_bothFakeCR_F3_p", T1.Phi(), M1.Pt(), M2.Pt(), puevWt_);

  }
}
void FakeRate::EleMuTauFakeWJets()
{
  if (selectTightElectronList.size() < 1) return;
  TLorentzVector E1;
  const Electron& tagElectron = selectTightElectronList[0];

  TLorentzVector M2;
  float probeIndx = -1;
  for (unsigned int indx = 0; indx < selectLooseMuonList.size(); ++indx) {
    const Muon& probe = selectLooseMuonList[indx];
    M2.SetPtEtaPhiE(probe.pt, probe.eta, probe.phi, probe.energy);
    if (AnaUtil::deltaR(E1, M2) < 0.5) continue;
    //if ((tagMuon.charge + probe.charge) == 0) continue;  //To avoid Z/Gamma jets//Do Not need for WJets MC
    if (probe.pt <= AnaUtil::cutValue(muonCutMap(), "ptProbe")) continue;
    probeIndx = indx;
    break;
  }
  if (probeIndx < 0) return;
  const Muon& probeMuon = selectLooseMuonList[probeIndx];

  const MET& mt = metColl()->at(0);
  float TagMet_MT = sqrt(2*tagElectron.pt*mt.met*(1-cos(AnaUtil::deltaPhi(tagElectron.phi, mt.metphi))));
  AnaUtil::fillHist1D("TagMet_MT_emt_WJets", TagMet_MT, puevWt_);
  if (TagMet_MT < 35.0) return;

  float ProbeMet_MT = sqrt(2*probeMuon.pt*mt.met*(1-cos(AnaUtil::deltaPhi(probeMuon.phi, mt.metphi))));
  if (ProbeMet_MT > 35.0) return;

  if (vetoMuon_wPt(10.0) != 1) return; //No Extra Muons
  if (vetoElectron_wPt(10.0) != 1) return; // No Extra Electrons
  if (bjetList.size() > 0) return; //No bJets
  if (selectTauList.size() > 0 ){
    if ((selectTauList[0].pt > 20.0) && (selectTauList[0].chargedIsoPtSum < 1.0)) return; //No Selected Taus
  }
  //Require 1 Jet to mimic the presence of a Tau object
  float nJet = 0;
  TLorentzVector J1;
  bool JetFlag = false;
  for (auto it = jetColl()->begin(); it != jetColl()->end(); ++it) {
    const Jet& jt = (*it);
    if (fabs(jt.eta) >= 2.4 || jt.pt <= 20.0) continue;
    nJet += 1;
    J1.SetPtEtaPhiE(jt.pt, jt.eta, jt.phi, jt.energy);
    if (AnaUtil::deltaR(E1, J1) < 0.5 || AnaUtil::deltaR(M2, J1) < 0.5) continue;
    JetFlag = true;
    //break; //just need 1 non-overlapping jet
  }
  if (!JetFlag) return;

  //Fill Here the isolation denominator and numerator
  AnaUtil::fillHist1D("JetToMuFake_WJets_emt_deno", probeMuon.pt, puevWt_);  

  //Check if Probe Muon Passes the Tight ID
  bool probeTID = true;
  if (fabs(probeMuon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
     (!probeMuon.isTrackerMuon) ||
     (!probeMuon.isGlobalMuonPromptTight) ||
     (probeMuon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) ||
     (probeMuon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches")) ||
     (probeMuon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) ||
     (probeMuon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) ||
     (probeMuon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) ||
     (probeMuon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) ||
     (fabs(probeMuon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0"))
  ) probeTID = false;

  if (fabs(probeMuon.eta) < 1.479){
    if ((probeMuon.chargedHadronIso/probeMuon.pt) < 0.15 && probeTID) AnaUtil::fillHist1D("JetToMuFake_WJets_emt_nume", probeMuon.pt, puevWt_);
  }
  else{
    if ((probeMuon.chargedHadronIso/probeMuon.pt) < 0.10 && probeTID) AnaUtil::fillHist1D("JetToMuFake_WJets_emt_nume", probeMuon.pt, puevWt_);
  }

  //Fill histogram for kNN training
  //skimForKNN(probeMuon, probeTID, nJet);
}

//MOre like EleEleTau scenario, just need this for kNN skim in Z->ee region
void FakeRate::EleMuTauFakeZJets()
{
  if (selectTightElectronList.size() < 2) return;
  TLorentzVector E1, E2, Mass;
  const Electron& tag1 = selectTightElectronList[0];
  const Electron& tag2 = selectTightElectronList[1];
  if ((tag1.charge + tag2.charge) != 0) return;
  E1.SetPtEtaPhiE(tag1.pt, tag1.eta, tag1.phi, tag1.energy);
  E2.SetPtEtaPhiE(tag2.pt, tag2.eta, tag2.phi, tag2.energy);

  Mass = E1+E2;
  if (Mass.M() < 70 || Mass.M() > 110) return;

  TLorentzVector M3;
  float probeIndx = -1;
  for (unsigned int indx = 0; indx < selectLooseMuonList.size(); ++indx) {
    const Muon& probe = selectLooseMuonList[indx];
    M3.SetPtEtaPhiE(probe.pt, probe.eta, probe.phi, probe.energy);
    if (AnaUtil::deltaR(E1, M3) < 0.5) continue;
    if (AnaUtil::deltaR(E2, M3) < 0.5) continue;
    if (probe.pt <= AnaUtil::cutValue(muonCutMap(), "ptProbe")) continue;
    probeIndx = indx;
    break;
  }
  if (probeIndx < 0) return;
  const Muon& probeMuon = selectLooseMuonList[probeIndx];

  const MET& mt = metColl()->at(0);
  if (mt.met >=20.0) return;

  if (vetoMuon_wPt(10.0) != 1) return; //No Extra Muons
  if (vetoElectron_wPt(10.0) != 2) return; // No Extra Electrons
  if (bjetList.size() > 0) return; //No bJets

  if (selectTauList.size() > 0 ){
    if ((selectTauList[0].pt > 20.0) && (selectTauList[0].chargedIsoPtSum < 1.0)) return; //No Selected Taus
  }

  float ProbeMet_MT = sqrt(2*probeMuon.pt*mt.met*(1-cos(AnaUtil::deltaPhi(probeMuon.phi, mt.metphi))));
  //if (ProbeMet_MT > 20.0) return;  // To supress WZ contamination

  //Fill Here the isolation denominator and numerator
  AnaUtil::fillHist1D("JetToMuFake_ZJets_emt_deno", probeMuon.pt, puevWt_);  

  //Check if Probe Muon Passes the Tight ID
  if (fabs(probeMuon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") ||
     (!probeMuon.isTrackerMuon) ||
     (!probeMuon.isGlobalMuonPromptTight) ||
     (probeMuon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) ||
     (probeMuon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches")) ||
     (probeMuon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) ||
     (probeMuon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) ||
     (probeMuon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) ||
     (probeMuon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) ||
     (fabs(probeMuon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0"))
  ) return;

  if (fabs(probeMuon.eta) < 1.479){
    if ((probeMuon.chargedHadronIso/probeMuon.pt) < 0.15) AnaUtil::fillHist1D("JetToMuFake_ZJets_emt_nume", probeMuon.pt, puevWt_);
  }
  else{
    if ((probeMuon.chargedHadronIso/probeMuon.pt) < 0.10) AnaUtil::fillHist1D("JetToMuFake_ZJets_emt_nume", probeMuon.pt, puevWt_);
  }
  //Fill histogram for kNN training
  //skimForKNN(probeMuon, probeTID, nJet);
}

void FakeRate::skimForKNN(const Muon& probeMuon, const bool& probeTID, const float& nJet)
{
  //Fill histogram for kNN training
  if (fabs(probeMuon.eta) < 1.479){
    if ((probeMuon.chargedHadronIso/probeMuon.pt) < 0.20 && probeTID){
      if (_skimSig) {
        SigVariables varList;
        varList.muEta      = probeMuon.eta;
        varList.muPt       = probeMuon.pt;
        varList.nJet       = nJet;
        _skimSig->fill(varList);
      }
    }
    else{
      if (_skimBkg) {
        BkgVariables varList;
        varList.muEta      = probeMuon.eta;
        varList.muPt       = probeMuon.pt;
        varList.nJet       = nJet;
        _skimBkg->fill(varList);
      }
    }
  }
  else{
    if ((probeMuon.chargedHadronIso/probeMuon.pt) < 0.15 && probeTID){
      if (_skimSig) {
        SigVariables varList;
        varList.muEta      = probeMuon.eta;
        varList.muPt       = probeMuon.pt;
        varList.nJet       = nJet;
        _skimSig->fill(varList);
      }
    }
    else{
      if (_skimBkg) {
        BkgVariables varList;
        varList.muEta      = probeMuon.eta;
        varList.muPt       = probeMuon.pt;
        varList.nJet       = nJet;
        _skimBkg->fill(varList);
      }
    }
  }
}

// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void FakeRate::endJob() 
{  
  histf()->cd();

  histf()->Write();
  histf()->Close();
  delete histf();

  fLog() << resetiosflags(ios::fixed);
  if (_skimSig) _skimSig->close();
  if (_skimBkg) _skimBkg->close();
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
bool FakeRate::readJob(const string& jobFile, int& nFiles)
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
    else if (key == "readMVAFK")
      _readMVAFK = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "readMVA")
      _readMVA = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "MVAdisFile")
      _MVAdisFile = value;
    else if (key == "MVAFKdisFile")
      _MVAFKdisFile = value;
    else if (key == "createMVATree")
      _createMVATree = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "mvaInputFileSigkNN")
      _mvaInputFileSigkNN = value;
    else if (key == "mvaInputFileBkgkNN")
      _mvaInputFileBkgkNN = value;


    tokens.clear();
  }
  // Close the file
  fin.close();

  printJob();

  return true;
}
void FakeRate::printJob(ostream& os) const
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("tau1CutList", _tau1CutMap));
  hmap.insert(pair<string, map<string, double> >("tau2CutList", _tau2CutMap));
  AnaUtil::showCuts(hmap, os);
}

