#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TAxis.h>
#include <TMath.h>
#include <TSystem.h>
#include <TMath.h>
#include <TString.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include<fstream>

void mc_closure()
{

TFile* f7 = new TFile("WJets_reweighted.root");
TH1F* W = (TH1F*)f7->Get("faketaupt_nume_WJetsAll_SSSS");
cout<<"Get no of bins"<<W->GetNbinsX()<<endl;
Int_t bins_no = W->GetNbinsX();
TFile* f10 = new TFile("WJets.root");
TH1F *Wiso = (TH1F*)f10->Get("tautauMass_CR");
//(61526.7 * 35900)/56871716.934994
Wiso->Scale(38.838435852);
//Wiso->SetFillColor(kGreen);
Wiso->Rebin(20);
//W1->Scale(1.0 / W1->Integral());




TH1F *Wfr = (TH1F*)f10->Get("tautauMass_AA");
Wfr->Scale(38.838435852);
Wfr->SetFillColor(kGreen);
Wfr->Rebin(20);
//W1->Scale(1.0 / W1->Integral());

for(int i=1 ; i <= Wfr->GetNbinsX() ; i++ )
{
Double_t weight_WJets = W->GetBinContent(i);
Double_t weight = ( weight_WJets/(1- weight_WJets));
Double_t entry_Wjets = Wfr->GetBinContent(i);
Wfr->SetBinContent(i, entry_Wjets * weight);
}
TCanvas* c1 = new TCanvas();
Wfr->Draw("hist");
Wiso->Draw("E1SAME");

TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
legend->SetTextFont(72);
legend->SetTextSize(0.02);
legend->AddEntry(Wfr,"both tau anti-iso","f");
legend->AddEntry(Wiso,"Tau OS to #mu AntiIso", "lep");
legend->Draw();
}
