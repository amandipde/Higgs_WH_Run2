#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"


Double_t weightfunZ ( Double_t x)
{  
return  ( (0.2221 - (0.0001765 * x)) );           
}

Double_t weightfunW ( Double_t x)
{

return  (0.1644 + (0.000748 * x)); 

}

Double_t weightfun(Double_t x)
{
return (0.4 * (weightfunZ(x)/(1-weightfunZ(x))) ) +  ( (1-0.4) *(weightfunW(x)/(1-weightfunW(x)) ) ) ;

}







void bkg_estimation()
{

THStack *hs = new THStack("hs"," M_{#tau#tau} in Signal region");


TFile* f = TFile::Open("ZZ.root");
TH1F *h = (TH1F*)f->Get("tautauMass_BL");
h->SetFillColor(40);
//(1.212 * 27260)/6671074.119412
// lumi 27.26
h->Scale(0.004952594);
h->Rebin(20);
hs->Add(h);


TFile* f1 = TFile::Open("WZ.root");
TH1F *h2 = (TH1F*)f1->Get("tautauMass_BL");
h2->SetFillColor(48);
//(5.595 * 27260)/26987319.630721
//19720 in pb inverse is luminosity for the data B,C,D,E,F
h2->Scale(0.005651532);
h2->Rebin(20);
hs->Add(h2);

TFile* f2 = TFile::Open("WW.root");
TH1F *h3 = (TH1F*)f2->Get("tautauMass_AA");
h3->Rebin(20);
h3->SetFillColor(20);
//(49.997 * 27260)/1998519.872133
h3->Scale(0.681963807);

 
/*
TFile* Q = new TFile("QCD.root");
TH1F *q1 = (TH1F*)Q->Get("tautauMass_CR");
//((720648000  * 27260)/22096073.028419) * 0.00042   for Sin mu B only 
q1->Scale(373.407667099);
q1->SetFillColor(kRed);
q1->Rebin(20);
//q1->Scale(1.0 / q1->Integral());
hs->Add(q1);

*/

TFile* g1 = new TFile("DY_all.root");
TH1F *D1 = (TH1F*)g1->Get("tautauMass_AA");
D1->Scale(3.361677);
D1->Rebin(20);
D1->SetFillColor(kYellow);



TFile* f10 = new TFile("WJets.root");
TH1F *W1 = (TH1F*)f10->Get("tautauMass_AA");
W1->Scale(29.491246834);
W1->Rebin(20);
W1->SetFillColor(kGreen);



Int_t binsDY = D1->GetNbinsX();
Int_t binsWjets  = W1->GetNbinsX();

for(int i=1 ; i <= binsDY ; i++ )
{
Double_t entry_DY = D1->GetBinContent(i);
Double_t entry_X_DY = D1->GetXaxis()->GetBinCenter(i);
Double_t weight_DY = weightfun(entry_X_DY);
//cout <<"the bincontent-> "<< entry_DY<< " entry->  " << entry_DY <<" weight->  "  << weight_DY <<"entry_DY * weight-  "<<entry_DY * weight_DY<<endl;
D1->SetBinContent(i, entry_DY * weight_DY);

Double_t entry_Wjets = W1->GetBinContent(i);
Double_t entry_X_Wjets = W1->GetXaxis()->GetBinCenter(i);
Double_t weight_Wjets = weightfun(entry_X_Wjets);
W1->SetBinContent(i, entry_Wjets * weight_Wjets);

Double_t entry_WW = h3->GetBinContent(i);
Double_t entry_X_WW = h3->GetXaxis()->GetBinCenter(i);
Double_t weight_WW = weightfun(entry_X_WW);
h3->SetBinContent(i, entry_WW * weight_WW);

}


TFile* f6 = TFile::Open("SingleMuon.root");
TH1F *h7 = (TH1F*)f6->Get("tautauMass_BL");
h7->SetMarkerColor(kBlack);
h7->SetMarkerStyle(20);
h7->SetMarkerSize(0.8);
h7->Rebin(20);



hs->Add(h3);
hs->Add(D1);
hs->Add(W1);

TCanvas* c2 = new TCanvas();

hs->SetMaximum(100);
	
hs->Draw("hist");
h7->Draw("E1SAME");


}
