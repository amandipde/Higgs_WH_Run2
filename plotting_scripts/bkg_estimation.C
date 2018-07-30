#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"


Double_t weightfunZ ( Double_t x)
{  
return  ( (0.2189 - (0.0001022 * x)) );           
}

Double_t weightfunW ( Double_t x)
{

return  (0.1932 - (0.0003106 * x)); 

}

Double_t weightfun(Double_t x)
{
return (0.2 * (weightfunZ(x)/(1-weightfunZ(x))) ) +  ( 0.8 *(weightfunW(x)/(1-weightfunW(x)) ) ) ;

}

void bkg_estimation()
{

THStack *hs = new THStack("hs"," M_{#tau#tau} in Signal region");


TFile* f = TFile::Open("ZZ.root");
TH1F *h = (TH1F*)f->Get("tautauMass_BL");
h->SetFillColor(40);
//(1.212 * 35900)/6671074.119412
// lumi 35.9
h->Scale(0.006522308);
h->Rebin(20);
//h->Scale(1.0 / h->Integral());
hs->Add(h);

TFile* f1 = TFile::Open("WZ.root");
TH1F *h2 = (TH1F*)f1->Get("tautauMass_BL");
h2->SetFillColor(48);
//(5.595 * 35900)/26987319.630721
//19720 in pb inverse is luminosity for the data B,C,D,E,F
h2->Scale(0.007442773);
h2->Rebin(20);
//h2->Scale(1.0 / h2->Integral());
hs->Add(h2);


TFile* f2 = TFile::Open("WW.root");
TH1F *h3 = (TH1F*)f2->Get("tautauMass_AA");
h3->SetFillColor(30);
//(49.997 * 35900)/1998519.872133
h3->Scale(0.898110809);
h3->Rebin(20);
//h3->Scale(1.0 / h3->Integral());


TFile* f11 = TFile::Open("TTbar.root");
TH1F *h11 = (TH1F*)f11->Get("tautauMass_AA");
//(831.76 * 35900)/74821215.094177	
h11->Scale(0.39908713);
h11->SetFillColor(kBlue-3);
h11->Rebin(20);
//h11->Scale(1.0 / h11->Integral());




 
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
//(5765.4 * 35900)/46751905.072229  fontr Sin mu B only 
D1->Scale(4.427153496);
D1->SetFillColor(kYellow);
D1->Rebin(20);
//D1->Scale(1.0 / D1->Integral());



TFile* f10 = new TFile("WJets.root");
TH1F *W1 = (TH1F*)f10->Get("tautauMass_AA");
//(61526.7 * 35900)/56871716.934994
W1->Scale(38.838435852);
W1->SetFillColor(kGreen);
W1->Rebin(20);
//W1->Scale(1.0 / W1->Integral());



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


Double_t entry_ttbar = h11->GetBinContent(i);
Double_t entry_X_ttbar = h11->GetXaxis()->GetBinCenter(i);
Double_t weight_ttbar = weightfun(entry_X_ttbar);
h11->SetBinContent(i, entry_ttbar * weight_ttbar);



}


TFile* f6 = TFile::Open("SingleMuon_all.root");
TH1F *h7 = (TH1F*)f6->Get("tautauMass_BL");
h7->SetMarkerColor(kBlack);
h7->SetMarkerStyle(20);
h7->SetMarkerSize(0.8);
h7->Rebin(20);



hs->Add(h3);
hs->Add(h11);
hs->Add(D1);
hs->Add(W1);

TCanvas* c2 = new TCanvas();

hs->SetMaximum(120);
	
hs->Draw("hist");
h7->Draw("E1SAME");




TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
legend->SetTextFont(72);
legend->SetTextSize(0.02);
legend->AddEntry(h,"ZZ","f");
legend->AddEntry(h2,"WZ","f");
legend->AddEntry(h3,"WW","f");
//legend->AddEntry(q1,"QCD","f");
legend->AddEntry(h11,"ttbar","f");
legend->AddEntry(D1,"DY+Jets","f");
legend->AddEntry(W1,"Wjets","f");

legend->AddEntry(h7,"Sin Muon Data","lep");
legend->Draw();


}
