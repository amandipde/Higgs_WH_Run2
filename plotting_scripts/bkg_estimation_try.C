#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"



void bkg_estimation_try()
{

TFile* ff = new TFile("DY_reweighted.root");
TH1F* D = (TH1F*)ff->Get("faketaupt_nume_DYAll");


TFile* f7 = new TFile("WJets_reweighted.root");
TH1F* W = (TH1F*)f7->Get("faketaupt_nume_WJetsAll_SSSS");

cout<<"Get no of bins"<<W->GetNbinsX()<<endl;
Int_t bins_no = W->GetNbinsX();




THStack *hs = new THStack("hs"," M_{#tau#tau} in Signal region");


TFile* f = TFile::Open("ZZ.root");
TH1F *h = (TH1F*)f->Get("tautauMass_BL");
h->SetFillColor(40);
h->Scale(0.006522308);
h->Rebin(20);
hs->Add(h);


TFile* f1 = TFile::Open("WZ.root");
TH1F *h2 = (TH1F*)f1->Get("tautauMass_BL");
h2->SetFillColor(48);
h2->Scale(0.007442773);
h2->Rebin(20);
hs->Add(h2);

TFile* f2 = TFile::Open("WW.root");
TH1F *h3 = (TH1F*)f2->Get("tautauMass_AA");
h3->Rebin(20);
h3->SetFillColor(20);
h3->Scale(0.898110809);


TFile* f11 = TFile::Open("TTbar.root");
TH1F *h11 = (TH1F*)f11->Get("tautauMass_BL");	
h11->Scale(0.39908713);
h11->SetFillColor(kBlue-3);
h11->Rebin(20);
//h11->Scale(1.0 / h11->Integral());




 

TFile* Q = new TFile("QCD.root");
TH1F *q1 = (TH1F*)Q->Get("tautauMass_BL"); 
q1->Scale(491.7584464);
q1->SetFillColor(kRed);
q1->Rebin(20);
//q1->Scale(1.0 / q1->Integral());
hs->Add(q1);



TFile* g1 = new TFile("DY_all.root");
TH1F *D1 = (TH1F*)g1->Get("tautauMass_AA");
D1->Scale(4.427153496);
D1->Rebin(20);
D1->SetFillColor(kYellow);



TFile* f10 = new TFile("WJets.root");
TH1F *W1 = (TH1F*)f10->Get("tautauMass_AA");
W1->Scale(38.838435852);
W1->Rebin(20);
W1->SetFillColor(kGreen);


for(int i=1 ; i <= bins_no ; i++ )
{

Double_t entry_DY = D1->GetBinContent(i);
Double_t weight_DY = D->GetBinContent(i);
Double_t weight_WJets = W->GetBinContent(i);

//Double_t entry_X_DY = D1->GetXaxis()->GetBinCenter(i);
Double_t weight = (0.2 * weight_DY/(1-weight_DY) ) + (0.8 * weight_WJets/(1- weight_WJets)) ;
cout <<"the bincontent-> "<< entry_DY<< " entry->  " << entry_DY <<" weight->  "  << weight_DY <<"entry_DY * weight-  "<<entry_DY * weight<<endl;
D1->SetBinContent(i, entry_DY * weight);


Double_t entry_Wjets = W1->GetBinContent(i);
W1->SetBinContent(i, entry_Wjets * weight);

Double_t entry_WW = h3->GetBinContent(i);
//Double_t entry_X_WW = h3->GetXaxis()->GetBinCenter(i);
//Double_t weight_WW = weightfun(entry_X_WW);
h3->SetBinContent(i, entry_WW * weight);


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

TCanvas* c1 = new TCanvas();

//hs->SetMaximum(300);
	
hs->Draw("hist");	
h7->Draw("E1SAME");




TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
legend->SetTextFont(72);
legend->SetTextSize(0.02);
legend->AddEntry(h,"ZZ","f");
legend->AddEntry(h2,"WZ","f");
legend->AddEntry(h3,"WW","f");
legend->AddEntry(q1,"QCD","f");
legend->AddEntry(h11,"ttbar","f");
legend->AddEntry(D1,"DY+Jets","f");
legend->AddEntry(W1,"Wjets","f");

legend->AddEntry(h7,"Sin Muon Data","lep");
legend->Draw();


}
