
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

void nPV()
{
THStack *hs = new THStack("hs"," m_{#tau#tau} in  OSAntiSSAnti ");	
TFile* f = TFile::Open("ZZ.root");
TH1F *h = (TH1F*)f->Get("tautauMass_AA");
h->SetFillColor(40);
//(1.212 * 27260)/6671074.119412
// lumi 27.26
h->Scale(0.004952594);
h->Rebin(20);
//h->Scale(1.0 / h->Integral());
hs->Add(h);

TFile* f1 = TFile::Open("WZ.root");
TH1F *h2 = (TH1F*)f1->Get("tautauMass_AA");
h2->SetFillColor(48);
//(5.595 * 27260)/26987319.630721
//19720 in pb inverse is luminosity for the data B,C,D,E,F
h2->Scale(0.005651532);
h2->Rebin(20);
//h2->Scale(1.0 / h2->Integral());
hs->Add(h2);

TFile* f2 = TFile::Open("WW.root");
TH1F *h3 = (TH1F*)f2->Get("tautauMass_AA");
h3->SetFillColor(30);
//(49.997 * 27260)/1998519.872133
h3->Scale(0.681963807);
h3->Rebin(20);
//h3->Scale(1.0 / h3->Integral());
hs->Add(h3);

/*
TFile* Q = new TFile("QCD.root");
TH1F *q1 = (TH1F*)Q->Get("tautauMass_AA");
//((720648000  * 27260)/22096073.028419) * 0.00042   for Sin mu B only 
q1->Scale(373.407667099);
q1->SetFillColor(kRed);
q1->Rebin(20);
//q1->Scale(1.0 / q1->Integral());
hs->Add(q1);

*/

TFile* f11 = TFile::Open("TTbar.root");
TH1F *h11 = (TH1F*)f11->Get("tautauMass_AA");
//(831.76 * 27260)/74821215.094177	
h11->Scale(0.30303942);
h11->SetFillColor(kBlue-3);
h11->Rebin(20);
//h11->Scale(1.0 / h11->Integral());
hs->Add(h11);



TFile* g1 = new TFile("DY_all.root");
TH1F *D1 = (TH1F*)g1->Get("tautauMass_AA");
//(5765.4 * 27260)/46751905.072229  fontr Sin mu B only 
D1->Scale(3.361677);
D1->SetFillColor(kYellow);
D1->Rebin(20);
//D1->Scale(1.0 / D1->Integral());
hs->Add(D1);



TFile* f10 = new TFile("WJets.root");
TH1F *W1 = (TH1F*)f10->Get("tautauMass_AA");
//(61526.7 * 27260)/56871716.934994
W1->Scale(29.491246834);
W1->SetFillColor(kGreen);
W1->Rebin(20);
//W1->Scale(1.0 / W1->Integral());
hs->Add(W1);



/*
TFile* f6 = TFile::Open("SingleMuon.root");
TH1F *h7 = (TH1F*)f6->Get("tautauMass_AA");
h7->SetMarkerColor(kBlack);
h7->SetMarkerStyle(20);
h7->SetMarkerSize(0.8);
//h7->Scale(1.0 / h7->Integral());
h7->Rebin(20);
*/

TCanvas* c2 = new TCanvas();
/*
TCanvas* c = new TCanvas("c","canvas",800,800);

TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
pad1->SetBottomMargin(0); 
//pad1->SetGridx();         
pad1->Draw();             
pad1->cd(); 

*/      
hs->Draw("hist");
hs->GetXaxis()->SetTitleOffset(1.2);
hs->GetXaxis()->SetTitle("m_{#tau#tau}");
hs->GetYaxis()->SetTitle("Events");


//hs->SetMaximum(50);


// Area Under the Curve 

Double_t area1 = D1->Integral(0,300);
Double_t area2 = W1->Integral(0,300);

cout<<"Area under the curve"<<" DY "<< area1 <<" under WJets "<<area2 <<endl;


//hs->GetXaxis()->SetRangeUser(0,200);
//hs->GetYaxis()->SetRangeUser(0,15);
//hs->GetYaxis()->SetLabelSize(0.);

//h7->Draw("E1SAME");

// draw the legend
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

//legend->AddEntry(h7,"Sin Muon Data","lep");
legend->Draw();




/*


c->cd();

TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);

pad2->SetTopMargin(0.1);
pad2->SetBottomMargin(0.2);
pad2->SetGridx(); // vertical grid
pad2->Draw();
pad2->cd(); 

// Define the ratio plot (data-mc)/mc
TH1F *F1 = (TH1F*)h7->Clone("F1");
F1->SetLineColor(kBlack);
F1->SetMinimum(0.5);  // Define Y ..
F1->SetMaximum(1.5); // .. range
F1->Sumw2();    
F1->SetStats(0);


h->Add(h2);
h->Add(h3);
h->Add(q1);
h->Add(h11);
h->Add(W1);
h->Add(D1);


F1->Divide(h);
F1->SetMarkerStyle(20);
F1->Draw("ep");
F1->SetTitle("");

// Y axis ratio plot setting 
   F1->GetYaxis()->SetTitle("Data/MC");
  // F1->GetYaxis()->SetNdivisions(505);
   F1->GetYaxis()->SetTitleSize(18);
   F1->GetYaxis()->SetTitleFont(43);
   //F1->GetYaxis()->SetTitleOffset(1.55);
   F1->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   F1->GetYaxis()->SetLabelSize(15);


// X axis ratio plot settings
   F1->GetXaxis()->SetTitleSize(20);
   F1->GetXaxis()->SetTitleFont(43);
   F1->GetXaxis()->SetTitleOffset(4.);
   F1->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   F1->GetXaxis()->SetLabelSize(15);

TLine *l=new TLine(0,1,200,1);
   l->Draw("SAME");
   pad2->Update();
   */

}