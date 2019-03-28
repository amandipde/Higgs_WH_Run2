#include <numeric>
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include  "math.h"
#include <vector>
void PuWeight()
{

std::vector<double> puBins_;
std::vector<double> puWeight_;
std::vector<double> mcpu;
TFile * f_data = new TFile("pudistributions_data_2017.root") ;
TH1D *pu_data = (TH1D*)f_data->Get("pileup");

int nbins_data = pu_data->GetXaxis()->GetNbins();
for (int i = 0; i < nbins_data; ++i) {

       double wt = pu_data->GetBinContent(i);
       puWeight_.push_back(wt);
}
TFile * f_mc = new TFile("pudistributions_mc_2017.root") ;
TH1D *pu_DY = (TH1D*)f_mc->Get("#DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8#RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1#MINIAODSIM");

int nbins_DY = pu_DY->GetXaxis()->GetNbins();
for (int i = 0; i < nbins_DY; ++i) {

       double wt = pu_DY->GetBinContent(i);
       mcpu.push_back(wt);
}
TH1D *reweight_DY = new TH1D("pileup","pileup distribution",100,0,100);
assert( puWeight_.size() == mcpu.size() );
 auto scl  = std::accumulate(mcpu.begin(),mcpu.end(),0.) / std::accumulate(puWeight_.begin(),puWeight_.end(),0.);  // rescale input distribs to unit area
for( size_t ib = 0; ib<puWeight_.size(); ++ib ) 
{ 
puWeight_[ib] *= scl / mcpu[ib];

//cout<< puWeight_[ib]<<endl;

reweight_DY->SetBinContent(ib,puWeight_[ib]);


}

TFile *Reweight_PU_DY = new TFile("Reweight_PU_DY_2017.root", "recreate");
reweight_DY->Write();
Reweight_PU_DY->Close();



}
