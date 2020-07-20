#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"

using namespace std;

TH2F* Multiply(TH2F* h1, TH2F* h2)
{
  TH2F* hproduct = 0;
  if( h1->GetNbinsX()>=h2->GetNbinsX() && h1->GetNbinsY()>=h2->GetNbinsY() )
    hproduct = (TH2F*)h1->Clone();
  else
    if( h2->GetNbinsX()>=h1->GetNbinsX() && h2->GetNbinsY()>=h1->GetNbinsY() )
      hproduct = (TH2F*)h2->Clone();
    else
    {
      cout<<"[ERROR]: cannot tell which is the denser-bin histogram"<<endl;
      return 0;
    }

  for(int ix=1; ix<hproduct->GetNbinsX()+1; ++ix)
  {
    float x = hproduct->GetXaxis()->GetBinCenter(ix);
    for(int iy=1; iy<hproduct->GetNbinsY()+1; ++iy)
    {
      float y = hproduct->GetYaxis()->GetBinCenter(iy);
      float v1 = h1->GetBinContent( h1->FindBin(x,y) );
      float v2 = h2->GetBinContent( h2->FindBin(x,y) );
      //cout<<"(x,y)=("<<x<<","<<y<<")\t v1="<<v1<<"\t v2="<<v2<<"product="<<v1*v2<<endl;
      hproduct->SetBinContent(ix,iy, v1*v2);
    }
  }
  return hproduct;
}

int main(int argc, char** argv)
{

  if(argc!= 4)
  {
    cout<<"usage: ReweightForBTagMTD.exe <oldfilename> <treename> <ReweightSetup>"<<endl;
    return -1;
  }
  TString oldfilename(argv[1]);
  TString treename(argv[2]);
  TString ReweightSetup(argv[3]);

  //Reweight for new btag efficiencies                                                                            
  TH2F* btaglooseeff_reweightmap = 0;
  float Ptmax_btag_reweightmap = -1;
  if(ReweightSetup=="DelphesMTD_DelphesNoMTD")
  {
    cout<<"loading btag reweight map from file /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/data/ScaleFactorsFor50ps.root"<<endl;
    TFile *btag_reweight_file = new TFile("/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/data/ScaleFactorsFor50ps.root");
    btaglooseeff_reweightmap = (TH2F*) btag_reweight_file->Get("DelphesFromMTDToNominalEffLoose");
    btaglooseeff_reweightmap->SetDirectory(0);
    btag_reweight_file->Close();
  }
  if(ReweightSetup=="DelphesMTD_PabloMTD")
  {
    cout<<"loading btag reweight map from file /afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/data/ScaleFactorsFor50ps.root"<<endl;
    TFile *btag_reweight_file = new TFile("/afs/cern.ch/user/f/fmonti/work/HHGGBBAnalysis/data/ScaleFactorsFor50ps.root");
    TH2F* btaglooseeff_reweightmap1 = (TH2F*) btag_reweight_file->Get("DelphesFromMTDToNominalEffLoose");
    TH2F* btaglooseeff_reweightmap2 = (TH2F*) btag_reweight_file->Get("eff_BTagSFLoose");
    btaglooseeff_reweightmap=Multiply(btaglooseeff_reweightmap1,btaglooseeff_reweightmap2);
    btaglooseeff_reweightmap->SetDirectory(0);
    btag_reweight_file->Close();
  }

  if(!btaglooseeff_reweightmap)
  {
    cout<<"[ERROR]: btag efficiency or fake rate reweight map not found"<<endl;
    return -1;
  }
  else
    Ptmax_btag_reweightmap = btaglooseeff_reweightmap->GetXaxis()->GetXmax();   

  cout<<"Reading file "<<oldfilename.Data()<<endl;
  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile(oldfilename.Data());
  TTree *oldtree = (TTree*)oldfile->Get(treename.Data());
  Long64_t nentries = oldtree->GetEntries();
  float evWeight;
  float dibjet_leadgenflav;
  float dibjet_subleadgenflav;
  float dibjet_leadbtaglevel;
  float dibjet_subleadbtaglevel;
  float dibjet_leadEta;
  float dibjet_subleadEta;
  float dibjet_leadPt;
  float dibjet_subleadPt;

  oldtree->SetBranchAddress("evWeight",&evWeight);
  oldtree->SetBranchAddress("dibjet_leadgenflav",&dibjet_leadgenflav);
  oldtree->SetBranchAddress("dibjet_subleadgenflav",&dibjet_subleadgenflav);
  oldtree->SetBranchAddress("dibjet_leadbtaglevel",&dibjet_leadbtaglevel);
  oldtree->SetBranchAddress("dibjet_subleadbtaglevel",&dibjet_subleadbtaglevel);
  oldtree->SetBranchAddress("dibjet_leadEta",&dibjet_leadEta);
  oldtree->SetBranchAddress("dibjet_subleadEta",&dibjet_subleadEta);
  oldtree->SetBranchAddress("dibjet_leadPt",&dibjet_leadPt);
  oldtree->SetBranchAddress("dibjet_subleadPt",&dibjet_subleadPt);

  //Create a new file + a clone of old tree in new file
  TString newfilename(oldfilename);
  if(newfilename.Contains("/"))
    newfilename.Insert(newfilename.Last('/')+1,ReweightSetup+"_reweightSECCO_"+treename+"_");
  else
    newfilename.Insert(0,ReweightSetup+"_reweightSECCO_"+treename+"_");
  TFile *newfile = new TFile(newfilename.Data(),"recreate");
  TTree *newtree = oldtree->CloneTree(0);

  //loop over events
  for (Long64_t i=0;i<nentries; i++) 
  {
    oldtree->GetEntry(i);
    /////////////////////////////////////////////////////////////////////////////////
    //temporary for validation
    //evWeight = 0.097 / 499949;    
    /////////////////////////////////////////////////////////////////////////////////

    if(dibjet_leadgenflav==5 && dibjet_leadbtaglevel>=4)
    {
      double pt_jet = std::min(Ptmax_btag_reweightmap-1,dibjet_leadPt);
      evWeight *= btaglooseeff_reweightmap->GetBinContent( btaglooseeff_reweightmap->FindBin(pt_jet,fabs(dibjet_leadEta)));
      //cout<<"SF="<<btaglooseeff_reweightmap->GetBinContent( btaglooseeff_reweightmap->FindBin(pt_jet,fabs(dibjet_leadEta)))<<endl;
    }

    if(dibjet_subleadgenflav==5 && dibjet_subleadbtaglevel>=4)
    {
      double pt_jet = std::min(Ptmax_btag_reweightmap-1,dibjet_subleadPt);
      evWeight *= btaglooseeff_reweightmap->GetBinContent( btaglooseeff_reweightmap->FindBin(pt_jet,fabs(dibjet_subleadEta)));
      //cout<<"SF="<<btaglooseeff_reweightmap->GetBinContent( btaglooseeff_reweightmap->FindBin(pt_jet,fabs(dibjet_subleadEta)))<<endl;
    }
    newtree->Fill();

  }
  newtree->AutoSave();

  delete oldfile;
  delete newfile;
  delete btaglooseeff_reweightmap;
  return 0;
}
