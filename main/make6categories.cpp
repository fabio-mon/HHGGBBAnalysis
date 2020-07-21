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

int main(int argc, char** argv)
{
  TString LTfolder(argv[1]);
  TString HTfolder(argv[2]);
  TString outfolder(argv[3]);

  vector<TString> filenames = {
    "LT_DoubleEG.root",
    "LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root",
    "LT_output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root",
    "LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root",
    "LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root",
    "LT_output_bbHToGG_M-125_13TeV_amcatnlo.root",
    "LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root"};

  /*
  vector<TString> filenames = {
    "intermediate_subset_plotTree_HHggbb_HT_klambda_0.000_withMVA.root",
    "intermediate_subset_plotTree_HHggbb_HT_klambda_2.500_withMVA.root",
    "intermediate_subset_plotTree_HHggbb_HT_klambda_6.000_withMVA.root"};
  */

  for( auto filename : filenames)
  {
    
    TString fullfilename_LT = LTfolder+"/"+filename;
    TString fullfilename_HT = HTfolder+"/"+filename;
    /*
    TString fullfilename_HT = HTfolder+"/"+filename;
    TString fullfilename_LT = LTfolder+"/"+filename.ReplaceAll("HT","LT");
    */
    cout<<"Merging files "<<endl;
    cout<<"\t"<<fullfilename_LT.Data()<<endl;
    cout<<"\t"<<fullfilename_HT.Data()<<endl;
    TChain* ch = new TChain("TCVARS","TCVARS");
    ch->Add(fullfilename_LT.Data());
    ch->Add(fullfilename_HT.Data());
    //branch ch
    float mtot;
    int cut_based_ct;
    ch->SetBranchAddress("mtot",&mtot);
    ch->SetBranchAddress("cut_based_ct",&cut_based_ct);
    
    //create outfile and clone tree in outfile
    TFile* outputFile = TFile::Open((outfolder+"/"+filename).Data(),"RECREATE");
    outputFile -> cd();
    TTree* newTree = ch -> CloneTree(0);
    int catID;
    newTree->Branch("catID",&catID);

    //loop over events
    Long64_t nentries = ch->GetEntries();
    for (Long64_t i=0;i<nentries; i++) 
    {
      ch->GetEntry(i);
      if(cut_based_ct<0)
	continue;
      catID = 
	0*(mtot<350)+
	2*(mtot>350 && mtot<480)+
	4*(mtot>480)+
	cut_based_ct;
      //cout<<catID<<endl;
      newTree->Fill();
    }
    newTree->AutoSave();
    outputFile->Close();
    delete ch;
  }
  return 0;
}
