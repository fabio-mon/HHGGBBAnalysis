#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/CMS_lumi.h"
#include "interface/SetTDRStyle.h"
#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"
#include "interface/RooFitUtils.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TMath.h"


using namespace std;


int main(int argc, char* argv[])
{
  
  vector<string> treename = {"ntuple/Event","ntuple/PhotonTight"};
  vector<string> filename = {"/eos/user/f/fmonti/HHGGBB/output/DelphesDump/p2ntuple_GG_HH_2B2Gamma_200PU_0.root","/eos/user/f/fmonti/HHGGBB/output/DelphesDump/p2ntuple_GG_HH_2B2Gamma_200PU_1.root"};

  //  map<string,TChain*> trees;
  //  trees["Event"] = new TChain("tree_0","");
  //  trees["TightPh"] = new TChain("tree_1","");
  

  //----------
  // get trees
  
  std::map<std::string,TChain*> trees;

  for(unsigned int ntree = 0; ntree < treename.size(); ++ntree)
  {
    trees[treename.at(ntree)] = new TChain(Form("tree_%i",ntree),"");
    for(unsigned int nfile = 0; nfile < filename.size(); ++nfile)
    {
      std::cout << ">>> Adding trees " << filename.at(nfile)+"/"+treename.at(ntree) << " to chain " << Form("tree_%i",ntree) << std::endl;
      trees[treename.at(ntree)] -> Add((filename.at(nfile)+"/"+treename.at(ntree)).c_str());
    }
  }

  //-----------
  // branch tree

  int event;
  trees["ntuple/Event"]->SetBranchStatus("Event",1);
  trees["ntuple/Event"]->SetBranchAddress("Event",&event);

  int nTightPh;
  trees["ntuple/PhotonTight"]->SetBranchStatus("PhotonTight_size", 1);
  trees["ntuple/PhotonTight"]->SetBranchAddress("PhotonTight_size",&nTightPh);

  float TightPh_pt[500];
  trees["ntuple/PhotonTight"]->SetBranchStatus("PT", 1);
  trees["ntuple/PhotonTight"]->SetBranchAddress("PT",TightPh_pt);


  //----------
  // add tree_i as friends to tree_0
  TChain* tree = trees[treename.at(0)];
  for(unsigned int ntree = 1; ntree < treename.size(); ++ntree)
  {
    std::cout << ">>> Adding chain " << Form("tree_%i",ntree) << " to chain tree_0" << std::endl;
    trees[treename.at(0)]->AddFriend(Form("tree_%i",ntree),"");
  }



  int nEntries = tree->GetEntries();
  int fCurrent = -1;
  std::cout << "Total entries = " << nEntries << std::endl;
  for(int i=0; i<nEntries; i++)
  {

    tree -> GetEntry(i);
    if (tree->GetTreeNumber() != fCurrent)
      fCurrent = tree->GetTreeNumber();
    if( i%1000==0 ) std::cout << "Processing entry "<< i << "\r" << std::flush;
    if(nTightPh>0)
      std::cout<<"N_Ph="<<nTightPh<<"\t\tpt[0]="<<TightPh_pt[0]<<"\t\tevent="<<event<<std::endl;
  }
  std::cout << std::endl;
  return 0;
}
