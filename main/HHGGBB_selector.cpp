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

#include "RooMsgService.h"
#include "RooRealVar.h"

using namespace std;
using namespace RooFit;

bool cutBased = false;

float lumi = 35.9;

int jetThreshold = 2;
int bjetThreshold = 1;
float jetPtThreshold = 0.;



int main(int argc, char* argv[])
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> HHGGBB_selector.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }
  
  
  //------------------
  // graphics settings
  
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = Form("%.1f fb^{-1} (13 TeV)",lumi);
  
  setTDRStyle();
  gStyle -> SetOptFit(0);
  gStyle -> SetOptStat(0);
  
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  
  
  //----------------------
  // parse the config file
  
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  //std::vector<std::string> input = opts.GetOpt<std::vector<std::string> >("Input.input");
  
  std::vector<std::string> filename = opts.GetOpt<std::vector<std::string> >("Input.filename");
  std::vector<std::string> treename = opts.GetOpt<std::vector<std::string> >("Input.treename");
  std::string label =  opts.GetOpt<std::string> ("Input.label");
  std::string type =  opts.GetOpt<std::string> ("Input.type");

   
  //------------------
  // define histograms

  
  
  //----------
  // get trees
  
  std::map<std::string,TChain*> trees;
  std::vector<std::string> onlytreename_str;
  for(unsigned int ntree = 0; ntree < treename.size(); ++ntree)
  {
    TString onlytreename(treename.at(ntree));//treename could be "ntuple/only_treename" but the method InitRawTreeVars wants only only_treename
    onlytreename.Remove(0,onlytreename.Last('/')+1);
    onlytreename_str.push_back(onlytreename.Data());
    trees[onlytreename_str.at(ntree)] = new TChain(onlytreename_str.at(ntree).c_str(),"");
    for(unsigned int nfile = 0; nfile < filename.size(); ++nfile)
    {
      std::cout << ">>> Adding trees " << filename.at(nfile)+"/"+treename.at(ntree) << " to chain " << onlytreename_str.at(ntree) << std::endl;
      trees[onlytreename_str.at(ntree)] -> Add((filename.at(nfile)+"/"+treename.at(ntree)).c_str());
    }
  }


  //---------------
  // tree variables
  RawTreeVars treeVars;
  
  //------------------
  // branch tree: the only functioning way consists in branching separately the trees and, only after, adding them as friends
  InitRawTreeVars(trees,treeVars);
  TChain *tree = trees[onlytreename_str.at(0)];


  //----------
  // add tree_i as friends to tree_0
  for(unsigned int ntree = 1; ntree < onlytreename_str.size(); ++ntree)
  {
    std::cout << ">>> Adding chain " << onlytreename_str.at(ntree) << " to chain " << onlytreename_str.at(0) << std::endl;
    tree->AddFriend(onlytreename_str.at(ntree).c_str(),"");
  }
  
  //------------------
  // loop over samples
  /*    
  //output trees for plots
  std::string outputPlotFolder = opts.GetOpt<std::string>("Output.outputPlotFolder");
    
  TFile* outFile = TFile::Open(Form("%s/plotTree_%s.root",outputPlotFolder.c_str(),label.c_str()),"RECREATE");
  outFile -> cd();
  TTree* outTree_2jets = new TTree("plotTree_2jets","plotTree_2jets");
  InitOutTreeVars(outTree_2jets,treeVars);
  */  
  int nEntries = tree->GetEntries();
  int fCurrent = -1;
  std::cout << "Total entries = " << nEntries << std::endl;
  for(int i=0; i<nEntries; i++)
  {

    tree -> GetEntry(i);
    if (tree->GetTreeNumber() != fCurrent)
      fCurrent = tree->GetTreeNumber();
    if( i%1000==0 ) std::cout << "Processing entry "<< i << "\r" << std::flush;
    if(treeVars.N_TightPh>0)
      std::cout<<treeVars.TightPh_pt[0]<<std::endl;
  }
  std::cout << std::endl;
    /*      
    // common cuts
    if( treeVars.dipho_mass < 100 || treeVars.dipho_mass > 180 ) continue;
    if( type == -1 && treeVars.dipho_mass > 115 && treeVars.dipho_mass < 135 ) continue;
    bool passCutBased = CutBasedSelection(treeVars,0.4,0.3,0.,-0.5,2.5,2.);
      
    if( treeVars.nJets < 2 ) continue;
      
    outTree_2jets -> Fill();
      
    // fill event counters - cut based
    if( type == 1 )
        nEvents_cutBased["bkg"][0.] += treeVars.weight;
    if( type == 2 && label == "ttH" )
        nEvents_cutBased["sig"][0.] += treeVars.weight;
    if( passCutBased )
    {
      if( type == 1 )
        nEvents_cutBased["bkg"][1.] += treeVars.weight;
      if( type == 2 && label == "ttH" )
        nEvents_cutBased["sig"][1.] += treeVars.weight;
    }
    if( cutBased && !passCutBased ) continue;
 
    outTree_2jets -> AutoSave();
    outFile -> Close();
    std::cout << "Processed tag " << label << ", " << nEntries << " events out of " << nEntries << std::endl;
  }
  
  
  MakePlot2(mass_oneCat_histo, "oneCategory");
   
  system("mv *.png ~/www/ttH/DiLeptonStudy_new/");
  */   
    return 0;
}
