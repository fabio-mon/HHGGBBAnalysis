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

  for(unsigned int ntree = 0; ntree < treename.size(); ++ntree)
  {
    trees[treename.at(ntree)] = new TChain(Form("tree_%i",ntree),"");
    for(unsigned int nfile = 0; nfile < filename.size(); ++nfile)
    {
      std::cout << ">>> Adding trees " << filename.at(nfile)+"/"+treename.at(ntree) << " to chain " << Form("tree_%i",ntree) << std::endl;
      trees[treename.at(ntree)] -> Add((filename.at(nfile)+"/"+treename.at(ntree)).c_str());
    }
  }

  //----------
  // add tree_i as friends to tree_0
  for(unsigned int ntree = 1; ntree < treename.size(); ++ntree)
  {
    std::cout << ">>> Adding chain " << Form("tree_%i",ntree) << " to chain tree_0" << std::endl;
    trees[treename.at(0)]->AddFriend(Form("tree_%i",ntree),"");
  }

  /*
  TCanvas cTEST;
  trees[treename.at(0)]->Draw("Lumi");
  cTEST.Print("plots/test_Lumi.pdf");

  trees[treename.at(0)]->Draw("PhotonTight_size");
  cTEST.Print("plots/test_nPh.pdf");

  trees[treename.at(0)]->Draw("PT");
  cTEST.Print("plots/genericPT.pdf");

  trees[treename.at(0)]->Draw("PhotonTight.PT");
  cTEST.Print("plots/PhPT.pdf");
  */


 
  //---------------
  // tree variables
  RawTreeVars treeVars;
  
  //------------------
  // branch tree
  InitTreeVars(treees[treename.at(0)],treeVars);
  
  //------------------
  // loop over samples

  /*    
  //output trees for plots
  std::string outputPlotFolder = opts.GetOpt<std::string>("Output.outputPlotFolder");
    
  TFile* outFile = TFile::Open(Form("%s/plotTree_%s.root",outputPlotFolder.c_str(),label.c_str()),"RECREATE");
  outFile -> cd();
  TTree* outTree_2jets = new TTree("plotTree_2jets","plotTree_2jets");
  InitOutTreeVars(outTree_2jets,treeVars);
    
  int nEntries = tree->GetEntries();
  for(int i=0; i<nEntries; i++)
  {
    tree -> GetEntry(i);
    if( i%1000==0 ) std::cout << "Processing tag " << label << ", event " << i << " out of " << nEntries << "\r" << std::flush;
      
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
