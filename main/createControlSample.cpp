#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h" 
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TLegend.h"



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> createControlSample.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }

  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  //--- open files and get trees
  std::vector<std::string> input = opts.GetOpt<std::vector<std::string> >("Input.input");
  
  TChain* t = new TChain("chain","the chain");
  for(unsigned int ii = 0; ii < input.size()/2; ++ii)
  {
    std::string inFileName = input.at(0+ii*2);
    std::string treeName = input.at(1+ii*2);
    t -> Add((inFileName+"/"+treeName).c_str());
  }
  long int nEntries = t->GetEntries();
  std::cout << "Added " << nEntries << " entries to the chain" << std::endl;
  
  TreeVars treeVars;
  InitTreeVars(t,treeVars);
  
  
  //------------------------------------
  // create histograms for normalization
  std::cout << ">>> creating histograms for normalization" << std::endl;
  
  TH1F* h1_data_oneCategory = new TH1F("h1_data_oneCategory","",80,100.,180.);
  TH1F* h1_cs_oneCategory   = new TH1F("h1_cs_oneCategory","",80,100.,180.);
  
  TH1F* h1_data_diLepton = new TH1F("h1_data_diLepton","",80,100.,180.);
  TH1F* h1_cs_diLepton   = new TH1F("h1_cs_diLepton","",80,100.,180.);
  
  TH1F* h1_data_singleLepton = new TH1F("h1_data_singleLepton","",80,100.,180.);
  TH1F* h1_cs_singleLepton   = new TH1F("h1_cs_singleLepton","",80,100.,180.);
  
  for(long int ii = 0; ii < nEntries; ++ii)
  {
    if( ii%1000 == 0 ) std::cout << ">>> Reading entry " << ii << " / " << nEntries << "\r" << std::flush;
    t -> GetEntry(ii);
    
    
    // common cuts                                                                                                                                                                                                                          
    if( treeVars.dipho_mass < 100 || treeVars.dipho_mass > 180 ) continue;
    if( treeVars.dipho_mass > 115 && treeVars.dipho_mass < 135 ) continue;
    
    
    //-------------
    // one category
    
    if( OneCategorySelection(treeVars,-1) == true )
    {
      h1_data_oneCategory -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    else if( OneCategorySelection(treeVars,-2) == true )
    {
      h1_cs_oneCategory -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    // one category
    
    
    //--------------------------
    // two categories - dilepton
    
    DiLeptonCategories cat;
    if( DiLeptonSelection(treeVars,-1,cat) == true )
    {
      h1_data_diLepton -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    else if( DiLeptonSelection(treeVars,-2,cat) == true )
    {
      h1_cs_diLepton -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    // two categories - dilepton
    
    
    //-------------------------------
    // two categories - single lepton
    
    if( SingleLeptonSelection(treeVars,-1) == true )
    {
      h1_data_singleLepton -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    else if( SingleLeptonSelection(treeVars,-2) == true)
    {
      h1_cs_singleLepton -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    // two categories - single lepton
  }
  
  float scaleFactor_oneCategory = h1_data_oneCategory->Integral() / h1_cs_oneCategory->Integral();
  float scaleFactor_diLepton = h1_data_diLepton->Integral() / h1_cs_diLepton->Integral();
  float scaleFactor_singleLepton = h1_data_singleLepton->Integral() / h1_cs_singleLepton->Integral();
  
  
  //-------------------------------------------------------
  //--- create a new file + a clone of old tree in new file
  std::string outputFileName = opts.GetOpt<std::string>("Output.outputFileName");
  std::string outputTreeName = opts.GetOpt<std::string>("Output.outputTreeName");
  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  TTree* outputTree_oneCategory = t->CloneTree(0);
  outputTree_oneCategory -> SetName((outputTreeName+"_oneCategory").c_str());
  TTree* outputTree_diLepton = t->CloneTree(0);
  outputTree_diLepton -> SetName((outputTreeName+"_diLepton").c_str());
  TTree* outputTree_singleLepton = t->CloneTree(0);
  outputTree_singleLepton -> SetName((outputTreeName+"_singleLepton").c_str());
  
  for(long int ii = 0; ii < nEntries; ++ii)
  {
    if( ii%1000 == 0 ) std::cout << ">>> Reading entry " << ii << " / " << nEntries << "\r" << std::flush;
    t -> GetEntry(ii);
    
    float oldWeight = treeVars.weight;
    
    if( OneCategorySelection(treeVars,-1) == false && OneCategorySelection(treeVars,-2) == true )
    {
      treeVars.weight = oldWeight * scaleFactor_oneCategory;
      outputTree_oneCategory -> Fill();
    }
    
    DiLeptonCategories cat;
    if( DiLeptonSelection(treeVars,-1,cat) == false && DiLeptonSelection(treeVars,-2,cat) == true )
    {
      treeVars.weight = oldWeight * scaleFactor_diLepton;
      outputTree_diLepton -> Fill();
    }
    
    if( SingleLeptonSelection(treeVars,-1) == false && SingleLeptonSelection(treeVars,-2) == true )
    {
      treeVars.weight = oldWeight * scaleFactor_singleLepton;
      outputTree_singleLepton -> Fill();
    }
  }
  
  outputTree_oneCategory -> AutoSave();
  outputTree_diLepton -> AutoSave();
  outputTree_singleLepton -> AutoSave();
  
  h1_data_oneCategory  -> Write();
  h1_cs_oneCategory    -> Write();
  h1_data_diLepton     -> Write();
  h1_cs_diLepton       -> Write();
  h1_data_singleLepton -> Write();
  h1_cs_singleLepton   -> Write();
     
  outputFile -> Close();
  
  
  return 0;
}
