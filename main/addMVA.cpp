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

#include "TMVA/Reader.h"



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> addMVA.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }

  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  
  //-------------------------
  // open files and get trees
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
  
  int nJobs = opts.GetOpt<int>("Input.nJobs");
  int jobId = opts.GetOpt<int>("Input.jobId");
  long int nEntriesPerJob = int(nEntries/nJobs);
  long int firstEntry = (jobId-1)*nEntriesPerJob;
  long int lastEntry = firstEntry + nEntriesPerJob;
  if( jobId == nJobs ) lastEntry = nEntries;

  
  TreeVars treeVars;
  InitTreeVars(t,treeVars);
  
  std::map<std::string,float*> varMap;
  varMap["mgg"] = &treeVars.mgg;
  varMap["mjj"] = &treeVars.mjj;
  varMap["evWeight"] = &treeVars.evWeight;
  varMap["MetPt"] = &treeVars.MetPt;
  varMap["DPhimin_met_bjet"] = &treeVars.DPhimin_met_bjet;
  varMap["DPhimax_met_bjet"] = &treeVars.DPhimax_met_bjet;
  varMap["DRmin_pho_bjet"] = &treeVars.DRmin_pho_bjet;
  varMap["costheta_HH"] = &treeVars.costheta_HH;
  varMap["costheta_bb"] = &treeVars.costheta_bb;
  varMap["event"] = &treeVars.event;
  varMap["nLep"] = &treeVars.nLep;
  varMap["nJets"] = &treeVars.nJets;
  
  
  //---------------
  // clone the tree
  std::string outputFileName = opts.GetOpt<std::string>("Output.outputFileName");
  TFile* outputFile = TFile::Open(outputFileName.c_str(),"RECREATE");
  outputFile -> cd();
  TTree* newTree = t -> CloneTree(0);
  
  std::vector<std::string> MVA_labels = opts.GetOpt<std::vector<std::string> >("Input.MVA_labels");
  std::map<std::string,float> mva;
  for(unsigned int ii = 0; ii < MVA_labels.size(); ++ii)
  {
    std::string MVA_label = MVA_labels.at(ii);
    newTree -> Branch(Form("mva_%s",MVA_label.c_str()),&mva[MVA_label]);
  }
  
  
  //-----
  // TMVA
  std::map<std::string,TMVA::Reader*> MVAReaders;
  std::map<std::string,std::string> MVA_methods;
  for(unsigned int ii = 0; ii < MVA_labels.size(); ++ii)
  {
    std::string MVA_label = MVA_labels.at(ii);
    
    MVA_methods[MVA_label] = opts.GetOpt<std::string>(Form("Input.%s.method",MVA_label.c_str()));
    std::string weightsFile = opts.GetOpt<std::string>(Form("Input.%s.weightsFile",MVA_label.c_str()));
    std::vector<std::string> inputVariables = opts.GetOpt<std::vector<std::string> >(Form("Input.%s.inputVariables",MVA_label.c_str()));
    std::vector<std::string> spectatorVariables = opts.GetOpt<std::vector<std::string> >(Form("Input.%s.spectatorVariables",MVA_label.c_str()));
    
    MVAReaders[MVA_label] = new TMVA::Reader( "!Color:!Silent" );
    
    for(unsigned int jj = 0; jj < inputVariables.size(); ++jj)
    {
      std::string inputVariable = inputVariables.at(jj);
      MVAReaders[MVA_label] -> AddVariable(inputVariable.c_str(),varMap[inputVariable.c_str()]);
      std::cout << "adding variable " << inputVariable << std::endl;
    }
    
    for(unsigned int jj = 0; jj < spectatorVariables.size(); ++jj)
    {
      std::string spectatorVariable = spectatorVariables.at(jj);
      MVAReaders[MVA_label] -> AddSpectator(spectatorVariable.c_str(),varMap[spectatorVariable.c_str()]);
    }
    
    MVAReaders[MVA_label] -> BookMVA( MVA_methods[MVA_label],weightsFile.c_str() );
  }
  
  
  //-----------------
  // loop over events
  std::cout << "jobId: " << jobId << " / " << nJobs << std::endl;
  std::cout << "tot entries: " << nEntries << std::endl;
  std::cout << "first entry: " << firstEntry << std::endl;
  std::cout << " last entry: " << lastEntry << std::endl;
  for(long int ii = firstEntry; ii < lastEntry; ++ii)
  {
    if( ii%1000 == 0 ) std::cout << ">>> Reading entry " << ii << " / " << nEntries << std::endl;
    t -> GetEntry(ii);
    
    // evaluate  MVA
    for(unsigned int ii = 0; ii < MVA_labels.size(); ++ii)
    {
      std::string MVA_label = MVA_labels.at(ii);
      mva[MVA_label] = MVAReaders[MVA_label] -> EvaluateMVA(MVA_methods[MVA_label].c_str());
    }
    
    newTree -> Fill();
  }
  
  
  newTree -> AutoSave();
  outputFile -> Close();
  
  
  return 0;
}
