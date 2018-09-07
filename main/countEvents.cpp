#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/SetTDRStyle.h"
#include "interface/Plotter.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> countEvents.cpp::usage:   " << argv[0] << " inputFileList   [default=0/debug=1]" << std::endl;
    return -1;
  }
  TDirectory* baseDir = gDirectory;
  
  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  std::string inputFileList(argv[1]);
  
  TH1F* h_eventCounter = NULL;
  std::ifstream list(inputFileList.c_str(),std::ios::in);
  std::string fileName;
  while(1)
  {
    getline(list,fileName,'\n');
    if( !list.good() ) break;
    
    if( debugMode ) std::cout << ">>> opening file " << fileName << "...";
    TFile* inFile = TFile::Open(fileName.c_str(),"READ");
    
    if( h_eventCounter == NULL )
    {
      h_eventCounter = (TH1F*)( inFile->Get("weightCounter/Event_weight"));
      if( h_eventCounter != NULL )
      {
        baseDir -> cd();
        h_eventCounter = (TH1F*)( h_eventCounter->Clone());
      }
    }
    else
    {
      TH1F* h_temp = (TH1F*)( inFile->Get("weightCounter/Event_weight"));
      if( h_temp != NULL )
      {
        h_eventCounter -> Add( h_temp );
      }
    } 
    
    inFile -> Close();
    if( debugMode ) std::cout << "<<< file closed" << std::endl;
  }
  
  std::cout << ">>> nTotEvents: " << std::fixed << std::setprecision(0) << h_eventCounter->GetEntries() << std::endl;
  
  return 0;
}
