#include <map>

void convertForMaxime()
{
  std::map<std::string,std::vector<std::string> > out_in_map;
  out_in_map["output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph_0.root"].push_back("plotTree_delphesntuples_loose.root");
  out_in_map["output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root"].push_back("plotTree_ggHGG_loose.root");
  out_in_map["output_VBFHToGG_M-125_13TeV_powheg_pythia8.root"].push_back("plotTree_VBFHgg_loose.root");
  out_in_map["output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root"].push_back("plotTree_ZHqqgg_loose.root");
  out_in_map["output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root"].push_back("plotTree_ttH_loose.root");
  out_in_map["output_bbHToGG_M-125_13TeV_amcatnlo.root"].push_back("plotTree_VBFHgg_loose.root");
  out_in_map["DoubleEG.root"].push_back("plotTree_DiPhotonJetsBox_loose.root");
  out_in_map["DoubleEG.root"].push_back("plotTree_ttGammaHadr_loose.root");
  
  std::vector<std::string> cats;
  cats.push_back("LT");
  cats.push_back("HT");
  
  for(auto cat: cats)
  {
    for(auto mapIt: out_in_map)
    {
      TFile* outFile = TFile::Open(Form("%s_%s",cat.c_str(),mapIt.first.c_str()),"RECREATE");
      std::vector<std::string> inFileNames = mapIt.second;
      
      TList* list = new TList();
      
      for(auto inFileName : inFileNames)
      {
        TFile* inFile = TFile::Open(Form("output/%s",inFileName.c_str()),"READ");
        TTree* tree;
        if( cat == "LT" ) tree = (TTree*)( inFile->Get("all_lowMx") );
        else              tree = (TTree*)( inFile->Get("all_highMx") );
        
        list -> Add(tree);
      }
      
      outFile -> cd();
      TTree* finalTree = TTree::MergeTrees(list);
      finalTree -> SetName("TCVARS");
      finalTree -> Write();
      
      outFile -> Close();
    }
  }
}
