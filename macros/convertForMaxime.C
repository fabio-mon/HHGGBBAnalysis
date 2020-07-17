#include <map>

//std::string COPY_mva_based_ct_TO_cut_based_ct(std::string inFileName,std::string treeName);
void convertForMaxime()
{
  std::map<std::string,std::vector<std::string> > out_in_map;
  out_in_map["output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph_0.root"].push_back("optimistic_subset_plotTree_HHggbb_withMVA.root");
  out_in_map["output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root"].push_back("optimistic_subset_plotTree_ggH_withMVA.root ");
  out_in_map["output_VBFHToGG_M-125_13TeV_powheg_pythia8.root"].push_back("optimistic_subset_plotTree_qqH_withMVA.root");
  out_in_map["output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root"].push_back("optimistic_subset_plotTree_VH_withMVA.root");
  out_in_map["output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root"].push_back("optimistic_subset_plotTree_ttH_withMVA.root");
  out_in_map["output_bbHToGG_M-125_13TeV_amcatnlo.root"].push_back("optimistic_subset_plotTree_bbH_withMVA.root");
  //out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_gg_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_gg_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_tt_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_ttghad_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_ttglep_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_ttgsemilepfromt_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_ttgsemilepfromtbar_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_ttgg_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_ggH_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_qqH_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_VH_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_ttH_withMVA.root");
  out_in_map["DoubleEG.root"].push_back("optimistic_subset_plotTree_bbH_withMVA.root");



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
	std::string treeName;
	if( cat == "LT" ) treeName="all_lowMx";
        else              treeName="all_highMx";
	std::string new_inFileName=inFileName;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//new_inFilename = COPY_mva_based_ct_TO_cut_based_ct(inFileName,treeName);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TFile* inFile = TFile::Open(Form("%s",new_inFileName.c_str()),"READ");
        TTree* tree;
        if( cat == "LT" ) tree = (TTree*)( inFile->Get(treeName.c_str()) );
        else              tree = (TTree*)( inFile->Get(treeName.c_str()) );
        
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

/*
std::string COPY_mva_based_ct_TO_cut_based_ct(std::string inFileName,std::string treeName)
{
  //cout<<"Reading file "<<oldfilename.Data()<<endl;
  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile(inFileName.c_str());
  TTree *oldtree = (TTree*)oldfile->Get(treeName.c_str());
  Long64_t nentries = oldtree->GetEntries();
  int cut_based_ct;
  int mva_based_ct;
  oldtree->SetBranchAddress("cut_based_ct",&cut_based_ct);
  oldtree->SetBranchAddress("mva_based_ct",&mva_based_ct);

  //Create a new file + a clone of old tree in new file
  TString newfilename(oldfilename);
  if(newfilename.Contains("/"))
    newfilename.Insert(newfilename.Last('/')+1,"subset_");
  else
    newfilename.Insert(0,"subset_");
  TFile *newfile = new TFile(newfilename.Data(),"recreate");
  TTree *newtree = oldtree->CloneTree(0);
  for (Long64_t i=0;i<nentries; i++) 
  {
    oldtree->GetEntry(i);
    cut_based_ct=mva_based_ct;
    newtree->Fill();
  }
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
  return newfilename;
}
*/
