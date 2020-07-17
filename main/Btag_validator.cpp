#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/CMS_lumi.h"
#include "interface/SetTDRStyle.h"
#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"

#include <iostream>

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
#include "TRandom3.h"

using namespace std;

TRandom3 rndm;
void generatenewBtag(RawTreeVars &treeVars, map<int,TH2F*> &eff_map_flav4, map<int,TH2F*> &eff_map_flav5);

int main(int argc, char* argv[])
{
  gStyle->SetOptStat(0);
  if( argc < 2 )
  {
    std::cerr << ">>>>> HHGGBB_selector.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }
  
  //----------------------
  // load btag maps
  map<int,TH2F*> eff_map_flav5;
  map<int,TH2F*> eff_map_flav4;
  TFile *effmapfile = new TFile("/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/neweffmap_NOMTD.root");
  //  TFile *effmapfile = new TFile("/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/neweffmap.root");
  eff_map_flav5[4]=(TH2F*)effmapfile->Get("h2_loose_flavor5_conditioned");
  eff_map_flav5[4]->SetDirectory(0);  
  eff_map_flav5[5]=(TH2F*)effmapfile->Get("h2_medium_flavor5_conditioned");
  eff_map_flav5[5]->SetDirectory(0);  
  eff_map_flav5[6]=(TH2F*)effmapfile->Get("h2_tight_flavor5_conditioned");
  eff_map_flav5[6]->SetDirectory(0);  
  eff_map_flav4[4]=(TH2F*)effmapfile->Get("h2_loose_flavor4_conditioned");
  eff_map_flav4[4]->SetDirectory(0);  
  eff_map_flav4[5]=(TH2F*)effmapfile->Get("h2_medium_flavor4_conditioned");
  eff_map_flav4[5]->SetDirectory(0);  
  eff_map_flav4[6]=(TH2F*)effmapfile->Get("h2_tight_flavor4_conditioned");
  eff_map_flav4[6]->SetDirectory(0);  
  effmapfile->Close();
  //----------------------
  // parse the config file
  
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  std::string filelist = opts.GetOpt<std::string>("Input.filelist");
  std::vector<std::string> treename = opts.GetOpt<std::vector<std::string> >("Input.treename");
  std::string label =  opts.GetOpt<std::string> ("Input.label");
  
  std::string Loose_Tight_Photon = "PhotonLoose";
  if(opts.OptExist("Input.Loose_Tight_Photon"))
    Loose_Tight_Photon = opts.GetOpt<std::string> ("Input.Loose_Tight_Photon");
  else
    cout<<"Option <Input.Loose_Tight_Photon> not found --> Analysis by default on "<<Loose_Tight_Photon<<endl;
    
  //------------------
  // define histograms
  
  map<TString,TProfile*> p;
  map<TString,TProfile*> p_new;
  p["p_effloose_flav5_pt"]      = new TProfile("p_effloose_flav5_pt","p_effloose_flav5_pt",15,0.,500.);
  p["p_effmedium_flav5_pt"]     = new TProfile("p_effmedium_flav5_pt","p_effmedium_flav5_pt",15,0.,500.);
  p["p_efftight_flav5_pt"]      = new TProfile("p_efftight_flav5_pt","p_efftight_flav5_pt",15,0.,500.);
  p_new["p_effloose_flav5_pt"]  = new TProfile("p_effloosenew_flav5_pt","p_effloosenew_flav5_pt",15,0.,500.);
  p_new["p_effmedium_flav5_pt"] = new TProfile("p_effmediumnew_flav5_pt","p_effmediumnew_flav5_pt",15,0.,500.);
  p_new["p_efftight_flav5_pt"]  = new TProfile("p_efftightnew_flav5_pt","p_efftightnew_flav5_pt",15,0.,500.);

  p["p_effloose_flav5_eta"]     = new TProfile("p_effloose_flav5_eta","p_effloose_flav5_eta",20,-3.5,3.5);
  p["p_effmedium_flav5_eta"]    = new TProfile("p_effmedium_flav5_eta","p_effmedium_flav5_eta",20,-3.5,3.5);
  p["p_efftight_flav5_eta"]     = new TProfile("p_efftight_flav5_eta","p_efftight_flav5_eta",20,-3.5,3.5);
  p_new["p_effloose_flav5_eta"] = new TProfile("p_effloosenew_flav5_eta","p_effloosenew_flav5_eta",20,-3.5,3.5);
  p_new["p_effmedium_flav5_eta"]= new TProfile("p_effmediumnew_flav5_eta","p_effmediumnew_flav5_eta",20,-3.5,3.5);
  p_new["p_efftight_flav5_eta"] = new TProfile("p_efftightnew_flav5_eta","p_efftightnew_flav5_eta",20,-3.5,3.5);

  p["p_effloose_flav4_pt"]      = new TProfile("p_effloose_flav4_pt","p_effloose_flav4_pt",15,0.,500.);
  p["p_effmedium_flav4_pt"]     = new TProfile("p_effmedium_flav4_pt","p_effmedium_flav4_pt",15,0.,500.);
  p["p_efftight_flav4_pt"]      = new TProfile("p_efftight_flav4_pt","p_efftight_flav4_pt",15,0.,500.);
  p_new["p_effloose_flav4_pt"]  = new TProfile("p_effloosenew_flav4_pt","p_effloosenew_flav4_pt",15,0.,500.);
  p_new["p_effmedium_flav4_pt"] = new TProfile("p_effmediumnew_flav4_pt","p_effmediumnew_flav4_pt",15,0.,500.);
  p_new["p_efftight_flav4_pt"]  = new TProfile("p_efftightnew_flav4_pt","p_efftightnew_flav4_pt",15,0.,500.);

  p["p_effloose_flav4_eta"]     = new TProfile("p_effloose_flav4_eta","p_effloose_flav4_eta",20,-3.5,3.5);
  p["p_effmedium_flav4_eta"]    = new TProfile("p_effmedium_flav4_eta","p_effmedium_flav4_eta",20,-3.5,3.5);
  p["p_efftight_flav4_eta"]     = new TProfile("p_efftight_flav4_eta","p_efftight_flav4_eta",20,-3.5,3.5);
  p_new["p_effloose_flav4_eta"] = new TProfile("p_effloosenew_flav4_eta","p_effloosenew_flav4_eta",20,-3.5,3.5);
  p_new["p_effmedium_flav4_eta"]= new TProfile("p_effmediumnew_flav4_eta","p_effmediumnew_flav4_eta",20,-3.5,3.5);
  p_new["p_efftight_flav4_eta"] = new TProfile("p_efftightnew_flav4_eta","p_efftightnew_flav4_eta",20,-3.5,3.5);

  p["p_effloose_flavlight_pt"]      = new TProfile("p_effloose_flavlight_pt","p_effloose_flavlight_pt",15,0.,500.);
  p["p_effmedium_flavlight_pt"]     = new TProfile("p_effmedium_flavlight_pt","p_effmedium_flavlight_pt",15,0.,500.);
  p["p_efftight_flavlight_pt"]      = new TProfile("p_efftight_flavlight_pt","p_efftight_flavlight_pt",15,0.,500.);
  p_new["p_effloose_flavlight_pt"]  = new TProfile("p_effloosenew_flavlight_pt","p_effloosenew_flavlight_pt",15,0.,500.);
  p_new["p_effmedium_flavlight_pt"] = new TProfile("p_effmediumnew_flavlight_pt","p_effmediumnew_flavlight_pt",15,0.,500.);
  p_new["p_efftight_flavlight_pt"]  = new TProfile("p_efftightnew_flavlight_pt","p_efftightnew_flavlight_pt",15,0.,500.);

  p["p_effloose_flavlight_eta"]     = new TProfile("p_effloose_flavlight_eta","p_effloose_flavlight_eta",20,-3.5,3.5);
  p["p_effmedium_flavlight_eta"]    = new TProfile("p_effmedium_flavlight_eta","p_effmedium_flavlight_eta",20,-3.5,3.5);
  p["p_efftight_flavlight_eta"]     = new TProfile("p_efftight_flavlight_eta","p_efftight_flavlight_eta",20,-3.5,3.5);
  p_new["p_effloose_flavlight_eta"] = new TProfile("p_effloosenew_flavlight_eta","p_effloosenew_flavlight_eta",20,-3.5,3.5);
  p_new["p_effmedium_flavlight_eta"]= new TProfile("p_effmediumnew_flavlight_eta","p_effmediumnew_flavlight_eta",20,-3.5,3.5);
  p_new["p_efftight_flavlight_eta"] = new TProfile("p_efftightnew_flavlight_eta","p_efftightnew_flavlight_eta",20,-3.5,3.5);


  //----------
  // get trees
  std::map<std::string,TChain*> trees;
  std::vector<std::string> onlytreename_str;
  
  rndm.SetSeed();

  for(unsigned int ntree = 0; ntree < treename.size(); ++ntree)
  {
    TString onlytreename(treename.at(ntree));//treename could be "ntuple/only_treename" but the method InitRawTreeVars wants only only_treename
    onlytreename.Remove(0,onlytreename.Last('/')+1);
    onlytreename_str.push_back(onlytreename.Data());
    trees[onlytreename_str.at(ntree)] = new TChain(onlytreename_str.at(ntree).c_str(),"");
  }
  
  std::ifstream list(filelist.c_str(),std::ios::in);
  std::string filename;
  while(1)
  {
    getline(list,filename,'\n');
    if( !list.good() ) break;
    if( filename.at(0) == '#' ) continue;
    std::cout << ">>> reading file " << filename << std::endl;
    
    for(unsigned int ntree = 0; ntree < treename.size(); ++ntree)
    {
      TString onlytreename(treename.at(ntree));//treename could be "ntuple/only_treename" but the method InitRawTreeVars wants only only_treename
      
      //std::cout << ">>> Adding trees " << treename.at(ntree) << " to chain " << onlytreename_str.at(ntree) << std::endl;
      trees[onlytreename_str.at(ntree)] -> Add((filename+"/"+treename.at(ntree)).c_str());
    }
  }
  
  //---------------
  // tree variables
  RawTreeVars treeVars;

  // branch tree: the only functioning way consists in branching separately the trees and, only after, adding them as friends
  InitRawTreeVars(trees,treeVars,Loose_Tight_Photon);
  TChain *tree = trees[onlytreename_str.at(0)];

  // add tree_i as friends to tree_0
  for(unsigned int ntree = 1; ntree < onlytreename_str.size(); ++ntree)
  {
    std::cout << ">>> Adding chain " << onlytreename_str.at(ntree) << " to chain " << onlytreename_str.at(0) << std::endl;
    tree->AddFriend(onlytreename_str.at(ntree).c_str(),"");
  }
  

  //-------------      
  // output trees
  std::string outputPlotFolder = opts.GetOpt<std::string>("Output.outputFolder");
  system(Form("mkdir -p %s",outputPlotFolder.c_str()));
  TFile* outFile = TFile::Open(Form("%s/plotTree_%s.root",outputPlotFolder.c_str(),label.c_str()),"RECREATE");
  outFile -> cd();  
  TTree* outTree  = new TTree("all", "all");
  TreeVars outtreeVars;
  InitOutTreeVars(outTree,outtreeVars);
  
  int nEntries = tree->GetEntries();
  std::cout << "Total entries = " << nEntries << std::endl;
  for(int i=0; i<nEntries; i++)
  {
    tree -> GetEntry(i);
    if( i%1000==0 ) std::cout << "Processing entry "<< i << "\r" << std::flush;

    /*    
    // count leptons
    outtreeVars.nLep = 0;
    for(int i=0;i<treeVars.N_TightEl;i++)
    {
      if( treeVars.TightEl_pt[i]<10 ) continue;
      if( fabs(treeVars.TightEl_eta[i])>2.5 ) continue;
      if( fabs(treeVars.TightEl_eta[i])>1.44 && fabs(treeVars.TightEl_eta[i])<1.57 ) continue;
      outtreeVars.nLep += 1;
    }
    for(int i=0;i<treeVars.N_TightMu;i++)
    {
      if( treeVars.TightMu_pt[i]<5 ) continue;
      if( fabs(treeVars.TightMu_eta[i])>2.4 ) continue;
      outtreeVars.nLep += 1;
    }
    if(outtreeVars.nLep>0)
      continue;
    */
    //Jets selections
    //int BTagLoose_mask  = 0b001000;
    //int BTagMedium_mask = 0b010000;
    //int BTagTight_mask  = 0b100000;
    //Jets selections
    int BTagLoose_mask  = 0b000001;
    int BTagMedium_mask = 0b000010;
    int BTagTight_mask  = 0b000100;

    generatenewBtag(treeVars,eff_map_flav4,eff_map_flav5);

    for(int i=0;i<treeVars.N_Jet;i++)
    {
      if(treeVars.Jet_pt[i]<25) continue;
      if(fabs(treeVars.Jet_eta[i])>3.5) continue;
      if(abs(treeVars.Jet_hadflav[i])==5)
      {
	int btag = treeVars.Jet_deepcsv[i];
	if(btag & BTagLoose_mask)
	{
	  p["p_effloose_flav5_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p["p_effloose_flav5_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p["p_effloose_flav5_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p["p_effloose_flav5_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag & BTagMedium_mask)
	{
	  p["p_effmedium_flav5_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p["p_effmedium_flav5_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p["p_effmedium_flav5_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p["p_effmedium_flav5_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag & BTagTight_mask)
	{
	  p["p_efftight_flav5_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p["p_efftight_flav5_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p["p_efftight_flav5_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p["p_efftight_flav5_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	btag = treeVars.Jet_mvav2[i];
	if(btag>=4 && btag<=6)
	{
	  p_new["p_effloose_flav5_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p_new["p_effloose_flav5_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p_new["p_effloose_flav5_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p_new["p_effloose_flav5_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag>=5 && btag<=6)
	{
	  p_new["p_effmedium_flav5_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p_new["p_effmedium_flav5_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p_new["p_effmedium_flav5_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p_new["p_effmedium_flav5_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag==6)
	{
	  p_new["p_efftight_flav5_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p_new["p_efftight_flav5_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p_new["p_efftight_flav5_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p_new["p_efftight_flav5_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}
      }
      else if(abs(treeVars.Jet_hadflav[i])==4)
      {
	int btag = treeVars.Jet_deepcsv[i];
	if(btag & BTagLoose_mask)
	{
	  p["p_effloose_flav4_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p["p_effloose_flav4_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p["p_effloose_flav4_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p["p_effloose_flav4_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag & BTagMedium_mask)
	{
	  p["p_effmedium_flav4_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p["p_effmedium_flav4_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p["p_effmedium_flav4_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p["p_effmedium_flav4_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag & BTagTight_mask)
	{
	  p["p_efftight_flav4_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p["p_efftight_flav4_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p["p_efftight_flav4_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p["p_efftight_flav4_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	btag = treeVars.Jet_mvav2[i];
	if(btag>=4 && btag<=6)
	{
	  p_new["p_effloose_flav4_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p_new["p_effloose_flav4_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p_new["p_effloose_flav4_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p_new["p_effloose_flav4_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag>=5 && btag<=6)
	{
	  p_new["p_effmedium_flav4_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p_new["p_effmedium_flav4_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p_new["p_effmedium_flav4_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p_new["p_effmedium_flav4_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag==6)
	{
	  p_new["p_efftight_flav4_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p_new["p_efftight_flav4_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p_new["p_efftight_flav4_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p_new["p_efftight_flav4_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}
      }
      else
      {
	int btag = treeVars.Jet_deepcsv[i];
	if(btag & BTagLoose_mask)
	{
	  p["p_effloose_flavlight_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p["p_effloose_flavlight_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p["p_effloose_flavlight_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p["p_effloose_flavlight_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag & BTagMedium_mask)
	{
	  p["p_effmedium_flavlight_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p["p_effmedium_flavlight_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p["p_effmedium_flavlight_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p["p_effmedium_flavlight_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag & BTagTight_mask)
	{
	  p["p_efftight_flavlight_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p["p_efftight_flavlight_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p["p_efftight_flavlight_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p["p_efftight_flavlight_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	btag = treeVars.Jet_mvav2[i];
	if(btag>=4 && btag<=6)
	{
	  p_new["p_effloose_flavlight_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p_new["p_effloose_flavlight_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p_new["p_effloose_flavlight_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p_new["p_effloose_flavlight_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag>=5 && btag<=6)
	{
	  p_new["p_effmedium_flavlight_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p_new["p_effmedium_flavlight_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p_new["p_effmedium_flavlight_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p_new["p_effmedium_flavlight_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}

	if(btag==6)
	{
	  p_new["p_efftight_flavlight_pt"]->Fill(treeVars.Jet_pt[i],1.);
	  p_new["p_efftight_flavlight_eta"]->Fill(treeVars.Jet_eta[i],1.);
	}
	else
	{
	  p_new["p_efftight_flavlight_pt"]->Fill(treeVars.Jet_pt[i],0.);
	  p_new["p_efftight_flavlight_eta"]->Fill(treeVars.Jet_eta[i],0.);
	}
      }
    }
    
    //outTree->Fill();
  }
  outTree -> AutoSave();

  TCanvas* c = new TCanvas();
  for(auto prof : p)
  {
    (prof.second)->SetTitle("std Delphes btag");
    (prof.second)->Draw();
    auto is_pt_plot = p.find("pt");
    if(is_pt_plot!=p.end())
      (prof.second)->GetXaxis()->SetTitle("p_{T} (GeV)");
    else
      (prof.second)->GetXaxis()->SetTitle("#eta");
    auto prof_new = p_new.find(prof.first);
    if(prof_new != p_new.end())
    {
      (prof_new->second)->SetLineColor(kRed);
      (prof_new->second)->SetTitle("new Delphes btag");
      (prof_new->second)->Draw("SAME");
    }
    c->BuildLegend();
    c->Print((prof.first+".png").Data());
    c->Print((prof.first+".pdf").Data());
  }
  system(Form("mv *.png %s",outputPlotFolder.c_str()));
  system(Form("mv *.pdf %s",outputPlotFolder.c_str()));

  for(auto prof : p)
    delete prof.second;

  for(auto prof_new : p_new)
    delete prof_new.second;

  return 0;
}

void generatenewBtag(RawTreeVars &treeVars, map<int,TH2F*> &eff_map_flav4, map<int,TH2F*> &eff_map_flav5)
{
  float eta,pt,P;
  for(int i=0;i<treeVars.N_Jet;i++)
  {
    eta=treeVars.Jet_eta[i];
    pt=treeVars.Jet_pt[i];
    if(pt<25 || fabs(eta)>3.5)
    {
      treeVars.Jet_mvav2[i]=0;
      continue;
    }

    switch(abs(treeVars.Jet_hadflav[i]))
    {
      case 5:
	P=eff_map_flav5[4]->GetBinContent(eff_map_flav5[4]->FindBin(pt,abs(eta)));
	//cout<<"Ploose"<<P<<endl;
	if(rndm.Uniform()<=P)//it's loose
        {
	  treeVars.Jet_mvav2[i]=4;
	  P=eff_map_flav5[5]->GetBinContent(eff_map_flav5[5]->FindBin(pt,abs(eta)));
	  //cout<<"Pmedium"<<P<<endl;
	  if(rndm.Uniform()<=P)//it's medium
	  {
	    treeVars.Jet_mvav2[i]=5;
	    P=eff_map_flav5[6]->GetBinContent(eff_map_flav5[6]->FindBin(pt,abs(eta)));
	    //cout<<"Ptight"<<P<<endl;
	    if(rndm.Uniform()<=P)//it's tight
	      treeVars.Jet_mvav2[i]=6;
	  }
	}
	else
	  treeVars.Jet_mvav2[i]=0;
	break;
      case 4:
	P=eff_map_flav4[4]->GetBinContent(eff_map_flav4[4]->FindBin(pt,abs(eta)));
	//cout<<"eta="<<eta<<"\tppt="<<pt<<"\tPloose"<<P<<endl;
	if(rndm.Uniform()<=P)//it's loose
        {
	  treeVars.Jet_mvav2[i]=4;
	  P=eff_map_flav4[5]->GetBinContent(eff_map_flav4[5]->FindBin(pt,abs(eta)));
	  if(rndm.Uniform()<=P)//it's medium
	  {
	    treeVars.Jet_mvav2[i]=5;
	    P=eff_map_flav4[6]->GetBinContent(eff_map_flav4[6]->FindBin(pt,abs(eta)));
	    if(rndm.Uniform()<=P)//it's tight
	      treeVars.Jet_mvav2[i]=6;
	  }
	}
	else
	  treeVars.Jet_mvav2[i]=0;
	break;
      default:
	P=0.1;
	if(rndm.Uniform()<=P)//it's loose
        {
	  treeVars.Jet_mvav2[i]=4;
	  P=0.1;
	  if(rndm.Uniform()<=P)//it's medium
	  {
	    treeVars.Jet_mvav2[i]=5;
	    P=0.1;
	    if(rndm.Uniform()<=P)//it's tight
	      treeVars.Jet_mvav2[i]=6;
	  }
	}
	else
	  treeVars.Jet_mvav2[i]=0;
    }
  }
}
