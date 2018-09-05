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



int main(int argc, char* argv[])
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> HHGGBB_selector.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }
  
  
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
  
  bool useMTD=false;
  if(opts.OptExist("Input.useMTD"))
    useMTD = opts.GetOpt<bool> ("Input.useMTD");
  else
    cout<<"Option <Input.useMTD> not found --> Analysis by default on useMTD="<<useMTD<<endl;
  
  float CrossSection = 1.;
  if(opts.OptExist("Input.CrossSection"))
    CrossSection = opts.GetOpt<float> ("Input.CrossSection");
  
  float NEventsMC = 1.;
  if(opts.OptExist("Input.NEventsMC"))
    NEventsMC = opts.GetOpt<float> ("Input.NEventsMC");
  float weightMC = 1./NEventsMC;
  
  int Nph_veto = -1;//discard events with N selected photon gen-matched > Nph_veto --> Useful to avoid double counting in GJet and QCD samples   
  if(opts.OptExist("Input.Nph_veto"))
    Nph_veto = opts.GetOpt<int> ("Input.Nph_veto");
  
  
  //------------------
  // define histograms
  
  std::map<std::string,TH1F*> h;
  h["mgg"] = new TH1F("mgg","mgg",100,90.,160.);
  h["dipho_sumpt"] = new TH1F("dipho_sumpt","dipho_sumpt",100,0.,300);
  h["dipho_deltaeta"] = new TH1F("dipho_deltaeta","dipho_deltaeta",100,0,8);
  h["dipho_deltaphi"] = new TH1F("dipho_deltaphi","dipho_deltaphi",100,0,3.14);
  h["dipho_leadPt"] = new TH1F("dipho_leadPt","dipho_leadPt",100,0,200);
  h["dipho_leadEta"] = new TH1F("dipho_leadEta","dipho_leadEta",100,-3.5,3.5);
  h["dipho_leadPhi"] = new TH1F("dipho_leadPhi","dipho_leadPhi",100,-3.14,3.14);
  h["dipho_leadptoM"] = new TH1F("dipho_leadptoM","dipho_leadptoM",100,0,3.5);
  h["dipho_subleadPt"] = new TH1F("dipho_subleadPt","dipho_subleadPt",100,0,200);
  h["dipho_subleadEta"] = new TH1F("dipho_subleadEta","dipho_subleadEta",100,-3.5,3.5);
  h["dipho_subleadPhi"] = new TH1F("dipho_subleadPhi","dipho_subleadPhi",100,-3.14,3.14);
  h["dipho_subleadptoM"] = new TH1F("dipho_subleadptoM","dipho_subleadptoM",100,0,3.5);
  //h["dipho_leaddeltaR_GenReco"] = new TH1F("dipho_leaddeltaR_GenReco","dipho_leaddeltaR_GenReco",300,0.,0.1);
  //h["dipho_subleaddeltaR_GenReco"] = new TH1F("dipho_subleaddeltaR_GenReco","dipho_subleaddeltaR_GenReco",300,0.,0.1);

  h["nJets"] = new TH1F("nJets","nJets",11,-0.5,10.5);
  h["nJets_bTagLoose"] = new TH1F("nJets_bTagLoose","nJets_bTagLoose",11,-0.5,10.5);
  h["nJets_bTagMedium"] = new TH1F("nJets_bTagMedium","nJets_bTagMedium",11,-0.5,10.5);
  h["nJets_bTagTight"] = new TH1F("nJets_bTagTight","nJets_bTagTight",11,-0.5,10.5);
  h["mjj"] = new TH1F("mjj","mjj",100,60,200);
  h["dibjet_sumpt"] = new TH1F("dibjet_sumpt","dibjet_sumpt",100,0,300);
  h["dibjet_deltaeta"] = new TH1F("dibjet_deltaeta","dibjet_deltaeta",100,0,8);
  h["dibjet_deltaphi"] = new TH1F("dibjet_deltaphi","dibjet_deltaphi",100,0,3.14);
  h["dibjet_leadPt"] = new TH1F("dibjet_leadPt","dibjet_leadPt",100,0,200);
  h["dibjet_leadbtagmedium"] = new TH1F("dibjet_leadbtagmedium","dibjet_leadbtagmedium",9,-1.5,7.5);
  h["dibjet_leadEta"] = new TH1F("dibjet_leadEta","dibjet_leadEta",100,-3.5,3.5);
  h["dibjet_leadPhi"] = new TH1F("dibjet_leadPhi","dibjet_leadPhi",100,-3.14,3.14);
  h["dibjet_leadptoM"] = new TH1F("dibjet_leadptoM","dibjet_leadptoM",100,0,3.5);
  h["dibjet_leadEnergy"] = new TH1F("dibjet_leadEnergy","dibjet_leadEnergy",100,0,300);
  h["dibjet_leadbtagscore"] = new TH1F("dibjet_leadbtagscore","dibjet_leadbtagscore",9,-1.5,7.5);
  h["dibjet_subleadbtagmedium"] = new TH1F("dibjet_subleadbtagmedium","dibjet_subleadbtagmedium",9,-1.5,7.5);
  h["dibjet_subleadPt"] = new TH1F("dibjet_subleadPt","dibjet_subleadPt",100,0,200);
  h["dibjet_subleadEta"] = new TH1F("dibjet_subleadEta","dibjet_subleadEta",100,-3.5,3.5);
  h["dibjet_subleadPhi"] = new TH1F("dibjet_subleadPhi","dibjet_subleadPhi",100,-3.14,3.14);
  h["dibjet_subleadptoM"] = new TH1F("dibjet_subleadptoM","dibjet_subleadptoM",100,0,3.5);
  h["dibjet_subleadEnergy"] = new TH1F("dibjet_subleadEnergy","dibjet_subleadEnergy",100,0,300);
  h["dibjet_subleadbtagscore"] = new TH1F("dibjet_subleadbtagscore","dibjet_subleadbtagscore",9,-1.5,7.5);

  h["mtot"] = new TH1F("M_x","M_x",100,0,1000);
  h["DRmin_sel_phots_sel_jets"] = new TH1F("DeltaRmin","DeltaRmin",300,0,4);
  h["costhetastar_HH"] = new TH1F("costhetastar_HH","costhetastar_HH",200,0,1);
  h["costhetastar_gg"] = new TH1F("costhetastar_gg","costhetastar_gg",200,0,1);
  h["costhetastar_bb"] = new TH1F("costhetastar_bb","costhetastar_bb",200,0,1);

  h["MET"] = new TH1F("MET","MET",100,0,1000);
  h["MET_phi"] = new TH1F("MET_phi","MET_phi",100,-3.14,3.14);
  //h["dipho_leadIso_o_E"] = new TH1F("dipho_leadIso_o_E","dipho_leadIso_o_E",300,0.,0.25);
  //h["dipho_subleadIso_o_E"] = new TH1F("dipho_subleadIso_o_E","dipho_subleadIso_o_E",300,0.,0.25);
  //h["dipho_leaddeltaEta_GenReco"] = new TH1F("dipho_leaddeltaEta_GenReco","dipho_leaddeltaEta_GenReco",300,0.,0.2);
  //h["dipho_subleaddeltaEta_GenReco"] = new TH1F("dipho_subleaddeltaEta_GenReco","dipho_subleaddeltaEta_GenReco",300,0.,0.2);
  //h["dipho_leaddeltaPhi_GenReco"] = new TH1F("dipho_leaddeltaPhi_GenReco","dipho_leaddeltaPhi_GenReco",300,0.,0.2);
  //h["dipho_subleaddeltaPhi_GenReco"] = new TH1F("dipho_subleaddeltaPhi_GenReco","dipho_subleaddeltaPhi_GenReco",300,0.,0.2);
  
  
  
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
  TTree* outTree_all_lowMx  = new TTree("all_lowMx", "all_lowMx");
  TTree* outTree_all_highMx = new TTree("all_highMx","all_highMx");
  TreeVars outtreeVars;
  InitOutTreeVars(outTree_all_highMx,outtreeVars);
  InitOutTreeVars(outTree_all_lowMx, outtreeVars);
  
  //Get photon resolution from external file (HHGGBBAnalysis/data/Delphes_Ph_resolution.root)
  TFile* Ph_resolution_file = TFile::Open("/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/Delphes_Ph_resolution.root","READ");
  TH2F* ph_resolution_map;
  Ph_resolution_file->GetObject("rel_resolution_plus1_in_quad",ph_resolution_map);
  ph_resolution_map->SetDirectory(0);
  Ph_resolution_file->Close();
  outFile -> cd();

  //------------------
  // loop over samples
  int Nev_highMx_JCR=0;
  int Nev_highMx_MPC=0;
  int Nev_highMx_HPC=0;
  int Nev_lowMx_JCR=0;
  int Nev_lowMx_MPC=0;
  int Nev_lowMx_HPC=0;
  int Nev_all_lowMx=0;
  int Nev_all_highMx=0;
  
  int Nev_preselected=0;
  int Nev_phselection=0;
  int Nev_jet_kin_preselection=0;
  int Nev_jetselection=0;
  int Nev_Phpromptprompt=0;
  int Nev_Phpromptfake=0;
  int Nev_Phfakefake=0;
  int Nev_bjetpromptprompt=0;
  int Nev_bjetpromptfake=0;
  int Nev_bjetfakefake=0;
  int Nev_bquarkpromptprompt=0;
  int Nev_bquarkpromptfake=0;
  int Nev_bquarkfakefake=0;
  
  int nEntries = tree->GetEntries();
  std::cout << "Total entries = " << nEntries << std::endl;
  for(int i=0; i<nEntries; i++)
  {
    tree -> GetEntry(i);
    if( i%1000==0 ) std::cout << "Processing entry "<< i << "\r" << std::flush;
    
    if(treeVars.N_SelectedPh<2) continue;
    
    ++Nev_preselected;
    //find higgs daugher gen-photons and flag their .isHdaug
    // if( ! FindGenPh_Hdaug(treeVars) ) continue;
    
    
    //find lead & sublead reco photons    
    //PrintRecoPhoton(treeVars);
    int pho_lead_i;
    int pho_sublead_i;
    FindLeadSublead_pho(treeVars,pho_lead_i,pho_sublead_i);
    TLorentzVector pho_lead,pho_sublead;
    pho_lead.SetPtEtaPhiE(treeVars.SelectedPh_pt[pho_lead_i],treeVars.SelectedPh_eta[pho_lead_i],treeVars.SelectedPh_phi[pho_lead_i],treeVars.SelectedPh_E[pho_lead_i]);
    pho_sublead.SetPtEtaPhiE(treeVars.SelectedPh_pt[pho_sublead_i],treeVars.SelectedPh_eta[pho_sublead_i],treeVars.SelectedPh_phi[pho_sublead_i],treeVars.SelectedPh_E[pho_sublead_i]);
    TLorentzVector dipho = pho_lead+pho_sublead;
    
    //Gen-matching
    // if(!PhoGenMatch(pho_lead,pho_sublead,treeVars,outtreeVars,0.1))//default DeltaRmax=0.03
    //   continue;
    
    //Cuts on photons
    if(!DiPhotonSelection(pho_lead,pho_sublead)) continue;
    
    ++Nev_phselection;
    
    if(Nph_veto>=0)
    {
      //cout<<"veto loop"<<endl;
      TLorentzVector gen_pho_match;
      int Nph_gen_match=0;
      if (PhoGenericGenMatch(pho_lead,treeVars,gen_pho_match,0.1)) Nph_gen_match++;
      if(Nph_gen_match>Nph_veto) continue;
      if (PhoGenericGenMatch(pho_sublead,treeVars,gen_pho_match,0.1)) Nph_gen_match++;
      if(Nph_gen_match>Nph_veto) continue;
      //cout<<"passed veto loop"<<endl;
    }
    
    
    //Clean jet collection (to do? maybe required for bkg study...)
    //Tag good jets with a bool, requirements are:
    // 1. Not matching with any reco photon or any reco electron
    // 2. ratio Echarge / Eneutral  
    //TagGoodJets(treeVars);
    
    //Jets selections
    int BTagMedium_mask;
    if(useMTD == false)
      BTagMedium_mask=0b000010;
    else
      BTagMedium_mask=0b010000;
    
    //PrintRecoJet(treeVars);
    outtreeVars.nJets=0;
    outtreeVars.nJets_bTagLoose=0;
    outtreeVars.nJets_bTagMedium=0;	 
    outtreeVars.nJets_bTagTight=0;
    for(int i=0;i<treeVars.N_Jet;i++)
    {
      if(treeVars.Jet_pt[i]<25) continue;
      if(fabs(treeVars.Jet_eta[i])>2.4) continue;
      if( DeltaR(treeVars.Jet_eta[i],treeVars.Jet_phi[i],pho_lead.Eta(),pho_lead.Phi()) < 0.4 ) continue;
      if( DeltaR(treeVars.Jet_eta[i],treeVars.Jet_phi[i],pho_sublead.Eta(),pho_sublead.Phi()) < 0.4 ) continue;
      
      outtreeVars.nJets++;
      outtreeVars.jet_pt[int(outtreeVars.nJets)-1] = treeVars.Jet_pt[i];
      outtreeVars.jet_eta[int(outtreeVars.nJets)-1] = treeVars.Jet_eta[i];                    
      outtreeVars.jet_phi[int(outtreeVars.nJets)-1] = treeVars.Jet_phi[i];
      outtreeVars.jet_mass[int(outtreeVars.nJets)-1] = treeVars.Jet_mass[i];
      outtreeVars.jet_mvav2[int(outtreeVars.nJets)-1] = treeVars.Jet_mvav2[i];
      outtreeVars.jet_hadflav[int(outtreeVars.nJets)-1] = treeVars.Jet_hadflav[i];//gen level info! handle with care!
      
      if(useMTD)
      {
        if(outtreeVars.jet_mvav2[int(outtreeVars.nJets)-1] & 0b001000)
          outtreeVars.nJets_bTagLoose++;
        
        if(outtreeVars.jet_mvav2[int(outtreeVars.nJets)-1] & 0b010000)
        {
          outtreeVars.nJets_bTagMedium++;
          outtreeVars.jet_BTagMedium[int(outtreeVars.nJets)-1] = 1;
        }
        else
          outtreeVars.jet_BTagMedium[int(outtreeVars.nJets)-1] = 0;
        
        if(outtreeVars.jet_mvav2[int(outtreeVars.nJets)-1] & 0b100000)
          outtreeVars.nJets_bTagTight++;
      }
      else
      {
        if(outtreeVars.jet_mvav2[int(outtreeVars.nJets)-1] & 0b000001)
          outtreeVars.nJets_bTagLoose++;
        
        if(outtreeVars.jet_mvav2[int(outtreeVars.nJets)-1] & 0b000010)
        {
          outtreeVars.nJets_bTagMedium++;
          outtreeVars.jet_BTagMedium[int(outtreeVars.nJets)-1] = 1;
        }
        else
          outtreeVars.jet_BTagMedium[int(outtreeVars.nJets)-1] = 0;
        
        if(outtreeVars.jet_mvav2[int(outtreeVars.nJets)-1] & 0b000100)
          outtreeVars.nJets_bTagTight++;
      }
    }
    
    
    if(outtreeVars.nJets<2) continue;
    
    Nev_jet_kin_preselection++;
    
    
    //Select the two jets with the higher BTag level, if they have the same value select the harder one
    int bjet_lead_i;
    int bjet_sublead_i;
    SelectBestScoreBJets2(outtreeVars,bjet_lead_i,bjet_sublead_i,useMTD);
    TLorentzVector bjet_lead,bjet_sublead;
    bjet_lead.SetPtEtaPhiM(outtreeVars.jet_pt[bjet_lead_i],outtreeVars.jet_eta[bjet_lead_i],outtreeVars.jet_phi[bjet_lead_i],outtreeVars.jet_mass[bjet_lead_i]);
    bjet_sublead.SetPtEtaPhiM(outtreeVars.jet_pt[bjet_sublead_i],outtreeVars.jet_eta[bjet_sublead_i],outtreeVars.jet_phi[bjet_sublead_i],outtreeVars.jet_mass[bjet_sublead_i]);
    
    TLorentzVector dibjet= bjet_lead+bjet_sublead;
    double mjj = dibjet.M();
    
    if(mjj<80 || mjj>200)
      continue;
    
    ++Nev_jetselection;
    
    
    // count leptons
    outtreeVars.nEle = 0;
    outtreeVars.nMu  = 0;
    outtreeVars.nLep = 0;
    for(int i=0;i<treeVars.N_TightEl;i++)
    {
      if( treeVars.TightEl_pt[i]<10 ) continue;
      if( fabs(treeVars.TightEl_eta[i])>2.5 ) continue;
      if( fabs(treeVars.TightEl_eta[i])>1.44 && fabs(treeVars.TightEl_eta[i])<1.57 ) continue;
      if( DeltaR(treeVars.TightEl_eta[i],treeVars.TightEl_phi[i],pho_lead.Eta(),pho_lead.Phi()) < 0.4 ) continue;
      if( DeltaR(treeVars.TightEl_eta[i],treeVars.TightEl_phi[i],pho_sublead.Eta(),pho_sublead.Phi()) < 0.4 ) continue;
      if( DeltaR(treeVars.TightEl_eta[i],treeVars.TightEl_phi[i],bjet_lead.Eta(),bjet_lead.Phi()) < 0.4 ) continue;
      if( DeltaR(treeVars.TightEl_eta[i],treeVars.TightEl_phi[i],bjet_sublead.Eta(),bjet_sublead.Phi()) < 0.4 ) continue;
      
      outtreeVars.nEle += 1;
      outtreeVars.nLep += 1;
    }
    for(int i=0;i<treeVars.N_TightMu;i++)
    {
      if( treeVars.TightMu_pt[i]<5 ) continue;
      if( fabs(treeVars.TightMu_eta[i])>2.4 ) continue;
      if( DeltaR(treeVars.TightMu_eta[i],treeVars.TightMu_phi[i],pho_lead.Eta(),pho_lead.Phi()) < 0.4 ) continue;
      if( DeltaR(treeVars.TightMu_eta[i],treeVars.TightMu_phi[i],pho_sublead.Eta(),pho_sublead.Phi()) < 0.4 ) continue;
      if( DeltaR(treeVars.TightMu_eta[i],treeVars.TightMu_phi[i],bjet_lead.Eta(),bjet_lead.Phi()) < 0.4 ) continue;
      if( DeltaR(treeVars.TightMu_eta[i],treeVars.TightMu_phi[i],bjet_sublead.Eta(),bjet_sublead.Phi()) < 0.4 ) continue;
      
      outtreeVars.nMu += 1;
      outtreeVars.nLep += 1;
    }
    
    
    //Find dR_min between selected gamma and jets
    vector<TLorentzVector> pho_selected;
    pho_selected.push_back(pho_lead);
    pho_selected.push_back(pho_sublead);
    vector<TLorentzVector> bjet_selected;
    bjet_selected.push_back(bjet_lead);
    bjet_selected.push_back(bjet_sublead);
    float DeltaRmin_bjet_pho = DeltaRmin(pho_selected,bjet_selected);
    float DeltaPhimin_met_bjet = std::min(DeltaPhi(treeVars.Met_phi[0],bjet_lead.Phi()),DeltaPhi(treeVars.Met_phi[0],bjet_sublead.Phi()));
    float DeltaPhimax_met_bjet = std::max(DeltaPhi(treeVars.Met_phi[0],bjet_lead.Phi()),DeltaPhi(treeVars.Met_phi[0],bjet_sublead.Phi()));
    
    
    //Find interesting angles in the frame of dihigggs
    TLorentzVector diHiggs= pho_lead+pho_sublead+bjet_lead+bjet_sublead;   
    // get angle between dipho and z-axis in diHiggs rest frame
    TLorentzVector boosted_dipho(dipho);
    boosted_dipho.Boost( -diHiggs.BoostVector() );
    float costheta_HH = fabs(boosted_dipho.CosTheta());
    // get angle between leading photon and z-axiz in dipho rest frame
    TLorentzVector boosted_pho_lead(pho_lead);
    boosted_pho_lead.Boost( -dipho.BoostVector() );
    float costheta_gg = fabs(boosted_pho_lead.CosTheta());
    // get angle between leading jet and z-axiz in dibjet rest frame
    TLorentzVector boosted_bjet_lead(bjet_lead);
    boosted_bjet_lead.Boost( -dibjet.BoostVector() );
    float costheta_bb = fabs(boosted_bjet_lead.CosTheta());
    
    
    //Fill outtree and histos
    outtreeVars.evWeight = CrossSection*weightMC;
    outtreeVars.cross_sec = CrossSection;
    outtreeVars.event = treeVars.event;
    outtreeVars.nvtx = treeVars.N_Vtx;
    outtreeVars.mgg = (dipho).M();
    outtreeVars.dipho_sumpt = (dipho).Pt();
    outtreeVars.dipho_deltaeta = DeltaEta( pho_lead.Eta() , pho_sublead.Eta() );
    outtreeVars.dipho_deltaphi = DeltaPhi( pho_lead.Phi() , pho_sublead.Phi() );
    outtreeVars.dipho_leadPt = pho_lead.Pt();
    outtreeVars.dipho_leadEta = pho_lead.Eta();
    outtreeVars.dipho_leadPhi = pho_lead.Phi();
    outtreeVars.dipho_leadptoM = pho_lead.Pt() / (dipho).M();
    outtreeVars.dipho_leadEnergy = pho_lead.E();
    outtreeVars.dipho_lead_sigmaEoE = ph_resolution_map->GetBinContent(ph_resolution_map->FindBin(fabs(pho_lead.Eta()),pho_lead.E()));
    //outtreeVars.dipho_leadDeltaRgenreco = DeltaR(pho_lead.Eta(),pho_lead.Phi(),outtreeVars.dipho_leadEta_gen,outtreeVars.dipho_leadPhi_gen);
    //outtreeVars.dipho_leadDeltaEtagenreco = DeltaEta(pho_lead.Eta(),outtreeVars.dipho_leadEta_gen);
    //outtreeVars.dipho_leadDeltaPhigenreco = DeltaPhi(pho_lead.Phi(),outtreeVars.dipho_leadPhi_gen);
    outtreeVars.dipho_subleadPt = pho_sublead.Pt();
    outtreeVars.dipho_subleadEta = pho_sublead.Eta();
    outtreeVars.dipho_subleadPhi = pho_sublead.Phi();
    outtreeVars.dipho_subleadptoM = pho_sublead.Pt() / (dipho).M();
    outtreeVars.dipho_subleadEnergy = pho_sublead.E();
    outtreeVars.dipho_sublead_sigmaEoE = ph_resolution_map->GetBinContent(ph_resolution_map->FindBin(fabs(pho_sublead.Eta()),pho_sublead.E()));
    //outtreeVars.dipho_subleadDeltaRgenreco = DeltaR(pho_sublead.Eta(),pho_sublead.Phi(),outtreeVars.dipho_subleadEta_gen,outtreeVars.dipho_subleadPhi_gen);
    //outtreeVars.dipho_subleadDeltaEtagenreco = DeltaEta(pho_sublead.Eta(),outtreeVars.dipho_subleadEta_gen);
    //outtreeVars.dipho_subleadDeltaPhigenreco = DeltaPhi(pho_sublead.Phi(),outtreeVars.dipho_subleadPhi_gen);
    outtreeVars.mjj = (dibjet).M();
    outtreeVars.dibjet_sumpt = (dibjet).Pt();
    outtreeVars.dibjet_deltaeta = DeltaEta( bjet_lead.Eta() , bjet_sublead.Eta() );
    outtreeVars.dibjet_deltaphi = DeltaPhi( bjet_lead.Phi() , bjet_sublead.Phi() );
    outtreeVars.dibjet_leadPt = bjet_lead.Pt();
    outtreeVars.dibjet_leadEta = bjet_lead.Eta();
    outtreeVars.dibjet_leadPhi = bjet_lead.Phi();
    outtreeVars.dibjet_leadptoM = bjet_lead.Pt() / (dibjet).M();
    outtreeVars.dibjet_leadEnergy = bjet_lead.E();
    outtreeVars.dibjet_leadbtagmedium = outtreeVars.jet_BTagMedium[bjet_lead_i];
    outtreeVars.dibjet_leadmvav2 = outtreeVars.jet_mvav2[bjet_lead_i];
    outtreeVars.dibjet_subleadPt = bjet_sublead.Pt();
    outtreeVars.dibjet_subleadEta = bjet_sublead.Eta();
    outtreeVars.dibjet_subleadPhi = bjet_sublead.Phi();
    outtreeVars.dibjet_subleadptoM = bjet_sublead.Pt() / (dibjet).M();
    outtreeVars.dibjet_subleadEnergy = bjet_sublead.E();
    outtreeVars.dibjet_subleadbtagmedium = outtreeVars.jet_BTagMedium[bjet_sublead_i];
    outtreeVars.dibjet_subleadmvav2 = outtreeVars.jet_mvav2[bjet_sublead_i];
    outtreeVars.mtot = diHiggs.M() - outtreeVars.mjj - outtreeVars.mgg + 250.;
    outtreeVars.DRmin_pho_bjet = DeltaRmin_bjet_pho; 
    outtreeVars.DPhimin_met_bjet = DeltaPhimin_met_bjet;
    outtreeVars.DPhimax_met_bjet = DeltaPhimax_met_bjet;
    outtreeVars.costheta_HH = costheta_HH; 
    outtreeVars.costheta_gg = costheta_gg; 
    outtreeVars.costheta_bb = costheta_bb; 
    outtreeVars.MetPt = treeVars.Met_pt[0];
    outtreeVars.MetPhi = treeVars.Met_phi[0];
    
    outtreeVars.ttHTagger = 0;
    
    h["mgg"] -> Fill(outtreeVars.mgg);
    h["dipho_sumpt"] -> Fill(outtreeVars.dipho_sumpt);
    h["dipho_deltaeta"] -> Fill(outtreeVars.dipho_deltaeta);
    h["dipho_deltaphi"] -> Fill(outtreeVars.dipho_deltaphi);
    h["dipho_leadPt"] -> Fill(outtreeVars.dipho_leadPt);
    h["dipho_leadEta"] -> Fill(outtreeVars.dipho_leadEta);
    h["dipho_leadPhi"] -> Fill(outtreeVars.dipho_leadPhi);
    h["dipho_leadptoM"] -> Fill(outtreeVars.dipho_leadptoM);
    h["dipho_subleadPt"] -> Fill(outtreeVars.dipho_subleadPt);
    h["dipho_subleadEta"] -> Fill(outtreeVars.dipho_subleadEta);
    h["dipho_subleadPhi"] -> Fill(outtreeVars.dipho_subleadPhi);
    h["dipho_subleadptoM"] -> Fill(outtreeVars.dipho_subleadptoM);
    
    h["nJets"] -> Fill(outtreeVars.nJets);
    h["nJets_bTagLoose"] -> Fill(outtreeVars.nJets_bTagLoose);
    h["nJets_bTagMedium"] -> Fill(outtreeVars.nJets_bTagMedium);
    h["nJets_bTagTight"] -> Fill(outtreeVars.nJets_bTagTight);
    h["mjj"] -> Fill(outtreeVars.mjj);
    h["dibjet_sumpt"] -> Fill(outtreeVars.dibjet_sumpt);
    h["dibjet_deltaeta"] -> Fill(outtreeVars.dibjet_deltaeta);
    h["dibjet_deltaphi"] -> Fill(outtreeVars.dibjet_deltaphi);

    h["dibjet_leadPt"] -> Fill(outtreeVars.dibjet_leadPt);
    h["dibjet_leadEta"] -> Fill(outtreeVars.dibjet_leadEta);
    h["dibjet_leadPhi"] -> Fill(outtreeVars.dibjet_leadPhi);
    h["dibjet_leadptoM"] -> Fill(outtreeVars.dibjet_leadptoM);
    h["dibjet_leadEnergy"] -> Fill(outtreeVars.dibjet_leadEnergy);
    h["dibjet_leadbtagmedium"] -> Fill(outtreeVars.dibjet_leadbtagmedium);
    
    h["dibjet_subleadPt"] -> Fill(outtreeVars.dibjet_subleadPt);
    h["dibjet_subleadEta"] -> Fill(outtreeVars.dibjet_subleadEta);
    h["dibjet_subleadPhi"] -> Fill(outtreeVars.dibjet_subleadPhi);
    h["dibjet_subleadptoM"] -> Fill(outtreeVars.dibjet_subleadptoM);
    h["dibjet_subleadEnergy"] -> Fill(outtreeVars.dibjet_subleadEnergy);
    h["dibjet_subleadbtagmedium"] -> Fill(outtreeVars.dibjet_subleadbtagmedium);
    
    h["mtot"] -> Fill(outtreeVars.mtot);
    h["DRmin_sel_phots_sel_jets"] -> Fill(outtreeVars.DRmin_pho_bjet); 
    h["costhetastar_HH"] -> Fill(outtreeVars.costheta_HH); 
    h["costhetastar_gg"] -> Fill(outtreeVars.costheta_gg); 
    h["costhetastar_bb"] -> Fill(outtreeVars.costheta_bb); 
    h["MET"] -> Fill(outtreeVars.MetPt); 
    h["MET_phi"] -> Fill(outtreeVars.MetPhi); 
    
    //JCR(jet control region): less than one jet with medium b-tag  
    //MPC: exactly one jet with medium b-tag
    //HPC: at least two jets with medium b-tag
    
    if( (outtreeVars.dibjet_leadmvav2 & BTagMedium_mask) && (outtreeVars.dibjet_subleadmvav2 & BTagMedium_mask) )
    {
      //cout<<"HPC"<<endl;
      outtreeVars.cut_based_ct = 0;
    }
    else if( ( (outtreeVars.dibjet_leadmvav2 & BTagMedium_mask)  && !(outtreeVars.dibjet_subleadmvav2 & BTagMedium_mask) ) ||
             ( !(outtreeVars.dibjet_leadmvav2 & BTagMedium_mask) &&  (outtreeVars.dibjet_subleadmvav2 & BTagMedium_mask) ) ) //must be after HPC!!!
    {
      //cout<<"MPC"<<endl;
      outtreeVars.cut_based_ct = 1;
    }
    else if( !(outtreeVars.dibjet_leadmvav2 & BTagMedium_mask) && !(outtreeVars.dibjet_subleadmvav2 & BTagMedium_mask) )
    {
      //cout<<"JCR"<<endl;
      outtreeVars.cut_based_ct = -1;
    }
    
    
    //fill low mass and high mass categories
    if(outtreeVars.mtot<=350)
    {
      ++Nev_all_lowMx;
      if(outtreeVars.cut_based_ct == -1) ++Nev_lowMx_JCR;
      if(outtreeVars.cut_based_ct ==  0) ++Nev_lowMx_HPC;
      if(outtreeVars.cut_based_ct ==  1) ++Nev_lowMx_MPC;
      outTree_all_lowMx->Fill();
    }
    else
    {
      ++Nev_all_highMx;
      if(outtreeVars.cut_based_ct == -1) ++Nev_highMx_JCR;
      if(outtreeVars.cut_based_ct ==  0) ++Nev_highMx_HPC;
      if(outtreeVars.cut_based_ct ==  1) ++Nev_highMx_MPC;
      outTree_all_highMx->Fill();
    }
  }

  std::cout << std::endl;

  cout<<"\n-----------------------------------------------------------------"<<endl;
  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"MC events = "<<NEventsMC<<endl;
  cout<<"2b2g skim = "<<nEntries<<" / "<<NEventsMC<<endl;
  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"\t\tpreselection acceptance = "<<Nev_preselected<<" / "<<NEventsMC<<" = "<<Nev_preselected/NEventsMC<<endl;
  cout<<"\t\tpass photon preselection = "<<Nev_phselection<<" / "<<nEntries<<endl;
  cout<<"\t\tpass kinematic jet preselection = "<< Nev_jet_kin_preselection <<" / "<<nEntries<<endl;
  cout<<"\t\tpass jet preselection = "<<Nev_jetselection<<" / "<<nEntries<<endl;
  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"\t\t\t-----------------------------------------------------------------"<<endl;  
  cout<<"\t\t\thigh mass jet control region acceptance = "<<Nev_highMx_JCR<<" / "<<NEventsMC<<" = "<<1.*Nev_highMx_JCR/NEventsMC<<endl;
  cout<<"\t\t\thigh mass medium purity cat. acceptance = "<<Nev_highMx_MPC<<" / "<<NEventsMC<<" = "<<1.*Nev_highMx_MPC/NEventsMC<<endl;
  cout<<"\t\t\thigh mass high purity cat. acceptance = "<<Nev_highMx_HPC<<" / "<<NEventsMC<<" = "<<1.*Nev_highMx_HPC/NEventsMC<<endl;
  cout<<"\t\t\t-----------------------------------------------------------------"<<endl;
  cout<<"\t\t\tlow mass jet control region acceptance = "<<Nev_lowMx_JCR<<" / "<<NEventsMC<<" = "<<1.*Nev_lowMx_JCR/NEventsMC<<endl;
  cout<<"\t\t\tlow mass medium purity cat. acceptance = "<<Nev_lowMx_MPC<<" / "<<NEventsMC<<" = "<<1.*Nev_lowMx_MPC/NEventsMC<<endl;
  cout<<"\t\t\tlow mass high purity cat. acceptance = "<<Nev_lowMx_HPC<<" / "<<NEventsMC<<" = "<<1.*Nev_lowMx_HPC/NEventsMC<<endl;
  cout<<"\t\t\t\tdiphoton promptprompt = "<<Nev_Phpromptprompt<<" / "<<NEventsMC<<" = "<<1.*Nev_Phpromptprompt/NEventsMC<<endl;
  cout<<"\t\t\t\tdiphoton promptfake = "<<Nev_Phpromptfake<<" / "<<NEventsMC<<" = "<<1.*Nev_Phpromptfake/NEventsMC<<endl;
  cout<<"\t\t\t\tdiphoton fakefake = "<<Nev_Phfakefake<<" / "<<NEventsMC<<" = "<<1.*Nev_Phfakefake/NEventsMC<<endl;
  cout<<"\t\t\t\tdibjet promptprompt = "<<Nev_bjetpromptprompt<<" / "<<NEventsMC<<" = "<<1.*Nev_bjetpromptprompt/NEventsMC<<endl;
  cout<<"\t\t\t\tdibjet promptfake = "<<Nev_bjetpromptfake<<" / "<<NEventsMC<<" = "<<1.*Nev_bjetpromptfake/NEventsMC<<endl;
  cout<<"\t\t\t\tdibjet fakefake = "<<Nev_bjetfakefake<<" / "<<NEventsMC<<" = "<<1.*Nev_bjetfakefake/NEventsMC<<endl;
  cout<<"\t\t\t\tdibquark promptprompt = "<<Nev_bquarkpromptprompt<<" / "<<NEventsMC<<" = "<<1.*Nev_bquarkpromptprompt/NEventsMC<<endl;
  cout<<"\t\t\t\tdibquark promptfake = "<<Nev_bquarkpromptfake<<" / "<<NEventsMC<<" = "<<1.*Nev_bquarkpromptfake/NEventsMC<<endl;
  cout<<"\t\t\t\tdibquark fakefake = "<<Nev_bquarkfakefake<<" / "<<NEventsMC<<" = "<<1.*Nev_bquarkfakefake/NEventsMC<<endl;



  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"N_preselected * XS /N_MC = "<< Nev_preselected*CrossSection/NEventsMC <<" fb"<<endl; 
  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"-----------------------------------------------------------------\n"<<endl;

  outTree_all_lowMx -> AutoSave();
  outTree_all_highMx -> AutoSave();
  outFile -> Close();
  
  // MakePlot3(h);
  // system(Form("mv *.png %s",outputPlotFolder.c_str()));
  // system(Form("mv *.pdf %s",outputPlotFolder.c_str()));
  
  
  
  return 0;
}
