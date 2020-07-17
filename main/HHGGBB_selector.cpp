#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/CMS_lumi.h"
#include "interface/SetTDRStyle.h"
#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"
#include "interface/RooFitUtils.h"
#include "interface/HH_reweight_components.h"

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
#include "TRandom3.h"

#include "RooMsgService.h"
#include "RooRealVar.h"

using namespace std;
TRandom3 rndm;
void generatenewBtag(RawTreeVars &treeVars, map<int,TH2F*> &eff_map_flav4, map<int,TH2F*> &eff_map_flav5);



int main(int argc, char* argv[])
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> HHGGBB_selector.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }
  
  rndm.SetSeed(15);
  //rndm.SetSeed();
  
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

  float default_CrossSection = 1.;
  if(opts.OptExist("Input.CrossSection"))
    default_CrossSection = opts.GetOpt<float> ("Input.CrossSection");
  
  float NEventsMC = 1.;
  if(opts.OptExist("Input.NEventsMC"))
    NEventsMC = opts.GetOpt<float> ("Input.NEventsMC");
  
  int Nph_veto = -1;//discard events with N selected photon gen-matched > Nph_veto --> Useful to avoid double counting in GJet and QCD samples   
  if(opts.OptExist("Input.Nph_veto"))
    Nph_veto = opts.GetOpt<int> ("Input.Nph_veto");

  bool do_klambda_reweight = false;
  HH_reweight_components* HHrew = NULL;
  float klambda=1;
  float BSM_CrossSection=1.;
  float weight_sum=1.;
  if(opts.OptExist("Input.klambda"))
  {
    do_klambda_reweight = true;
    klambda = opts.GetOpt<float> ("Input.klambda");
    string reweight_filename = opts.GetOpt<string> ("Input.klambda_reweightfile");
    HHrew = new HH_reweight_components(reweight_filename.c_str());
    BSM_CrossSection = opts.GetOpt<float> ("Input.BSM_CrossSection");
    weight_sum = opts.GetOpt<float> ("Input.weight_sum");
  }

  //----------------------
  // load btag maps
  map<int,TH2F*> eff_map_flav5;
  map<int,TH2F*> eff_map_flav4;
  TFile *effmapfile=0;

  if(useMTD)
  {
    cout<<"loading bTag map from file /afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/neweffmap.root"<<endl;
    effmapfile= new TFile("/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/neweffmap.root","READ");
  }
  else
  {
    cout<<"loading bTag map from file /afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/neweffmap_NOMTD.root"<<endl;
    effmapfile= new TFile("/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/neweffmap_NOMTD.root","READ");
  }

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


  //Reweight for new photonID efficiency                                                                                                                    
  TH2F* photonIDeff_reweightmap = 0;
  float Ptmax_eff_reweightmap = 0;
  if(opts.OptExist("Input.photonIDeff_reweightfile"))
  {
    string photonIDeff_reweight_filename = opts.GetOpt<string> ("Input.photonIDeff_reweightfile");
    TFile *photonIDeff_reweight_file = new TFile(photonIDeff_reweight_filename.c_str());
    cout<<"loading photonIDeff reweight map from file "<<photonIDeff_reweight_filename<<endl;
    photonIDeff_reweightmap = (TH2F*) photonIDeff_reweight_file->Get("eff_PhotonIsoSF");
    if(!photonIDeff_reweightmap)
    {
      cout<<"[ERROR]: photonIDeff reweight map not found"<<endl;
      return -1;
    }
    else
    {
      photonIDeff_reweightmap->SetDirectory(0);
      photonIDeff_reweight_file->Close();
      Ptmax_eff_reweightmap = photonIDeff_reweightmap->GetXaxis()->GetXmax();
    }
  }

  //Reweight for new photonID fake rate                                                                            
  TH2F* photonIDfake_reweightmap = 0;
  float Ptmax_fake_reweightmap = 0;
  if(opts.OptExist("Input.photonIDfake_reweightfile"))
  {
    string photonIDfake_reweight_filename = opts.GetOpt<string> ("Input.photonIDfake_reweightfile");
    TFile *photonIDfake_reweight_file = new TFile(photonIDfake_reweight_filename.c_str());
    cout<<"loading photonIDfake reweight map from file "<<photonIDfake_reweight_filename<<endl;
    photonIDfake_reweightmap = (TH2F*) photonIDfake_reweight_file->Get("frate_PhotonIsoSF");
    if(!photonIDfake_reweightmap)
    {
      cout<<"[ERROR]: photonIDfake reweight map not found"<<endl;
      return -1;
    }
    else
    {
      photonIDfake_reweightmap->SetDirectory(0);
      photonIDfake_reweight_file->Close();
      Ptmax_fake_reweightmap = photonIDfake_reweightmap->GetXaxis()->GetXmax();
    }
  }

  //Reweight for new btag efficiencies                                                                            
  TH2F* btaglooseeff_reweightmap = 0;
  TH2F* btagmediumeff_reweightmap = 0;
  TH2F* btagtighteff_reweightmap = 0;
  TH2F* btagloosefake_reweightmap = 0;
  TH2F* btagmediumfake_reweightmap = 0;
  TH2F* btagtightfake_reweightmap = 0;
  if(opts.OptExist("Input.btag_reweightfile"))
  {
    string btag_reweight_filename = opts.GetOpt<string> ("Input.btag_reweightfile");
    TFile *btag_reweight_file = new TFile(btag_reweight_filename.c_str());
    cout<<"loading btag reweight map from file "<<btag_reweight_filename<<endl;
    btaglooseeff_reweightmap = (TH2F*) btag_reweight_file->Get("eff_BTagSFLoose");
    btagmediumeff_reweightmap = (TH2F*) btag_reweight_file->Get("eff_BTagSFMedium");
    btagtighteff_reweightmap = (TH2F*) btag_reweight_file->Get("eff_BTagSFTight");
    btagloosefake_reweightmap = (TH2F*) btag_reweight_file->Get("frate_BTagSFLoose");
    btagmediumfake_reweightmap = (TH2F*) btag_reweight_file->Get("frate_BTagSFMedium");
    btagtightfake_reweightmap = (TH2F*) btag_reweight_file->Get("frate_BTagSFTight");
    if(!btaglooseeff_reweightmap   || 
       !btagmediumeff_reweightmap  || 
       !btagtighteff_reweightmap   || 
       !btagloosefake_reweightmap  || 
       !btagmediumfake_reweightmap || 
       !btagtightfake_reweightmap    )
    {
      cout<<"[ERROR]: some btag efficiency or fake rate reweight map not found"<<endl;
      return -1;
    }
    else
    {
      btaglooseeff_reweightmap->SetDirectory(0);
      btagmediumeff_reweightmap->SetDirectory(0);  
      btagtighteff_reweightmap->SetDirectory(0);   
      btagloosefake_reweightmap->SetDirectory(0);  
      btagmediumfake_reweightmap->SetDirectory(0); 
      btagtightfake_reweightmap->SetDirectory(0);    
      btag_reweight_file->Close();
    }
  }
    
  
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

  int Nbadph=0;
  int Nfakeph=0;
  
  int Nev_preselected=0;
  int Nev_phselection=0;
  int Nev_jet_kin_preselection=0;
  int Nev_jetselection=0;
  int Nev_Phpromptprompt=0;
  int Nev_Phpromptfake=0;
  int Nev_Phfakefake=0;
  int Nev_bjetpromptprompt=0;
  int Nev_bjetpromptfake=0;
  int Nev_cjetpromptprompt=0;
  int Nev_cjetpromptfake=0;
  int Nev_jetfakefake=0;
  int Nev_bcjet=0;
  int Nev_bquarkpromptprompt=0;
  int Nev_bquarkpromptfake=0;
  int Nev_bquarkfakefake=0;

  int Nev_withbquarkHdaugh=0;
  int bquarkmatch=0;

  int nEntries = tree->GetEntries();
  std::cout << "Total entries = " << nEntries << std::endl;
  for(int i=0; i<nEntries; i++)
  {
    tree -> GetEntry(i);
    if( i%1000==0 ) std::cout << "Processing entry "<< i << "\r" << std::flush;

    float CrossSection = default_CrossSection;
    float weightMC = 1./NEventsMC;

    if(treeVars.N_SelectedPh<2) continue;
    
    ++Nev_preselected;
    //find higgs daugher gen-photons and flag their .isHdaug
    if(do_klambda_reweight)    
      if( ! FindGenPh_Hdaug(treeVars) ) continue;
    
    
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

    //Reweight gen-matched photons for MTD efficiency and reweight NOT-gen-matched photons for MTD fake-rate                                      
    int Nph_gen_match=0;
    TLorentzVector gen_pho_match;

    //cout<<"LEAD PHOTON\n\tPt="<< pho_lead.Pt()<<"\tEta="<< pho_lead.Eta();
    if (PhoGenericGenMatch(pho_lead,treeVars,gen_pho_match,0.1)) 
    {
      Nph_gen_match++;
      //cout<<"\tphoton GEN-MATCH "<<endl;
      //cout<<"\teta"<<pho_lead.Eta()-gen_pho_match.Eta()<<"\tphi"<<pho_lead.Phi()-gen_pho_match.Phi()<<"\tPt"<<pho_lead.Pt()-gen_pho_match.Pt()<<endl;
      if(photonIDeff_reweightmap)
      {
	float pho_pt = std::min((float)pho_lead.Pt(),Ptmax_eff_reweightmap-1);
	weightMC *= photonIDeff_reweightmap->GetBinContent( photonIDeff_reweightmap->FindBin(pho_pt,fabs(pho_lead.Eta())));
	//cout<<"\tSF="<<photonIDeff_reweightmap->GetBinContent( photonIDeff_reweightmap->FindBin(pho_pt,fabs(pho_lead.Eta())))<<endl;
      }
    }
    else
    {
      //cout<<"\tNOT photon GEN-MATCH"<<endl;
      if(photonIDfake_reweightmap)
      {
	float pho_pt = std::min((float)pho_lead.Pt(),Ptmax_fake_reweightmap-1);
	weightMC *= photonIDfake_reweightmap->GetBinContent( photonIDfake_reweightmap->FindBin(pho_pt,fabs(pho_lead.Eta())));
	//cout<<"\tSF="<<photonIDfake_reweightmap->GetBinContent( photonIDfake_reweightmap->FindBin(pho_pt,fabs(pho_lead.Eta())))<<endl;
      }
    }

    //cout<<"SUBLEAD PHOTON\n\tPt="<< pho_sublead.Pt()<<"\tEta="<< pho_sublead.Eta();
    if (PhoGenericGenMatch(pho_sublead,treeVars,gen_pho_match,0.1)) 
    {
      Nph_gen_match++;
      //cout<<"\tphoton GEN-MATCH "<<endl;
      //cout<<"\teta"<<pho_sublead.Eta()-gen_pho_match.Eta()<<"\tphi"<<pho_sublead.Phi()-gen_pho_match.Phi()<<"\tPt"<<pho_sublead.Pt()-gen_pho_match.Pt()<<endl;
      if(photonIDeff_reweightmap)
      {
	float pho_pt = std::min((float)pho_sublead.Pt(),Ptmax_eff_reweightmap-1);
	weightMC *= photonIDeff_reweightmap->GetBinContent( photonIDeff_reweightmap->FindBin(pho_pt,fabs(pho_sublead.Eta())));
        //cout<<"\tSF="<<photonIDeff_reweightmap->GetBinContent( photonIDeff_reweightmap->FindBin(pho_pt,fabs(pho_sublead.Eta())))<<endl;
      }
    }
    else
    {
      //cout<<"\tNOT photon GEN-MATCH"<<endl;
      if(photonIDfake_reweightmap)
      {
	float pho_pt = std::min((float)pho_sublead.Pt(),Ptmax_fake_reweightmap-1);
	weightMC *= photonIDfake_reweightmap->GetBinContent( photonIDfake_reweightmap->FindBin(pho_pt,fabs(pho_sublead.Eta())));
        //cout<<"\tSF="<<photonIDfake_reweightmap->GetBinContent( photonIDfake_reweightmap->FindBin(pho_pt,fabs(pho_sublead.Eta())))<<endl;
      }
    }


    if(Nph_veto>=0 && Nph_gen_match>Nph_veto) continue;

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

    //reparametrize the btag
    generatenewBtag(treeVars,eff_map_flav4,eff_map_flav5);
    
    //PrintRecoJet(treeVars);
    outtreeVars.nJets=0;
    outtreeVars.nJets_bTagLoose=0;
    outtreeVars.nJets_bTagMedium=0;	 
    outtreeVars.nJets_bTagTight=0;
    for(int i=0;i<treeVars.N_Jet;i++)
    {
      if(treeVars.Jet_pt[i]<25) continue;
      if(fabs(treeVars.Jet_eta[i])>2.4) continue;
      //if(treeVars.Jet_pt[i]<20) continue;
      //if(fabs(treeVars.Jet_eta[i])>3) continue;
      if( DeltaR(treeVars.Jet_eta[i],treeVars.Jet_phi[i],pho_lead.Eta(),pho_lead.Phi()) < 0.4 ) continue;
      if( DeltaR(treeVars.Jet_eta[i],treeVars.Jet_phi[i],pho_sublead.Eta(),pho_sublead.Phi()) < 0.4 ) continue;
      
      outtreeVars.nJets++;
      outtreeVars.jet_pt[int(outtreeVars.nJets)-1] = treeVars.Jet_pt[i];
      outtreeVars.jet_eta[int(outtreeVars.nJets)-1] = treeVars.Jet_eta[i];                    
      outtreeVars.jet_phi[int(outtreeVars.nJets)-1] = treeVars.Jet_phi[i];
      outtreeVars.jet_mass[int(outtreeVars.nJets)-1] = treeVars.Jet_mass[i];
      outtreeVars.jet_mvav2[int(outtreeVars.nJets)-1] = treeVars.Jet_mvav2[i];
      outtreeVars.jet_hadflav[int(outtreeVars.nJets)-1] = treeVars.Jet_hadflav[i];//gen level info! handle with care!
      
      if(outtreeVars.jet_mvav2[int(outtreeVars.nJets)-1] >= 4)
      {
	outtreeVars.nJets_bTagLoose++;
	outtreeVars.jet_BTagLoose[int(outtreeVars.nJets)-1] = 1;
      }
      else
	outtreeVars.jet_BTagLoose[int(outtreeVars.nJets)-1] = 0;

      if(outtreeVars.jet_mvav2[int(outtreeVars.nJets)-1] >= 5)
      {
	outtreeVars.nJets_bTagMedium++;
	outtreeVars.jet_BTagMedium[int(outtreeVars.nJets)-1] = 1;
      }
      else
	outtreeVars.jet_BTagMedium[int(outtreeVars.nJets)-1] = 0;
      
      if(outtreeVars.jet_mvav2[int(outtreeVars.nJets)-1] == 6)
      {
	outtreeVars.nJets_bTagTight++;
	outtreeVars.jet_BTagTight[int(outtreeVars.nJets)-1] = 1;
      }
      else
	outtreeVars.jet_BTagTight[int(outtreeVars.nJets)-1] = 1;
    }

    if(outtreeVars.nJets<2) continue;
    
    Nev_jet_kin_preselection++;
    
    
    //Select the two jets with the higher BTag level, if they have the same value select the harder one
    int bjet_lead_i;
    int bjet_sublead_i;
    SelectBestScoreBJets3(outtreeVars,bjet_lead_i,bjet_sublead_i);
    TLorentzVector bjet_lead,bjet_sublead;
    bjet_lead.SetPtEtaPhiM(outtreeVars.jet_pt[bjet_lead_i],outtreeVars.jet_eta[bjet_lead_i],outtreeVars.jet_phi[bjet_lead_i],outtreeVars.jet_mass[bjet_lead_i]);
    bjet_sublead.SetPtEtaPhiM(outtreeVars.jet_pt[bjet_sublead_i],outtreeVars.jet_eta[bjet_sublead_i],outtreeVars.jet_phi[bjet_sublead_i],outtreeVars.jet_mass[bjet_sublead_i]);
    
    TLorentzVector dibjet= bjet_lead+bjet_sublead;
    double mjj = dibjet.M();
    
    if(mjj<80 || mjj>200)
      continue;
    
    ++Nev_jetselection;
    
    //reweight for MTD btag efficiencies
    outtreeVars.dibjet_leadPt = bjet_lead.Pt();
    outtreeVars.dibjet_leadEta = bjet_lead.Eta();
    outtreeVars.dibjet_leadmvav2 = outtreeVars.jet_mvav2[bjet_lead_i];
    outtreeVars.dibjet_leadbtaglevel = outtreeVars.dibjet_leadmvav2;
    outtreeVars.dibjet_leadgenflav = outtreeVars.jet_hadflav[bjet_lead_i];
    outtreeVars.dibjet_subleadPt = bjet_sublead.Pt();
    outtreeVars.dibjet_subleadEta = bjet_sublead.Eta();
    outtreeVars.dibjet_subleadmvav2 = outtreeVars.jet_mvav2[bjet_sublead_i];
    outtreeVars.dibjet_subleadbtaglevel = outtreeVars.dibjet_subleadmvav2;
    outtreeVars.dibjet_subleadgenflav = outtreeVars.jet_hadflav[bjet_sublead_i];

    btagReweight(outtreeVars,weightMC,
		 btaglooseeff_reweightmap, btagmediumeff_reweightmap, btagtighteff_reweightmap,
		 btagloosefake_reweightmap,btagmediumfake_reweightmap,btagtightfake_reweightmap);

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

    //find mHH gen and costhetaHH* gen for klambda scan
    if(do_klambda_reweight)
    {
      if(!FindHHGen(treeVars,outtreeVars))
	cout<<"[WARNING]:gen Higgs not found"<<endl;
      //cout<<"reweight("<<outtreeVars.mHH_gen<<","<<outtreeVars.costhetaHH_gen<<")="<<HHrew->get_weight(klambda, outtreeVars.mHH_gen, outtreeVars.costhetaHH_gen)<<endl;
      //cout<<"weightSUM="<< weight_sum<<endl;
      //cout<<"BSM/SM XS = "<<BSM_CrossSection<<endl;
      CrossSection *= HHrew->get_weight(klambda, outtreeVars.mHH_gen, outtreeVars.costhetaHH_gen) / (weight_sum/nEntries) * BSM_CrossSection;
      //cout<<"CrossSection="<<CrossSection<<endl;
      //cout<<"--------------"<<endl;
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
    //cout<<"SF_tot="<<weightMC*NEventsMC<<endl;
    //getchar();
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
    outtreeVars.dibjet_leadPhi = bjet_lead.Phi();
    outtreeVars.dibjet_leadptoM = bjet_lead.Pt() / (dibjet).M();
    outtreeVars.dibjet_leadEnergy = bjet_lead.E();
    outtreeVars.dibjet_leadbtagloose = outtreeVars.jet_BTagLoose[bjet_lead_i];
    outtreeVars.dibjet_leadbtagmedium = outtreeVars.jet_BTagMedium[bjet_lead_i];
    outtreeVars.dibjet_leadbtagtight = outtreeVars.jet_BTagTight[bjet_lead_i];

    outtreeVars.dibjet_subleadPt = bjet_sublead.Pt();
    outtreeVars.dibjet_subleadEta = bjet_sublead.Eta();
    outtreeVars.dibjet_subleadPhi = bjet_sublead.Phi();
    outtreeVars.dibjet_subleadptoM = bjet_sublead.Pt() / (dibjet).M();
    outtreeVars.dibjet_subleadEnergy = bjet_sublead.E();
    outtreeVars.dibjet_subleadbtagloose = outtreeVars.jet_BTagLoose[bjet_sublead_i];
    outtreeVars.dibjet_subleadbtagmedium = outtreeVars.jet_BTagMedium[bjet_sublead_i];
    outtreeVars.dibjet_subleadbtagtight = outtreeVars.jet_BTagTight[bjet_sublead_i];

    outtreeVars.diHiggs_Eta      = diHiggs.Eta();  
    outtreeVars.diHiggs_Pt       = diHiggs.Pt();
    outtreeVars.diHiggs_Energy   = diHiggs.Energy();
    outtreeVars.diHiggs_Rapidity = diHiggs.Rapidity();

    outtreeVars.mtot = diHiggs.M() - outtreeVars.mjj - outtreeVars.mgg + 250.;
    outtreeVars.DRmin_pho_bjet = DeltaRmin_bjet_pho; 
    outtreeVars.DPhimin_met_bjet = DeltaPhimin_met_bjet;
    outtreeVars.DPhimax_met_bjet = DeltaPhimax_met_bjet;
    outtreeVars.costheta_HH = costheta_HH; 
    outtreeVars.costheta_gg = costheta_gg; 
    outtreeVars.costheta_bb = costheta_bb; 
    outtreeVars.MetPt = treeVars.Met_pt[0];
    outtreeVars.MetPhi = treeVars.Met_phi[0];
    
    
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

    /*
    /////////////////////////////////////////////////////////////////////////
    //temporary for bjet selection efficiency
    if( ! Findbquark_Hdaug(treeVars) ) continue;
    Nev_withbquarkHdaugh++;
    //store b-quarks higgs daugher
    vector<TLorentzVector> bquark_gen;
    for(int i=0;i<treeVars.N_GenPart;i++)
      if(treeVars.GenPart_isHdaug[i]==true && fabs(treeVars.GenPart_pid[i])==5)
      {
	float gen_E=sqrt(treeVars.GenPart_pt[i]*treeVars.GenPart_pt[i] + treeVars.GenPart_pz[i]*treeVars.GenPart_pz[i]+125.*125);    
	TLorentzVector bquark_gen_aux;
	bquark_gen_aux.SetPtEtaPhiE(treeVars.GenPart_pt[i],treeVars.GenPart_eta[i],treeVars.GenPart_phi[i],gen_E);
	bquark_gen.push_back(bquark_gen_aux);
      }
    outtreeVars.bquark_eta1=bquark_gen.at(0).Eta();
    outtreeVars.bquark_phi1=bquark_gen.at(0).Phi();
    outtreeVars.bquark_E1=bquark_gen.at(0).E();
    outtreeVars.bquark_pt1=bquark_gen.at(0).Pt();
    outtreeVars.bquark_eta2=bquark_gen.at(1).Eta();
    outtreeVars.bquark_phi2=bquark_gen.at(1).Phi();
    outtreeVars.bquark_E2=bquark_gen.at(1).E();
    outtreeVars.bquark_pt2=bquark_gen.at(1).Pt();

    //try to match b quarks with the selected recojets
    if( DeltaR( bquark_gen.at(0).Eta() , bquark_gen.at(0).Phi() , bjet_lead.Eta()    , bjet_lead.Phi()    ) < 0.4 && 
	DeltaR( bquark_gen.at(1).Eta() , bquark_gen.at(1).Phi() , bjet_sublead.Eta() , bjet_sublead.Phi() ) < 0.4 )
      bquarkmatch++;
    else
      if( DeltaR( bquark_gen.at(1).Eta() , bquark_gen.at(1).Phi() , bjet_lead.Eta()    , bjet_lead.Phi()    ) < 0.4 && 
	  DeltaR( bquark_gen.at(0).Eta() , bquark_gen.at(0).Phi() , bjet_sublead.Eta() , bjet_sublead.Phi() ) < 0.4 )
	bquarkmatch++;
    /////////////////////////////////////////////////////////////////////////
    */

    
    

    
    //JCR(jet control region): less than one jet with medium b-tag  
    //MPC: exactly one jet with medium b-tag
    //HPC: at least two jets with medium b-tag
    
    if( outtreeVars.dibjet_leadmvav2>=5 && outtreeVars.dibjet_subleadmvav2 >=5 )
    {
      //cout<<"HPC"<<endl;
      outtreeVars.cut_based_ct = 0;
      int NbjetPrompt=0;
    }
    else if(( outtreeVars.dibjet_leadmvav2>=5 && outtreeVars.dibjet_subleadmvav2<5  ) ||
            ( outtreeVars.dibjet_leadmvav2<5  && outtreeVars.dibjet_subleadmvav2>=5 )    )
    {
      //cout<<"MPC"<<endl;
      outtreeVars.cut_based_ct = 1;
    }
    else if( outtreeVars.dibjet_leadmvav2<5 &&  outtreeVars.dibjet_subleadmvav2<5 )
    {
      //cout<<"JCR"<<endl;
      outtreeVars.cut_based_ct = -1;
    }
    
    /////////////////////////////////////////////////////////////////////////
    //temporary for fake photons study
    if(outtreeVars.cut_based_ct>=0)
      if(!PhoGenMatch(pho_lead , pho_sublead , treeVars , outtreeVars, 0.07) )
      {
	Nbadph++;
	TLorentzVector aux;
	if(RecoJetGenericMatch (pho_lead , treeVars , aux, 0.07) || RecoJetGenericMatch (pho_sublead , treeVars , aux, 0.07) )
	  Nfakeph++;
      }
    /////////////////////////////////////////////////////////////////////////
    
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
      if(outtreeVars.cut_based_ct ==  0)
      {
	++Nev_highMx_HPC;
	/*
      if (abs(outtreeVars.jet_hadflav[bjet_lead_i]) == 5 && abs(outtreeVars.jet_hadflav[bjet_sublead_i]) == 5) 	Nev_bjetpromptprompt++;
      if (abs(outtreeVars.jet_hadflav[bjet_lead_i]) == 4 && abs(outtreeVars.jet_hadflav[bjet_sublead_i]) == 4) 	Nev_cjetpromptprompt++;
      if (abs(outtreeVars.jet_hadflav[bjet_lead_i]) == 5 && abs(outtreeVars.jet_hadflav[bjet_sublead_i]) == 4 ||
	  abs(outtreeVars.jet_hadflav[bjet_lead_i]) == 4 && abs(outtreeVars.jet_hadflav[bjet_sublead_i]) == 5) 	Nev_bcjet++;
      if (abs(outtreeVars.jet_hadflav[bjet_lead_i]) == 5 && 
	  abs(outtreeVars.jet_hadflav[bjet_sublead_i]) != 5 && abs(outtreeVars.jet_hadflav[bjet_sublead_i]) != 4 || 
	  abs(outtreeVars.jet_hadflav[bjet_sublead_i]) == 5 && 
	  abs(outtreeVars.jet_hadflav[bjet_lead_i]) != 5 && abs(outtreeVars.jet_hadflav[bjet_lead_i]) != 4) 	Nev_bjetpromptfake++;
      if (abs(outtreeVars.jet_hadflav[bjet_lead_i]) == 4 && 
	  abs(outtreeVars.jet_hadflav[bjet_sublead_i]) != 5 && abs(outtreeVars.jet_hadflav[bjet_sublead_i]) != 4 || 
	  abs(outtreeVars.jet_hadflav[bjet_sublead_i]) == 4 && 
	  abs(outtreeVars.jet_hadflav[bjet_lead_i]) != 5 && abs(outtreeVars.jet_hadflav[bjet_lead_i]) != 4) 	Nev_cjetpromptfake++;
      if (abs(outtreeVars.jet_hadflav[bjet_lead_i]) != 5 && abs(outtreeVars.jet_hadflav[bjet_lead_i]) != 4 &&
	  abs(outtreeVars.jet_hadflav[bjet_sublead_i]) != 5 && abs(outtreeVars.jet_hadflav[bjet_sublead_i]) != 4) Nev_jetfakefake++;
	*/
      }
      if(outtreeVars.cut_based_ct ==  1) 
      {
	++Nev_highMx_MPC;
      }
      outTree_all_highMx->Fill();
    }
  }

  std::cout << std::endl;

  cout<<"Nbadph="<<Nbadph<<endl;
  cout<<"Nfakeph="<<Nfakeph<<endl;
  cout<<"Njetgoodmatched/Nev="<<bquarkmatch<<"/"<<Nev_withbquarkHdaugh<<"="<<1.*bquarkmatch/Nev_withbquarkHdaugh<<endl;
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
  cout<<"\t\t\t\tbb = "<<Nev_bjetpromptprompt<<" / "<<NEventsMC<<" = "<<1.*Nev_bjetpromptprompt/NEventsMC<<endl;
  cout<<"\t\t\t\tcc = "<<Nev_cjetpromptprompt<<" / "<<NEventsMC<<" = "<<1.*Nev_cjetpromptprompt/NEventsMC<<endl;
  cout<<"\t\t\t\tbc = "<<Nev_bcjet<<" / "<<NEventsMC<<" = "<<1.*Nev_bcjet/NEventsMC<<endl;
  cout<<"\t\t\t\tbj = "<<Nev_bjetpromptfake<<" / "<<NEventsMC<<" = "<<1.*Nev_bjetpromptfake/NEventsMC<<endl;
  cout<<"\t\t\t\tcj = "<<Nev_cjetpromptfake<<" / "<<NEventsMC<<" = "<<1.*Nev_cjetpromptfake/NEventsMC<<endl;
  cout<<"\t\t\t\tjj = "<<Nev_jetfakefake<<" / "<<NEventsMC<<" = "<<1.*Nev_jetfakefake/NEventsMC<<endl;
  cout<<"\t\t\t-----------------------------------------------------------------"<<endl;
  cout<<"\t\t\tlow mass jet control region acceptance = "<<Nev_lowMx_JCR<<" / "<<NEventsMC<<" = "<<1.*Nev_lowMx_JCR/NEventsMC<<endl;
  cout<<"\t\t\tlow mass medium purity cat. acceptance = "<<Nev_lowMx_MPC<<" / "<<NEventsMC<<" = "<<1.*Nev_lowMx_MPC/NEventsMC<<endl;
  cout<<"\t\t\tlow mass high purity cat. acceptance = "<<Nev_lowMx_HPC<<" / "<<NEventsMC<<" = "<<1.*Nev_lowMx_HPC/NEventsMC<<endl;
  /*  cout<<"\t\t\t\tdiphoton promptprompt = "<<Nev_Phpromptprompt<<" / "<<NEventsMC<<" = "<<1.*Nev_Phpromptprompt/NEventsMC<<endl;
  cout<<"\t\t\t\tdiphoton promptfake = "<<Nev_Phpromptfake<<" / "<<NEventsMC<<" = "<<1.*Nev_Phpromptfake/NEventsMC<<endl;
  cout<<"\t\t\t\tdiphoton fakefake = "<<Nev_Phfakefake<<" / "<<NEventsMC<<" = "<<1.*Nev_Phfakefake/NEventsMC<<endl;
  cout<<"\t\t\t\tdibquark promptprompt = "<<Nev_bquarkpromptprompt<<" / "<<NEventsMC<<" = "<<1.*Nev_bquarkpromptprompt/NEventsMC<<endl;
  cout<<"\t\t\t\tdibquark promptfake = "<<Nev_bquarkpromptfake<<" / "<<NEventsMC<<" = "<<1.*Nev_bquarkpromptfake/NEventsMC<<endl;
  cout<<"\t\t\t\tdibquark fakefake = "<<Nev_bquarkfakefake<<" / "<<NEventsMC<<" = "<<1.*Nev_bquarkfakefake/NEventsMC<<endl;
  */


  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"N_preselected * XS /N_MC = "<< Nev_preselected*default_CrossSection/NEventsMC <<" fb"<<endl; 
  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"-----------------------------------------------------------------\n"<<endl;
  cout<<"WANRNING!!!: remember to re-adjust the photon pt and jet pt,eta selection!!!!"<<endl;
  outTree_all_lowMx -> AutoSave();
  outTree_all_highMx -> AutoSave();
  outFile -> Close();
  
  // MakePlot3(h);
  // system(Form("mv *.png %s",outputPlotFolder.c_str()));
  // system(Form("mv *.pdf %s",outputPlotFolder.c_str()));
  
  for(auto h2 : eff_map_flav4)
    delete h2.second;
  for(auto h2 : eff_map_flav5)
    delete h2.second;
  
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
