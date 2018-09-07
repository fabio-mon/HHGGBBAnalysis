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
  lumi_sqrtS = "";//Form("%.1f fb^{-1} (13 TeV)",lumi);
  
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

  std::string Loose_Tight_Photon = "PhotonTight";
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

  int Nph_veto = -1;//discard events with N selected photon gen-matched > Nph_veto --> Useful to avoid double counting in GJet and QCD samples   
  if(opts.OptExist("Input.Nph_veto"))
    Nph_veto = opts.GetOpt<int> ("Input.Nph_veto");
  
  //------------------
  // define histograms
  std::map<std::string,TH1F*> h;
  h["dipho_mass"] = new TH1F("dipho_mass","dipho_mass",100,90.,160.);
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
  h["dibjet_mass"] = new TH1F("dibjet_mass","dibjet_mass",100,60,200);
  h["dibjet_sumpt"] = new TH1F("dibjet_sumpt","dibjet_sumpt",100,0,300);
  h["dibjet_deltaeta"] = new TH1F("dibjet_deltaeta","dibjet_deltaeta",100,0,8);
  h["dibjet_deltaphi"] = new TH1F("dibjet_deltaphi","dibjet_deltaphi",100,0,3.14);
  h["dibjet_leadPt"] = new TH1F("dibjet_leadPt","dibjet_leadPt",100,0,200);
  h["dibjet_leadEta"] = new TH1F("dibjet_leadEta","dibjet_leadEta",100,-3.5,3.5);
  h["dibjet_leadPhi"] = new TH1F("dibjet_leadPhi","dibjet_leadPhi",100,-3.14,3.14);
  h["dibjet_leadptoM"] = new TH1F("dibjet_leadptoM","dibjet_leadptoM",100,0,3.5);
  h["dibjet_leadEnergy"] = new TH1F("dibjet_leadEnergy","dibjet_leadEnergy",100,0,300);
  h["dibjet_leadbtagscore"] = new TH1F("dibjet_leadbtagscore","dibjet_leadbtagscore",9,-1.5,7.5);
  h["dibjet_subleadPt"] = new TH1F("dibjet_subleadPt","dibjet_subleadPt",100,0,200);
  h["dibjet_subleadEta"] = new TH1F("dibjet_subleadEta","dibjet_subleadEta",100,-3.5,3.5);
  h["dibjet_subleadPhi"] = new TH1F("dibjet_subleadPhi","dibjet_subleadPhi",100,-3.14,3.14);
  h["dibjet_subleadptoM"] = new TH1F("dibjet_subleadptoM","dibjet_subleadptoM",100,0,3.5);
  h["dibjet_subleadEnergy"] = new TH1F("dibjet_subleadEnergy","dibjet_subleadEnergy",100,0,300);
  h["dibjet_subleadbtagscore"] = new TH1F("dibjet_subleadbtagscore","dibjet_subleadbtagscore",9,-1.5,7.5);

  h["Mx"] = new TH1F("M_x","M_x",100,0,1000);
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
  InitRawTreeVars(trees,treeVars,Loose_Tight_Photon);
  TChain *tree = trees[onlytreename_str.at(0)];
  

  //----------
  // add tree_i as friends to tree_0
  for(unsigned int ntree = 1; ntree < onlytreename_str.size(); ++ntree)
  {
    std::cout << ">>> Adding chain " << onlytreename_str.at(ntree) << " to chain " << onlytreename_str.at(0) << std::endl;
    tree->AddFriend(onlytreename_str.at(ntree).c_str(),"");
  }
  
  //------------------      
  //output trees
  std::string outputPlotFolder = opts.GetOpt<std::string>("Output.outputFolder");
  system(Form("mkdir -p %s",outputPlotFolder.c_str()));
  TFile* outFile = TFile::Open(Form("%s/plotTree_%s.root",outputPlotFolder.c_str(),label.c_str()),"RECREATE");
  outFile -> cd();  
  TTree* outTree_preselected = new TTree("preselected","preselected");
  TTree* outTree_preselected_lowMx = new TTree("preselected_lowMx","preselected_lowMx");
  TTree* outTree_preselected_highMx = new TTree("preselected_highMx","preselected_highMx");
  TTree* outTree_all_lowMx  = new TTree("all_lowMx", "all_lowMx");
  TTree* outTree_all_highMx = new TTree("all_highMx","all_highMx");
  TreeVars outtreeVars;
  InitOutTreeVars(outTree_preselected,outtreeVars);
  InitOutTreeVars(outTree_preselected_lowMx,outtreeVars);
  InitOutTreeVars(outTree_preselected_highMx,outtreeVars);
  InitOutTreeVars(outTree_all_highMx,outtreeVars);
  InitOutTreeVars(outTree_all_lowMx, outtreeVars);
  
  
  //----------
  // get number of input events to include in the weight
  std::cout << ">>> Merging Event_weight histograms and temporary storing in /tmp/$USER/"<<endl;
  float weightMC=1.;
  TH1F* h_weight;
  TString only_filename_tstr(filename.at(0));
  only_filename_tstr.Remove(0,only_filename_tstr.Last('/')+1);
  only_filename_tstr.ReplaceAll("*","all");
  TString hadd_command="hadd -T -f -k /tmp/$USER/"+only_filename_tstr;
  int Nev_MC = 1;
  for(unsigned i=0; i<filename.size(); i++)
  {
    //TString filename_tstr(filename);
    hadd_command+=' ';
    hadd_command+=filename.at(i).c_str();
  }
    
  system(hadd_command.Data()); // -T=merge only histos -f=overwrite target file -k=do not exit in case of corrupted files
  TFile merged_input(("/tmp/$USER/"+only_filename_tstr).Data(),"READ");
  if(merged_input.IsOpen())
  {
    merged_input.GetObject("weightCounter/Event_weight",h_weight);
    outFile -> cd();
    if(h_weight && h_weight->GetEntries()>0)
    {
      weightMC = 1./h_weight->GetEntries();
      Nev_MC = h_weight->GetEntries();
    }
    else
    {
      cout<<"[WARNING]: histogram /tmp/$USER/"<<only_filename_tstr.Data()<<"/weightCounter/Event_weight not found or empty --> Set by default weight 1"<<endl;
      h_weight = new TH1F("Event_weight","Event_weight",1,0,1);
      h_weight ->Fill(0.);
      weightMC = 1.;
      Nev_MC = 1;
    }

    TDirectory *weightfolder = outFile->mkdir("weightCounter");
    weightfolder->cd();
    h_weight->Write("Event_weight");
    merged_input.Close();
  }
  outFile -> cd();
  

  

  //------------------
  // loop over samples

  int Nev_preselected=0;
  int Nev_preselected_lowMx=0;
  int Nev_preselected_highMx=0;
  int Nev_highMx_JCR=0;
  int Nev_highMx_MPC=0;
  int Nev_highMx_HPC=0;
  int Nev_all_lowMx=0;
  int Nev_all_highMx=0;
  
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

    //find higgs daugher gen-photons and flag their .isHdaug
    //if( ! FindGenPh_Hdaug(treeVars) ) continue;

    //find lead & sublead reco photons    
    //PrintRecoPhoton(treeVars);
    int pho_lead_i;
    int pho_sublead_i;
    FindLeadSublead_pho(treeVars,pho_lead_i,pho_sublead_i);
    //cout<<"lead="<<pho_lead_i<<"\t\tsublead="<<pho_sublead_i<<endl;
    TLorentzVector pho_lead,pho_sublead;
    pho_lead.SetPtEtaPhiE(treeVars.SelectedPh_pt[pho_lead_i],treeVars.SelectedPh_eta[pho_lead_i],treeVars.SelectedPh_phi[pho_lead_i],treeVars.SelectedPh_E[pho_lead_i]);
    pho_sublead.SetPtEtaPhiE(treeVars.SelectedPh_pt[pho_sublead_i],treeVars.SelectedPh_eta[pho_sublead_i],treeVars.SelectedPh_phi[pho_sublead_i],treeVars.SelectedPh_E[pho_sublead_i]);
    TLorentzVector dipho = pho_lead+pho_sublead;

    //Gen-matching
    //if(!PhoGenMatch(pho_lead,pho_sublead,treeVars,outtreeVars,0.1))//default DeltaRmax=0.03
    //   continue;
    
    //Cuts on photons
    if(!DiPhotonSelection(pho_lead,pho_sublead))
    {
      //cout<<"photonselection_NOTpassed"<<endl;
      continue;
    }
    ++Nev_phselection;
    //cout<<"photonselection_passed"<<endl;

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
    int BTagOffset;//--->See AnalysisUtils::GetBTagLevel
    if(useMTD == false)
      BTagOffset=0;
    else
      BTagOffset=3;

    //PrintRecoJet(treeVars);
    outtreeVars.nJets=0;
    outtreeVars.nJets_bTagLoose=0;
    outtreeVars.nJets_bTagMedium=0;	 
    outtreeVars.nJets_bTagTight=0;
    for(int i=0;i<treeVars.N_Jet;i++)
    {
      if(treeVars.Jet_pt[i]<25) continue;
      if(fabs(treeVars.Jet_eta[i])>3) continue;
      if( DeltaR(treeVars.Jet_eta[i],treeVars.Jet_phi[i],pho_lead.Eta(),pho_lead.Phi()) < 0.4 ) continue;
      if( DeltaR(treeVars.Jet_eta[i],treeVars.Jet_phi[i],pho_sublead.Eta(),pho_sublead.Phi()) < 0.4 ) continue;
      int BTag = GetBTagLevel(treeVars.Jet_mvav2[i],useMTD);
      //cout<<i<<"\tbtagvalue="<<treeVars.Jet_mvav2[i]<<"\tbtaglevel="<<BTag<<"\tbtagoffset="<<BTagOffset<<endl;
      //if(BTag>BTagOffset && BTag<4+BTagOffset)
      {
	outtreeVars.nJets++;
	outtreeVars.jet_pt[outtreeVars.nJets-1] = treeVars.Jet_pt[i];
	outtreeVars.jet_eta[outtreeVars.nJets-1] = treeVars.Jet_eta[i];                    
	outtreeVars.jet_phi[outtreeVars.nJets-1] = treeVars.Jet_phi[i];
	outtreeVars.jet_mass[outtreeVars.nJets-1] = treeVars.Jet_mass[i];
	outtreeVars.jet_BTagLevel[outtreeVars.nJets-1] = BTag;
	outtreeVars.jet_mvav2[outtreeVars.nJets-1] = treeVars.Jet_mvav2[i];
	outtreeVars.jet_hadflav[outtreeVars.nJets-1] = treeVars.Jet_hadflav[i];//gen level info! handle with care!

	if(BTag-BTagOffset==1)
	  outtreeVars.nJets_bTagLoose++;
	else
	  if(BTag-BTagOffset==2)
	    outtreeVars.nJets_bTagMedium++;
	  else
	    if(BTag-BTagOffset==3)
	      outtreeVars.nJets_bTagTight++;
      }
    }
    


    if(outtreeVars.nJets<2) 
    {
      //cout<<"NOT pass jet selection"<<endl;
      continue;
    }
    Nev_jet_kin_preselection++;
    //cout<<"pass jet selection"<<endl;

    //Select the two jets with the higher BTag level, if they have the same value select the harder one
    int bjet_lead_i;
    int bjet_sublead_i;
    SelectBestScoreBJets2(outtreeVars,bjet_lead_i,bjet_sublead_i,useMTD);
    TLorentzVector bjet_lead,bjet_sublead;
    bjet_lead.SetPtEtaPhiM(outtreeVars.jet_pt[bjet_lead_i],outtreeVars.jet_eta[bjet_lead_i],outtreeVars.jet_phi[bjet_lead_i],outtreeVars.jet_mass[bjet_lead_i]);
    bjet_sublead.SetPtEtaPhiM(outtreeVars.jet_pt[bjet_sublead_i],outtreeVars.jet_eta[bjet_sublead_i],outtreeVars.jet_phi[bjet_sublead_i],outtreeVars.jet_mass[bjet_sublead_i]);

    TLorentzVector dibjet= bjet_lead+bjet_sublead;
    double dibjet_mass = dibjet.M();
    //cout<<"Mjj="<<dibjet_mass<<endl;
    if(dibjet_mass<70 || dibjet_mass>200)
      continue;
    ++Nev_jetselection;
    //cout<<"pass jets invariant mass selection"<<endl;
    //cout<<"\n\n\n\n\n\n"<<endl;

    //Find dR_min between selected gamma and jets
    vector<TLorentzVector> pho_selected;
    pho_selected.push_back(pho_lead);
    pho_selected.push_back(pho_sublead);
    vector<TLorentzVector> bjet_selected;
    bjet_selected.push_back(bjet_lead);
    bjet_selected.push_back(bjet_sublead);
    float DeltaRmin_bjet_pho = DeltaRmin(pho_selected,bjet_selected);

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
    outtreeVars.weight = CrossSection*weightMC;
    outtreeVars.cross_sec = CrossSection;
    outtreeVars.nvtx = treeVars.N_Vtx;
    outtreeVars.dipho_mass = (dipho).M();
    outtreeVars.dipho_sumpt = (dipho).Pt();
    outtreeVars.dipho_deltaeta = DeltaEta( pho_lead.Eta() , pho_sublead.Eta() );
    outtreeVars.dipho_deltaphi = DeltaPhi( pho_lead.Phi() , pho_sublead.Phi() );
    outtreeVars.dipho_leadPt = pho_lead.Pt();
    outtreeVars.dipho_leadEta = pho_lead.Eta();
    outtreeVars.dipho_leadPhi = pho_lead.Phi();
    outtreeVars.dipho_leadptoM = pho_lead.Pt() / (dipho).M();
    outtreeVars.dipho_leadEnergy = pho_lead.E();
    //outtreeVars.dipho_leadDeltaRgenreco = DeltaR(pho_lead.Eta(),pho_lead.Phi(),outtreeVars.dipho_leadEta_gen,outtreeVars.dipho_leadPhi_gen);
    //outtreeVars.dipho_leadDeltaEtagenreco = DeltaEta(pho_lead.Eta(),outtreeVars.dipho_leadEta_gen);
    //outtreeVars.dipho_leadDeltaPhigenreco = DeltaPhi(pho_lead.Phi(),outtreeVars.dipho_leadPhi_gen);
    outtreeVars.dipho_subleadPt = pho_sublead.Pt();
    outtreeVars.dipho_subleadEta = pho_sublead.Eta();
    outtreeVars.dipho_subleadPhi = pho_sublead.Phi();
    outtreeVars.dipho_subleadptoM = pho_sublead.Pt() / (dipho).M();
    outtreeVars.dipho_subleadEnergy = pho_sublead.E();
    //outtreeVars.dipho_subleadDeltaRgenreco = DeltaR(pho_sublead.Eta(),pho_sublead.Phi(),outtreeVars.dipho_subleadEta_gen,outtreeVars.dipho_subleadPhi_gen);
    //outtreeVars.dipho_subleadDeltaEtagenreco = DeltaEta(pho_sublead.Eta(),outtreeVars.dipho_subleadEta_gen);
    //outtreeVars.dipho_subleadDeltaPhigenreco = DeltaPhi(pho_sublead.Phi(),outtreeVars.dipho_subleadPhi_gen);
    outtreeVars.dibjet_mass = (dibjet).M();
    outtreeVars.dibjet_sumpt = (dibjet).Pt();
    outtreeVars.dibjet_deltaeta = DeltaEta( bjet_lead.Eta() , bjet_sublead.Eta() );
    outtreeVars.dibjet_deltaphi = DeltaPhi( bjet_lead.Phi() , bjet_sublead.Phi() );
    outtreeVars.dibjet_leadPt = bjet_lead.Pt();
    outtreeVars.dibjet_leadEta = bjet_lead.Eta();
    outtreeVars.dibjet_leadPhi = bjet_lead.Phi();
    outtreeVars.dibjet_leadptoM = bjet_lead.Pt() / (dibjet).M();
    outtreeVars.dibjet_leadEnergy = bjet_lead.E();
    outtreeVars.dibjet_leadbtagscore = outtreeVars.jet_BTagLevel[bjet_lead_i];
    outtreeVars.dibjet_leadmvav2 = outtreeVars.jet_mvav2[bjet_lead_i];
    outtreeVars.dibjet_subleadPt = bjet_sublead.Pt();
    outtreeVars.dibjet_subleadEta = bjet_sublead.Eta();
    outtreeVars.dibjet_subleadPhi = bjet_sublead.Phi();
    outtreeVars.dibjet_subleadptoM = bjet_sublead.Pt() / (dibjet).M();
    outtreeVars.dibjet_subleadEnergy = bjet_sublead.E();
    outtreeVars.dibjet_subleadbtagscore = outtreeVars.jet_BTagLevel[bjet_sublead_i];
    outtreeVars.dibjet_subleadmvav2 = outtreeVars.jet_mvav2[bjet_sublead_i];
    outtreeVars.Mx = diHiggs.M() - outtreeVars.dibjet_mass - outtreeVars.dipho_mass + 250.;
    outtreeVars.DRmin_pho_bjet = DeltaRmin_bjet_pho; 
    outtreeVars.costheta_HH = costheta_HH; 
    outtreeVars.costheta_gg = costheta_gg; 
    outtreeVars.costheta_bb = costheta_bb; 
    outtreeVars.MetPt = treeVars.Met_pt[0];
    outtreeVars.MetPhi = treeVars.Met_phi[0];
    
    outtreeVars.ttHTagger = 0;
    
    h["dipho_mass"] -> Fill(outtreeVars.dipho_mass);
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
    h["dibjet_mass"] -> Fill(outtreeVars.dibjet_mass);
    h["dibjet_sumpt"] -> Fill(outtreeVars.dibjet_sumpt);
    h["dibjet_deltaeta"] -> Fill(outtreeVars.dibjet_deltaeta);
    h["dibjet_deltaphi"] -> Fill(outtreeVars.dibjet_deltaphi);

    h["dibjet_leadPt"] -> Fill(outtreeVars.dibjet_leadPt);
    h["dibjet_leadEta"] -> Fill(outtreeVars.dibjet_leadEta);
    h["dibjet_leadPhi"] -> Fill(outtreeVars.dibjet_leadPhi);
    h["dibjet_leadptoM"] -> Fill(outtreeVars.dibjet_leadptoM);
    h["dibjet_leadEnergy"] -> Fill(outtreeVars.dibjet_leadEnergy);
    h["dibjet_leadbtagscore"] -> Fill(outtreeVars.dibjet_leadbtagscore);

    h["dibjet_subleadPt"] -> Fill(outtreeVars.dibjet_subleadPt);
    h["dibjet_subleadEta"] -> Fill(outtreeVars.dibjet_subleadEta);
    h["dibjet_subleadPhi"] -> Fill(outtreeVars.dibjet_subleadPhi);
    h["dibjet_subleadptoM"] -> Fill(outtreeVars.dibjet_subleadptoM);
    h["dibjet_subleadEnergy"] -> Fill(outtreeVars.dibjet_subleadEnergy);
    h["dibjet_subleadbtagscore"] -> Fill(outtreeVars.dibjet_subleadbtagscore);
    h["Mx"] -> Fill(outtreeVars.Mx);
    h["DRmin_sel_phots_sel_jets"] -> Fill(outtreeVars.DRmin_pho_bjet); 
    h["costhetastar_HH"] -> Fill(outtreeVars.costheta_HH); 
    h["costhetastar_gg"] -> Fill(outtreeVars.costheta_gg); 
    h["costhetastar_bb"] -> Fill(outtreeVars.costheta_bb); 
    h["MET"] -> Fill(outtreeVars.MetPt); 
    h["MET_phi"] -> Fill(outtreeVars.MetPhi); 
    
    
    //h["dipho_leaddeltaR_GenReco"] -> Fill( DeltaR(pho_lead.Eta(),pho_lead.Phi(),outtreeVars.dipho_leadEta_gen,outtreeVars.dipho_leadPhi_gen) );
    //h["dipho_subleaddeltaR_GenReco"] -> Fill( DeltaR(pho_sublead.Eta(),pho_sublead.Phi(),outtreeVars.dipho_subleadEta_gen,outtreeVars.dipho_subleadPhi_gen) );
    //h["dipho_leaddeltaEta_GenReco"] -> Fill( DeltaEta(pho_lead.Eta(),outtreeVars.dipho_leadEta_gen) );
    //h["dipho_subleaddeltaEta_GenReco"] -> Fill( DeltaEta(pho_sublead.Eta(),outtreeVars.dipho_subleadEta_gen) );
    //h["dipho_leaddeltaPhi_GenReco"] -> Fill( DeltaPhi(pho_lead.Phi(),outtreeVars.dipho_leadPhi_gen) );
    //h["dipho_subleaddeltaPhi_GenReco"] -> Fill( DeltaPhi(pho_sublead.Phi(),outtreeVars.dipho_subleadPhi_gen) );
    ++Nev_preselected;
    outTree_preselected->Fill();
    
    //fill low mass and high mass categories
    if(outtreeVars.Mx<=350)
    {
      ++Nev_preselected_lowMx;
      outTree_preselected_lowMx->Fill();
    }
    else
    {
      ++Nev_preselected_highMx;
      outTree_preselected_highMx->Fill();
    }
    
    //cout<<"high mass"<<endl;
    //high mass events: categorization based on b-tag level of the two selected jets
    //JCR(jet control region): less than one jet with medium b-tag  
    //MPC: exactly one jet with medium b-tag
    //HPC: at least two jets with medium b-tag

    int BTagMedium_mask;
    if(useMTD == false)
      BTagMedium_mask=0b000010;
    else
      BTagMedium_mask=0b010000;

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
    if(outtreeVars.Mx<=350)
    {
      ++Nev_all_lowMx;
      outTree_all_lowMx->Fill();
    }
    else
    {
      ++Nev_all_highMx;
      outTree_all_highMx->Fill();
    }
    
    // if( (outtreeVars.dibjet_leadmvav2 & BTagMedium_mask) && (outtreeVars.dibjet_subleadmvav2 & BTagMedium_mask) )
    // {
    //   //cout<<"HPC"<<endl;
    //   ++Nev_highMx_HPC;
    //   outTree_highMx_HPC->Fill();
    //   TLorentzVector v;
    //   int NphPrompt=0;
    //   if (PhoGenericGenMatch(pho_lead,treeVars,v,0.1)) NphPrompt++;
    //   if (PhoGenericGenMatch(pho_sublead,treeVars,v,0.1)) NphPrompt++;
    //   if(NphPrompt==2) 
    //     Nev_Phpromptprompt++;
    //   else
    //     if(NphPrompt==1)
    //       Nev_Phpromptfake++;
    //     else
    //       Nev_Phfakefake++;
    //   int NbjetPrompt=0;
    //   if (outtreeVars.jet_hadflav[bjet_lead_i]==5 || outtreeVars.jet_hadflav[bjet_lead_i]==4) NbjetPrompt++; 
    //   if (outtreeVars.jet_hadflav[bjet_sublead_i]==5 || outtreeVars.jet_hadflav[bjet_lead_i]==4) NbjetPrompt++; 
    //   if(NbjetPrompt==2) 
    //     Nev_bjetpromptprompt++;
    //   else
    //     if(NbjetPrompt==1)
    //       Nev_bjetpromptfake++;
    //     else
    //       Nev_bjetfakefake++;

    //   int NbquarkPrompt=0;
    //   if (bquarkGenericGenMatch(bjet_lead,treeVars,v,0.4)) 	NbquarkPrompt++;
    //   if (bquarkGenericGenMatch(bjet_sublead,treeVars,v,0.4))	NbquarkPrompt++;
    //   if(NbquarkPrompt==2) 
    //     Nev_bquarkpromptprompt++;
    //   else
    //     if(NbquarkPrompt==1)
    //       Nev_bquarkpromptfake++;
    //     else
    //       Nev_bquarkfakefake++;
    // }
    // else
    //   if( (outtreeVars.dibjet_leadmvav2 & BTagMedium_mask) && !(outtreeVars.dibjet_subleadmvav2 & BTagMedium_mask) ||
    //      !(outtreeVars.dibjet_leadmvav2 & BTagMedium_mask) &&  (outtreeVars.dibjet_subleadmvav2 & BTagMedium_mask) )//must be after HPC!!!
    //   {
    //     //cout<<"MPC"<<endl;
    //     ++Nev_highMx_MPC;
    //     outTree_highMx_MPC->Fill();
    //   }
    //   else
    //     if ( !(outtreeVars.dibjet_leadmvav2 & BTagMedium_mask) && !(outtreeVars.dibjet_subleadmvav2 & BTagMedium_mask) )  
    //     {
    //       //cout<<"JCR"<<endl;
    //       ++Nev_highMx_JCR;
    //       outTree_highMx_JCR->Fill();
    //     }
	  
    //additional possible selections for high purity category
    //if(outtreeVars.dipho_mass<115 || outtreeVars.dipho_mass>135) continue;
    //if(outtreeVars.dibjet_sumpt<20 || outtreeVars.dibjet_sumpt>550) continue;
    //if(outtreeVars.dipho_sumpt<20 || outtreeVars.dipho_sumpt>550) continue;
    //if(outtreeVars.Mx<300 || outtreeVars.Mx>950) continue;
    //if(outtreeVars.costheta_HH>0.8) continue;
    //if(outtreeVars.dipho_deltaphi>2.6) continue;

  }

  std::cout << std::endl;

  cout<<"\n-----------------------------------------------------------------"<<endl;
  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"MC events = "<<Nev_MC<<endl;
  cout<<"2b2g skim = "<<nEntries<<" / "<<Nev_MC<<endl;
  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"preselected events = "<<Nev_preselected<<" / "<<nEntries<<endl;
  cout<<"\t\tpass photon preselection = "<<Nev_phselection<<" / "<<nEntries<<endl;
  cout<<"\t\tpass kinematic jet preselection = "<< Nev_jet_kin_preselection <<" / "<<nEntries<<endl;
  cout<<"\t\tpass jet preselection = "<<Nev_jetselection<<" / "<<nEntries<<endl;
  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"\t\tpreselection acceptance = "<<Nev_preselected<<" / "<<Nev_MC<<" = "<<Nev_preselected/Nev_MC<<endl;
  cout<<"\t\tpreselection low mass acceptance = "<<Nev_preselected_lowMx<<" / "<<Nev_MC<<" = "<<1.*Nev_preselected_lowMx/Nev_MC<<endl;
  cout<<"\t\tpreselection high mass acceptance = "<<Nev_preselected_highMx<<" / "<<Nev_MC<<" = "<<1.*Nev_preselected_highMx/Nev_MC<<endl;
  cout<<"\t\t\thigh mass jet control region acceptance = "<<Nev_highMx_JCR<<" / "<<Nev_MC<<" = "<<1.*Nev_highMx_JCR/Nev_MC<<endl;
  cout<<"\t\t\thigh mass medium purity cat. acceptance = "<<Nev_highMx_MPC<<" / "<<Nev_MC<<" = "<<1.*Nev_highMx_MPC/Nev_MC<<endl;
  cout<<"\t\t\thigh mass high purity cat. acceptance = "<<Nev_highMx_HPC<<" / "<<Nev_MC<<" = "<<1.*Nev_highMx_HPC/Nev_MC<<endl;
  cout<<"\t\t\t\tdiphoton promptprompt = "<<Nev_Phpromptprompt<<" / "<<Nev_MC<<" = "<<1.*Nev_Phpromptprompt/Nev_MC<<endl;
  cout<<"\t\t\t\tdiphoton promptfake = "<<Nev_Phpromptfake<<" / "<<Nev_MC<<" = "<<1.*Nev_Phpromptfake/Nev_MC<<endl;
  cout<<"\t\t\t\tdiphoton fakefake = "<<Nev_Phfakefake<<" / "<<Nev_MC<<" = "<<1.*Nev_Phfakefake/Nev_MC<<endl;
  cout<<"\t\t\t\tdibjet promptprompt = "<<Nev_bjetpromptprompt<<" / "<<Nev_MC<<" = "<<1.*Nev_bjetpromptprompt/Nev_MC<<endl;
  cout<<"\t\t\t\tdibjet promptfake = "<<Nev_bjetpromptfake<<" / "<<Nev_MC<<" = "<<1.*Nev_bjetpromptfake/Nev_MC<<endl;
  cout<<"\t\t\t\tdibjet fakefake = "<<Nev_bjetfakefake<<" / "<<Nev_MC<<" = "<<1.*Nev_bjetfakefake/Nev_MC<<endl;
  cout<<"\t\t\t\tdibquark promptprompt = "<<Nev_bquarkpromptprompt<<" / "<<Nev_MC<<" = "<<1.*Nev_bquarkpromptprompt/Nev_MC<<endl;
  cout<<"\t\t\t\tdibquark promptfake = "<<Nev_bquarkpromptfake<<" / "<<Nev_MC<<" = "<<1.*Nev_bquarkpromptfake/Nev_MC<<endl;
  cout<<"\t\t\t\tdibquark fakefake = "<<Nev_bquarkfakefake<<" / "<<Nev_MC<<" = "<<1.*Nev_bquarkfakefake/Nev_MC<<endl;



  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"N_preselected * XS /N_MC = "<< Nev_preselected*CrossSection/Nev_MC <<" fb"<<endl; 
  cout<<"-----------------------------------------------------------------"<<endl;
  cout<<"-----------------------------------------------------------------\n"<<endl;

  outTree_preselected -> AutoSave();
  outTree_preselected_lowMx -> AutoSave();
  outTree_preselected_highMx -> AutoSave();
  // outTree_highMx_JCR -> AutoSave();
  // outTree_highMx_MPC -> AutoSave();
  // outTree_highMx_HPC -> AutoSave();
  outTree_all_lowMx -> AutoSave();
  outTree_all_highMx -> AutoSave();
  outFile -> Close();
  
  MakePlot3(h);

  system(Form("mv *.png %s",outputPlotFolder.c_str()));
  system(Form("mv *.pdf %s",outputPlotFolder.c_str()));  



  return 0;
}


