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
  TTree* outTree = new TTree("plotTree","plotTree");
  TreeVars outtreeVars;
  InitOutTreeVars(outTree,outtreeVars);

  
  //----------
  // get number of input events to include in the weight
  std::cout << ">>> Merging Event_weight histograms and temporary storing in /tmp/$USER/"<<endl;
  float weightMC=1.;
  TH1F* h_weight;
  TString only_filename_tstr(filename.at(0));
  only_filename_tstr.Remove(0,only_filename_tstr.Last('/')+1);
  only_filename_tstr.ReplaceAll("*","all");
  TString hadd_command="hadd -T -f -k /tmp/$USER/"+only_filename_tstr;
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
      weightMC = 1./h_weight->GetEntries();
    else
    {
      cout<<"[WARNING]: histogram /tmp/$USER/"<<only_filename_tstr.Data()<<"/weightCounter/Event_weight not found or empty --> Set by default weight 1"<<endl;
      h_weight = new TH1F("Event_weight","Event_weight",1,0,1);
      h_weight ->Fill(0.);
      weightMC = 1.;
    }

    TDirectory *weightfolder = outFile->mkdir("weightCounter");
    weightfolder->cd();
    h_weight->Write("Event_weight");
    merged_input.Close();
  }
  outFile -> cd();
  

  

  //------------------
  // loop over samples

  int Nev_selected=0;
  int Nev_phselection=0;
  int Nev_jet_kin_preselection=0;
  int Nev_jetselection=0;
  int nEntries = tree->GetEntries();
  std::cout << "Total entries = " << nEntries << std::endl;
  for(int i=0; i<nEntries; i++)
  {

    tree -> GetEntry(i);
    if( i%1000==0 ) std::cout << "Processing entry "<< i << "\r" << std::flush;
    cout<<"Nphselected="<<treeVars.N_SelectedPh<<endl;
    cout<<"Nphloose="<<treeVars.N_LoosePh<<endl;
    cout<<"Nphtight="<<treeVars.N_TightPh<<endl;

    //    std::cout<<"\n\nEVENT "<<i<<"\n #GEN="<<treeVars.N_GenPart;
    //std::cout<<""<<treeVars.N_GenPart;
    //for(int i=0;i<treeVars.N_GenPart;i++)
    //  std::cout<<"\t"<<treeVars.GenPart_pid[i]<<","<<treeVars.GenPart_st[i]<<","<<treeVars.GenPart_pt[i]<<","<<treeVars.GenPart_mass[i];
    //std::cout<<"\n #GENPH="<<treeVars.N_GenPh;
    //for(int i=0;i<treeVars.N_GenPh;i++)
    //  std::cout<<"\t"<<treeVars.GenPh_st[i]<<","<<treeVars.GenPh_pt[i];
    
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
    
    //Gen-matching
    //if(!PhoGenMatch(pho_lead,pho_sublead,treeVars,outtreeVars,0.1))//default DeltaRmax=0.03
    //   continue;
    
    //Cuts on photons
    if(!DiPhotonSelection(pho_lead,pho_sublead))
    {
      cout<<"photonselection_NOTpassed"<<endl;
      continue;
    }
    ++Nev_phselection;
    cout<<"photonselection_passed"<<endl;

    //Clean jet collection (to do? maybe required for bkg study...)
    //Tag good jets with a bool, requirements are:
    // 1. Not matching with any reco photon or any reco electron
    // 2. ratio Echarge / Eneutral  
    //TagGoodJets(treeVars);

    
    if( ! FindGenJet_Hdaug(treeVars) ) continue;
    //cout<<"GenJet_Hdaug found"<<endl;

    //Jets selections
    int BTagOffset;
    //BTag level:
    // 0 no BTag
    // 1 BTag loose
    // 2 BTag medium
    // 3 BTag tight
    // 4 BTag loose with MTD
    // 5 BTag medium with MTD
    // 6 BTag tight with MTD
    if(useMTD == false)
      BTagOffset=0;
    else
      BTagOffset=3;

    //PrintRecoJet(treeVars);
    outtreeVars.nJets=0;
    outtreeVars.nJets_bTagLoose=0;
    outtreeVars.nJets_bTagMedium=0;	 
    outtreeVars.nJets_bTagTight=0;
    //select in output only jets with a certain minimum pt and b-tagged
    for(int i=0;i<treeVars.N_Jet;i++)
    {
      if(treeVars.Jet_pt[i]<25) continue;
      if(fabs(treeVars.Jet_eta[i])>3) continue;
      if( DeltaR(treeVars.Jet_eta[i],treeVars.Jet_phi[i],pho_lead.Eta(),pho_lead.Phi()) < 0.4 ) continue;
      if( DeltaR(treeVars.Jet_eta[i],treeVars.Jet_phi[i],pho_sublead.Eta(),pho_sublead.Phi()) < 0.4 ) continue;
      int BTag = GetBTagLevel(treeVars.Jet_mvav2[i],useMTD);
      //cout<<i<<"\tbtagvalue="<<treeVars.Jet_mvav2[i]<<"\tbtaglevel="<<BTag<<endl;
      if(BTag>BTagOffset && BTag<4+BTagOffset)
      {
	outtreeVars.nJets++;
	outtreeVars.jet_pt[outtreeVars.nJets-1] = treeVars.Jet_pt[i];
	outtreeVars.jet_eta[outtreeVars.nJets-1] = treeVars.Jet_eta[i];                    
	outtreeVars.jet_phi[outtreeVars.nJets-1] = treeVars.Jet_phi[i];
	outtreeVars.jet_mass[outtreeVars.nJets-1] = treeVars.Jet_mass[i];
	outtreeVars.jet_BTagLevel[outtreeVars.nJets-1] = BTag;
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
    SelectBestScoreBJets(outtreeVars,bjet_lead_i,bjet_sublead_i,useMTD);
    TLorentzVector bjet_lead,bjet_sublead;
    bjet_lead.SetPtEtaPhiM(outtreeVars.jet_pt[bjet_lead_i],outtreeVars.jet_eta[bjet_lead_i],outtreeVars.jet_phi[bjet_lead_i],outtreeVars.jet_mass[bjet_lead_i]);
    bjet_sublead.SetPtEtaPhiM(outtreeVars.jet_pt[bjet_sublead_i],outtreeVars.jet_eta[bjet_sublead_i],outtreeVars.jet_phi[bjet_sublead_i],outtreeVars.jet_mass[bjet_sublead_i]);

    double dibjet_mass = (bjet_lead+bjet_sublead).M();
    //cout<<"Mjj="<<dibjet_mass<<endl;
    if(dibjet_mass<70 || dibjet_mass>180)
      continue;
    ++Nev_jetselection;
    //cout<<"pass jets invariant mass selection"<<endl;
    //cout<<"\n\n\n\n\n\n"<<endl;

    ++Nev_selected;

    //Fill outtree and histos
    outtreeVars.weight = CrossSection*weightMC;
    outtreeVars.cross_sec = CrossSection;

    outtreeVars.nvtx = treeVars.N_Vtx;
    outtreeVars.dipho_mass = (pho_lead+pho_sublead).M();
    outtreeVars.dipho_sumpt = (pho_lead+pho_sublead).Pt();
    outtreeVars.dipho_deltaeta = DeltaEta( pho_lead.Eta() , pho_sublead.Eta() );
    outtreeVars.dipho_deltaphi = DeltaPhi( pho_lead.Phi() , pho_sublead.Phi() );

    outtreeVars.dipho_leadPt = pho_lead.Pt();
    outtreeVars.dipho_leadEta = pho_lead.Eta();
    outtreeVars.dipho_leadPhi = pho_lead.Phi();
    outtreeVars.dipho_leadptoM = pho_lead.Pt() / (pho_lead+pho_sublead).M();
    outtreeVars.dipho_leadEnergy = pho_lead.E();
    //outtreeVars.dipho_leadIso = treeVars.TightPh_iso[pho_lead_i];
    //outtreeVars.dipho_leadDeltaRgenreco = DeltaR(pho_lead.Eta(),pho_lead.Phi(),outtreeVars.dipho_leadEta_gen,outtreeVars.dipho_leadPhi_gen);
    //outtreeVars.dipho_leadDeltaEtagenreco = DeltaEta(pho_lead.Eta(),outtreeVars.dipho_leadEta_gen);
    //outtreeVars.dipho_leadDeltaPhigenreco = DeltaPhi(pho_lead.Phi(),outtreeVars.dipho_leadPhi_gen);

    outtreeVars.dipho_subleadPt = pho_sublead.Pt();
    outtreeVars.dipho_subleadEta = pho_sublead.Eta();
    outtreeVars.dipho_subleadPhi = pho_sublead.Phi();
    outtreeVars.dipho_subleadptoM = pho_sublead.Pt() / (pho_lead+pho_sublead).M();
    outtreeVars.dipho_subleadEnergy = pho_sublead.E();
    //outtreeVars.dipho_subleadIso = treeVars.TightPh_iso[pho_sublead_i];
    //outtreeVars.dipho_subleadDeltaRgenreco = DeltaR(pho_sublead.Eta(),pho_sublead.Phi(),outtreeVars.dipho_subleadEta_gen,outtreeVars.dipho_subleadPhi_gen);
    //outtreeVars.dipho_subleadDeltaEtagenreco = DeltaEta(pho_sublead.Eta(),outtreeVars.dipho_subleadEta_gen);
    //outtreeVars.dipho_subleadDeltaPhigenreco = DeltaPhi(pho_sublead.Phi(),outtreeVars.dipho_subleadPhi_gen);

    outtreeVars.dibjet_mass = (bjet_lead+bjet_sublead).M();
    outtreeVars.dibjet_sumpt = (bjet_lead+bjet_sublead).Pt();
    outtreeVars.dibjet_deltaeta = DeltaEta( bjet_lead.Eta() , bjet_sublead.Eta() );
    outtreeVars.dibjet_deltaphi = DeltaPhi( bjet_lead.Phi() , bjet_sublead.Phi() );

    outtreeVars.dibjet_leadPt = bjet_lead.Pt();
    outtreeVars.dibjet_leadEta = bjet_lead.Eta();
    outtreeVars.dibjet_leadPhi = bjet_lead.Phi();
    outtreeVars.dibjet_leadptoM = bjet_lead.Pt() / (bjet_lead+bjet_sublead).M();
    outtreeVars.dibjet_leadEnergy = bjet_lead.E();
    outtreeVars.dibjet_leadbtagscore = outtreeVars.jet_BTagLevel[bjet_lead_i];

    outtreeVars.dibjet_subleadPt = bjet_sublead.Pt();
    outtreeVars.dibjet_subleadEta = bjet_sublead.Eta();
    outtreeVars.dibjet_subleadPhi = bjet_sublead.Phi();
    outtreeVars.dibjet_subleadptoM = bjet_sublead.Pt() / (bjet_lead+bjet_sublead).M();
    outtreeVars.dibjet_subleadEnergy = bjet_sublead.E();
    outtreeVars.dibjet_subleadbtagscore = outtreeVars.jet_BTagLevel[bjet_sublead_i];


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
    
    //h["dipho_leaddeltaR_GenReco"] -> Fill( DeltaR(pho_lead.Eta(),pho_lead.Phi(),outtreeVars.dipho_leadEta_gen,outtreeVars.dipho_leadPhi_gen) );
    //h["dipho_subleaddeltaR_GenReco"] -> Fill( DeltaR(pho_sublead.Eta(),pho_sublead.Phi(),outtreeVars.dipho_subleadEta_gen,outtreeVars.dipho_subleadPhi_gen) );
    //h["dipho_leaddeltaEta_GenReco"] -> Fill( DeltaEta(pho_lead.Eta(),outtreeVars.dipho_leadEta_gen) );
    //h["dipho_subleaddeltaEta_GenReco"] -> Fill( DeltaEta(pho_sublead.Eta(),outtreeVars.dipho_subleadEta_gen) );
    //h["dipho_leaddeltaPhi_GenReco"] -> Fill( DeltaPhi(pho_lead.Phi(),outtreeVars.dipho_leadPhi_gen) );
    //h["dipho_subleaddeltaPhi_GenReco"] -> Fill( DeltaPhi(pho_sublead.Phi(),outtreeVars.dipho_subleadPhi_gen) );

    outTree->Fill();

  }
  std::cout << std::endl;

  cout<<"\n\n\n\n\nselected events = "<<Nev_selected<<" / "<<nEntries<<endl;
  cout<<"\t\tpass photon selection = "<<Nev_phselection<<" / "<<nEntries<<endl;
  cout<<"\t\tpass kinematic bjet preselection = "<< Nev_jet_kin_preselection <<" / "<<nEntries<<endl;
  cout<<"\t\tpass jet selection = "<<Nev_jetselection<<" / "<<nEntries<<endl;
  cout<<"\n\n\n\n\n"<<endl;

  outTree -> AutoSave();
  outFile -> Close();

  MakePlot3(h);

  system(Form("mv *.png %s",outputPlotFolder.c_str()));
  system(Form("mv *.pdf %s",outputPlotFolder.c_str()));  



  return 0;
}


