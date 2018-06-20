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
  h["nJets"] = new TH1F("nJets","nJets",11,-0.5,10.5);
  h["dipho_leaddeltaR_GenReco"] = new TH1F("dipho_leaddeltaR_GenReco","dipho_leaddeltaR_GenReco",300,0.,0.1);
  h["dipho_subleaddeltaR_GenReco"] = new TH1F("dipho_subleaddeltaR_GenReco","dipho_subleaddeltaR_GenReco",300,0.,0.1);
  h["dipho_leadIso_o_E"] = new TH1F("dipho_leadIso_o_E","dipho_leadIso_o_E",300,0.,0.25);
  h["dipho_subleadIso_o_E"] = new TH1F("dipho_subleadIso_o_E","dipho_subleadIso_o_E",300,0.,0.25);
  h["dipho_leaddeltaEta_GenReco"] = new TH1F("dipho_leaddeltaEta_GenReco","dipho_leaddeltaEta_GenReco",300,0.,0.2);
  h["dipho_subleaddeltaEta_GenReco"] = new TH1F("dipho_subleaddeltaEta_GenReco","dipho_subleaddeltaEta_GenReco",300,0.,0.2);
  h["dipho_leaddeltaPhi_GenReco"] = new TH1F("dipho_leaddeltaPhi_GenReco","dipho_leaddeltaPhi_GenReco",300,0.,0.2);
  h["dipho_subleaddeltaPhi_GenReco"] = new TH1F("dipho_subleaddeltaPhi_GenReco","dipho_subleaddeltaPhi_GenReco",300,0.,0.2);

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

  //------------------
  // loop over samples
  int nEntries = tree->GetEntries();
  std::cout << "Total entries = " << nEntries << std::endl;
  for(int i=0; i<nEntries; i++)
  {

    tree -> GetEntry(i);
    if( i%1000==0 ) std::cout << "Processing entry "<< i << "\r" << std::flush;
    //    std::cout<<"\n\nEVENT "<<i<<"\n #GEN="<<treeVars.N_GenPart;
    //std::cout<<""<<treeVars.N_GenPart;
    //for(int i=0;i<treeVars.N_GenPart;i++)
    //  std::cout<<"\t"<<treeVars.GenPart_pid[i]<<","<<treeVars.GenPart_st[i]<<","<<treeVars.GenPart_pt[i]<<","<<treeVars.GenPart_mass[i];
    //std::cout<<"\n #GENPH="<<treeVars.N_GenPh;
    //for(int i=0;i<treeVars.N_GenPh;i++)
    //  std::cout<<"\t"<<treeVars.GenPh_st[i]<<","<<treeVars.GenPh_pt[i];
    
    if(treeVars.N_TightPh<2) continue;

    //find higgs daugher gen-photons and flag their .isHdaug
    if( ! FindGenPh_Hdaug(treeVars) ) continue;

    //find lead & sublead reco photons
    int pho_lead_i;
    int pho_sublead_i;
    FindLeadSublead_pho(treeVars,pho_lead_i,pho_sublead_i);
    TLorentzVector pho_lead,pho_sublead;
    pho_lead.SetPtEtaPhiE(treeVars.TightPh_pt[pho_lead_i],treeVars.TightPh_eta[pho_lead_i],treeVars.TightPh_phi[pho_lead_i],treeVars.TightPh_E[pho_lead_i]);
    pho_sublead.SetPtEtaPhiE(treeVars.TightPh_pt[pho_sublead_i],treeVars.TightPh_eta[pho_sublead_i],treeVars.TightPh_phi[pho_sublead_i],treeVars.TightPh_E[pho_sublead_i]);

    //Gen-matching
    if(!PhoGenMatch(pho_lead,pho_sublead,treeVars,outtreeVars,0.1))//default DeltaRmax=0.03
       continue;
    
    //Cuts on photons
    if(!DiPhotonSelection(pho_lead,pho_sublead))
       continue;

    //Jets selections
    if(!JetSelection(treeVars,outtreeVars))//require njetmin with above a minimum pt value, fill also outtreeVars.nJets with the number of hard jets
      continue;
    //cout<<"-----------"<<endl;
    //cout<<treeVars.N_Vtx<<endl;
    //Fill outtree and histos
    outtreeVars.weight = 1.;
    outtreeVars.dipho_mass = (pho_lead+pho_sublead).M();
    outtreeVars.dipho_sumpt = (pho_lead+pho_sublead).Pt();
    outtreeVars.dipho_deltaeta = DeltaEta( pho_lead.Eta() , pho_sublead.Eta() );
    outtreeVars.dipho_deltaphi = DeltaPhi( pho_lead.Phi() , pho_sublead.Phi() );

    outtreeVars.dipho_leadPt = pho_lead.Pt();
    outtreeVars.dipho_leadEta = pho_lead.Eta();
    outtreeVars.dipho_leadPhi = pho_lead.Phi();
    outtreeVars.dipho_leadptoM = pho_lead.Pt() / (pho_lead+pho_sublead).M();
    outtreeVars.dipho_leadEnergy = treeVars.TightPh_E[pho_lead_i];
    outtreeVars.dipho_leadIso = treeVars.TightPh_iso[pho_lead_i];
    outtreeVars.dipho_leadDeltaRgenreco = DeltaR(pho_lead.Eta(),pho_lead.Phi(),outtreeVars.dipho_leadEta_gen,outtreeVars.dipho_leadPhi_gen);
    outtreeVars.dipho_leadDeltaEtagenreco = DeltaEta(pho_lead.Eta(),outtreeVars.dipho_leadEta_gen);
    outtreeVars.dipho_leadDeltaPhigenreco = DeltaPhi(pho_lead.Phi(),outtreeVars.dipho_leadPhi_gen);

    outtreeVars.dipho_subleadPt = pho_sublead.Pt();
    outtreeVars.dipho_subleadEta = pho_sublead.Eta();
    outtreeVars.dipho_subleadPhi = pho_sublead.Phi();
    outtreeVars.dipho_subleadptoM = pho_sublead.Pt() / (pho_lead+pho_sublead).M();
    outtreeVars.dipho_subleadEnergy = treeVars.TightPh_E[pho_sublead_i];
    outtreeVars.dipho_subleadIso = treeVars.TightPh_iso[pho_sublead_i];
    outtreeVars.dipho_subleadDeltaRgenreco = DeltaR(pho_sublead.Eta(),pho_sublead.Phi(),outtreeVars.dipho_subleadEta_gen,outtreeVars.dipho_subleadPhi_gen);
    outtreeVars.dipho_subleadDeltaEtagenreco = DeltaEta(pho_sublead.Eta(),outtreeVars.dipho_subleadEta_gen);
    outtreeVars.dipho_subleadDeltaPhigenreco = DeltaPhi(pho_sublead.Phi(),outtreeVars.dipho_subleadPhi_gen);

    outtreeVars.nvtx = treeVars.N_Vtx;
    outtreeVars.nJets = treeVars.N_Jet;

    h["dipho_mass"] -> Fill(outtreeVars.dipho_mass);
    h["dipho_sumpt"] -> Fill(outtreeVars.dipho_sumpt);
    h["dipho_deltaeta"] -> Fill(outtreeVars.dipho_deltaeta);
    h["dipho_deltaphi"] -> Fill(outtreeVars.dipho_deltaphi);
    h["dipho_leadPt"] -> Fill(outtreeVars.dipho_leadPt);
    h["dipho_leadEta"] -> Fill(outtreeVars.dipho_leadEta);
    h["dipho_leadPhi"] -> Fill(outtreeVars.dipho_leadPhi);
    h["dipho_leadptoM"] -> Fill(outtreeVars.dipho_leadptoM);
    h["dipho_leadIso_o_E"] -> Fill(outtreeVars.dipho_leadIso/outtreeVars.dipho_leadEnergy);
    h["dipho_subleadPt"] -> Fill(outtreeVars.dipho_subleadPt);
    h["dipho_subleadEta"] -> Fill(outtreeVars.dipho_subleadEta);
    h["dipho_subleadPhi"] -> Fill(outtreeVars.dipho_subleadPhi);
    h["dipho_subleadptoM"] -> Fill(outtreeVars.dipho_subleadptoM);
    h["dipho_subleadIso_o_E"] -> Fill(outtreeVars.dipho_subleadIso/outtreeVars.dipho_subleadEnergy);
    h["nJets"] -> Fill(outtreeVars.nJets);
    h["dipho_leaddeltaR_GenReco"] -> Fill( DeltaR(pho_lead.Eta(),pho_lead.Phi(),outtreeVars.dipho_leadEta_gen,outtreeVars.dipho_leadPhi_gen) );
    h["dipho_subleaddeltaR_GenReco"] -> Fill( DeltaR(pho_sublead.Eta(),pho_sublead.Phi(),outtreeVars.dipho_subleadEta_gen,outtreeVars.dipho_subleadPhi_gen) );
    h["dipho_leaddeltaEta_GenReco"] -> Fill( DeltaEta(pho_lead.Eta(),outtreeVars.dipho_leadEta_gen) );
    h["dipho_subleaddeltaEta_GenReco"] -> Fill( DeltaEta(pho_sublead.Eta(),outtreeVars.dipho_subleadEta_gen) );
    h["dipho_leaddeltaPhi_GenReco"] -> Fill( DeltaPhi(pho_lead.Phi(),outtreeVars.dipho_leadPhi_gen) );
    h["dipho_subleaddeltaPhi_GenReco"] -> Fill( DeltaPhi(pho_sublead.Phi(),outtreeVars.dipho_subleadPhi_gen) );

    outTree->Fill();

  }
  std::cout << std::endl;
  outTree -> AutoSave();
  outFile -> Close();

  MakePlot3(h);

  system(Form("mv *.png %s",outputPlotFolder.c_str()));
  system(Form("mv *.pdf %s",outputPlotFolder.c_str()));  

    /*      
    // common cuts
    if( treeVars.dipho_mass < 100 || treeVars.dipho_mass > 180 ) continue;
    if( type == -1 && treeVars.dipho_mass > 115 && treeVars.dipho_mass < 135 ) continue;
    bool passCutBased = CutBasedSelection(treeVars,0.4,0.3,0.,-0.5,2.5,2.);
      
    if( treeVars.nJets < 2 ) continue;
      
    outTree_2jets -> Fill();
      
    // fill event counters - cut based
    if( type == 1 )
        nEvents_cutBased["bkg"][0.] += treeVars.weight;
    if( type == 2 && label == "ttH" )
        nEvents_cutBased["sig"][0.] += treeVars.weight;
    if( passCutBased )
    {
      if( type == 1 )
        nEvents_cutBased["bkg"][1.] += treeVars.weight;
      if( type == 2 && label == "ttH" )
        nEvents_cutBased["sig"][1.] += treeVars.weight;
    }
    if( cutBased && !passCutBased ) continue;

    outTree_2jets -> AutoSave();
    outFile -> Close();
 
    std::cout << "Processed tag " << label << ", " << nEntries << " events out of " << nEntries << std::endl;
  }
  
  
  */   
    return 0;
}


