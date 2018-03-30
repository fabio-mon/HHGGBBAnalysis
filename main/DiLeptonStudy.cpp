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
    std::cerr << ">>>>> DiLeptonStudy.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }
  
  
  //------------------
  // graphics settings
  
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = Form("%.1f fb^{-1} (13 TeV)",lumi);
  
  setTDRStyle();
  gStyle -> SetOptFit(0);
  gStyle -> SetOptStat(0);
  
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  
  
  //----------------------
  // parse the config file
  
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  std::vector<std::string> input = opts.GetOpt<std::vector<std::string> >("Input.input");
  
  
  //------------------------
  // define roofit variables
  
  RooRealVar mass_("mass_", "m_{#gamma#gamma}", 100, 180);
  RooRealVar weight_("weight_", "weight", -100, 100);
  RooRealVar mva_("mva_", "mva", -1., 1.);  
  
  RooDataSet dataOneCat("dataOneCat","dataOneCat", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet csOneCat("csOneCat","csOneCat", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet tthOneCat("tthOneCat","tthOneCat", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  
  RooDataSet dataDiLepton("dataDiLepton","dataDiLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet csDiLepton("csDiLepton","csDiLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet tthDiLepton("tthDiLepton","tthDiLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  
  RooDataSet dataSingleLepton("dataSingleLepton","dataSingleLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet csSingleLepton("csSingleLepton","csSingleLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet tthSingleLepton("tthSingleLepton","tthSingleLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  
  
  //------------------
  // define histograms
  
  std::map<std::string,TH1F*> mass_oneCat_histo;
  std::map<std::string,TH1F*> mass_diMu_histo;
  std::map<std::string,TH1F*> mass_diEle_histo;
  std::map<std::string,TH1F*> mass_Mixed_histo;
  std::map<std::string,TH1F*> mass_diLepton_histo;
  std::map<std::string,TH1F*> mass_singleLepton_histo;
  
  
  //----------
  // get trees
  
  std::map<std::string,TChain*> trees;
  std::map<std::string,int> types; // -1 = data; -2 = control sample; 1 = MC bkg; 2 = MC signal;

  for(unsigned int n = 0; n < input.size()/4; ++n)
  {
    std::string inFileName = input.at(0+n*4);
    std::string treeName   = input.at(1+n*4);
    std::string label      = input.at(2+n*4);
    int type               = atoi(input.at(3+n*4).c_str());
    
    if( trees[label] == 0 )
    {
      trees[label] = new TChain(Form("tree_%s",label.c_str()),"");
      
      mass_oneCat_histo[label]    = new TH1F(("mass_oneCat_histo"   +label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
      mass_diMu_histo[label]      = new TH1F(("mass_diMu_histo"     +label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
      mass_diEle_histo[label]     = new TH1F(("mass_diEle_histo"    +label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
      mass_Mixed_histo[label]     = new TH1F(("mass_Mixed_histo"    +label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
      mass_diLepton_histo[label]  = new TH1F(("mass_diLepton_histo" +label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
      mass_singleLepton_histo[label] = new TH1F(("mass_singleLepton_histo"+label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
    }
    
    std::cout << ">>> Adding trees " << inFileName+"/"+treeName << " to chain " << "tree_"+label << std::endl;
    trees[label] -> Add((inFileName+"/"+treeName).c_str());
    types[label] = type;
  }
  
  
  //---------------
  // tree variables
  TreeVars treeVars;
  
  
  
  //------------------
  // loop over samples
  std::map<std::string,std::map<float,float> > nEvents_cutBased;
  std::map<std::string,std::map<float,float> > nEvents_mvaCut;
  std::map<std::string,std::map<float,float> > nEvents_mvaCut_new;
  for(std::map<std::string,TChain*>::const_iterator treeIt = trees.begin(); treeIt != trees.end(); ++treeIt)
  {
    std::string label = treeIt -> first;
    TChain* tree = treeIt -> second;
    int type = types[label];
    
    InitTreeVars(tree,treeVars);
    float lumiFactor = type > 0 ? lumi : 1.;
    
    
    //output trees for plots
    std::string outputPlotFolder = opts.GetOpt<std::string>("Output.outputPlotFolder");
    
    TFile* outFile = TFile::Open(Form("%s/plotTree_%s.root",outputPlotFolder.c_str(),label.c_str()),"RECREATE");
    outFile -> cd();
    TTree* outTree_2jets = new TTree("plotTree_2jets","plotTree_2jets");
    TTree* outTree_2jets_1bTagLoose = new TTree("plotTree_2jets_1bTagLoose","plotTree_2jets_1bTagLoose");
    TTree* outTree_diLepton = new TTree("plotTree_diLepton","plotTree_diLepton");
    TTree* outTree_singleLepton = new TTree("plotTree_singleLepton","plotTree_singleLepton");
    TTree* outTree_oneCategory = new TTree("plotTree_oneCategory","plotTree_oneCategory");
    InitOutTreeVars(outTree_2jets,treeVars);
    InitOutTreeVars(outTree_2jets_1bTagLoose,treeVars);
    InitOutTreeVars(outTree_diLepton,treeVars);
    InitOutTreeVars(outTree_singleLepton,treeVars);
    InitOutTreeVars(outTree_oneCategory,treeVars);
    
    
    int nEntries = tree->GetEntries();
    for(int i=0; i<nEntries; i++)
    {
      tree -> GetEntry(i);
      if( i%1000==0 ) std::cout << "Processing tag " << label << ", event " << i << " out of " << nEntries << "\r" << std::flush;
      
      // common cuts
      if( treeVars.dipho_mass < 100 || treeVars.dipho_mass > 180 ) continue;
      if( type == -1 && treeVars.dipho_mass > 115 && treeVars.dipho_mass < 135 ) continue;
      bool passCutBased = CutBasedSelection(treeVars,0.4,0.3,0.,-0.5,2.5,2.);
      
      if( treeVars.nJets < 2 ) continue;
      
      outTree_2jets -> Fill();
      
      if( treeVars.nJets_bTagLoose < 1 ) continue;
      
      outTree_2jets_1bTagLoose -> Fill();
      
      
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
      
      // fill event counters - mva
      for(float mvaCut = -1.; mvaCut < 1.; mvaCut+=0.025)
      {
        if( treeVars.dipho_mva > mvaCut )
        {
          if( type == 1 )
            nEvents_mvaCut["bkg"][mvaCut] += treeVars.weight;
          if( type == 2 && label == "ttH" )
            nEvents_mvaCut["sig"][mvaCut] += treeVars.weight;
        }
        
        // if( treeVars.dipho_mva_new > mvaCut )
        if( treeVars.dipho_mva > mvaCut )
        {
          if( type == 1 )
            nEvents_mvaCut_new["bkg"][mvaCut] += treeVars.weight;
          if( type == 2 && label == "ttH" )
            nEvents_mvaCut_new["sig"][mvaCut] += treeVars.weight;
        }
      }
      
      
      if( cutBased && !passCutBased ) continue;
      
      
      //-------------
      // one category
      
      if( OneCategorySelection(treeVars,type) == true )
      {
        mass_oneCat_histo[label] -> Fill(treeVars.dipho_mass, treeVars.weight*lumiFactor);
        
        mass_ = treeVars.dipho_mass;
        // mva_ = treeVars.dipho_mva_new;
        mva_ = treeVars.dipho_mva;
        if( label == "ttH"            ) tthOneCat.add(RooArgSet(mass_,mva_),  treeVars.weight*lumiFactor);
        if( label == "CS_oneCategory" ) csOneCat.add(RooArgSet(mass_,mva_),   treeVars.weight*lumiFactor);
        if( label == "data"           ) dataOneCat.add(RooArgSet(mass_,mva_), treeVars.weight*lumiFactor);
        
        outTree_oneCategory -> Fill();
      } // one category
      
      
      //--------------------------
      // two categories - dilepton
      
      DiLeptonCategories cat = None;
      if( DiLeptonSelection(treeVars,type,cat) )
      {
        mass_diLepton_histo[label] -> Fill(treeVars.dipho_mass, treeVars.weight*lumiFactor);
        
        mass_ = treeVars.dipho_mass;
        // mva_ = treeVars.dipho_mva_new;
        mva_ = treeVars.dipho_mva;
        if( label == "ttH"         ) tthDiLepton.add(RooArgSet(mass_,mva_),  treeVars.weight*lumiFactor);
        if( label == "CS_diLepton" ) csDiLepton.add(RooArgSet(mass_,mva_),   treeVars.weight*lumiFactor);
        if( label == "data"        ) dataDiLepton.add(RooArgSet(mass_,mva_), treeVars.weight*lumiFactor);
        
        outTree_diLepton -> Fill();
        
        if( type != -2 ) continue;
      } // two categories - dilepton
      
      
      //-------------------------------
      // two categories - single lepton
      
      if( SingleLeptonSelection(treeVars,type) )
      {
        mass_singleLepton_histo[label] -> Fill(treeVars.dipho_mass, treeVars.weight*lumiFactor);
        
        mass_ = treeVars.dipho_mass;
        // mva_ = treeVars.dipho_mva_new;
        mva_ = treeVars.dipho_mva;
        if( label == "ttH"             ) tthSingleLepton.add(RooArgSet(mass_,mva_),  treeVars.weight*lumiFactor);
        if( label == "CS_singleLepton" ) csSingleLepton.add(RooArgSet(mass_,mva_),   treeVars.weight*lumiFactor);
        if( label == "data"            ) dataSingleLepton.add(RooArgSet(mass_,mva_), treeVars.weight*lumiFactor);
        
        outTree_singleLepton -> Fill();
      } // two categories - single lepton
      
      
    } // loop over events

    outTree_2jets -> AutoSave();
    outTree_2jets_1bTagLoose -> AutoSave();
    outTree_diLepton -> AutoSave();
    outTree_singleLepton -> AutoSave();
    outTree_oneCategory -> AutoSave();
    outFile -> Close();
    
    std::cout << "Processed tag " << label << ", " << nEntries << " events out of " << nEntries << std::endl;
  }
  
  
  MakePlot2(mass_oneCat_histo, "oneCategory");
  // MakePlot2(mass_diMu_histo, "diMuon");
  // MakePlot2(mass_diEle_histo, "diElectron");
  // MakePlot2(mass_Mixed_histo, "mixed");
  MakePlot2(mass_diLepton_histo, "diLepton");
  MakePlot2(mass_singleLepton_histo, "singleLepton");
  
  
  system("mv *.png ~/www/ttH/DiLeptonStudy_new/");
  system("mv *.pdf ~/www/ttH/DiLeptonStudy_new/");
  
  
  
  std::cout << "#########################################################################" << std::endl;
  std::cout << "################################## ROC ##################################" << std::endl;
  
  TFile* outFile_ROC = new TFile("ROC.root","RECREATE");
  outFile_ROC -> cd();
  
  TGraph* g_ROC_cutBased = new TGraph();
  g_ROC_cutBased -> SetPoint(g_ROC_cutBased->GetN(),nEvents_cutBased["sig"][1.]/nEvents_cutBased["sig"][0.],1.-nEvents_cutBased["bkg"][1.]/nEvents_cutBased["bkg"][0.]);
  
  TGraph* g_ROC_mva = new TGraph();
  TGraph* g_ROC_mva_new = new TGraph();
  for(float mvaCut = -1.; mvaCut < 1.; mvaCut+=0.025)
  {
    g_ROC_mva -> SetPoint(g_ROC_mva->GetN(),nEvents_mvaCut["sig"][mvaCut]/nEvents_mvaCut["sig"][-1.],1.-nEvents_mvaCut["bkg"][mvaCut]/nEvents_mvaCut["bkg"][-1.]);
    g_ROC_mva_new -> SetPoint(g_ROC_mva_new->GetN(),nEvents_mvaCut_new["sig"][mvaCut]/nEvents_mvaCut_new["sig"][-1.],1.-nEvents_mvaCut_new["bkg"][mvaCut]/nEvents_mvaCut_new["bkg"][-1.]);
  }
  
  g_ROC_cutBased -> Write("g_ROC_cutBased");
  g_ROC_mva -> Write("g_ROC_mva");
  g_ROC_mva_new -> Write("g_ROC_mva_new");
  outFile_ROC -> Close();
  
  
  
  bool doFit = 1;
  bool diSimultaneous = 0;
  
  if(doFit)
  {
    float oneCatSignificance = makeFits(&tthOneCat, &csOneCat, -1., "OneCategory", 0, 0, 0);
    float diLeptonSignificance = makeFits(&tthDiLepton, &csDiLepton, -1., "DiLepton", 0, 0, 0);
    float singleLeptonSignificance = makeFits(&tthSingleLepton, &csSingleLepton, -1., "SingleLepton", 0, 0, 0);
    
    float purityOneCat = (mass_oneCat_histo["ttH"] -> Integral() / (mass_oneCat_histo["ttH"]->Integral() + mass_oneCat_histo["ggH"]->Integral() + mass_oneCat_histo["VBF"]->Integral() + mass_oneCat_histo["VH"]->Integral() + mass_oneCat_histo["bbH"]->Integral() + mass_oneCat_histo["tHq"]->Integral() + mass_oneCat_histo["tHW"]->Integral() ) )*100.;
    float purityDiLepton = (mass_diLepton_histo["ttH"] -> Integral() / (mass_diLepton_histo["ttH"]->Integral() + mass_diLepton_histo["ggH"]->Integral() + mass_diLepton_histo["VBF"]->Integral() + mass_diLepton_histo["VH"]->Integral() + mass_diLepton_histo["bbH"]->Integral() + mass_diLepton_histo["tHq"]->Integral() + mass_diLepton_histo["tHW"]->Integral() ) )*100.;
    float puritySingleLepton = (mass_singleLepton_histo["ttH"] -> Integral() / (mass_singleLepton_histo["ttH"]->Integral() + mass_singleLepton_histo["ggH"]->Integral() + mass_singleLepton_histo["VBF"]->Integral() + mass_singleLepton_histo["VH"]->Integral() + mass_singleLepton_histo["bbH"]->Integral() + mass_singleLepton_histo["tHq"]->Integral() + mass_singleLepton_histo["tHW"]->Integral() ) )*100.;
    
    std::cout << "#################################################################################" << std::endl;
    std::cout << "################################## FIT RESULTS ##################################" << std::endl;
    std::cout << "Significance with one leptonic category: " << oneCatSignificance << std::endl;
    std::cout << "Significance with two leptonic categories: " << sqrt(diLeptonSignificance*diLeptonSignificance + singleLeptonSignificance*singleLeptonSignificance) << " (" << diLeptonSignificance <<  " from dileptonic category, " << singleLeptonSignificance << " from one lepton category)" << std::endl << std::endl;
    std::cout << "Signal events diLepton:        " <<  mass_diLepton_histo["ttH"] -> Integral() << " events, tag purity: " << purityDiLepton     << std::endl; 
    std::cout << "Signal events one Lepton:      " << mass_singleLepton_histo["ttH"] -> Integral() << " events, tag purity: " << puritySingleLepton << std::endl; 
    std::cout << "Signal events single category: " <<    mass_oneCat_histo["ttH"] -> Integral() << " events, tag purity: " << purityOneCat       << std::endl; 
    
    if(diSimultaneous)
    {
      float significanceCombined = makeFitSimulataneous(&tthDiLepton, &csDiLepton, &tthSingleLepton, &csSingleLepton, "ChiLoSa");
      std::cout << std::endl << "Significance from the simultanoeus fit of the two categories " << significanceCombined << std::endl;
    }
    std::cout << "#################################################################################" << std::endl;
    
    
    bool doMVAScan = true;
    if( doMVAScan )
    {
      std::cout << std::endl;
      std::cout << "####################################################################################" << std::endl;
      std::cout << "################################## DOING MVA SCAN ##################################" << std::endl;
      
      TGraph* g_oneCatSignificance = new TGraph();
      TGraph* g_diLeptonSignificance = new TGraph();
      TGraph* g_singleLeptonSignificance = new TGraph();
      
      TCut massCut = Form("mass_ > 110");
      RooDataSet* tthOneCat_red = (RooDataSet*)tthOneCat.reduce(massCut);
      RooDataSet* csOneCat_red  = (RooDataSet*)csOneCat.reduce(massCut);
      RooDataSet* tthDiLepton_red = (RooDataSet*)tthDiLepton.reduce(massCut);
      RooDataSet* csDiLepton_red  = (RooDataSet*)csDiLepton.reduce(massCut);
      RooDataSet* tthSingleLepton_red = (RooDataSet*)tthSingleLepton.reduce(massCut);
      RooDataSet* csSingleLepton_red  = (RooDataSet*)csSingleLepton.reduce(massCut);
      
      float sumEntries_sig_oneCat = tthOneCat_red -> sumEntries();
      float sumEntries_bkg_oneCat = csOneCat_red -> sumEntries();
      float sumEntries_sig_diLepton = tthDiLepton_red -> sumEntries();
      float sumEntries_bkg_diLepton = csDiLepton_red -> sumEntries();
      float sumEntries_sig_singleLepton = tthSingleLepton_red -> sumEntries();
      float sumEntries_bkg_singleLepton = csSingleLepton_red -> sumEntries();
      
      TGraph* g_sigEff_oneCat = new TGraph();
      TGraph* g_bkgEff_oneCat = new TGraph();
      TGraph* g_sigEff_diLepton = new TGraph();
      TGraph* g_bkgEff_diLepton = new TGraph();
      TGraph* g_sigEff_singleLepton = new TGraph();
      TGraph* g_bkgEff_singleLepton = new TGraph();
      
      for(float mvaMin = -1.; mvaMin < 1.; mvaMin+=0.1)
      {
        oneCatSignificance = makeFits(&tthOneCat, &csOneCat, mvaMin, "OneCategory", 0, 0, 0);
        diLeptonSignificance = makeFits(&tthDiLepton, &csDiLepton, mvaMin, "DiLepton", 0, 0, 0);
        singleLeptonSignificance = makeFits(&tthSingleLepton, &csSingleLepton, mvaMin, "SingleLepton", 0, 0, 0);    
        
        g_oneCatSignificance -> SetPoint(g_oneCatSignificance->GetN(),mvaMin,oneCatSignificance);
        g_diLeptonSignificance -> SetPoint(g_diLeptonSignificance->GetN(),mvaMin,diLeptonSignificance);
        g_singleLeptonSignificance -> SetPoint(g_singleLeptonSignificance->GetN(),mvaMin,singleLeptonSignificance);
        
        TCut mvaCut = Form("mva_ > %f && mass_ > 110",mvaMin);
        tthOneCat_red = (RooDataSet*)tthOneCat.reduce(mvaCut);
        csOneCat_red  = (RooDataSet*)csOneCat.reduce(mvaCut);
        tthDiLepton_red = (RooDataSet*)tthDiLepton.reduce(mvaCut);
        csDiLepton_red  = (RooDataSet*)csDiLepton.reduce(mvaCut);
        tthSingleLepton_red = (RooDataSet*)tthSingleLepton.reduce(mvaCut);
        csSingleLepton_red  = (RooDataSet*)csSingleLepton.reduce(mvaCut);
        
        g_sigEff_oneCat -> SetPoint(g_sigEff_oneCat->GetN(),mvaMin,tthOneCat_red->sumEntries()/sumEntries_sig_oneCat);
        g_bkgEff_oneCat -> SetPoint(g_bkgEff_oneCat->GetN(),mvaMin,csOneCat_red->sumEntries()/sumEntries_bkg_oneCat);
        g_sigEff_diLepton -> SetPoint(g_sigEff_diLepton->GetN(),mvaMin,tthDiLepton_red->sumEntries()/sumEntries_sig_diLepton);
        g_bkgEff_diLepton -> SetPoint(g_bkgEff_diLepton->GetN(),mvaMin,csDiLepton_red->sumEntries()/sumEntries_bkg_diLepton);
        g_sigEff_singleLepton -> SetPoint(g_sigEff_singleLepton->GetN(),mvaMin,tthSingleLepton_red->sumEntries()/sumEntries_sig_singleLepton);
        g_bkgEff_singleLepton -> SetPoint(g_bkgEff_singleLepton->GetN(),mvaMin,csSingleLepton_red->sumEntries()/sumEntries_bkg_singleLepton);
      }
      
      TCanvas* c = new TCanvas();
      c -> cd();
      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,1.,2.) );
      hPad -> SetTitle(";diphoton MVA cut;significance");
      
      g_oneCatSignificance -> SetLineColor(kBlack);
      g_oneCatSignificance -> SetLineStyle(2);
      g_oneCatSignificance -> Draw("L,same");
      g_diLeptonSignificance -> SetLineColor(kRed);
      g_diLeptonSignificance -> Draw("L,same");
      g_singleLeptonSignificance -> SetLineColor(kMagenta);
      g_singleLeptonSignificance -> Draw("L,same");
      
      c -> SaveAs("c_significance_vs_mva.png");
      c -> SaveAs("c_significance_vs_mva.pdf");
      
      
      c = new TCanvas();
      c -> cd();
      c -> SetLogy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,0.01,1.,1.) );
      hPad -> SetTitle(";diphoton MVA cut;efficiency");
      
      g_sigEff_oneCat -> SetLineColor(kRed);
      g_sigEff_oneCat -> Draw("L,same");
      g_bkgEff_oneCat -> SetLineColor(kBlack);
      g_bkgEff_oneCat -> Draw("L,same");
      
      c -> SaveAs("c_eff_vs_mva_oneCategory.png");
      c -> SaveAs("c_eff_vs_mva_oneCategory.pdf");
      
      
      c = new TCanvas();
      c -> cd();
      c -> SetLogy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,0.01,1.,1.) );
      hPad -> SetTitle(";diphoton MVA cut;efficiency");
      
      g_sigEff_diLepton -> SetLineColor(kRed);
      g_sigEff_diLepton -> Draw("L,same");
      g_bkgEff_diLepton -> SetLineColor(kBlack);
      g_bkgEff_diLepton -> Draw("L,same");
      
      c -> SaveAs("c_eff_vs_mva_diLepton.png");
      c -> SaveAs("c_eff_vs_mva_diLepton.pdf");
      
      
      c = new TCanvas();
      c -> cd();
      c -> SetLogy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,0.01,1.,1.) );
      hPad -> SetTitle(";diphoton MVA cut;efficiency");
      
      g_sigEff_singleLepton -> SetLineColor(kRed);
      g_sigEff_singleLepton -> Draw("L,same");
      g_bkgEff_singleLepton -> SetLineColor(kBlack);
      g_bkgEff_singleLepton -> Draw("L,same");
      
      c -> SaveAs("c_eff_vs_mva_singleLepton.png");
      c -> SaveAs("c_eff_vs_mva_singleLepton.pdf");
    }
    
    system("mv *.png ~/www/ttH/DiLeptonStudy_new/Optimization/");
    system("mv *.pdf ~/www/ttH/DiLeptonStudy_new/Optimization/");
  }
}
