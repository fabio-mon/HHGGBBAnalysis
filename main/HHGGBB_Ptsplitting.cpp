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

//program to take a Plottree ntuple, produce and output file containing several Plottree, one for each Pt interval 

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
  std::vector<float> Ptlimit = opts.GetOpt<std::vector<float>>("Input.Ptlimit");
  std::vector<std::pair<float,float>> Ptrange;
  for(unsigned i=0;i<Ptlimit.size()-1;i++)
    Ptrange.push_back(std::make_pair(Ptlimit[i],Ptlimit[i+1]));
    

  //------------------
  // define histograms
  
  std::map<std::pair<float,float>,TH1F*> h_mass_Ptrange_EB;
  std::map<std::pair<float,float>,TH1F*> h_mass_Ptrange_Hg;
  std::map<std::pair<float,float>,TH1F*> h_Ereco_Egen_ratio_Ptrange_EB;
  std::map<std::pair<float,float>,TH1F*> h_Ereco_Egen_ratio_Ptrange_Hg;
  for(unsigned i=0;i<Ptrange.size();++i)
  {
    float Pt_min = Ptrange.at(i).first;
    float Pt_max = Ptrange.at(i).second;
    h_mass_Ptrange_EB[Ptrange.at(i)] = new TH1F(Form("dipho_mass_EB_Pt%.0f%.0f",Pt_min,Pt_max),Form("dipho_mass_EB_Pt%.0f%.0f",Pt_min,Pt_max),300,100.,150);
    h_mass_Ptrange_Hg[Ptrange.at(i)] = new TH1F(Form("dipho_mass_Hgcal_Pt%.0f%.0f",Pt_min,Pt_max),Form("dipho_mass_Hgcal_Pt%.0f%.0f",Pt_min,Pt_max),300,100.,150);
    h_Ereco_Egen_ratio_Ptrange_EB[Ptrange.at(i)] = new TH1F(Form("Ereco_Egen_diff_EB_Pt%.0f%.0f",Pt_min,Pt_max),Form("Ereco_Egen_diff_EB_Pt%.0f%.0f",Pt_min,Pt_max),200,0.9,1.1);
    h_Ereco_Egen_ratio_Ptrange_Hg[Ptrange.at(i)] = new TH1F(Form("Ereco_Egen_diff_Hgcal_Pt%.0f%.0f",Pt_min,Pt_max),Form("Ereco_Egen_diff_Hgcal_Pt%.0f%.0f",Pt_min,Pt_max),200,0.9,1.1);
  }

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
    }
    
    std::cout << ">>> Adding trees " << inFileName+"/"+treeName << " to chain " << "tree_"+label << std::endl;
    trees[label] -> Add((inFileName+"/"+treeName).c_str());
    types[label] = type;
  }
  
  
  //---------------
  // tree variables
  TreeVars treeVars;
  
  //output trees for plots
  std::string outputPlotFolder = opts.GetOpt<std::string>("Output.outputPlotFolder");
  
  
  //------------------
  // loop over samples
  std::map<std::pair<float,float>,int> nEvents_Ptrange;

  for(std::map<std::string,TChain*>::const_iterator treeIt = trees.begin(); treeIt != trees.end(); ++treeIt)
  {
    std::string label = treeIt -> first;
    TChain* tree = treeIt -> second;
    int type = types[label];
    
    InitTreeVars(tree,treeVars);
    float lumiFactor = type > 0 ? lumi : 1.;
    
    
    
    TFile* outFile = TFile::Open(Form("%s/plotTree_%s.root",outputPlotFolder.c_str(),label.c_str()),"RECREATE");
    outFile -> cd();
    std::map<std::pair<float,float>,TTree*> outTree_Ptrange;
    for(unsigned i=0;i<Ptrange.size();++i)
    {
      float Pt_min = Ptrange.at(i).first;
      float Pt_max = Ptrange.at(i).second;
      outTree_Ptrange[Ptrange.at(i)] = new TTree(Form("plotTree_Pt%.0f%.0f",Pt_min,Pt_max),Form("plotTree_Pt%.0f%.0f",Pt_min,Pt_max));
      InitOutTreeVars(outTree_Ptrange[Ptrange.at(i)],treeVars);
    }
    
    
    int nEntries = tree->GetEntries();
    for(int i=0; i<nEntries; i++)
    {
      tree -> GetEntry(i);
      if( i%1000==0 ) std::cout << "Processing tag " << label << ", event " << i << " out of " << nEntries << "\r" << std::flush;
      unsigned j;
      for(j=0;j<Ptrange.size();++j)
	if(treeVars.dipho_leadPt>=Ptrange.at(j).first && treeVars.dipho_leadPt<Ptrange.at(j).second)
	  break;
      if(treeVars.dipho_leadPt > Ptrange.at(0).first && treeVars.dipho_leadPt < Ptrange.back().second)
      {
	outTree_Ptrange[Ptrange.at(j)] -> Fill();
	if(fabs(treeVars.dipho_leadEta)<1.4)
	{
	  h_mass_Ptrange_EB[Ptrange.at(j)] -> Fill (treeVars.dipho_mass);
	  h_Ereco_Egen_ratio_Ptrange_EB[Ptrange.at(j)] -> Fill(treeVars.dipho_leadEnergy/treeVars.dipho_leadEnergy_gen);
	}
	else
	  if(fabs(treeVars.dipho_leadEta)>1.6)
	  {
	    h_mass_Ptrange_Hg[Ptrange.at(j)] -> Fill (treeVars.dipho_mass);
	    h_Ereco_Egen_ratio_Ptrange_Hg[Ptrange.at(j)] -> Fill(treeVars.dipho_leadEnergy/treeVars.dipho_leadEnergy_gen);
	  }
      }

      for(j=0;j<Ptrange.size();++j)
	if(treeVars.dipho_subleadPt>=Ptrange.at(j).first && treeVars.dipho_subleadPt<Ptrange.at(j).second)
	  break;
      if(treeVars.dipho_subleadPt > Ptrange.at(0).first && treeVars.dipho_subleadPt < Ptrange.back().second)
      {
	if(fabs(treeVars.dipho_leadEta)<1.4)
	  h_Ereco_Egen_ratio_Ptrange_EB[Ptrange.at(j)] -> Fill(treeVars.dipho_subleadEnergy/treeVars.dipho_subleadEnergy_gen);
	else
	  if(fabs(treeVars.dipho_leadEta)>1.6)
	  {
	    h_Ereco_Egen_ratio_Ptrange_Hg[Ptrange.at(j)] -> Fill(treeVars.dipho_subleadEnergy/treeVars.dipho_subleadEnergy_gen);
	  }
      }

    }
    for(unsigned i=0;i<Ptrange.size();++i)
      outTree_Ptrange[Ptrange.at(i)] -> AutoSave();
    outFile -> Close();
    
    std::cout << "Processed tag " << label << ", " << nEntries << " events out of " << nEntries << std::endl;
  }

  //Get dipho_mass_resolution as a func of lead photon Pt
  TGraphErrors *gr_mass_res_Pt_EB = new TGraphErrors();
  TGraphErrors *gr_mass_res_Pt_Hg = new TGraphErrors();
  gr_mass_res_Pt_EB->SetTitle("#sigma(M^{#gamma#gamma}) vs lead photon Pt EB");
  gr_mass_res_Pt_Hg->SetTitle("#sigma(M^{#gamma#gamma}) vs lead photon Pt Hgcal");
  for(unsigned i=0;i<Ptrange.size();++i)
  {
    float Pt_min = Ptrange.at(i).first;
    float Pt_max = Ptrange.at(i).second;
    gr_mass_res_Pt_EB->SetPoint(i,(Pt_min+Pt_max)/2,h_mass_Ptrange_EB[Ptrange.at(i)]->GetRMS());
    gr_mass_res_Pt_EB->SetPointError(i,(Pt_max-Pt_min)/2,h_mass_Ptrange_EB[Ptrange.at(i)]->GetRMSError());
    gr_mass_res_Pt_Hg->SetPoint(i,(Pt_min+Pt_max)/2,h_mass_Ptrange_Hg[Ptrange.at(i)]->GetRMS());
    gr_mass_res_Pt_Hg->SetPointError(i,(Pt_max-Pt_min)/2,h_mass_Ptrange_Hg[Ptrange.at(i)]->GetRMSError());
  }

  //Get Energy_resolution as a func of lead photon Pt
  TGraphErrors *gr_E_res_Pt_EB = new TGraphErrors();
  TGraphErrors *gr_E_res_Pt_Hg = new TGraphErrors();
  gr_E_res_Pt_EB->SetTitle("#sigma_{E}/E vs photon Pt EB");
  gr_E_res_Pt_Hg->SetTitle("#sigma_{E}/E vs photon Pt Hgcal");
  for(unsigned i=0;i<Ptrange.size();++i)
  {
    float Pt_min = Ptrange.at(i).first;
    float Pt_max = Ptrange.at(i).second;
    gr_E_res_Pt_EB->SetPoint(i,(Pt_min+Pt_max)/2,h_Ereco_Egen_ratio_Ptrange_EB[Ptrange.at(i)]->GetRMS());
    gr_E_res_Pt_EB->SetPointError(i,(Pt_max-Pt_min)/2,h_Ereco_Egen_ratio_Ptrange_EB[Ptrange.at(i)]->GetRMSError());
    gr_E_res_Pt_Hg->SetPoint(i,(Pt_min+Pt_max)/2,h_Ereco_Egen_ratio_Ptrange_Hg[Ptrange.at(i)]->GetRMS());
    gr_E_res_Pt_Hg->SetPointError(i,(Pt_max-Pt_min)/2,h_Ereco_Egen_ratio_Ptrange_Hg[Ptrange.at(i)]->GetRMSError());
  }

  //Draw plots
  TCanvas *c1 = new TCanvas();
  c1->cd();
  for(unsigned i=0;i<Ptrange.size();++i)
  {
    h_mass_Ptrange_EB[Ptrange.at(i)] -> Draw();
    h_mass_Ptrange_EB[Ptrange.at(i)] -> GetXaxis() -> SetTitle("diphoton mass (GeV)");
    c1->Print(Form("c_dipho_mass_Pt_EB%.0f_%.0f.png",Ptrange.at(i).first,Ptrange.at(i).second));
    c1->Print(Form("c_dipho_mass_Pt_EB%.0f_%.0f.pdf",Ptrange.at(i).first,Ptrange.at(i).second));

    h_mass_Ptrange_Hg[Ptrange.at(i)] -> Draw();
    h_mass_Ptrange_Hg[Ptrange.at(i)] -> GetXaxis() -> SetTitle("diphoton mass (GeV)");
    c1->Print(Form("c_dipho_mass_Pt_Hg%.0f_%.0f.png",Ptrange.at(i).first,Ptrange.at(i).second));
    c1->Print(Form("c_dipho_mass_Pt_Hg%.0f_%.0f.pdf",Ptrange.at(i).first,Ptrange.at(i).second));
  }

  gr_mass_res_Pt_EB->SetMarkerStyle(20);
  gr_mass_res_Pt_EB->Draw("AP");
  gr_mass_res_Pt_EB-> GetXaxis() -> SetTitle("lead phot P_{T} (GeV)");
  gr_mass_res_Pt_EB-> GetYaxis() -> SetTitle("#sigma(M_{#gamma#gamma}) (GeV)");
  gr_mass_res_Pt_EB-> GetYaxis() -> SetRangeUser(1.5,4.5);
  gr_mass_res_Pt_EB-> Draw("AP");
  c1->Update();
  c1->Print("c_dipho_mass_res_Pt_EB.png");
  c1->Print("c_dipho_mass_res_Pt_EB.pdf");


  gr_mass_res_Pt_Hg->SetMarkerStyle(20);
  gr_mass_res_Pt_Hg->Draw("AP");
  gr_mass_res_Pt_Hg-> GetXaxis() -> SetTitle("lead phot P_{T} (GeV)");
  gr_mass_res_Pt_Hg-> GetYaxis() -> SetTitle("#sigma(M_{#gamma#gamma}) (GeV)");
  gr_mass_res_Pt_Hg-> GetYaxis() -> SetRangeUser(1.5,4.5);
  gr_mass_res_Pt_Hg-> Draw("AP");
  c1->Update();
  c1->Print("c_dipho_mass_res_Pt_Hg.png");
  c1->Print("c_dipho_mass_res_Pt_Hg.pdf");

  for(unsigned i=0;i<Ptrange.size();++i)
  {
    h_Ereco_Egen_ratio_Ptrange_EB[Ptrange.at(i)] -> Draw();
    h_Ereco_Egen_ratio_Ptrange_EB[Ptrange.at(i)] -> GetXaxis() -> SetTitle("E/E_{TRUE}");
    c1->Print(Form("c_Ereco_Egen_ratio_Pt_EB%.0f_%.0f.png",Ptrange.at(i).first,Ptrange.at(i).second));
    c1->Print(Form("c_Ereco_Egen_ratio_Pt_EB%.0f_%.0f.pdf",Ptrange.at(i).first,Ptrange.at(i).second));
  }

  for(unsigned i=0;i<Ptrange.size();++i)
  {
    h_Ereco_Egen_ratio_Ptrange_Hg[Ptrange.at(i)] -> Draw();
    h_Ereco_Egen_ratio_Ptrange_Hg[Ptrange.at(i)] -> GetXaxis() -> SetTitle("E/E_{TRUE}");
    c1->Print(Form("c_Ereco_Egen_ratio_Pt_Hg%.0f_%.0f.png",Ptrange.at(i).first,Ptrange.at(i).second));
    c1->Print(Form("c_Ereco_Egen_ratio_Pt_Hg%.0f_%.0f.pdf",Ptrange.at(i).first,Ptrange.at(i).second));
  }

  c1->SetLogx();
  gr_E_res_Pt_EB->SetMarkerStyle(20);
  gr_E_res_Pt_EB->Draw("AP");
  gr_E_res_Pt_EB-> GetXaxis() -> SetTitle("phot P_{T} (GeV)");
  gr_E_res_Pt_EB-> GetXaxis() -> SetRangeUser(10,2000);
  gr_E_res_Pt_EB-> GetYaxis() -> SetTitle("#sigma_{E}/E");
  gr_E_res_Pt_EB-> GetYaxis() -> SetRangeUser(0,0.04);
  gr_E_res_Pt_EB->Draw("AP");
  c1->Update();
  c1->Print("c_dipho_E_res_Pt_EB.png");
  c1->Print("c_dipho_E_res_Pt_EB.pdf");

  c1->SetLogx();
  gr_E_res_Pt_Hg->SetMarkerStyle(20);
  gr_E_res_Pt_Hg->Draw("AP");
  gr_E_res_Pt_Hg-> GetXaxis() -> SetTitle("phot P_{T} (GeV)");
  gr_E_res_Pt_Hg-> GetXaxis() -> SetRangeUser(10,2000);
  gr_E_res_Pt_Hg-> GetYaxis() -> SetTitle("#sigma_{E}/E");
  gr_E_res_Pt_Hg-> GetYaxis() -> SetRangeUser(0,0.04);
  gr_E_res_Pt_Hg->Draw("AP");
  c1->Update();
  c1->Print("c_dipho_E_res_Pt_Hg.png");
  c1->Print("c_dipho_E_res_Pt_Hg.pdf");
  TFile *outgraph = new TFile("outgraph.root","RECREATE");
  outgraph->cd();
  gr_E_res_Pt_EB->Write("gr_E_res_Pt_EB");
  gr_E_res_Pt_Hg->Write("gr_E_res_Pt_Hg");
  outgraph->Close();

  system(Form("mv *.png %s",outputPlotFolder.c_str()));
  system(Form("mv *.pdf %s",outputPlotFolder.c_str()));
 
}
