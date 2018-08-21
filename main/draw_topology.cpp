#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"

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
#include "TLine.h"
#include "TEllipse.h"

using namespace std;

TLine* draw_transv(const vector<TLorentzVector> &vec, int color, int linestyle, TPad *c, float Ptscale);
TLine* draw_long(const vector<TLorentzVector> &vec, int color, int linestyle, TPad *c, float Pscale);

int main(int argc, char* argv[])
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> draw_topology.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }
    
  gStyle -> SetOptFit(0);
  gStyle -> SetOptStat(0);
  
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

  int Nplots=100;
  if(opts.OptExist("Input.Nplots"))
    Nplots = opts.GetOpt<int> ("Input.Nplots");
  else
    cout<<"Option <Input.Nplots> not found --> Analysis by default on Nplots="<<Nplots<<endl;


  //Jets selections
  int BTagOffset;
  if(useMTD == false)
    BTagOffset=0;
  else
    BTagOffset=3;

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
  //output file
  std::string outputPlotFolder = opts.GetOpt<std::string>("Output.outputFolder");
  system(Form("mkdir -p %s",outputPlotFolder.c_str()));
  //TFile* outFile = TFile::Open(Form("%s/plotToplogy_%s.root",outputPlotFolder.c_str(),label.c_str()),"RECREATE");
  //outFile -> cd();
  //TTree* outTree = new TTree("plotTree","plotTree");
  //InitOutTreeVars(outTree,outtreeVars);
  //TCanvas *c_transv = new TCanvas("transverse plane","transverse plane",700,700);  
  //TCanvas *c_long = new TCanvas("longitudinal plane","longitudinal plane",700,700);  
  TreeVars outtreeVars;
  TCanvas *c = new TCanvas ("plot","plot",1460,440);// 730,220);
  TPad* TransvPlanePad = new TPad("TransvPlanePad", "TransvPlanePad", 10./730  , 10./220., 210./730 , 210./220);
  TPad* LongPlanePad =   new TPad("LongPlanePad", "LongPlanePad",     220./730 , 10./220., 720./730 , 210./220);
  TransvPlanePad->Draw();        
  LongPlanePad->Draw();

  //------------------
  // loop over samples

  int nEntries = tree->GetEntries();
  std::cout << "Total entries = " << nEntries << std::endl;
  int Nplotsdone=0;
  for(int ientry=0; ientry<nEntries; ientry++)
  {

    tree -> GetEntry(ientry);
    if( ientry%100==0 ) std::cout << "Processing entry "<< ientry << "\r" << std::flush;
    
    if(treeVars.N_SelectedPh<2) continue;

    
    //find higgs daugher gen-photons and flag their .isHdaug
    //if( ! FindGenPh_Hdaug(treeVars) ) continue;
    //store phgen higgs daughter
    //vector<TLorentzVector> pho_gen;
    //for(int i=0;i<treeVars.N_GenPh;i++)
    //  if(treeVars.GenPh_isHdaug[i]==true)
    //  {
    //  TLorentzVector pho_gen_aux;
    //	pho_gen_aux.SetPtEtaPhiE(treeVars.GenPh_pt[i],treeVars.GenPh_eta[i],treeVars.GenPh_phi[i],treeVars.GenPh_E[i]);
    //	pho_gen.push_back(pho_gen_aux);
    //  }
    

    //find lead & sublead reco photons
    //PrintRecoPhoton(treeVars);    
    int pho_lead_i;
    int pho_sublead_i;
    FindLeadSublead_pho(treeVars,pho_lead_i,pho_sublead_i);
    //store lead & sublead reco photons
    TLorentzVector pho_lead,pho_sublead;
    pho_lead.SetPtEtaPhiE(treeVars.SelectedPh_pt[pho_lead_i],treeVars.SelectedPh_eta[pho_lead_i],treeVars.SelectedPh_phi[pho_lead_i],treeVars.SelectedPh_E[pho_lead_i]);
    pho_sublead.SetPtEtaPhiE(treeVars.SelectedPh_pt[pho_sublead_i],treeVars.SelectedPh_eta[pho_sublead_i],treeVars.SelectedPh_phi[pho_sublead_i],treeVars.SelectedPh_E[pho_sublead_i]);
    vector<TLorentzVector> lead_recoph;
    lead_recoph.push_back(pho_lead);
    lead_recoph.push_back(pho_sublead);
    
    //Cuts on photons
    if(!DiPhotonSelection(pho_lead,pho_sublead))
      continue;
    //cout<<"pass photon selection"<<endl;
    
    //find higgs daugher b-quarks and flag their .isHdaug                                                                                        
    //if( ! Findbquark_Hdaug(treeVars) ) continue;
    //store b-quarks higgs daugher
    //vector<TLorentzVector> bquark_gen;
    //for(int i=0;i<treeVars.N_GenPart;i++)
    //  if(treeVars.GenPart_isHdaug[i]==true)
    //  {
    //	TLorentzVector bquark_gen_aux;
    //	bquark_gen_aux.SetPtEtaPhiE(treeVars.GenPart_pt[i],treeVars.GenPart_eta[i],treeVars.GenPart_phi[i],treeVars.GenPart_E[i]);
    //	bquark_gen.push_back(bquark_gen_aux);
    //  }
    

    
    //find and store recobjets
    vector<TLorentzVector> recobjet;
    for(int i=0;i<treeVars.N_Jet;i++)
    {
      if(treeVars.Jet_hadflav[i]==5)
      {
	TLorentzVector recobjet_aux;
	recobjet_aux.SetPtEtaPhiM(treeVars.Jet_pt[i],treeVars.Jet_eta[i],treeVars.Jet_phi[i],treeVars.Jet_mass[i]);
	recobjet.push_back(recobjet_aux);
      }
    }


    
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
    //select in output only jets with a certain minimum pt and b-tagged
    for(int i=0;i<treeVars.N_Jet;i++)
    {
      if(treeVars.Jet_pt[i]<25) continue;
      if(fabs(treeVars.Jet_eta[i])>3) continue;
      if( DeltaR(treeVars.Jet_eta[i],treeVars.Jet_phi[i],pho_lead.Eta(),pho_lead.Phi()) < 0.4 ) continue;
      if( DeltaR(treeVars.Jet_eta[i],treeVars.Jet_phi[i],pho_sublead.Eta(),pho_sublead.Phi()) < 0.4 ) continue;
      //cout<<i<<"\tbtagvalue="<<treeVars.Jet_mvav2[i]<<"\tbtaglevel="<<BTag<<"\tbtagoffset="<<BTagOffset<<endl;
      //if(BTag>BTagOffset && BTag<4+BTagOffset)
      {
	outtreeVars.nJets++;
	outtreeVars.jet_pt[outtreeVars.nJets-1] = treeVars.Jet_pt[i];
	outtreeVars.jet_eta[outtreeVars.nJets-1] = treeVars.Jet_eta[i];                    
	outtreeVars.jet_phi[outtreeVars.nJets-1] = treeVars.Jet_phi[i];
	outtreeVars.jet_mass[outtreeVars.nJets-1] = treeVars.Jet_mass[i];
	outtreeVars.jet_mvav2[outtreeVars.nJets-1] = treeVars.Jet_mvav2[i];

      }
    }
    


    if(outtreeVars.nJets<2) 
    {
      //cout<<"NOT pass jet selection"<<endl;
      continue;
    }
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
    //cout<<"pass jets invariant mass selection"<<endl;
    //cout<<"\n\n\n\n\n\n"<<endl;

    int BTagMedium_mask;
    if(useMTD == false)
      BTagMedium_mask=0b000010;
    else
      BTagMedium_mask=0b010000;
    //if( outtreeVars.dibjet_leadbtagscore-BTagOffset>=2 && outtreeVars.dibjet_subleadbtagscore-BTagOffset>=2 )
    if( !((outtreeVars.jet_mvav2[bjet_lead_i] & BTagMedium_mask) && (outtreeVars.jet_mvav2[bjet_sublead_i] & BTagMedium_mask)) ) continue;
    //cout<<"HPC"<<endl;
    //getchar();
    vector<TLorentzVector> selected_bjet;
    selected_bjet.push_back(bjet_lead);
    selected_bjet.push_back(bjet_sublead);
    
    LongPlanePad->cd();
    LongPlanePad->DrawFrame(-5,-2,5,2);
    //draw_long(pho_gen,      kCyan,        1,  LongPlanePad,100);
    draw_long(lead_recoph,  kBlue,        2,  LongPlanePad,100);
    //draw_long(bquark_gen,   kMagenta-10,  1,  LongPlanePad,100);
    draw_long(recobjet,     kRed,         4,  LongPlanePad,100);
    draw_long(selected_bjet,kGreen,       5,  LongPlanePad,100);
    //if(ientry<100)
    //  LongPlanePad->Print(Form("longplane_%i.png",ientry));

		      
    TransvPlanePad->cd();
    TransvPlanePad->DrawFrame(-2,-2,2,2);
    TLegend leg(0.1,0.7,0.48,0.9);
    //leg.AddEntry( draw_transv(pho_gen,      kCyan,        1,  TransvPlanePad,100) , "gen photons H daughter","l");
    leg.AddEntry( draw_transv(lead_recoph,  kBlue,        2,  TransvPlanePad,100) , "selected photons","l");
    //leg.AddEntry( draw_transv(bquark_gen,   kMagenta-10,  1,  TransvPlanePad,100) , "b quarks H daughter","l");
    leg.AddEntry( draw_transv(recobjet,     kRed,         4,  TransvPlanePad,100) , "reco b jets","l");
    leg.AddEntry( draw_transv(selected_bjet,kGreen,       5,  TransvPlanePad,100) , "selected b jets","l");
    leg.Draw();
    TEllipse el(0,0,1.5);
    el.SetFillStyle(0);
    el.SetLineColor(1);
    el.Draw();

    if(Nplotsdone<Nplots)
    {
      c->Print(Form("topology_%i.png",ientry));
      ++Nplotsdone;
    }
  }
  std::cout << std::endl;
  /*
  cout<<"\n\n\n\n\nselected events = "<<Nev_selected<<" / "<<nEntries<<endl;
  cout<<"\t\tpass photon selection = "<<Nev_phselection<<" / "<<nEntries<<endl;
  cout<<"\t\tpass kinematic bjet preselection = "<< Nev_jet_kin_preselection <<" / "<<nEntries<<endl;
  cout<<"\t\tpass jet selection = "<<Nev_jetselection<<" / "<<nEntries<<endl;
  cout<<"\n\n\n\n\n"<<endl;
  */

  //outFile -> Close();

  system(Form("mv *.png %s",outputPlotFolder.c_str()));
  system(Form("mv *.pdf %s",outputPlotFolder.c_str()));  
  
  return 0;
}


TLine* draw_transv(const vector<TLorentzVector> &vec, int color, int linestyle, TPad *c, float Ptscale)
{
  TLine* line = new TLine();
  line->SetLineColor(color);
  line->SetLineWidth(3);
  line->SetLineStyle(linestyle);
  for (int i=0;i<vec.size();i++)
  {
    float module = vec[i].Pt()/Ptscale;
    float phi = vec[i].Phi();
    c->cd();
    line->DrawLine(0.,0.,module*cos(phi),module*sin(phi));
  }
  return line;
  
} 


TLine* draw_long(const vector<TLorentzVector> &vec, int color, int linestyle, TPad *c, float Pscale)
{
  TLine* line = new TLine();
  line->SetLineColor(color);
  line->SetLineWidth(3);
  line->SetLineStyle(linestyle);
  for (int i=0;i<vec.size();i++)
  {
    c->cd();
    line->DrawLine(0. , 0. , vec[i].Pz()/Pscale , vec[i].Py()/Pscale );
  }
  return line;
  
} 
