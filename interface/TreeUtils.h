#ifndef TREE_UTILS_H
#define TREE_UTILS_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "TChain.h"

/*** tree variables ***/
struct TreeVars
{

  int nvtx;
  float weight;
  float cross_sec;
  int event;
  
  float dipho_sumpt;
  float dipho_mass;
  float dipho_mass_gen;
  float dipho_deltaeta;
  float dipho_deltaphi;
  float dipho_leadPt;
  float dipho_leadEta;
  float dipho_leadEta_gen;
  float dipho_leadPhi;
  float dipho_leadPhi_gen;
  float dipho_leadEnergy;
  float dipho_leadEnergy_gen;
  float dipho_leadptoM;
  float dipho_lead_sigmaEoE;
  float dipho_leadDeltaRgenreco;
  float dipho_leadDeltaEtagenreco;
  float dipho_leadDeltaPhigenreco;
  float dipho_subleadPt;
  float dipho_subleadEta;
  float dipho_subleadEta_gen;
  float dipho_subleadPhi;
  float dipho_subleadPhi_gen;
  float dipho_subleadEnergy;
  float dipho_subleadEnergy_gen;
  float dipho_subleadptoM;
  float dipho_sublead_sigmaEoE;
  float dipho_subleadDeltaRgenreco;
  float dipho_subleadDeltaEtagenreco;
  float dipho_subleadDeltaPhigenreco;

  int nJets;
  int nJets_bTagLoose;
  int nJets_bTagMedium;
  int nJets_bTagTight;

  float dibjet_mass;
  float dibjet_sumpt;
  float dibjet_deltaeta;
  float dibjet_deltaphi;
  float dibjet_leadPt;
  float dibjet_leadEta;
  float dibjet_leadPhi;
  float dibjet_leadptoM;
  float dibjet_leadEnergy;
  int   dibjet_leadbtagmedium;
  int   dibjet_leadmvav2;
  float dibjet_subleadPt;
  float dibjet_subleadEta;
  float dibjet_subleadPhi;
  float dibjet_subleadptoM;
  float dibjet_subleadEnergy;
  int   dibjet_subleadbtagmedium;
  int   dibjet_subleadmvav2;

  float Mx;
  float DRmin_pho_bjet; 
  float costheta_HH; 
  float costheta_gg; 
  float costheta_bb; 

  float jet_pt[200];
  float jet_eta[200];
  float jet_phi[200];
  int   jet_bdiscriminant[200];
  int   jet_BTagMedium[200];
  int   jet_mvav2[200];
  float jet_mass[200];
  int   jet_hadflav[200];
  
  float MetPt;
  float MetPhi;
  
  int cut_based_ct;
  int ttHTagger;
};

struct RawTreeVars
{
  static constexpr int maxpart=50;
  static constexpr int maxjets=200;
  static constexpr int maxweights=500;

   Int_t run,event,lumi;

   //gen level event
   int   N_Weight;
   float Weight[maxweights];

   int   N_GenPart;
   int   GenPart_pid[maxpart];
   int   GenPart_ch[maxpart];
   int   GenPart_st[maxpart];
   float GenPart_p[maxpart];
   float GenPart_px[maxpart];
   float GenPart_py[maxpart];
   float GenPart_pz[maxpart];
   float GenPart_E[maxpart];
   float GenPart_pt[maxpart];
   float GenPart_eta[maxpart];
   float GenPart_phi[maxpart];
   float GenPart_mass[maxpart];
   float GenPart_relIso[maxpart];
   bool  GenPart_isHdaug[maxpart];

   int   N_GenJet;
   bool  GenJet_isHdaug[maxpart];
   float GenJet_pt[maxjets];
   float GenJet_eta[maxjets];
   float GenJet_phi[maxjets];
   float GenJet_mass[maxjets];
  
   int   N_GenPh;
   int   GenPh_st[maxpart];
   bool  GenPh_isHdaug[maxpart];
   float GenPh_p[maxpart];
   float GenPh_px[maxpart];
   float GenPh_py[maxpart];
   float GenPh_pz[maxpart];
   float GenPh_E[maxpart];
   float GenPh_pt[maxpart];
   float GenPh_eta[maxpart];
   float GenPh_phi[maxpart];

   //reco level event

   int   N_Vtx;
   float Vtx_pt2[maxjets];

   int   N_LooseEl;
   int   LooseEl_ch[maxpart];
   int   LooseEl_g[maxpart];
   float LooseEl_pt[maxpart];
   float LooseEl_eta[maxpart];
   float LooseEl_phi[maxpart];
   float LooseEl_mass[maxpart];
   float LooseEl_relIso[maxpart];
   float LooseEl_sf[maxpart];

   int   N_TightEl;
   int   TightEl_ch[maxpart];
   int   TightEl_g[maxpart];
   float TightEl_pt[maxpart];
   float TightEl_eta[maxpart];
   float TightEl_phi[maxpart];
   float TightEl_mass[maxpart];
   float TightEl_relIso[maxpart];
   float TightEl_sf[maxpart];
  
   int   N_MedEl;
   int   MedEl_ch[maxpart];
   int   MedEl_g[maxpart];
   float MedEl_pt[maxpart];
   float MedEl_eta[maxpart];
   float MedEl_phi[maxpart];
   float MedEl_mass[maxpart];
   float MedEl_relIso[maxpart];
   float MedEl_sf[maxpart];

   int   N_LooseMu;
   int   LooseMu_ch[maxpart];
   int   LooseMu_g[maxpart];
   float LooseMu_pt[maxpart];
   float LooseMu_eta[maxpart];
   float LooseMu_phi[maxpart];
   float LooseMu_mass[maxpart];
   float LooseMu_relIso[maxpart];
   float LooseMu_sf[maxpart];

   int   N_TightMu;
   int   TightMu_ch[maxpart];
   int   TightMu_g[maxpart];
   float TightMu_pt[maxpart];
   float TightMu_eta[maxpart];
   float TightMu_phi[maxpart];
   float TightMu_mass[maxpart];
   float TightMu_relIso[maxpart];
   float TightMu_sf[maxpart];

   int   N_Tau;
   int   Tau_ch[maxpart];
   int   Tau_g[maxpart];
   float Tau_pt[maxpart];
   float Tau_eta[maxpart];
   float Tau_phi[maxpart];
   float Tau_mass[maxpart];
   float Tau_dm[maxpart];
   float Tau_chargedIso[maxpart];
   float Tau_sf[maxpart];

   int   N_Jet;
   int   Jet_id[maxjets];
   int   Jet_g[maxjets];
   float Jet_pt[maxjets];
   float Jet_eta[maxjets];
   float Jet_phi[maxjets];
   float Jet_mass[maxjets];
   int   Jet_mvav2[maxjets];
   int   Jet_deepcsv[maxjets];
   int   Jet_flav[maxjets];
   int   Jet_hadflav[maxjets];
   int   Jet_pid[maxjets];
   float Jet_sf[maxjets];

   int   N_Met;
   float Met_pt[maxpart];
   float Met_phi[maxpart];
   float Met_eta[maxpart];
   float Met_sf[maxpart];

   int   N_LoosePh;
   int   LoosePh_g[maxpart];
   int   LoosePh_isEB[maxpart];
   float LoosePh_pt[maxpart];
   float LoosePh_eta[maxpart];
   float LoosePh_phi[maxpart];
   float LoosePh_E[maxpart];
   float LoosePh_pt_multi[maxpart];
   float LoosePh_eta_multi[maxpart];
   float LoosePh_phi_multi[maxpart];
   float LoosePh_E_multi[maxpart];
   float LoosePh_sf[maxpart];

   int   N_TightPh;
   int   TightPh_g[maxpart];
   int   TightPh_isEB[maxpart];
   float TightPh_pt[maxpart];
   float TightPh_eta[maxpart];
   float TightPh_phi[maxpart];
   float TightPh_E[maxpart];
   float TightPh_pt_multi[maxpart];
   float TightPh_eta_multi[maxpart];
   float TightPh_phi_multi[maxpart];
   float TightPh_E_multi[maxpart];
   float TightPh_sf[maxpart];

   int   N_SelectedPh;
   int   SelectedPh_g[maxpart];
   int   SelectedPh_isEB[maxpart];
   float SelectedPh_pt[maxpart];
   float SelectedPh_eta[maxpart];
   float SelectedPh_phi[maxpart];
   float SelectedPh_E[maxpart];
   float SelectedPh_pt_multi[maxpart];
   float SelectedPh_eta_multi[maxpart];
   float SelectedPh_phi_multi[maxpart];
   float SelectedPh_E_multi[maxpart];
   float SelectedPh_sf[maxpart];

};

  
void InitTreeVars(TChain* chain, TreeVars& treeVars);
void InitOutTreeVars(TTree* tree, TreeVars& treeVars);

void InitRawTreeVars(std::map<std::string,TChain*> &chain, RawTreeVars &treeVars, std::string Loose_Tight_Photon="PhotonTight");

#endif
