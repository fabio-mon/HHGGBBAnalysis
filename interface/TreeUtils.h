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
  
  float dipho_sumpt;
  float dipho_mass;
  float dipho_vtxProb;
  float dipho_sigmaRV;
  float dipho_sigmaWV;
  float dipho_deltaphi;
  float dipho_cosDeltaphi;
  float dipho_leadPt;
  float dipho_leadEta;
  float dipho_leadPhi;
  float dipho_leadEnergy;
  float dipho_leadR9;
  float dipho_lead_ptoM;
  float dipho_lead_sigmaEoE;
  float dipho_leadIDMVA;
  float dipho_subleadPt;
  float dipho_subleadEta;
  float dipho_subleadPhi;
  float dipho_subleadEnergy;
  float dipho_subleadR9;
  float dipho_sublead_ptoM;
  float dipho_sublead_sigmaEoE;
  float dipho_subleadIDMVA;
  float dipho_mva;
  
  float MetPt;
  float MetPhi;
  float ttHMVA;
  
  float nJets;
  float nJets_bTagLoose;
  float nJets_bTagMedium;
  float nJets_bTagTight;
  
  float jet_pt[9];
  float jet_eta[9];
  float jet_phi[9];
  float jet_bdiscriminant[9];
  
  float mu_pt[2];
  float mu_eta[2];
  float mu_phi[2];
  float mu_energy[2];
  float mu_isLoose[2];
  float mu_isMedium[2];
  float mu_miniIso[2];
  float mu_trackIso[2];
  float mu_charge[2];
  float mu_sumChargedHadronPt[2];
  float mu_sumNeutralHadronEt[2];
  float mu_sumPhotonEt[2];
  float mu_sumPUPt[2];
  float ele_pt[2];
  float ele_eta[2];
  float ele_phi[2];
  float ele_energy[2];
  float ele_passVetoId[2];
  float ele_passLooseId[2];
  float ele_passMediumId[2];
  float ele_passTightId[2];
  float ele_MVAMediumId[2];
  float ele_MVATightId[2];
  float ele_miniIso[2];
  float ele_ecalEnergy[2];
  float ele_SCx[2];
  float ele_SCy[2];
  float ele_SCz[2];
  float ele_charge[2];
  float ele_SCeta[2];
  float ele_SCphi[2];
  float ele_dEtaTrk[2];
  float ele_dPhiTrk[2];
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
   float GenPart_pid[maxpart];
   float GenPart_ch[maxpart];
   float GenPart_st[maxpart];
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

   int   N_GenJet;
   float GenJet_pt[maxjets];
   float GenJet_eta[maxjets];
   float GenJet_phi[maxjets];
   float GenJet_mass[maxjets];

   int   N_GenPh;
   float GenPh_st[maxpart];
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
   float LooseEl_ch[maxpart];
   float LooseEl_g[maxpart];
   float LooseEl_pt[maxpart];
   float LooseEl_eta[maxpart];
   float LooseEl_phi[maxpart];
   float LooseEl_mass[maxpart];
   float LooseEl_relIso[maxpart];
   float LooseEl_sf[maxpart];

   int   N_TightEl;
   float TightEl_ch[maxpart];
   float TightEl_g[maxpart];
   float TightEl_pt[maxpart];
   float TightEl_eta[maxpart];
   float TightEl_phi[maxpart];
   float TightEl_mass[maxpart];
   float TightEl_relIso[maxpart];
   float TightEl_sf[maxpart];
  
   int   N_MedEl;
   float MedEl_ch[maxpart];
   float MedEl_g[maxpart];
   float MedEl_pt[maxpart];
   float MedEl_eta[maxpart];
   float MedEl_phi[maxpart];
   float MedEl_mass[maxpart];
   float MedEl_relIso[maxpart];
   float MedEl_sf[maxpart];

   int   N_LooseMu;
   float LooseMu_ch[maxpart];
   float LooseMu_g[maxpart];
   float LooseMu_pt[maxpart];
   float LooseMu_eta[maxpart];
   float LooseMu_phi[maxpart];
   float LooseMu_mass[maxpart];
   float LooseMu_relIso[maxpart];
   float LooseMu_sf[maxpart];

   int   N_TightMu;
   float TightMu_ch[maxpart];
   float TightMu_g[maxpart];
   float TightMu_pt[maxpart];
   float TightMu_eta[maxpart];
   float TightMu_phi[maxpart];
   float TightMu_mass[maxpart];
   float TightMu_relIso[maxpart];
   float TightMu_sf[maxpart];

   int   N_Tau;
   float Tau_ch[maxpart];
   float Tau_g[maxpart];
   float Tau_pt[maxpart];
   float Tau_eta[maxpart];
   float Tau_phi[maxpart];
   float Tau_mass[maxpart];
   float Tau_dm[maxpart];
   float Tau_chargedIso[maxpart];
   float Tau_sf[maxpart];

   int   N_Jet;
   float Jet_id[maxjets];
   float Jet_g[maxjets];
   float Jet_pt[maxjets];
   float Jet_eta[maxjets];
   float Jet_phi[maxjets];
   float Jet_mass[maxjets];
   float Jet_mvav2[maxjets];
   float Jet_deepcsv[maxjets];
   float Jet_flav[maxjets];
   float Jet_hadflav[maxjets];
   float Jet_pid[maxjets];
   float Jet_sf[maxjets];

   int   N_Met;
   float Met_pt[maxpart];
   float Met_phi[maxpart];
   float Met_eta[maxpart];
   float Met_sf[maxpart];

   int   N_LoosePh;
   float LoosePh_g[maxpart];
   float LoosePh_isEB[maxpart];
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
   float TightPh_g[maxpart];
   float TightPh_isEB[maxpart];
   float TightPh_pt[maxpart];
   float TightPh_eta[maxpart];
   float TightPh_phi[maxpart];
   float TightPh_E[maxpart];
   float TightPh_pt_multi[maxpart];
   float TightPh_eta_multi[maxpart];
   float TightPh_phi_multi[maxpart];
   float TightPh_E_multi[maxpart];
   float TightPh_sf[maxpart];

};

  
void InitTreeVars(TChain* chain, TreeVars& treeVars);
void InitOutTreeVars(TTree* tree, TreeVars& treeVars);

void InitRawTreeVars(std::map<std::string,TChain*> &chain, RawTreeVars &treeVars);

#endif
