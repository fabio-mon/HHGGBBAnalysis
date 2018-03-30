#ifndef TREE_UTILS_H
#define TREE_UTILS_H

#include <iostream>
#include <vector>
#include <map>

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
  
void InitTreeVars(TChain* chain, TreeVars& treeVars);
void InitOutTreeVars(TTree* tree, TreeVars& treeVars);

#endif
