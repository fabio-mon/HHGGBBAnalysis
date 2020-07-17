#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

#include "interface/CMS_lumi.h"
#include "interface/TreeUtils.h"

#include <iostream>
#include <iomanip>
#include <map>

#include "TLorentzVector.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2F.h"

#define MZ 91.187
#define PI 3.14159265359

enum DiLeptonCategories { None, DiMuon, DiElectron, Mixed };

float DeltaEta(const float& eta1, const float& eta2);

float DeltaPhi(const float& phi1, const float& phi2);

float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2);

void MakePlot3(std::map<std::string,TH1F*> &h);

bool DiPhotonSelection(const TLorentzVector &pho_lead ,const TLorentzVector &pho_sublead);
bool FindGenPh_Hdaug(RawTreeVars &treeVars, float deltaMthr=10.);
bool Findbquark_Hdaug(RawTreeVars &treeVars);
bool FindGenJet_Hdaug(RawTreeVars &treeVars, float deltaMthr=30.);
bool JetSelection(const RawTreeVars &treeVars, TreeVars &outtreeVars);//const TLorentzVector &bjet_lead, const TLorentzVector &bjet_sublead, const TLorentzVector &pho_lead, const TLorentzVector &pho_sublead)
void FindLeadSublead_pho(const RawTreeVars &treeVars, int &pho_lead_i, int &pho_sublead_i);
bool FindLeadSublead_bjet(const RawTreeVars &treeVars, int &bjet_lead_i, int &bjet_sublead_i);
bool RecoJetGenericMatch(const TLorentzVector &reco_pho , const RawTreeVars& treeVars , TLorentzVector &reco_jet_match, float DeltaRmax=0.03);
bool PhoGenericGenMatch(const TLorentzVector &reco_pho , const RawTreeVars& treeVars , TLorentzVector &gen_pho_match, float DeltaRmax=0.03);
bool bquarkGenericGenMatch(const TLorentzVector &reco_jet , const RawTreeVars& treeVars , TLorentzVector &gen_b_match, float DeltaRmax=0.03);
bool PhoGenMatch(const TLorentzVector &pho_lead , const TLorentzVector &pho_sublead , const RawTreeVars& treeVars , TreeVars &outtreeVars, float DeltaRmax=0.03);
void ReducebtagEfficiency(TreeVars &outtreeVars, const bool &useMTD, const float &prob );
void IncreasebtagEfficiency(TreeVars &outtreeVars, const bool &useMTD, const float &prob );
bool SelectBestScoreBJets2(TreeVars outtreeVars,int &bjet_lead_i,int &bjet_sublead_i, const bool &useMTD);
bool SelectBestScoreBJets3(const TreeVars &outtreeVars,int &bjet_lead_i,int &bjet_sublead_i);
float DeltaRmin_phoRECO_phoGEN(const TLorentzVector &reco_pho , const RawTreeVars& treeVars);
float DeltaRmin_phoRECO_jetRECO(const TLorentzVector &reco_pho , const RawTreeVars& treeVars);
float DeltaRmin(const vector<TLorentzVector> &coll1 , const vector<TLorentzVector> &coll2);
void PrintRecoPhoton(const RawTreeVars& treeVars);
void PrintRecoJet(const RawTreeVars& treeVars);
bool FindHHGen(const RawTreeVars &treeVars, TreeVars &outtreeVars);
void btagReweight(TreeVars &outtreeVars, float &weightMC,
		  TH2F* btaglooseeff_reweightmap, TH2F* btagmediumeff_reweightmap, TH2F* btagtighteff_reweightmap,
		  TH2F* btagloosefake_reweightmap,TH2F* btagmediumfake_reweightmap,TH2F* btagtightfake_reweightmap);

#endif
