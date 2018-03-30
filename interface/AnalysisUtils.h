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

#define MZ 91.187
#define PI 3.14159265359

#define bDiscriminantThresholdLoose  0.5426
#define bDiscriminantThresholdMedium 0.8484
#define bDiscriminantThresholdTight  0.9535

enum DiLeptonCategories { None, DiMuon, DiElectron, Mixed };

float DeltaEta(const float& eta1, const float& eta2);

float DeltaPhi(const float& phi1, const float& phi2);

float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2);


bool DiMuonSelections(TLorentzVector mu1, TLorentzVector mu2, float charge1, float charge2, TLorentzVector ph1, TLorentzVector ph2, float Iso1, float Iso2);
bool DiEleSelections(TLorentzVector ele1, TLorentzVector ele2, float charge1, float charge2, TLorentzVector ph1, TLorentzVector ph2, float miniIso1, float miniIso2, float dTrk1, float dTrk2);
bool MixedSelections(TLorentzVector mu, TLorentzVector ele, float charge1, float charge2, TLorentzVector ph1, TLorentzVector ph2);
bool SingleMuSelections(TLorentzVector mu1,  TLorentzVector ph1, TLorentzVector ph2, float miniIso);
bool SingleEleSelections(TLorentzVector ele1, TLorentzVector ph1, TLorentzVector ph2, float miniIso, float drTrk);
bool SingleMuSelectionsStandard(TLorentzVector mu1, TLorentzVector ph1, TLorentzVector ph2, float miniIso);
bool SingleEleSelectionsStandard(TLorentzVector ele1, TLorentzVector ph1, TLorentzVector ph2, float miniIso, float drTrk);

bool OneCategorySelection(const TreeVars& treeVars, const int& type);
bool DiLeptonSelection(const TreeVars& treeVars, const int& type, DiLeptonCategories& cat);
bool SingleLeptonSelection(const TreeVars& treeVars, const int& type);

bool CutBasedSelection(const TreeVars& treeVars,
                       const float& min_lead_ptoM, const float& min_sublead_ptoM,
                       const float& min_leadIDMVA, const float& min_subleadIDMVA,
                       const float& max_deltaphi, const float& max_deltaeta);

void MakePlot(TH1F**, TString title);

void MakePlot2(std::map<std::string,TH1F*>& histos, TString title);

#endif
