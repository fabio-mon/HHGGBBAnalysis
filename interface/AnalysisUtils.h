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

void MakePlot3(std::map<std::string,TH1F*> &h);

bool DiPhotonSelection(const TLorentzVector &pho_lead ,const TLorentzVector &pho_sublead);
bool FindGenPh_Hdaug(RawTreeVars &treeVars, float deltaMthr=10.);
bool JetSelection(const RawTreeVars &treeVars, TreeVars &outtreeVars);
void  FindLeadSublead_pho(const RawTreeVars &treeVars, int &pho_lead_i, int &pho_sublead_i);
bool PhoGenMatch(const TLorentzVector &pho_lead , const TLorentzVector &pho_sublead , const RawTreeVars& treeVars , TreeVars &outtreeVars, float DeltaRmax=0.03);
float DeltaRmin(const TLorentzVector &reco_pho , const RawTreeVars& treeVars);

#endif
