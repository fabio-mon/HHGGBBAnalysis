#include "interface/TreeUtils.h"

void InitTreeVars(TChain* chain, TreeVars& treeVars)
{
  chain -> SetBranchAddress("weight",&treeVars.weight);
  chain -> SetBranchAddress("event",            &treeVars.event);
  chain -> SetBranchAddress("nvtx",            &treeVars.nvtx);
  chain -> SetBranchAddress("cross_sec",            &treeVars.cross_sec);

  chain -> SetBranchAddress("dipho_sumpt",   &treeVars.dipho_sumpt);
  chain -> SetBranchAddress("dipho_mass",    &treeVars.dipho_mass);
  //chain -> SetBranchAddress("dipho_mass_gen",    &treeVars.dipho_mass_gen);
  //  chain -> SetBranchAddress("dipho_sigmaRV", &treeVars.dipho_sigmaRV);
  chain -> SetBranchAddress("dipho_deltaeta",&treeVars.dipho_deltaeta);
  chain -> SetBranchAddress("dipho_deltaphi",&treeVars.dipho_deltaphi);
  
  chain -> SetBranchAddress("dipho_leadPt",                  &treeVars.dipho_leadPt);
  chain -> SetBranchAddress("dipho_leadEta",                 &treeVars.dipho_leadEta);
  chain -> SetBranchAddress("dipho_leadPhi",                 &treeVars.dipho_leadPhi);  
  //chain -> SetBranchAddress("dipho_leadEta_gen",                 &treeVars.dipho_leadEta_gen);
  //chain -> SetBranchAddress("dipho_leadPhi_gen",                 &treeVars.dipho_leadPhi_gen);
  chain -> SetBranchAddress("dipho_leadptoM",                &treeVars.dipho_leadptoM);
  chain -> SetBranchAddress("dipho_leadEnergy",              &treeVars.dipho_leadEnergy);
  //chain -> SetBranchAddress("dipho_leadEnergy_gen",          &treeVars.dipho_leadEnergy_gen);
  //chain -> SetBranchAddress("dipho_leadIsoVarRhoCorr",      &treeVars.dipho_leadIso);
  //chain -> SetBranchAddress("dipho_leadDeltaRgenreco",      &treeVars.dipho_leadDeltaRgenreco);
  //chain -> SetBranchAddress("dipho_leadDeltaEtagenreco",      &treeVars.dipho_leadDeltaEtagenreco);
  //chain -> SetBranchAddress("dipho_leadDeltaPhigenreco",      &treeVars.dipho_leadDeltaPhigenreco);
  //  chain -> SetBranchAddress("dipho_lead_sigmaEoE",&treeVars.dipho_lead_sigmaEoE);
  
  chain -> SetBranchAddress("dipho_subleadPt",               &treeVars.dipho_subleadPt);
  chain -> SetBranchAddress("dipho_subleadEta",              &treeVars.dipho_subleadEta);
  chain -> SetBranchAddress("dipho_subleadPhi",              &treeVars.dipho_subleadPhi);
  //chain -> SetBranchAddress("dipho_subleadEta_gen",              &treeVars.dipho_subleadEta_gen);
  //chain -> SetBranchAddress("dipho_subleadPhi_gen",              &treeVars.dipho_subleadPhi_gen);
  chain -> SetBranchAddress("dipho_subleadptoM",             &treeVars.dipho_subleadptoM);
  chain -> SetBranchAddress("dipho_subleadEnergy",           &treeVars.dipho_subleadEnergy);
  //chain -> SetBranchAddress("dipho_subleadEnergy_gen",       &treeVars.dipho_subleadEnergy_gen);
  //chain -> SetBranchAddress("dipho_subleadIsoVarRhoCorr",      &treeVars.dipho_subleadIso);
  //chain -> SetBranchAddress("dipho_subleadDeltaRgenreco",      &treeVars.dipho_subleadDeltaRgenreco);
  //chain -> SetBranchAddress("dipho_subleadDeltaEtagenreco",      &treeVars.dipho_subleadDeltaEtagenreco);
  //chain -> SetBranchAddress("dipho_subleadDeltaPhigenreco",      &treeVars.dipho_subleadDeltaPhigenreco);
  //  chain -> SetBranchAddress("dipho_sublead_sigmaEoE",&treeVars.dipho_sublead_sigmaEoE);
  
  chain -> SetBranchAddress("nJets",                               &treeVars.nJets);
  chain -> SetBranchAddress("nJets_bTagLoose",			&treeVars.nJets_bTagLoose);
  chain -> SetBranchAddress("nJets_bTagMedium",			&treeVars.nJets_bTagMedium);
  chain -> SetBranchAddress("nJets_bTagTight",			&treeVars.nJets_bTagTight);
  //chain -> SetBranchAddress("jet_pt",			        treeVars.jet_pt);
  //chain -> SetBranchAddress("jet_eta",			        treeVars.jet_eta);
  //chain -> SetBranchAddress("jet_phi",			        treeVars.jet_phi);
  //chain -> SetBranchAddress("jet_BTagLevel",			treeVars.jet_BTagLevel);

  chain -> SetBranchAddress("dibjet_mass",                      &treeVars.dibjet_mass);
  chain -> SetBranchAddress("dibjet_sumpt",                     &treeVars.dibjet_sumpt);
  chain -> SetBranchAddress("dibjet_deltaeta",                  &treeVars.dibjet_deltaeta);
  chain -> SetBranchAddress("dibjet_deltaphi",                  &treeVars.dibjet_deltaphi);
  
  chain -> SetBranchAddress("dibjet_leadPt",                    &treeVars.dibjet_leadPt);
  chain -> SetBranchAddress("dibjet_leadEta",                   &treeVars.dibjet_leadEta);
  chain -> SetBranchAddress("dibjet_leadPhi",                   &treeVars.dibjet_leadPhi);
  chain -> SetBranchAddress("dibjet_leadptoM",                  &treeVars.dibjet_leadptoM);
  chain -> SetBranchAddress("dibjet_leadEnergy",                &treeVars.dibjet_leadEnergy);
  chain -> SetBranchAddress("dibjet_leadbtagmedium",             &treeVars.dibjet_leadbtagmedium);
  
  chain -> SetBranchAddress("dibjet_subleadPt",                 &treeVars.dibjet_subleadPt);
  chain -> SetBranchAddress("dibjet_subleadEta",                &treeVars.dibjet_subleadEta);
  chain -> SetBranchAddress("dibjet_subleadPhi",                &treeVars.dibjet_subleadPhi);
  chain -> SetBranchAddress("dibjet_subleadptoM",               &treeVars.dibjet_subleadptoM);
  chain -> SetBranchAddress("dibjet_subleadEnergy",             &treeVars.dibjet_subleadEnergy);
  chain -> SetBranchAddress("dibjet_subleadbtagmedium",          &treeVars.dibjet_subleadbtagmedium);

  chain -> SetBranchAddress("Mx",                               &treeVars.Mx);
  chain -> SetBranchAddress("DRmin_pho_bjet",                   &treeVars.DRmin_pho_bjet); 
  chain -> SetBranchAddress("costheta_HH",                      &treeVars.costheta_HH); 
  chain -> SetBranchAddress("costheta_gg",                      &treeVars.costheta_gg); 
  chain -> SetBranchAddress("costheta_bb",                      &treeVars.costheta_bb); 
  chain -> SetBranchAddress("MetPt",                            &treeVars.MetPt); 
  chain -> SetBranchAddress("MetPhi",                           &treeVars.MetPhi); 

}



void InitOutTreeVars(TTree* tree, TreeVars& treeVars)
{
  tree -> Branch("evWeight",  &treeVars.weight);
  tree -> Branch("event",     &treeVars.event);
  tree -> Branch("nvtx",      &treeVars.nvtx);
  tree -> Branch("cross_sec", &treeVars.cross_sec);

  tree -> Branch("dipho_sumpt",   &treeVars.dipho_sumpt);
  tree -> Branch("mgg",           &treeVars.dipho_mass);
  //tree -> Branch("dipho_mass_gen",    &treeVars.dipho_mass_gen);
  tree -> Branch("dipho_deltaeta",&treeVars.dipho_deltaeta);
  tree -> Branch("dipho_deltaphi",&treeVars.dipho_deltaphi);
  
  tree -> Branch("dipho_leadPt",                  &treeVars.dipho_leadPt);
  tree -> Branch("dipho_leadEta",                 &treeVars.dipho_leadEta);
  tree -> Branch("dipho_leadPhi",                 &treeVars.dipho_leadPhi);
  //tree -> Branch("dipho_leadEta_gen",             &treeVars.dipho_leadEta_gen);
  //tree -> Branch("dipho_leadPhi_gen",             &treeVars.dipho_leadPhi_gen);
  tree -> Branch("dipho_leadptoM",                &treeVars.dipho_leadptoM);
  tree -> Branch("dipho_leadEnergy",              &treeVars.dipho_leadEnergy);
  //tree -> Branch("dipho_leadEnergy_gen",          &treeVars.dipho_leadEnergy_gen);
  //tree -> Branch("dipho_leadIsoVarRhoCorr",      &treeVars.dipho_leadIso);
  //tree -> Branch("dipho_leadDeltaRgenreco",      &treeVars.dipho_leadDeltaRgenreco);
  //tree -> Branch("dipho_leadDeltaEtagenreco",      &treeVars.dipho_leadDeltaEtagenreco);
  //tree -> Branch("dipho_leadDeltaPhigenreco",      &treeVars.dipho_leadDeltaPhigenreco);
  //  tree -> Branch("dipho_lead_sigmaEoE",&treeVars.dipho_lead_sigmaEoE);
  
  tree -> Branch("dipho_subleadPt",               &treeVars.dipho_subleadPt);
  tree -> Branch("dipho_subleadEta",              &treeVars.dipho_subleadEta);
  tree -> Branch("dipho_subleadPhi",              &treeVars.dipho_subleadPhi);
  //tree -> Branch("dipho_subleadEta_gen",             &treeVars.dipho_subleadEta_gen);
  //tree -> Branch("dipho_subleadPhi_gen",             &treeVars.dipho_subleadPhi_gen);
  tree -> Branch("dipho_subleadptoM",             &treeVars.dipho_subleadptoM);
  tree -> Branch("dipho_subleadEnergy",           &treeVars.dipho_subleadEnergy);
  //tree -> Branch("dipho_subleadEnergy_gen",       &treeVars.dipho_subleadEnergy_gen);
  //tree -> Branch("dipho_subleadIsoVarRhoCorr",      &treeVars.dipho_subleadIso);
  //tree -> Branch("dipho_subleadDeltaRgenreco",      &treeVars.dipho_subleadDeltaRgenreco);
  //tree -> Branch("dipho_subleadDeltaEtagenreco",      &treeVars.dipho_subleadDeltaEtagenreco);
  //tree -> Branch("dipho_subleadDeltaPhigenreco",      &treeVars.dipho_subleadDeltaPhigenreco);
  //  tree -> Branch("dipho_sublead_sigmaEoE",&treeVars.dipho_sublead_sigmaEoE);
  
  tree -> Branch("nJets",                               &treeVars.nJets);
  tree -> Branch("nJets_bTagLoose",			&treeVars.nJets_bTagLoose);
  tree -> Branch("nJets_bTagMedium",			&treeVars.nJets_bTagMedium);
  tree -> Branch("nJets_bTagTight",			&treeVars.nJets_bTagTight);
  //tree -> Branch("jet_pt",			        treeVars.jet_pt);
  //tree -> Branch("jet_eta",			        treeVars.jet_eta);
  //tree -> Branch("jet_phi",			        treeVars.jet_phi);
  //tree -> Branch("jet_BTagLevel",			treeVars.jet_BTagLevel);

  tree -> Branch("mjj",                              &treeVars.dibjet_mass);
  tree -> Branch("dibjet_sumpt",                     &treeVars.dibjet_sumpt);
  tree -> Branch("dibjet_deltaeta",                  &treeVars.dibjet_deltaeta);
  tree -> Branch("dibjet_deltaphi",                  &treeVars.dibjet_deltaphi);
  
  tree -> Branch("dibjet_leadPt",                    &treeVars.dibjet_leadPt);
  tree -> Branch("dibjet_leadEta",                   &treeVars.dibjet_leadEta);
  tree -> Branch("dibjet_leadPhi",                   &treeVars.dibjet_leadPhi);
  tree -> Branch("dibjet_leadptoM",                  &treeVars.dibjet_leadptoM);
  tree -> Branch("dibjet_leadEnergy",                &treeVars.dibjet_leadEnergy);
  tree -> Branch("dibjet_leadbtagmedium",             &treeVars.dibjet_leadbtagmedium);
  
  tree -> Branch("dibjet_subleadPt",                 &treeVars.dibjet_subleadPt);
  tree -> Branch("dibjet_subleadEta",                &treeVars.dibjet_subleadEta);
  tree -> Branch("dibjet_subleadPhi",                &treeVars.dibjet_subleadPhi);
  tree -> Branch("dibjet_subleadptoM",               &treeVars.dibjet_subleadptoM);
  tree -> Branch("dibjet_subleadEnergy",             &treeVars.dibjet_subleadEnergy);
  tree -> Branch("dibjet_subleadbtagmedium",          &treeVars.dibjet_subleadbtagmedium);
  
  tree -> Branch("mtot",                               &treeVars.Mx);
  tree -> Branch("DRmin_pho_bjet",                   &treeVars.DRmin_pho_bjet); 
  tree -> Branch("costheta_HH",                      &treeVars.costheta_HH); 
  tree -> Branch("costheta_gg",                      &treeVars.costheta_gg); 
  tree -> Branch("costheta_bb",                      &treeVars.costheta_bb); 
  tree -> Branch("MetPt",                            &treeVars.MetPt); 
  tree -> Branch("MetPhi",                           &treeVars.MetPhi); 
  
  tree -> Branch("cut_based_ct", &treeVars.cut_based_ct); 
  tree -> Branch("ttHTagger",    &treeVars.ttHTagger); 
}


void InitRawTreeVars(std::map<std::string,TChain*> &chain, RawTreeVars& treeVars, std::string Loose_Tight_Photon)
{
  // can be "PhotonLoose" or "PhotonTight"
  
  //chain->SetMakeClass(1);
  //for(std::map<std::string,TChain*>::iterator it=chain.begin(); it!=chain.end(); ++it)
  //   it->second->SetBranchStatus("*", 0);

  //event header
  chain["Event"]->SetBranchAddress("Run",               &treeVars.run);
  chain["Event"]->SetBranchAddress("Event",             &treeVars.event);
  chain["Event"]->SetBranchAddress("Lumi",              &treeVars.lumi);

  //gen level event
  chain["Event"]->SetBranchAddress("Weight_size",       &treeVars.N_Weight);
  chain["Event"]->SetBranchAddress("Weight",            treeVars.Weight);

  chain["Particle"]->SetBranchAddress("Particle_size",  &treeVars.N_GenPart);
  chain["Particle"]->SetBranchAddress("PID",            treeVars.GenPart_pid);
  chain["Particle"]->SetBranchAddress("Charge",         treeVars.GenPart_ch);
  chain["Particle"]->SetBranchAddress("Status",         treeVars.GenPart_st);
  chain["Particle"]->SetBranchAddress("P",              treeVars.GenPart_p);
  chain["Particle"]->SetBranchAddress("Px",             treeVars.GenPart_px);
  chain["Particle"]->SetBranchAddress("Py",             treeVars.GenPart_py);
  chain["Particle"]->SetBranchAddress("Pz",             treeVars.GenPart_pz);
  chain["Particle"]->SetBranchAddress("E",              treeVars.GenPart_E);
  chain["Particle"]->SetBranchAddress("PT",             treeVars.GenPart_pt);
  chain["Particle"]->SetBranchAddress("Eta",            treeVars.GenPart_eta);
  chain["Particle"]->SetBranchAddress("Phi",            treeVars.GenPart_phi);
  chain["Particle"]->SetBranchAddress("Mass",           treeVars.GenPart_mass);
  chain["Particle"]->SetBranchAddress("IsolationVar",   treeVars.GenPart_relIso);

  chain["GenJet"]->SetBranchAddress("GenJet_size",     &treeVars.N_GenJet);
  chain["GenJet"]->SetBranchAddress("PT",              treeVars.GenJet_pt);
  chain["GenJet"]->SetBranchAddress("Eta",             treeVars.GenJet_eta);
  chain["GenJet"]->SetBranchAddress("Phi",             treeVars.GenJet_phi);
  chain["GenJet"]->SetBranchAddress("Mass",            treeVars.GenJet_mass);

  chain["GenPhoton"]->SetBranchAddress("GenPhoton_size",  &treeVars.N_GenPh);
  chain["GenPhoton"]->SetBranchAddress("Status",          treeVars.GenPh_st);
  chain["GenPhoton"]->SetBranchAddress("P",               treeVars.GenPh_p);
  chain["GenPhoton"]->SetBranchAddress("Px",              treeVars.GenPh_px);
  chain["GenPhoton"]->SetBranchAddress("Py",              treeVars.GenPh_py);
  chain["GenPhoton"]->SetBranchAddress("Pz",              treeVars.GenPh_pz);
  chain["GenPhoton"]->SetBranchAddress("E",               treeVars.GenPh_E);
  chain["GenPhoton"]->SetBranchAddress("PT",              treeVars.GenPh_pt);
  chain["GenPhoton"]->SetBranchAddress("Eta",             treeVars.GenPh_eta);
  chain["GenPhoton"]->SetBranchAddress("Phi",             treeVars.GenPh_phi);
  //chain["GenPhoton"]->SetBranchAddress("isHdaug",         treeVars.GenPh_isHdaug); 

  //reco level event
  chain["Vertex"]->SetBranchAddress("Vertex_size",    &treeVars.N_Vtx);
  chain["Vertex"]->SetBranchAddress("SumPT2",         &treeVars.Vtx_pt2);

  chain["ElectronLoose"]->SetBranchAddress("ElectronLoose_size", &treeVars.N_LooseEl);
  chain["ElectronLoose"]->SetBranchAddress("Charge",       treeVars.LooseEl_ch);
  chain["ElectronLoose"]->SetBranchAddress("Particle",     treeVars.LooseEl_g);
  chain["ElectronLoose"]->SetBranchAddress("PT",           treeVars.LooseEl_pt);
  chain["ElectronLoose"]->SetBranchAddress("Eta",          treeVars.LooseEl_eta);
  chain["ElectronLoose"]->SetBranchAddress("Phi",          treeVars.LooseEl_phi);
  chain["ElectronLoose"]->SetBranchAddress("Mass",         treeVars.LooseEl_mass);
  chain["ElectronLoose"]->SetBranchAddress("IsolationVar", treeVars.LooseEl_relIso);
  chain["ElectronLoose"]->SetBranchAddress("SF", treeVars.LooseEl_sf);

  chain["ElectronTight"]->SetBranchAddress("ElectronTight_size", &treeVars.N_TightEl);
  chain["ElectronTight"]->SetBranchAddress("Charge",       treeVars.TightEl_ch);
  chain["ElectronTight"]->SetBranchAddress("Particle",     treeVars.TightEl_g);
  chain["ElectronTight"]->SetBranchAddress("PT",           treeVars.TightEl_pt);
  chain["ElectronTight"]->SetBranchAddress("Eta",          treeVars.TightEl_eta);
  chain["ElectronTight"]->SetBranchAddress("Phi",          treeVars.TightEl_phi);
  chain["ElectronTight"]->SetBranchAddress("Mass",         treeVars.TightEl_mass);
  chain["ElectronTight"]->SetBranchAddress("IsolationVar", treeVars.TightEl_relIso);
  chain["ElectronTight"]->SetBranchAddress("SF",           treeVars.TightEl_sf);
  
  chain["ElectronMedium"]->SetBranchAddress("ElectronMedium_size", &treeVars.N_MedEl);
  chain["ElectronMedium"]->SetBranchAddress("Charge",               treeVars.MedEl_ch);
  chain["ElectronMedium"]->SetBranchAddress("Particle",             treeVars.MedEl_g);
  chain["ElectronMedium"]->SetBranchAddress("PT",                   treeVars.MedEl_pt);
  chain["ElectronMedium"]->SetBranchAddress("Eta",                  treeVars.MedEl_eta);
  chain["ElectronMedium"]->SetBranchAddress("Phi",                  treeVars.MedEl_phi);
  chain["ElectronMedium"]->SetBranchAddress("Mass",                 treeVars.MedEl_mass);
  chain["ElectronMedium"]->SetBranchAddress("IsolationVar",         treeVars.MedEl_relIso);
  chain["ElectronMedium"]->SetBranchAddress("SF",           treeVars.MedEl_sf);

  chain["MuonLoose"]->SetBranchAddress("MuonLoose_size", &treeVars.N_LooseMu);
  chain["MuonLoose"]->SetBranchAddress("Charge",       treeVars.LooseMu_ch);
  chain["MuonLoose"]->SetBranchAddress("Particle",     treeVars.LooseMu_g);
  chain["MuonLoose"]->SetBranchAddress("PT",           treeVars.LooseMu_pt);
  chain["MuonLoose"]->SetBranchAddress("Eta",          treeVars.LooseMu_eta);
  chain["MuonLoose"]->SetBranchAddress("Phi",          treeVars.LooseMu_phi);
  chain["MuonLoose"]->SetBranchAddress("Mass",         treeVars.LooseMu_mass);
  chain["MuonLoose"]->SetBranchAddress("IsolationVar", treeVars.LooseMu_relIso);
  chain["MuonLoose"]->SetBranchAddress("SF",           treeVars.LooseMu_sf);

  chain["MuonTight"]->SetBranchAddress("MuonTight_size", &treeVars.N_TightMu);
  chain["MuonTight"]->SetBranchAddress("Charge",       treeVars.TightMu_ch);
  chain["MuonTight"]->SetBranchAddress("Particle",     treeVars.TightMu_g);
  chain["MuonTight"]->SetBranchAddress("PT",           treeVars.TightMu_pt);
  chain["MuonTight"]->SetBranchAddress("Eta",          treeVars.TightMu_eta);
  chain["MuonTight"]->SetBranchAddress("Phi",          treeVars.TightMu_phi);
  chain["MuonTight"]->SetBranchAddress("Mass",         treeVars.TightMu_mass);
  chain["MuonTight"]->SetBranchAddress("IsolationVar", treeVars.TightMu_relIso);
  chain["MuonTight"]->SetBranchAddress("SF",           treeVars.TightMu_sf);

  chain["TauAll"]->SetBranchAddress("TauAll_size", &treeVars.N_Tau);
  chain["TauAll"]->SetBranchAddress("Charge",       treeVars.Tau_ch);
  chain["TauAll"]->SetBranchAddress("Particle",     treeVars.Tau_g);
  chain["TauAll"]->SetBranchAddress("PT",           treeVars.Tau_pt);
  chain["TauAll"]->SetBranchAddress("Eta",          treeVars.Tau_eta);
  chain["TauAll"]->SetBranchAddress("Phi",          treeVars.Tau_phi);
  chain["TauAll"]->SetBranchAddress("Mass",         treeVars.Tau_mass);
  chain["TauAll"]->SetBranchAddress("DM",         treeVars.Tau_dm);
  chain["TauAll"]->SetBranchAddress("IsolationVar", treeVars.Tau_chargedIso);
  chain["TauAll"]->SetBranchAddress("SF",           treeVars.Tau_sf);

  chain["JetPUPPI"]->SetBranchAddress("JetPUPPI_size",             &treeVars.N_Jet);
  chain["JetPUPPI"]->SetBranchAddress("ID",                        treeVars.Jet_id);
  chain["JetPUPPI"]->SetBranchAddress("GenJet",                    treeVars.Jet_g);
  chain["JetPUPPI"]->SetBranchAddress("PT",                        treeVars.Jet_pt);
  chain["JetPUPPI"]->SetBranchAddress("Eta",                       treeVars.Jet_eta);
  chain["JetPUPPI"]->SetBranchAddress("Phi",                       treeVars.Jet_phi);
  chain["JetPUPPI"]->SetBranchAddress("Mass",                      treeVars.Jet_mass);
  chain["JetPUPPI"]->SetBranchAddress("MVAv2",                     treeVars.Jet_mvav2);
  chain["JetPUPPI"]->SetBranchAddress("DeepCSV",                   treeVars.Jet_deepcsv);
  chain["JetPUPPI"]->SetBranchAddress("PartonFlavor",              treeVars.Jet_flav);
  chain["JetPUPPI"]->SetBranchAddress("HadronFlavor",              treeVars.Jet_hadflav);
  chain["JetPUPPI"]->SetBranchAddress("GenPartonPID",              treeVars.Jet_pid);
  chain["JetPUPPI"]->SetBranchAddress("SF",                        treeVars.Jet_sf);

  chain["PuppiMissingET"]->SetBranchAddress("PuppiMissingET_size", &treeVars.N_Met);
  chain["PuppiMissingET"]->SetBranchAddress("MET",                 treeVars.Met_pt);
  chain["PuppiMissingET"]->SetBranchAddress("Phi",                 treeVars.Met_phi);
  chain["PuppiMissingET"]->SetBranchAddress("Eta",                 treeVars.Met_eta);
  chain["PuppiMissingET"]->SetBranchAddress("SF",                  treeVars.Met_sf);

  if(Loose_Tight_Photon == "PhotonTight")
  {
    chain["PhotonLoose"]->SetBranchAddress("PhotonLoose_size",       &treeVars.N_LoosePh);
    chain["PhotonLoose"]->SetBranchAddress("Particle",               treeVars.LoosePh_g);
    chain["PhotonLoose"]->SetBranchAddress("IsEB",                   treeVars.LoosePh_isEB);
    chain["PhotonLoose"]->SetBranchAddress("PT",                     treeVars.LoosePh_pt);
    chain["PhotonLoose"]->SetBranchAddress("Eta",                    treeVars.LoosePh_eta);
    chain["PhotonLoose"]->SetBranchAddress("Phi",                    treeVars.LoosePh_phi);
    chain["PhotonLoose"]->SetBranchAddress("E",                      treeVars.LoosePh_E);
    chain["PhotonLoose"]->SetBranchAddress("PT_multi",               treeVars.LoosePh_pt_multi);
    chain["PhotonLoose"]->SetBranchAddress("Eta_multi",              treeVars.LoosePh_eta_multi);
    chain["PhotonLoose"]->SetBranchAddress("Phi_multi",              treeVars.LoosePh_phi_multi);
    chain["PhotonLoose"]->SetBranchAddress("E_multi",                treeVars.LoosePh_E_multi);
    chain["PhotonLoose"]->SetBranchAddress("SF",                     treeVars.LoosePh_sf);
    
    chain["PhotonTight"]->SetBranchAddress("PhotonTight_size",       &treeVars.N_SelectedPh);
    chain["PhotonTight"]->SetBranchAddress("Particle",               treeVars.SelectedPh_g);
    chain["PhotonTight"]->SetBranchAddress("IsEB",                   treeVars.SelectedPh_isEB);
    chain["PhotonTight"]->SetBranchAddress("PT",                     treeVars.SelectedPh_pt);
    chain["PhotonTight"]->SetBranchAddress("Eta",                    treeVars.SelectedPh_eta);
    chain["PhotonTight"]->SetBranchAddress("Phi",                    treeVars.SelectedPh_phi);
    chain["PhotonTight"]->SetBranchAddress("E",                      treeVars.SelectedPh_E);
    chain["PhotonTight"]->SetBranchAddress("PT_multi",               treeVars.SelectedPh_pt_multi);
    chain["PhotonTight"]->SetBranchAddress("Eta_multi",              treeVars.SelectedPh_eta_multi);
    chain["PhotonTight"]->SetBranchAddress("Phi_multi",              treeVars.SelectedPh_phi_multi);
    chain["PhotonTight"]->SetBranchAddress("E_multi",                treeVars.SelectedPh_E_multi);
    chain["PhotonTight"]->SetBranchAddress("SF",                     treeVars.SelectedPh_sf);    
  }

  else
  {
    chain["PhotonLoose"]->SetBranchAddress("PhotonLoose_size",       &treeVars.N_SelectedPh);
    chain["PhotonLoose"]->SetBranchAddress("Particle",               treeVars.SelectedPh_g);
    chain["PhotonLoose"]->SetBranchAddress("IsEB",                   treeVars.SelectedPh_isEB);
    chain["PhotonLoose"]->SetBranchAddress("PT",                     treeVars.SelectedPh_pt);
    chain["PhotonLoose"]->SetBranchAddress("Eta",                    treeVars.SelectedPh_eta);
    chain["PhotonLoose"]->SetBranchAddress("Phi",                    treeVars.SelectedPh_phi);
    chain["PhotonLoose"]->SetBranchAddress("E",                      treeVars.SelectedPh_E);
    chain["PhotonLoose"]->SetBranchAddress("PT_multi",               treeVars.SelectedPh_pt_multi);
    chain["PhotonLoose"]->SetBranchAddress("Eta_multi",              treeVars.SelectedPh_eta_multi);
    chain["PhotonLoose"]->SetBranchAddress("Phi_multi",              treeVars.SelectedPh_phi_multi);
    chain["PhotonLoose"]->SetBranchAddress("E_multi",                treeVars.SelectedPh_E_multi);
    chain["PhotonLoose"]->SetBranchAddress("SF",                     treeVars.SelectedPh_sf);
    
    chain["PhotonTight"]->SetBranchAddress("PhotonTight_size",       &treeVars.N_TightPh);
    chain["PhotonTight"]->SetBranchAddress("Particle",               treeVars.TightPh_g);
    chain["PhotonTight"]->SetBranchAddress("IsEB",                   treeVars.TightPh_isEB);
    chain["PhotonTight"]->SetBranchAddress("PT",                     treeVars.TightPh_pt);
    chain["PhotonTight"]->SetBranchAddress("Eta",                    treeVars.TightPh_eta);
    chain["PhotonTight"]->SetBranchAddress("Phi",                    treeVars.TightPh_phi);
    chain["PhotonTight"]->SetBranchAddress("E",                      treeVars.TightPh_E);
    chain["PhotonTight"]->SetBranchAddress("PT_multi",               treeVars.TightPh_pt_multi);
    chain["PhotonTight"]->SetBranchAddress("Eta_multi",              treeVars.TightPh_eta_multi);
    chain["PhotonTight"]->SetBranchAddress("Phi_multi",              treeVars.TightPh_phi_multi);
    chain["PhotonTight"]->SetBranchAddress("E_multi",                treeVars.TightPh_E_multi);
    chain["PhotonTight"]->SetBranchAddress("SF",                     treeVars.TightPh_sf);    
  }

  /*
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress((Loose_Tight_Photon+"_size").c_str(),       &treeVars.N_TightPh);
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress("Particle",               treeVars.TightPh_g);
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress("IsEB",                   treeVars.TightPh_isEB);
  chain[Loose_Tight_Photon.c_str()]->SetBranchStatus("PT", 1);
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress("PT",                     treeVars.TightPh_pt);
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress("Eta",                    treeVars.TightPh_eta);
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress("Phi",                    treeVars.TightPh_phi);
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress("E",                      treeVars.TightPh_E);
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress("PT_multi",               treeVars.TightPh_pt_multi);
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress("Eta_multi",              treeVars.TightPh_eta_multi);
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress("Phi_multi",              treeVars.TightPh_phi_multi);
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress("E_multi",                treeVars.TightPh_E_multi);
  chain[Loose_Tight_Photon.c_str()]->SetBranchAddress("SF",                     treeVars.TightPh_sf);
  */

}
