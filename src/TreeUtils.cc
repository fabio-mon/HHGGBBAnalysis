#include "interface/TreeUtils.h"



void InitTreeVars(TChain* chain, TreeVars& treeVars)
{
  chain -> SetBranchAddress("nvtx", &treeVars.nvtx);
  chain -> SetBranchAddress("weight", &treeVars.weight);
  
  chain -> SetBranchAddress("dipho_sumpt", &treeVars.dipho_sumpt);
  chain -> SetBranchAddress("dipho_mass", &treeVars.dipho_mass);
  chain -> SetBranchAddress("dipho_vtxProb", &treeVars.dipho_vtxProb);
  chain -> SetBranchAddress("dipho_sigmaRV", &treeVars.dipho_sigmaRV);
  chain -> SetBranchAddress("dipho_sigmaWV", &treeVars.dipho_sigmaWV);
  chain -> SetBranchAddress("dipho_deltaphi", &treeVars.dipho_deltaphi);
  chain -> SetBranchAddress("dipho_cosDeltaphi", &treeVars.dipho_cosDeltaphi);
  chain -> SetBranchAddress("dipho_leadPt", &treeVars.dipho_leadPt);
  chain -> SetBranchAddress("dipho_leadEta", &treeVars.dipho_leadEta);
  chain -> SetBranchAddress("dipho_leadPhi", &treeVars.dipho_leadPhi);
  chain -> SetBranchAddress("dipho_leadEnergy", &treeVars.dipho_leadEnergy);
  chain -> SetBranchAddress("dipho_leadR9", &treeVars.dipho_leadR9);
  chain -> SetBranchAddress("dipho_lead_ptoM", &treeVars.dipho_lead_ptoM);
  chain -> SetBranchAddress("dipho_lead_sigmaEoE", &treeVars.dipho_lead_sigmaEoE);
  chain -> SetBranchAddress("dipho_leadIDMVA", &treeVars.dipho_leadIDMVA);
  chain -> SetBranchAddress("dipho_subleadPt", &treeVars.dipho_subleadPt);
  chain -> SetBranchAddress("dipho_subleadEta", &treeVars.dipho_subleadEta);
  chain -> SetBranchAddress("dipho_subleadPhi", &treeVars.dipho_subleadPhi);
  chain -> SetBranchAddress("dipho_subleadEnergy", &treeVars.dipho_subleadEnergy);
  chain -> SetBranchAddress("dipho_subleadR9", &treeVars.dipho_subleadR9);
  chain -> SetBranchAddress("dipho_sublead_ptoM", &treeVars.dipho_sublead_ptoM);
  chain -> SetBranchAddress("dipho_sublead_sigmaEoE", &treeVars.dipho_sublead_sigmaEoE);
  chain -> SetBranchAddress("dipho_subleadIDMVA", &treeVars.dipho_subleadIDMVA);
  chain -> SetBranchAddress("dipho_mva", &treeVars.dipho_mva);
  
  chain -> SetBranchAddress("nJets", &treeVars.nJets);
  chain -> SetBranchAddress("nJets_bTagLoose", &treeVars.nJets_bTagLoose);
  chain -> SetBranchAddress("nJets_bTagMedium", &treeVars.nJets_bTagMedium);
  chain -> SetBranchAddress("nJets_bTagTight", &treeVars.nJets_bTagTight);
  
  chain -> SetBranchAddress("MetPt", &treeVars.MetPt);
  chain -> SetBranchAddress("MetPhi", &treeVars.MetPhi);
  
  for(int i=1; i<6; i++)
  {
    chain -> SetBranchAddress(("jet_pt"+ std::to_string(i)).c_str(), &treeVars.jet_pt[i-1]);
    chain -> SetBranchAddress(("jet_eta"+ std::to_string(i)).c_str(), &treeVars.jet_eta[i-1]);
    chain -> SetBranchAddress(("jet_phi"+ std::to_string(i)).c_str(), &treeVars.jet_phi[i-1]);
    chain -> SetBranchAddress(("jet_bdiscriminant"+ std::to_string(i)).c_str(), &treeVars.jet_bdiscriminant[i-1]);
  }
  
  for(int i=1; i<3; i++)
  {
    chain -> SetBranchAddress(("mu_pt"+ std::to_string(i)).c_str(), &treeVars.mu_pt[i-1]);
    chain -> SetBranchAddress(("mu_eta"+ std::to_string(i)).c_str(), &treeVars.mu_eta[i-1]);
    chain -> SetBranchAddress(("mu_phi"+ std::to_string(i)).c_str(), &treeVars.mu_phi[i-1]);
    chain -> SetBranchAddress(("mu_energy"+ std::to_string(i)).c_str(), &treeVars.mu_energy[i-1]);
    chain -> SetBranchAddress(("mu_isLoose"+ std::to_string(i)).c_str(), &treeVars.mu_isLoose[i-1]);
    chain -> SetBranchAddress(("mu_isMedium"+ std::to_string(i)).c_str(), &treeVars.mu_isMedium[i-1]);
    chain -> SetBranchAddress(("mu_MiniIso"+ std::to_string(i)).c_str(), &treeVars.mu_miniIso[i-1]);
    chain -> SetBranchAddress(("mu_charge"+ std::to_string(i)).c_str(), &treeVars.mu_charge[i-1]);
    chain -> SetBranchAddress(("mu_trackIso"+ std::to_string(i)).c_str(), &treeVars.mu_trackIso[i-1]);
    chain -> SetBranchAddress(("mu_sumChargedHadronPt"+ std::to_string(i)).c_str(), &treeVars.mu_sumChargedHadronPt[i-1]);
    chain -> SetBranchAddress(("mu_sumNeutralHadronEt"+ std::to_string(i)).c_str(), &treeVars.mu_sumNeutralHadronEt[i-1]);
    chain -> SetBranchAddress(("mu_sumPhotonEt"+ std::to_string(i)).c_str(), &treeVars.mu_sumPhotonEt[i-1]);
    chain -> SetBranchAddress(("mu_sumPUPt"+ std::to_string(i)).c_str(), &treeVars.mu_sumPUPt[i-1]);
    
    chain -> SetBranchAddress(("ele_pt"+ std::to_string(i)).c_str(), &treeVars.ele_pt[i-1]);
    chain -> SetBranchAddress(("ele_eta"+ std::to_string(i)).c_str(), &treeVars.ele_eta[i-1]);
    chain -> SetBranchAddress(("ele_phi"+ std::to_string(i)).c_str(), &treeVars.ele_phi[i-1]);
    chain -> SetBranchAddress(("ele_energy"+ std::to_string(i)).c_str(), &treeVars.ele_energy[i-1]);
    chain -> SetBranchAddress(("ele_SCeta"+ std::to_string(i)).c_str(), &treeVars.ele_SCeta[i-1]);
    chain -> SetBranchAddress(("ele_SCphi"+ std::to_string(i)).c_str(), &treeVars.ele_SCphi[i-1]);
    chain -> SetBranchAddress(("ele_passLooseId"+ std::to_string(i)).c_str(), &treeVars.ele_passLooseId[i-1]);
    chain -> SetBranchAddress(("ele_passVetoId"+ std::to_string(i)).c_str(), &treeVars.ele_passVetoId[i-1]);
    chain -> SetBranchAddress(("ele_passMediumId"+ std::to_string(i)).c_str(), &treeVars.ele_passMediumId[i-1]);
    chain -> SetBranchAddress(("ele_passTightId"+ std::to_string(i)).c_str(), &treeVars.ele_passTightId[i-1]);
    chain -> SetBranchAddress(("ele_MVAMediumId"+ std::to_string(i)).c_str(), &treeVars.ele_MVAMediumId[i-1]);
    chain -> SetBranchAddress(("ele_MVATightId"+ std::to_string(i)).c_str(), &treeVars.ele_MVATightId[i-1]);
    chain -> SetBranchAddress(("ele_MiniIso"+ std::to_string(i)).c_str(), &treeVars.ele_miniIso[i-1]);
    chain -> SetBranchAddress(("ele_ecalEnergy"+ std::to_string(i)).c_str(), &treeVars.ele_ecalEnergy[i-1]);
    chain -> SetBranchAddress(("ele_SCx"+ std::to_string(i)).c_str(), &treeVars.ele_SCx[i-1]);
    chain -> SetBranchAddress(("ele_SCy"+ std::to_string(i)).c_str(), &treeVars.ele_SCy[i-1]);
    chain -> SetBranchAddress(("ele_SCz"+ std::to_string(i)).c_str(), &treeVars.ele_SCz[i-1]);
    chain -> SetBranchAddress(("ele_charge"+ std::to_string(i)).c_str(), &treeVars.ele_charge[i-1]);
    chain -> SetBranchAddress(("ele_dEtaSCTrackAtVtx"+ std::to_string(i)).c_str(), &treeVars.ele_dEtaTrk[i-1]);
    chain -> SetBranchAddress(("ele_dPhiSCTrackAtVtx"+ std::to_string(i)).c_str(), &treeVars.ele_dPhiTrk[i-1]);
  }
}



void InitOutTreeVars(TTree* tree, TreeVars& treeVars)
{
  tree -> Branch("weight",&treeVars.weight);
  
  tree -> Branch("dipho_mass",    &treeVars.dipho_mass);
  tree -> Branch("dipho_sigmaRV", &treeVars.dipho_sigmaRV);
  tree -> Branch("dipho_deltaphi",&treeVars.dipho_deltaphi);
  tree -> Branch("dipho_mva",     &treeVars.dipho_mva);
  
  tree -> Branch("dipho_leadEta",      &treeVars.dipho_leadEta);
  tree -> Branch("dipho_leadPhi",      &treeVars.dipho_leadPhi);
  tree -> Branch("dipho_lead_ptoM",    &treeVars.dipho_lead_ptoM);
  tree -> Branch("dipho_lead_sigmaEoE",&treeVars.dipho_lead_sigmaEoE);
  tree -> Branch("dipho_leadIDMVA",    &treeVars.dipho_leadIDMVA);
  
  tree -> Branch("dipho_subleadEta",      &treeVars.dipho_subleadEta);
  tree -> Branch("dipho_subleadPhi",      &treeVars.dipho_subleadPhi);
  tree -> Branch("dipho_sublead_ptoM",    &treeVars.dipho_sublead_ptoM);
  tree -> Branch("dipho_sublead_sigmaEoE",&treeVars.dipho_sublead_sigmaEoE);
  tree -> Branch("dipho_subleadIDMVA",    &treeVars.dipho_subleadIDMVA);
  
  tree -> Branch("nJets",           &treeVars.nJets);
  tree -> Branch("nJets_bTagLoose", &treeVars.nJets_bTagLoose);
  tree -> Branch("nJets_bTagMedium",&treeVars.nJets_bTagMedium);
  tree -> Branch("nJets_bTagTight", &treeVars.nJets_bTagTight);
}


void InitRawTreeVars(std::map<std::string,TChain*> &chain, RawTreeVars& treeVars)
{

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

  chain["PhotonTight"]->SetBranchStatus("PhotonTight_size", 1);
  chain["PhotonTight"]->SetBranchAddress("PhotonTight_size",       &treeVars.N_TightPh);
  chain["PhotonTight"]->SetBranchStatus("Particle", 1);
  chain["PhotonTight"]->SetBranchAddress("Particle",               treeVars.TightPh_g);
  chain["PhotonTight"]->SetBranchAddress("IsEB",                   treeVars.TightPh_isEB);
  chain["PhotonTight"]->SetBranchStatus("PT", 1);
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
