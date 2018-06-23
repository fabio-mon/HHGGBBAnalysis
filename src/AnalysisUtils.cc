#include "interface/AnalysisUtils.h"
#include "TLorentzVector.h"

float DeltaEta(const float& eta1, const float& eta2)
{
  return fabs( eta1 - eta2 );
}

float DeltaPhi(const float& phi1, const float& phi2)
{
  float dphi = fabs( phi1 - phi2 );
  if( dphi > PI ) dphi = 2*PI - dphi;
  return dphi;
}

float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2)
{
  return sqrt( DeltaEta(eta1,eta2)*DeltaEta(eta1,eta2) + 
               DeltaPhi(phi1,phi2)*DeltaPhi(phi1,phi2) );
}

void MakePlot3(std::map<std::string,TH1F*> &h)
{
  TCanvas* c = new TCanvas();
  c -> cd();

  for(std::map<std::string,TH1F*>::iterator it=h.begin(); it!=h.end(); ++it)
  {
    it->second -> SetMarkerStyle(20);	       
    it->second -> SetMarkerSize(1);
    it->second -> SetMarkerColor(kBlack);
    it->second -> SetFillStyle(0);
    it->second -> Draw("E1");

    if(it->first == "dipho_mass")
      it->second -> GetXaxis() -> SetTitle("diphoton mass (GeV/c^{2})");
    if(it->first == "dipho_sumpt")
      it->second -> GetXaxis() -> SetTitle("diphoton P_{T} sum (GeV/c)");
    if(it->first == "dipho_deltaeta")
      it->second -> GetXaxis() -> SetTitle("diphoton #Delta#eta");
    if(it->first == "dipho_deltaphi")
      it->second -> GetXaxis() -> SetTitle("diphoton #Delta#Phi");
    if(it->first == "dipho_leadPt")
      it->second -> GetXaxis() -> SetTitle("lead phot P_{T} (GeV/c)");
    if(it->first == "dipho_leadEta")
      it->second -> GetXaxis() -> SetTitle("lead phot #eta");
    if(it->first == "dipho_leadPhi")
      it->second -> GetXaxis() -> SetTitle("lead phot #Phi");
    if(it->first == "dipho_leadptoM")
      it->second -> GetXaxis() -> SetTitle("lead phot P_{T} / diphoton_mass (c)");
    if(it->first == "dipho_subleadPt")
      it->second -> GetXaxis() -> SetTitle("sublead phot P_{T} (GeV/c)");
    if(it->first == "dipho_subleadEta")
      it->second -> GetXaxis() -> SetTitle("sublead phot #eta");
    if(it->first == "dipho_subleadPhi")
      it->second -> GetXaxis() -> SetTitle("sublead phot #Phi");
    if(it->first == "dipho_subleadptoM")
      it->second -> GetXaxis() -> SetTitle("sublead phot P_{T} / diphoton_mass (c)");
    if(it->first == "nJets")
      it->second -> GetXaxis() -> SetTitle("# Jets");
    CMS_lumi(c, 0, 0);
    c -> SaveAs(("c_" +it->first + ".png").c_str());
    c -> SaveAs(("c_" +it->first + ".pdf").c_str());

  }

}

bool DiPhotonSelection(const TLorentzVector &pho_lead ,const TLorentzVector &pho_sublead)
{
  double lead_pt=pho_lead.Pt();
  double sublead_pt=pho_sublead.Pt();
  if(lead_pt<20 || sublead_pt<20) return false;
  if(fabs(pho_lead.Eta())>3 || fabs(pho_sublead.Eta())>3) return false;
  double dipho_mass=(pho_lead+pho_sublead).M();
  if(dipho_mass<100. || dipho_mass>180.) return false;
  if(lead_pt/dipho_mass<0.33) return false;
  if(sublead_pt/dipho_mass<0.25) return false;
  
  return true;
}

bool JetSelection(const RawTreeVars &treeVars, TreeVars &outtreeVars)//const TLorentzVector &bjet_lead, const TLorentzVector &bjet_sublead, const TLorentzVector &pho_lead, const TLorentzVector &pho_sublead)
{
  //check that selected jets are far from any reco photons
  //int N_hardJet=0;
  //for(int i=0;i<treeVars.N_TightPh;i++)
  //{
  //  if(treeVars.Jet_pt[i]<25) continue;
  //  if( DeltaRmin_phoRECO_jetRECO(pho_lead,treeVars)<0.4 ) continue;
  //  if( DeltaRmin_phoRECO_jetRECO(pho_sublead,treeVars)<0.4 ) continue;
  //  if(treeVars.Jet_mvav2[i]>0)
  //    N_hardJet++;
  //}      
  
  //if(N_hardJet<2) return false;
  return true;
}

bool SelectBestScoreBJets(const TreeVars &outtreeVars,int &bjet_lead_i,int &bjet_sublead_i,const bool &useMTD)
{
  int BTagOffset;
  if(useMTD == false)
    BTagOffset=0;
  else
    BTagOffset=3;

  //find jet with higher btag score
  int ijet1=0;
  int btag1=outtreeVars.jet_BTagLevel[0];
  float pt1 = outtreeVars.jet_pt[0];
  //---------------------------------------------------------------
  cout<<"0\tbtag"<<btag1<<"\tpt1"<<pt1<<endl;
  //---------------------------------------------------------------
  for(int i=1; i<outtreeVars.nJets; i++)
  {
    //---------------------------------------------------------------
    cout<<i<<"\tbtag"<<outtreeVars.jet_BTagLevel[i]<<"\tpt1"<<outtreeVars.jet_pt[i]<<endl;
    //---------------------------------------------------------------
    //If(BTag>BTagOffset && BTag<4+BTagOffset) <-- already required to fill outtreeVars
    if(outtreeVars.jet_BTagLevel[i] > btag1)
    {
      ijet1=i;
      btag1=outtreeVars.jet_BTagLevel[i];
      pt1 = outtreeVars.jet_pt[i];
    }
    else
      if(outtreeVars.jet_BTagLevel[i] == btag1)
	if(outtreeVars.jet_pt[i] > pt1)
	{
	  ijet1=i;
	  btag1=outtreeVars.jet_BTagLevel[i];
	  pt1 = outtreeVars.jet_pt[i];
	}
  }

  //find second jet with higher btag score
  int ijet2=-1;
  int btag2=-1;
  float pt2=-1;
  for(int i=0; i<outtreeVars.nJets; i++)
  {
    if(i==ijet1) continue;
    if(ijet2==-1)
    {
      ijet2=i;
      btag2=outtreeVars.jet_BTagLevel[i];
      pt2 = outtreeVars.jet_pt[i];
    }
    else
      if(outtreeVars.jet_BTagLevel[i] > btag2)
      {
	ijet2=i;
	btag2=outtreeVars.jet_BTagLevel[i];
	pt2 = outtreeVars.jet_pt[i];
      }
      else
	if(outtreeVars.jet_BTagLevel[i] == btag2)
	  if(outtreeVars.jet_pt[i] > pt2)
	  {
	    ijet2=i;
	    btag2=outtreeVars.jet_BTagLevel[i];
	    pt2 = outtreeVars.jet_pt[i];
	  }
  }

  //---------------------------------------------------------------
  cout<<"best score i "<<ijet1<<" "<<ijet2<<endl;
  //---------------------------------------------------------------
  if(ijet2==-1)
    return false;

  if(pt1>pt2)
  {
    bjet_lead_i = ijet1;
    bjet_sublead_i = ijet2;
  }
  else
  {
    bjet_lead_i = ijet2;
    bjet_sublead_i = ijet1;
  }

  return true;
}


int GetBTagLevel(int BTag)
{
  if(BTag & 0b100000) return 6;
  if(BTag & 0b010000) return 5;
  if(BTag & 0b001000) return 4;
  if(BTag & 0b000100) return 3;
  if(BTag & 0b000010) return 2;
  if(BTag & 0b000001) return 1;
  return 0;
}


bool FindGenPh_Hdaug(RawTreeVars &treeVars, float deltaMthr)
//function that find the pair of photons that minimize M_diphogen - M_higgs
//if min(M_diphogen - M_higgs)>deltaM_max it returns false, otherwise true
{
  if(treeVars.N_GenPh<2) return false;
  //initialize treeVars.GenPh_isHdaug
  for(int i=0;i<treeVars.N_GenPh;i++)
    treeVars.GenPh_isHdaug[i]=false;

  double deltaMaux,deltaMmin;
  int imin,jmin;
  TLorentzVector Vi,Vj;
  for(int i=0;i<treeVars.N_GenPh;i++)
  {
    Vi.SetPtEtaPhiE(treeVars.GenPh_pt[i],treeVars.GenPh_eta[i],treeVars.GenPh_phi[i],treeVars.GenPh_E[i]);
    for(int j=i+1;j<treeVars.N_GenPh;j++) 
    {  
      //cout<<"("<<i<<","<<j<<")\t";
      Vj.SetPtEtaPhiE(treeVars.GenPh_pt[j],treeVars.GenPh_eta[j],treeVars.GenPh_phi[j],treeVars.GenPh_E[j]);
      deltaMaux=fabs((Vi+Vj).M()-125.);
      //cout<<deltaMaux<<endl;
      if(i==0 && j==1)//first loop
      {
	imin=0;
	jmin=1;
	deltaMmin=deltaMaux;
      }
      else
	if(deltaMaux<deltaMmin)
	{
	  imin=i;
	  jmin=j;
	  deltaMmin=deltaMaux;
	}
    }
  }
  //cout<<"MIN "<<"("<<imin<<","<<jmin<<")"<<endl;
  treeVars.GenPh_isHdaug[imin]=true;
  treeVars.GenPh_isHdaug[jmin]=true;
  if(deltaMmin<deltaMthr) return true;
  return false;

}


void  FindLeadSublead_pho(const RawTreeVars &treeVars, int &pho_lead_i, int &pho_sublead_i)
{
  if(treeVars.N_TightPh<2) return;
  float ptmax1, ptmax2;
  if(treeVars.TightPh_pt[0]>treeVars.TightPh_pt[1])
  {
    pho_lead_i=0;
    ptmax1=treeVars.TightPh_pt[0];
    pho_sublead_i=1;
    ptmax2=treeVars.TightPh_pt[1];
  }
  else
  {
    pho_lead_i=1;
    ptmax1=treeVars.TightPh_pt[1];
    pho_sublead_i=0;
    ptmax2=treeVars.TightPh_pt[0];
  }
  //cout<<"\n0 "<<treeVars.TightPh_pt[0]<<endl;
  //cout<<"1 "<<treeVars.TightPh_pt[1]<<endl;
  for(int i=2;i<treeVars.N_TightPh;i++)
  {
    //cout<<i<<" "<<treeVars.TightPh_pt[i]<<endl;
    if(treeVars.TightPh_pt[i]>ptmax1)
    {
      ptmax2=ptmax1;
      pho_sublead_i=pho_lead_i;
      ptmax1=treeVars.TightPh_pt[i];
      pho_lead_i=i;
    }
    else
      if(treeVars.TightPh_pt[i]>ptmax2)
      {
        ptmax2=treeVars.TightPh_E[i];
        pho_sublead_i=i;
      }
  }
  //cout<<"i_LEAD="<<pho_lead_i<<"\ti_SUBLEAD="<<pho_sublead_i<<endl;
}


bool  FindLeadSublead_bjet(const RawTreeVars &treeVars, int &bjet_lead_i, int &bjet_sublead_i)
{  
  if(treeVars.N_Jet<2) return false;

  bjet_lead_i=-1;
  bjet_sublead_i=-1;
  float ptmax1, ptmax2;
  for(int i=0;i<treeVars.N_Jet;i++)
  {
    if(treeVars.Jet_mvav2[i]==0) continue;//NOT a bjet
    cout<<i<<" "<<treeVars.Jet_pt[i]<<endl;

    if(bjet_lead_i==-1)
    {
      bjet_lead_i=i;
      ptmax1=treeVars.Jet_pt[i];
    }
    else
      if(bjet_sublead_i==-1)
      {
	if(treeVars.Jet_pt[i]>ptmax1)
	{
	  bjet_sublead_i=bjet_lead_i;
	  ptmax2=ptmax1;
	  bjet_lead_i=i;
	  ptmax1=treeVars.Jet_pt[i];
	}
	else
	{
	  bjet_sublead_i=i;
	  ptmax2=treeVars.Jet_pt[i];
	}
      }
      else
	if(treeVars.Jet_pt[i]>ptmax1)
        {
	  ptmax2=ptmax1;
	  bjet_sublead_i=bjet_lead_i;
	  ptmax1=treeVars.Jet_pt[i];
	  bjet_lead_i=i;
	}
	else
	  if(treeVars.Jet_pt[i]>ptmax2)
	  {
	    ptmax2=treeVars.Jet_pt[i];
	    bjet_sublead_i=i;
	  }
  }
  cout<<"i_LEAD="<<bjet_lead_i<<"\ti_SUBLEAD="<<bjet_sublead_i<<endl;
  if(bjet_sublead_i==-1)
    return false;
  return true;
}


bool RecoJetGenericMatch (const TLorentzVector &reco_pho , const RawTreeVars& treeVars , TLorentzVector& reco_jet_match, float DeltaRmax)
{

  float reco_pho_eta = reco_pho.Eta();
  float reco_pho_phi = reco_pho.Phi();
  float reco_jet_eta;
  float reco_jet_phi;

  for(int i=0;i<treeVars.N_Jet;++i)
  {
    reco_jet_eta=treeVars.Jet_eta[i];
    reco_jet_phi=treeVars.Jet_phi[i];
    if( DeltaR(reco_pho_eta,reco_pho_phi,reco_jet_eta,reco_jet_phi) < DeltaRmax )
    {
      reco_jet_match.SetPtEtaPhiM(treeVars.Jet_pt[i],reco_jet_eta,reco_jet_phi,treeVars.Jet_mass[i]);
      return true;
    }
  }

  return false;
}


bool PhoGenericGenMatch (const TLorentzVector &reco_pho , const RawTreeVars& treeVars ,  TLorentzVector& gen_pho_match, float DeltaRmax)
{

  float reco_pho_eta = reco_pho.Eta();
  float reco_pho_phi = reco_pho.Phi();
  float gen_pho_eta;
  float gen_pho_phi;

  for(int i=0;i<treeVars.N_GenPh;++i)
    {
      gen_pho_eta=treeVars.GenPh_eta[i];
      gen_pho_phi=treeVars.GenPh_phi[i];
      if(fabs(treeVars.GenPh_E[i]-reco_pho.E()) > 5) continue;
      if( DeltaR(reco_pho_eta,reco_pho_phi,gen_pho_eta,gen_pho_phi) < DeltaRmax )
      {
	gen_pho_match.SetPtEtaPhiE(treeVars.GenPh_pt[i],gen_pho_eta,gen_pho_phi,treeVars.GenPh_E[i]);
	//cout<<"PtGEN-PtRECO="<<treeVars.GenPh_pt[i]-reco_pho.Pt()<<endl;
	return true;
      }
    }

  return false;
}


bool PhoGenMatch(const TLorentzVector &pho_lead , const TLorentzVector &pho_sublead , const RawTreeVars& treeVars , TreeVars &outtreeVars, float DeltaRmax)
//gen-matching with photons that are higgs daughter
{

  float pho_lead_eta = pho_lead.Eta();
  float pho_lead_phi = pho_lead.Phi();
  float pho_sublead_eta = pho_sublead.Eta();
  float pho_sublead_phi = pho_sublead.Phi();
  //cout<<"LEAD (ETA,PHI)="<<pho_lead.Eta()<<","<<pho_lead.Phi()<<")"<<endl;
  //cout<<"SUBLEAD (ETA,PHI)="<<pho_sublead.Eta()<<","<<pho_sublead.Phi()<<")"<<endl;
  float pho_gen1_eta;
  float pho_gen1_phi=-999;
  float pho_gen2_eta;
  float pho_gen2_phi=-999;

  outtreeVars.dipho_leadEnergy_gen = 0;
  outtreeVars.dipho_leadEta_gen = 0;
  outtreeVars.dipho_leadPhi_gen = -999;
    
  outtreeVars.dipho_subleadEnergy_gen = 0;
  outtreeVars.dipho_subleadEta_gen = 0;
  outtreeVars.dipho_subleadPhi_gen = -999;

  //find genphotons flagged as higgs daughter
  int i=0;
  int i_pho_gen1;
  int i_pho_gen2;
  while(pho_gen1_phi==-999 && i<treeVars.N_GenPh)
  {
    //cout<<i<<" GEN (ETA,PHI)="<<treeVars.GenPh_eta[i]<<","<<treeVars.GenPh_phi[i]<<")"<<endl;
    if(treeVars.GenPh_isHdaug[i]==true)
    {
      pho_gen1_eta=treeVars.GenPh_eta[i];
      pho_gen1_phi=treeVars.GenPh_phi[i];
      i_pho_gen1=i;
    }
    ++i;
  }
  //cout<<"-----"<<endl;
  if(pho_gen1_phi==-999) return false;//-->higgs daughter 1 not found

  while(pho_gen2_phi==-999 && i<treeVars.N_GenPh)
  {
    //cout<<i<<" GEN (ETA,PHI)="<<treeVars.GenPh_eta[i]<<","<<treeVars.GenPh_phi[i]<<")"<<endl;
    if(treeVars.GenPh_isHdaug[i]==true)
      {
	pho_gen2_eta=treeVars.GenPh_eta[i];
	pho_gen2_phi=treeVars.GenPh_phi[i];
	i_pho_gen2=i;
      }
    ++i;
  }
  //cout<<"-----"<<endl;
  if(pho_gen2_phi==-999) return false;//-->higgs daughter 2 not found

  //eta-phi gen-matching
  //cout<<"DeltaR1lead = \t"<<DeltaR(pho_gen1_eta,pho_gen1_phi,   pho_lead_eta,   pho_lead_phi)<<endl;
  //cout<<"DeltaR2sublead = \t"<<DeltaR(pho_gen2_eta,pho_gen2_phi,pho_sublead_eta,pho_sublead_phi)<<endl;
  TLorentzVector phlead_gen;
  TLorentzVector phsublead_gen;
  
  if(DeltaR(pho_gen1_eta,pho_gen1_phi,pho_lead_eta,pho_lead_phi)<DeltaRmax)
  {
    outtreeVars.dipho_leadEnergy_gen = treeVars.GenPh_E[i_pho_gen1];
    outtreeVars.dipho_leadEta_gen = treeVars.GenPh_eta[i_pho_gen1];
    outtreeVars.dipho_leadPhi_gen = treeVars.GenPh_phi[i_pho_gen1];
    phlead_gen.SetPtEtaPhiE(treeVars.GenPh_pt[i_pho_gen1],outtreeVars.dipho_leadEta_gen,outtreeVars.dipho_leadPhi_gen,outtreeVars.dipho_leadEnergy_gen);
  }

  if(DeltaR(pho_gen2_eta,pho_gen2_phi,pho_sublead_eta,pho_sublead_phi)<DeltaRmax)
  {
    outtreeVars.dipho_subleadEnergy_gen = treeVars.GenPh_E[i_pho_gen2];
    outtreeVars.dipho_subleadEta_gen = treeVars.GenPh_eta[i_pho_gen2];
    outtreeVars.dipho_subleadPhi_gen = treeVars.GenPh_phi[i_pho_gen2];
    phsublead_gen.SetPtEtaPhiE(treeVars.GenPh_pt[i_pho_gen2],outtreeVars.dipho_subleadEta_gen,outtreeVars.dipho_subleadPhi_gen,outtreeVars.dipho_subleadEnergy_gen);
  }
  
  if(DeltaR(pho_gen1_eta,pho_gen1_phi,   pho_lead_eta,   pho_lead_phi)<DeltaRmax &&
     DeltaR(pho_gen2_eta,pho_gen2_phi,pho_sublead_eta,pho_sublead_phi)<DeltaRmax)
  {
    outtreeVars.dipho_mass_gen = (phlead_gen+phsublead_gen).M();
    return true;
  }


  if( DeltaR(pho_gen2_eta,pho_gen2_phi, pho_lead_eta, pho_lead_phi)<DeltaRmax)
  {
    outtreeVars.dipho_leadEnergy_gen = treeVars.GenPh_E[i_pho_gen2];
    outtreeVars.dipho_leadEta_gen = treeVars.GenPh_eta[i_pho_gen2];
    outtreeVars.dipho_leadPhi_gen = treeVars.GenPh_phi[i_pho_gen2];
    phlead_gen.SetPtEtaPhiE(treeVars.GenPh_pt[i_pho_gen2],outtreeVars.dipho_leadEta_gen,outtreeVars.dipho_leadPhi_gen,outtreeVars.dipho_leadEnergy_gen);
  }

  if( DeltaR(pho_gen1_eta,pho_gen1_phi,pho_sublead_eta,pho_sublead_phi)<DeltaRmax )
  {
    outtreeVars.dipho_subleadEnergy_gen = treeVars.GenPh_E[i_pho_gen1];
    outtreeVars.dipho_subleadEta_gen = treeVars.GenPh_eta[i_pho_gen1];
    outtreeVars.dipho_subleadPhi_gen = treeVars.GenPh_phi[i_pho_gen1];
    phsublead_gen.SetPtEtaPhiE(treeVars.GenPh_pt[i_pho_gen1],outtreeVars.dipho_subleadEta_gen,outtreeVars.dipho_subleadPhi_gen,outtreeVars.dipho_subleadEnergy_gen);
  }

  if( DeltaR(pho_gen2_eta,pho_gen2_phi,   pho_lead_eta,   pho_lead_phi)<DeltaRmax && 
      DeltaR(pho_gen1_eta,pho_gen1_phi,pho_sublead_eta,pho_sublead_phi)<DeltaRmax) 
  {
    outtreeVars.dipho_mass_gen = (phlead_gen+phsublead_gen).M();
    return true;
  }

  return false;

}

float DeltaRmin_phoRECO_phoGEN(const TLorentzVector &reco_pho , const RawTreeVars& treeVars)
{
  float reco_pho_eta = reco_pho.Eta();
  float reco_pho_phi = reco_pho.Phi();
  float gen_pho_eta;
  float gen_pho_phi;
  float deltaR;
  float deltaRmin=-1;

  for(int i=0;i<treeVars.N_GenPh;++i)
  {
    gen_pho_eta=treeVars.GenPh_eta[i];
    gen_pho_phi=treeVars.GenPh_phi[i];
    if(fabs(reco_pho.E()-treeVars.GenPh_E[i])>5) continue;
    deltaR = DeltaR(reco_pho_eta,reco_pho_phi,gen_pho_eta,gen_pho_phi);
    if(deltaRmin==-1)
      deltaRmin=deltaR;
    else
      if(deltaR<deltaRmin)
	deltaRmin=deltaR;
  }
  return deltaRmin;
}

float DeltaRmin_phoRECO_jetRECO(const TLorentzVector &reco_pho , const RawTreeVars& treeVars)
{
  float reco_pho_eta = reco_pho.Eta();
  float reco_pho_phi = reco_pho.Phi();
  float reco_jet_eta;
  float reco_jet_phi;
  float deltaR;
  float deltaRmin=-1;

  for(int i=0;i<treeVars.N_Jet;++i)
  {
    reco_jet_eta=treeVars.Jet_eta[i];
    reco_jet_phi=treeVars.Jet_phi[i];
    deltaR = DeltaR(reco_pho_eta,reco_pho_phi,reco_jet_eta,reco_jet_phi);
    if(deltaRmin==-1)
      deltaRmin=deltaR;
    else
      if(deltaR<deltaRmin)
	deltaRmin=deltaR;
  }
  return deltaRmin;
}


void PrintRecoPhoton(const RawTreeVars& treeVars)
{
  cout<<"Event "<<treeVars.event<<"\tRECO Tight photon collection -> "<<treeVars.N_TightPh<<" entries"<<endl;
  for(int i=0;i<treeVars.N_TightPh;i++)
    cout<<i<<"\tEta="<<treeVars.TightPh_eta[i]<< "\tPhi="<<treeVars.TightPh_phi[i]<< "\tPt="<<treeVars.TightPh_pt[i]<<endl;
  cout<<endl;
}

void PrintRecoJet(const RawTreeVars& treeVars)
{
  cout<<"Event "<<treeVars.event<<"\tRECO PUPPI jet collection -> "<<treeVars.N_Jet<<" entries"<<endl;
  for(int i=0;i<treeVars.N_Jet;i++)
    cout<<i<<"\tEta="<<treeVars.Jet_eta[i]<< "\tPhi="<<treeVars.Jet_phi[i]<< "\tPt="<<treeVars.Jet_pt[i]<<endl;
  cout<<endl;
}
