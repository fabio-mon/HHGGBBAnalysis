#include "interface/AnalysisUtils.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
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

    if(it->first == "mgg")
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
  if(lead_pt<30 || sublead_pt<20) return false;
  if(fabs(pho_lead.Eta())>2.5 || fabs(pho_sublead.Eta())>2.5) return false;
  if(fabs(pho_lead.Eta())>1.44 && fabs(pho_lead.Eta())<1.57) return false;
  if(fabs(pho_sublead.Eta())>1.44 && fabs(pho_sublead.Eta())<1.57) return false;
  double mgg=(pho_lead+pho_sublead).M();
  //---------------------------------------------------------------
  //cout<<"Mgg="<<mgg<<endl;
  //cout<<"lead_pt/Mgg="<<lead_pt/mgg<<endl;
  //cout<<"sublead_pt/Mgg="<<sublead_pt/mgg<<endl;
  
  //-------------------------------------------
  if(mgg<100. || mgg>180.) return false;
  if(lead_pt/mgg<0.33) return false;
  if(sublead_pt/mgg<0.25) return false;
  
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

//reduce of a fraction "prob" the medium btagged jets from selected jets collection
void ReducebtagEfficiency(TreeVars &outtreeVars, const bool &useMTD, const float &prob )
{
  TRandom3 random;
  random.SetSeed(0);
  int BTagMedium_mask;
  if(useMTD == false)
    BTagMedium_mask=0b000010;
  else
    BTagMedium_mask=0b010000;

  //random deletion
  for(int i=0; i<outtreeVars.nJets; i++)
  {
    //---------------------------------------------------------------
    //cout<<i<<"\tbtag="<<(outtreeVars.jet_mvav2[i] & BTagMedium_mask)<<"\tflavor="<<outtreeVars.jet_hadflav[i];
    //---------------------------------------------------------------
    if(outtreeVars.jet_mvav2[i] & BTagMedium_mask)
      if(fabs(outtreeVars.jet_hadflav[i])==4 || fabs(outtreeVars.jet_hadflav[i])==5)
	if(random.Uniform()<prob)
	{
	  //cout<<"deleted"<<endl;
	  if(useMTD)
	    outtreeVars.jet_mvav2[i] -= 16;
	  else
	    outtreeVars.jet_mvav2[i] -= 2;	    
	  outtreeVars.nJets_bTagMedium--;
	  outtreeVars.jet_BTagMedium[i] = 0;
	}
    //cout<<endl;
  }
}

//increase of a fraction "prob" the medium btagged jets from selected jets collection
void IncreasebtagEfficiency(TreeVars &outtreeVars, const bool &useMTD, const float &prob )
{
  TRandom3 random;
  random.SetSeed(0);
  int BTagMedium_mask;
  if(useMTD == false)
    BTagMedium_mask=0b000010;
  else
    BTagMedium_mask=0b010000;

  //random deletion
  for(int i=0; i<outtreeVars.nJets; i++)
  {
    //---------------------------------------------------------------
    //cout<<i<<"\tbtag="<<(outtreeVars.jet_mvav2[i] & BTagMedium_mask)<<"\tflavor="<<outtreeVars.jet_hadflav[i];
    //---------------------------------------------------------------
    if(!(outtreeVars.jet_mvav2[i] & BTagMedium_mask))
      if(fabs(outtreeVars.jet_hadflav[i])==4 || fabs(outtreeVars.jet_hadflav[i])==5)
	if(random.Uniform()<prob)
	{
	  //cout<<"added"<<endl;
	  if(useMTD)
	    outtreeVars.jet_mvav2[i] += 16;
	  else
	    outtreeVars.jet_mvav2[i] += 2;	    
	  outtreeVars.nJets_bTagMedium++;
	  outtreeVars.jet_BTagMedium[i] = 1;
	}
    //cout<<endl;
  }
}


bool SelectBestScoreBJets2(const TreeVars &outtreeVars,int &bjet_lead_i,int &bjet_sublead_i,const bool &useMTD )
{
  int BTagMedium_mask;
  if(useMTD == false)
    BTagMedium_mask=0b000010;
  else
    BTagMedium_mask=0b010000;

  //select the 2 jets with medium btag score if present , or take the harder ones
  int ijet1=0;
  int btag1=outtreeVars.jet_mvav2[0] & BTagMedium_mask;
  float pt1 = outtreeVars.jet_pt[0];
  //---------------------------------------------------------------
  //cout<<"0\tbtag="<<btag1<<"\tpt1="<<pt1<<endl;
  //---------------------------------------------------------------
  for(int i=1; i<outtreeVars.nJets; i++)
  {
    //---------------------------------------------------------------
    //cout<<i<<"\tbtag="<<(outtreeVars.jet_mvav2[i] & BTagMedium_mask)<<"\tpt1="<<outtreeVars.jet_pt[i]<<endl;
    //---------------------------------------------------------------
    //If(BTag>BTagOffset && BTag<4+BTagOffset) <-- already required to fill outtreeVars
    if((outtreeVars.jet_mvav2[i] & BTagMedium_mask) > btag1)
    {
      ijet1=i;
      btag1=(outtreeVars.jet_mvav2[i] & BTagMedium_mask);
      pt1 = outtreeVars.jet_pt[i];
    }
    else
      if((outtreeVars.jet_mvav2[i] & BTagMedium_mask) == btag1)
	if(outtreeVars.jet_pt[i] > pt1)
	{
	  ijet1=i;
	  btag1=(outtreeVars.jet_mvav2[i] & BTagMedium_mask);
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
      btag2=(outtreeVars.jet_mvav2[i] & BTagMedium_mask);
      pt2 = outtreeVars.jet_pt[i];
    }
    else
      if((outtreeVars.jet_mvav2[i] & BTagMedium_mask) > btag2)
      {
	ijet2=i;
	btag2=(outtreeVars.jet_mvav2[i] & BTagMedium_mask);
	pt2 = outtreeVars.jet_pt[i];
      }
      else
	if((outtreeVars.jet_mvav2[i] & BTagMedium_mask) == btag2)
	  if(outtreeVars.jet_pt[i] > pt2)
	  {
	    ijet2=i;
	    btag2=(outtreeVars.jet_mvav2[i] & BTagMedium_mask);
	    pt2 = outtreeVars.jet_pt[i];
	  }
  }

  //---------------------------------------------------------------
  //cout<<"best score i "<<ijet1<<" "<<ijet2<<endl;
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
 
bool Findbquark_Hdaug(RawTreeVars &treeVars)
{
  int nflag=0;
  for(int i=0;i<treeVars.N_GenPart;i++)
  {
    treeVars.GenPart_isHdaug[i]=false;
    if(fabs(treeVars.GenPart_pid[i])==5 && treeVars.GenPart_st[i]==23)
    {
      treeVars.GenPart_isHdaug[i]=true;
      nflag++;
    }
  }
  if(nflag!=2)
  {
    cout<<"n bquark h daug="<<nflag<<" --->skip event"<<endl;
    for(int i=0;i<treeVars.N_GenPart;i++)
      if(fabs(treeVars.GenPart_pid[i])==5)
	cout<<"status="<<treeVars.GenPart_st[i]<<"\tpt="<<treeVars.GenPart_pt[i]<<"\teta="<<treeVars.GenPart_eta[i]<<"\tphi="<<treeVars.GenPart_phi[i]<<endl;
    return false;
  }
  return true;
}

bool FindGenJet_Hdaug(RawTreeVars &treeVars, float deltaMthr)
//function that find the pair of jenjets that minimize M_digenjet - M_higgs
//if min(M_digenjet - M_higgs)>deltaM_max it returns false, otherwise true
{
  //cout<<"in function  FindGenJet_Hdaug"<<endl;
  if(treeVars.N_GenJet<2) return false;
  //initialize treeVars.GenJet_isHdaug
  for(int i=0;i<treeVars.N_GenJet;i++)
    treeVars.GenJet_isHdaug[i]=false;

  double deltaMaux,deltaMmin;
  int imin=-1,jmin=-1;
  TLorentzVector Vi,Vj;
  for(int i=0;i<treeVars.N_GenJet;i++)
  {
    if(treeVars.GenJet_pt[i]<20) continue;
    Vi.SetPtEtaPhiM(treeVars.GenJet_pt[i],treeVars.GenJet_eta[i],treeVars.GenJet_phi[i],treeVars.GenJet_mass[i]);
    for(int j=i+1;j<treeVars.N_GenJet;j++) 
    {  
      //cout<<"("<<i<<","<<j<<")\t";
      if(treeVars.GenJet_pt[j]<20) continue;
      Vj.SetPtEtaPhiM(treeVars.GenJet_pt[j],treeVars.GenJet_eta[j],treeVars.GenJet_phi[j],treeVars.GenJet_mass[j]);
      deltaMaux=fabs((Vi+Vj).M()-125.);
      //cout<<deltaMaux<<endl;
      if(imin==-1 && jmin==-1)//first loop
      {
	imin=i;
	jmin=j;
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
  treeVars.GenJet_isHdaug[imin]=true;
  treeVars.GenJet_isHdaug[jmin]=true;
  if(deltaMmin<deltaMthr) return true;
  return false;

}


void  FindLeadSublead_pho(const RawTreeVars &treeVars, int &pho_lead_i, int &pho_sublead_i)
{
  if(treeVars.N_SelectedPh<2) return;
  float ptmax1, ptmax2;
  if(treeVars.SelectedPh_pt[0]>treeVars.SelectedPh_pt[1])
  {
    pho_lead_i=0;
    ptmax1=treeVars.SelectedPh_pt[0];
    pho_sublead_i=1;
    ptmax2=treeVars.SelectedPh_pt[1];
  }
  else
  {
    pho_lead_i=1;
    ptmax1=treeVars.SelectedPh_pt[1];
    pho_sublead_i=0;
    ptmax2=treeVars.SelectedPh_pt[0];
  }

  //cout<<"\n0 "<<treeVars.SelectedPh_pt[0]<<endl;
  //cout<<"1 "<<treeVars.SelectedPh_pt[1]<<endl;
  for(int i=2;i<treeVars.N_SelectedPh;i++)
  {
    //cout<<i<<" "<<treeVars.SelectedPh_pt[i]<<endl;
    if(treeVars.SelectedPh_pt[i]>ptmax1)
    {
      ptmax2=ptmax1;
      pho_sublead_i=pho_lead_i;
      ptmax1=treeVars.SelectedPh_pt[i];
      pho_lead_i=i;
    }
    else
      if(treeVars.SelectedPh_pt[i]>ptmax2)
      {
        ptmax2=treeVars.SelectedPh_E[i];
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
    //cout<<i<<" "<<treeVars.Jet_pt[i]<<endl;

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
  //cout<<"i_LEAD="<<bjet_lead_i<<"\ti_SUBLEAD="<<bjet_sublead_i<<endl;
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

  //printf("(eta,phi)reco=(%.2f,%.2f)\tE=%.2f\n\n",reco_pho_eta,reco_pho_phi,reco_pho.E());
  for(int i=0;i<treeVars.N_GenPh;++i)
    {
      if(treeVars.GenPh_E[i]<10) continue;
      gen_pho_eta=treeVars.GenPh_eta[i];
      gen_pho_phi=treeVars.GenPh_phi[i];
      //printf("(eta,phi)gen=(%.2f,%.2f)\tE=%.2f\n",gen_pho_eta,gen_pho_phi,treeVars.GenPh_E[i]);
      //cout<<"E_GEN-E_RECO="<<treeVars.GenPh_E[i]-reco_pho.E()<<endl;
      if(fabs(treeVars.GenPh_E[i]-reco_pho.E()) > 18) continue;
      if( DeltaR(reco_pho_eta,reco_pho_phi,gen_pho_eta,gen_pho_phi) < DeltaRmax )
      {
	gen_pho_match.SetPtEtaPhiE(treeVars.GenPh_pt[i],gen_pho_eta,gen_pho_phi,treeVars.GenPh_E[i]);
	return true;
      }
    }

  return false;
}

bool bquarkGenericGenMatch(const TLorentzVector &reco_jet , const RawTreeVars& treeVars , TLorentzVector &gen_b_match, float DeltaRmax)
{

  float reco_jet_eta = reco_jet.Eta();
  float reco_jet_phi = reco_jet.Phi();
  float gen_b_eta;
  float gen_b_phi;
  float gen_b_E;


  //printf("########\n(eta,phi)reco=(%.2f,%.2f)\tPt=%.2fE=%.2f\n",reco_jet_eta,reco_jet_phi,reco_jet.Pt(),reco_jet.E());
  for(int i=0;i<treeVars.N_GenPart;++i)
  {
     if(treeVars.GenPart_pt[i]<10) continue;
     if(abs(treeVars.GenPart_pid[i])!=5 && abs(treeVars.GenPart_pid[i])!=4) continue;
     //cout<<"pid "<<treeVars.GenPart_pid[i]<<endl;
     gen_b_eta=treeVars.GenPart_eta[i];
     gen_b_phi=treeVars.GenPart_phi[i];
     gen_b_E=sqrt(treeVars.GenPart_pt[i]*treeVars.GenPart_pt[i] + treeVars.GenPart_pz[i]*treeVars.GenPart_pz[i]);
     //printf("(eta,phi)gen=(%.2f,%.2f)\tpt=%.2f\tE=%.2f\t",gen_b_eta,gen_b_phi,treeVars.GenPart_pt[i],gen_b_E);
     //cout<<"E_GEN-E_RECO="<<gen_b_E-reco_jet.E()<<endl;
     if(fabs(gen_b_E-reco_jet.E()) > 100) continue;
     if( DeltaR(reco_jet_eta,reco_jet_phi,gen_b_eta,gen_b_phi) < DeltaRmax )
     {
       gen_b_match.SetPtEtaPhiM(treeVars.GenPart_pt[i],gen_b_eta,gen_b_phi,treeVars.GenPart_mass[i]);
       return true;
     }
  }
  //printf("########\n");
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

float DeltaRmin(const vector<TLorentzVector> &coll1 , const vector<TLorentzVector> &coll2)
{
  unsigned n1 = coll1.size();
  unsigned n2 = coll2.size();
  if(n1==0)
  {
    cerr<<"> DeltaRmin error: collection 1 is empty"<<endl;
    return -1;
  }
  if(n2==0)
  {
    cerr<<"> DeltaRmin error: collection 2 is empty"<<endl;
    return -1;
  }
  float DR;
  float DRmin = DeltaR( coll1.at(0).Eta() , coll1.at(0).Phi() , coll2.at(0).Eta() , coll2.at(0).Phi() );
  for(unsigned i1=1; i1<n1; ++i1)
    for(unsigned i2=1; i2<n2; ++i2)
    {
      DR=DeltaR( coll1.at(i1).Eta() , coll1.at(i1).Phi() , coll2.at(i2).Eta() , coll2.at(i2).Phi() );
      if(DR<DRmin)
	DRmin=DR;
    }
  return DRmin;
}

/*
int count_bc_quarks(const RawTreeVars& treeVars, int &b_pair, int &b_single, int &c_pair, int &c_single)
{
  int b=0,bbar=0,c=0,cbar=0;
  map<int,vector<TLorentzVector>> v;
  

}
*/
void PrintRecoPhoton(const RawTreeVars& treeVars)
{
  cout<<"Event "<<treeVars.event<<endl;
  /*
  cout<<"\tRECO Tight photon collection -> "<<treeVars.N_TightPh<<" entries"<<endl;
  for(int i=0;i<treeVars.N_TightPh;i++)
    cout<<i<<"\tEta="<<treeVars.TightPh_eta[i]<< "\tPhi="<<treeVars.TightPh_phi[i]<< "\tPt="<<treeVars.TightPh_pt[i]<<endl;
  cout<<endl;

  cout<<"\tRECO Loose photon collection -> "<<treeVars.N_LoosePh<<" entries"<<endl;
  for(int i=0;i<treeVars.N_LoosePh;i++)
    cout<<i<<"\tEta="<<treeVars.LoosePh_eta[i]<< "\tPhi="<<treeVars.LoosePh_phi[i]<< "\tPt="<<treeVars.LoosePh_pt[i]<<endl;
  cout<<endl;
  */
  cout<<"\tRECO Selected photon collection -> "<<treeVars.N_SelectedPh<<" entries"<<endl;
  for(int i=0;i<treeVars.N_SelectedPh;i++)
    cout<<i<<"\tEta="<<treeVars.SelectedPh_eta[i]<< "\tPhi="<<treeVars.SelectedPh_phi[i]<< "\tPt="<<treeVars.SelectedPh_pt[i]<<endl;
  cout<<endl;


}

void PrintRecoJet(const RawTreeVars& treeVars)
{
  cout<<"Event "<<treeVars.event<<"\tRECO PUPPI jet collection -> "<<treeVars.N_Jet<<" entries"<<endl;
  for(int i=0;i<treeVars.N_Jet;i++)
    cout<<i<<"\tEta="<<treeVars.Jet_eta[i]<< "\tPhi="<<treeVars.Jet_phi[i]<< "\tPt="<<treeVars.Jet_pt[i]<<"\tbtagvalue="<<treeVars.Jet_mvav2[i]<<endl;
  cout<<endl;
}
