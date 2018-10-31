#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////
// remove training for new mass categories:  highMx --> mtot>=480 ; low_Mx --> 350<mtot<480
///////////////////////////////////////////////////////////////////////////////////////////////                                                       

int main(int argc, char** argv)
{
  TString oldfilename(argv[1]);
  TString treename(argv[2]);
  TString scenario(argv[3]);

  cout<<"Reading file "<<oldfilename.Data()<<endl;
  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile(oldfilename.Data());
  TTree *oldtree = (TTree*)oldfile->Get(treename.Data());
  Long64_t nentries = oldtree->GetEntries();
  float event;
  float evWeight;
  float mtot;
  float HHTagger_v2;
  float HHTagger_v18;
  float ttHTagger;
  float ttHTagger_v4;
  int cut_based_ct;
  float dibjet_leadbtaglevel;
  float dibjet_subleadbtaglevel;
  float dibjet_leadEta;
  float dibjet_subleadEta;

  oldtree->SetBranchAddress("event",&event);
  oldtree->SetBranchAddress("evWeight",&evWeight);
  oldtree->SetBranchAddress("mtot",&mtot);
  oldtree->SetBranchAddress("cut_based_ct",&cut_based_ct);
  oldtree->SetBranchAddress("dibjet_leadbtaglevel",&dibjet_leadbtaglevel);
  oldtree->SetBranchAddress("dibjet_subleadbtaglevel",&dibjet_subleadbtaglevel);
  oldtree->SetBranchAddress("dibjet_leadEta",&dibjet_leadEta);
  oldtree->SetBranchAddress("dibjet_subleadEta",&dibjet_subleadEta);

  if(treename=="all_highMx")
  {
    oldtree->SetBranchAddress("ttHTagger_v4",&ttHTagger_v4);
    oldtree->SetBranchAddress("HHTagger_v18",&HHTagger_v18);
  }
  else 
    if(treename=="all_lowMx")
    {
      oldtree->SetBranchAddress("ttHTagger",&ttHTagger);
      oldtree->SetBranchAddress("HHTagger_v2",&HHTagger_v2);
    }

  //Create a new file + a clone of old tree in new file
  TString newfilename(oldfilename);
  if(newfilename.Contains("/"))
    newfilename.Insert(newfilename.Last('/')+1,scenario+"_subset_");
  else
    newfilename.Insert(0,scenario+"_subset_");
  TFile *newfile = new TFile(newfilename.Data(),"recreate");
  TTree *newtree = oldtree->CloneTree(0);

  //reweight cross sections
  float HHreweight = 1;//1.078223;
  float ttHreweight = 1./0.75;
  float ggreweight = 0;
  if(scenario=="pessimistic")
    ggreweight = 0.916667; // (weight*1.1)/(1.1/1.2) = (weight*1.1)/0.916667 
  else
    if(scenario=="intermediate")
      ggreweight = 1.;
    else
      if(scenario=="optimistic")
	ggreweight = 1.1;

  //check if i have to reweight
  bool is_gg  = oldfilename.Contains("_gg_");
  bool is_HH  = oldfilename.Contains("_HHggbb_");
  bool is_ttH = oldfilename.Contains("_ttH_");

  //loop over events
  for (Long64_t i=0;i<nentries; i++) 
  {
    oldtree->GetEntry(i);
    evWeight *= 2.;
    if(is_gg)
      evWeight /= ggreweight;
    if(is_HH)
      evWeight /= HHreweight;
    if(is_ttH)
      evWeight /= ttHreweight;

    if(((int)event)%2==0)  // filtering unwanted entries.
    {
      if(treename=="all_highMx")
      {
	cut_based_ct=-1;
	if(!(dibjet_subleadbtaglevel>=4 && dibjet_leadbtaglevel>=4 && dibjet_subleadEta<=2.4 && dibjet_leadEta<=2.4))
	  cut_based_ct=-1;
	else
        {
	  if(ttHTagger_v4<-0.3 || HHTagger_v18<0.65) 
	    cut_based_ct=-1;
	  else
	    if(HHTagger_v18>0.65 && HHTagger_v18<0.93)
	      cut_based_ct=1;
	    else
	      if(HHTagger_v18>0.93)
		cut_based_ct=0;
	}
      }
      else

	if(treename=="all_lowMx")
        {
	  cut_based_ct=-1;
	  if(!((dibjet_subleadbtaglevel>=5 || dibjet_leadbtaglevel>=5) && dibjet_subleadEta<=2.4 && dibjet_leadEta<=2.4))
	    cut_based_ct=-1;
	  else
          {
	    if(ttHTagger_v4<-0.3 || HHTagger_v2<0.5) 
	      cut_based_ct=-1;
	    else
	      if(HHTagger_v2>0.5 && HHTagger_v2<0.9)
		cut_based_ct=1;
	      else
		if(HHTagger_v2>0.9)
		  cut_based_ct=0;
	  }
	}
      newtree->Fill();
    }

  }
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
}
