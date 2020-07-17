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

double functionGF(double kl, double kt, double c2, double cg, double c2g, std::array<double, 15> const &A)
{
  // this can be extended to 5D coefficients; currently c2, cg, c2g are unused
  // return ( A1*pow(kt,4) + A3*pow(kt,2)*pow(kl,2) + A7*kl*pow(kt,3) );
  return ( A[0]*pow(kt,4) + A[1]*pow(c2,2) + (A[2]*pow(kt,2) + A[3]*pow(cg,2))*pow(kl,2) + A[4]*pow(c2g,2) + ( A[5]*c2 + A[6]*kt*kl )*pow(kt,2) + (A[7]*kt*kl + A[8]*cg*kl )*c2 + A[9]*c2*c2g + (A[10]*cg*kl + A[11]*c2g)*pow(kt,2)+ (A[12]*kl*cg + A[13]*c2g )*kt*kl + A[14]*cg*c2g*kl );
}

///////////////////////////////////////////////////////////////////////////////////////////////
// remove training for new mass categories:  highMx --> mtot>=480 ; low_Mx --> 350<mtot<480
///////////////////////////////////////////////////////////////////////////////////////////////                                                       

int main(int argc, char** argv)
{
  TString oldfilename(argv[1]);
  TString treename(argv[2]);
  TString scenario(argv[3]);
  float kl=1;
  if(argc>4)
    float kl = atof(argv[4]);

  std::array<double, 15> A_integralXS = 
  {
    2.100318379,
    10.2,
    0.287259045,
      0.098882779,
    1.321736614,
    -8.42431259,
    -1.388017366,
    2.8,
    0.518124457,
    -2.163473227,
    -0.550668596,
    5.871490593,
    0.296671491,
    -1.172793054,
    0.653429812
  };

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
  float HHTagger_v18LT;
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
      oldtree->SetBranchAddress("HHTagger_v18LT",&HHTagger_v18LT);
    }

  //Create a new file + a clone of old tree in new file
  TString newfilename(oldfilename);
  if(newfilename.Contains("/"))
    newfilename.Insert(newfilename.Last('/')+1,scenario+"_subset_");
  else
    newfilename.Insert(0,scenario+"_subset_");
  TFile *newfile = new TFile(newfilename.Data(),"recreate");
  TTree *newtree = oldtree->CloneTree(0);
  newtree->SetName("TCVARS");
  newtree->SetTitle("TCVARS");

  //reweight cross sections --> NO: ALREADY DONE IN HHGGBBselector.exe
  //if(argc>4)
  //  cout << "sigma_BSM / sigma_SM = " << functionGF(kl, 1, 0, 0, 0, A_integralXS) << endl;

  float HHreweight = 1.;
  //if(argc>4)
  //  HHreweight = 1./functionGF(kl, 1, 0, 0, 0, A_integralXS);//1.078223;
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
	  //cout<<ttHTagger_v4<<"\t"<<HHTagger_v18LT<<endl;
	  cut_based_ct=-1;
	  if(!(dibjet_subleadbtaglevel>=4 && dibjet_leadbtaglevel>=4 && dibjet_subleadEta<=2.4 && dibjet_leadEta<=2.4))
	    cut_based_ct=-1;
	  else
          {
	    if(ttHTagger<-0.3 || HHTagger_v18LT<0.5) 
	      cut_based_ct=-1;
	    else
	      if(HHTagger_v18LT>0.5 && HHTagger_v18LT<0.9)
		cut_based_ct=1;
	      else
		if(HHTagger_v18LT>0.9)
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
