void RemoveTrainingEvents(TString oldfilename,TString treename)
{
  cout<<"Reading file "<<oldfilename.Data()<<endl;
  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile(oldfilename.Data());
  TTree *oldtree = (TTree*)oldfile->Get(treename.Data());
  Long64_t nentries = oldtree->GetEntries();
  float event;
  float evWeight;
  float HHTagger;
  float ttHTagger;
  float HHTagger_v7bis;
  float ttHTagger_v4;
  float HHTagger_v7bisreducedtrain;
  float ttHTagger_v4reducedtrain;
  int cut_based_ct;

  oldtree->SetBranchAddress("event",&event);
  oldtree->SetBranchAddress("evWeight",&evWeight);
  oldtree->SetBranchAddress("HHTagger",&HHTagger);
  oldtree->SetBranchAddress("ttHTagger",&ttHTagger);
  oldtree->SetBranchAddress("HHTagger_v7bis",&HHTagger_v7bis);
  oldtree->SetBranchAddress("ttHTagger_v4",&ttHTagger_v4);
  oldtree->SetBranchAddress("HHTagger_v7bisreducedtrain",&HHTagger_v7bisreducedtrain);
  oldtree->SetBranchAddress("ttHTagger_v4reducedtrain",&ttHTagger_v4reducedtrain);
  oldtree->SetBranchAddress("cut_based_ct",&cut_based_ct);

  //Create a new file + a clone of old tree in new file
  TString newfilename(oldfilename);
  if(newfilename.Contains("/"))
    newfilename.Insert(newfilename.Last('/')+1,"subset2_");
  else
    newfilename.Insert(0,"subset2_");
  TFile *newfile = new TFile(newfilename.Data(),"recreate");
  TTree *newtree = oldtree->CloneTree(0);
  for (Long64_t i=0;i<nentries; i++) 
  {
    oldtree->GetEntry(i);
    /*
    evWeight *= 1.5;
    if(((int)event)%3!=0)  // filtering unwanted entries.
    {
      if(cut_based_ct!=0)
      {
	cut_based_ct=-1;
      	newtree->Fill();
      }
      else
      {
	if(ttHTagger_v4reducedtrain<-0.2 || HHTagger_v7bisreducedtrain<0.87) 
	  cut_based_ct=-1;
	else
	  if(HHTagger_v7bisreducedtrain>0.87 && HHTagger_v7bisreducedtrain<0.935)
	    cut_based_ct=1;
	  else
	    if(HHTagger_v7bisreducedtrain>0.935)
	      cut_based_ct=0;
	newtree->Fill(); // For those entry passing the filter, store the modified version of the data (and the untouched part too).    
      }
    }
    */
    
    evWeight *= 2.;
    if(((int)event)%2==0)  // filtering unwanted entries.
    {
      if(cut_based_ct==-1)
      	newtree->Fill();
      else
      {
	if(ttHTagger_v4<-0.2 || HHTagger_v7bis<0.5) 
	  cut_based_ct=-1;
	else
	  if(HHTagger_v7bis>0.5 && HHTagger_v7bis<0.93)
	    cut_based_ct=1;
	  else
	    if(HHTagger_v7bis>0.93)
	      cut_based_ct=0;
	newtree->Fill(); // For those entry passing the filter, store the modified version of the data (and the untouched part too).    
      }
    }

  }
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
}

void RemoveTrainingEvents2(TString oldfilename,TString treename)//apply also 1.25 k-factor
{
  cout<<"Reading file "<<oldfilename.Data()<<endl;
  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile(oldfilename.Data());
  TTree *oldtree = (TTree*)oldfile->Get(treename.Data());
  Long64_t nentries = oldtree->GetEntries();
  float event;
  float evWeight;
  float HHTagger;
  float ttHTagger;
  int cut_based_ct;
  oldtree->SetBranchAddress("event",&event);
  oldtree->SetBranchAddress("evWeight",&evWeight);
  oldtree->SetBranchAddress("HHTagger",&HHTagger);
  oldtree->SetBranchAddress("ttHTagger",&ttHTagger);
  oldtree->SetBranchAddress("cut_based_ct",&cut_based_ct);

  //Create a new file + a clone of old tree in new file
  TString newfilename(oldfilename);
  if(newfilename.Contains("/"))
    newfilename.Insert(newfilename.Last('/')+1,"subset_withkfactor_");
  else
    newfilename.Insert(0,"subset_");
  TFile *newfile = new TFile(newfilename.Data(),"recreate");
  TTree *newtree = oldtree->CloneTree(0);
  for (Long64_t i=0;i<nentries; i++) 
  {
    oldtree->GetEntry(i);
    // Do some data modifications. 
    evWeight *= 2.*1.25;

    if(((int)event)%2==0)  // filtering unwanted entries.
    {
      //      if(cut_based_ct==-1)
      //	newtree->Fill();
      //      else
      //      {
	if(ttHTagger<-0.5 || HHTagger<0.5) 
	  cut_based_ct=-1;
	else
	  if(HHTagger>0.5 && HHTagger<0.9)
	    cut_based_ct=1;
	  else
	    if(HHTagger>0.9)
	      cut_based_ct=0;
	newtree->Fill(); // For those entry passing the filter, store the modified version of the data (and the untouched part too).
	//      }
    }
  }
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
}


void iterate_RemoveTrainingEvents()
{
  vector<TString> oldfilename_list = 
  {
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_HHggbb_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_VH_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_bbH_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_ggH_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_gg_withMVA.root",
    //"/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_gjet_withMVA.root",
    //"/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_qcd_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_qqH_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_ttH_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_tt_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_ttgg_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_ttghad_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_ttglep_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_ttgsemilepfromt_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v5/plotTree_ttgsemilepfromtbar_withMVA.root"
  };

  vector<TString> oldfilename_list_LT =
    {
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_HHggbb_withMVA_LT.root",
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_VH_withMVA_LT.root",
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_bbH_withMVA_LT.root",
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_ggH_withMVA_LT.root",
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_gg_withMVA_LT.root",
      //"/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_gjet_withMVA_LT.root",
      //"/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_qcd_withMVA_LT.root",
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_qqH_withMVA_LT.root",
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_ttH_withMVA_LT.root",
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_tt_withMVA_LT.root",
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_ttgg_withMVA_LT.root",
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_ttghad_withMVA_LT.root",
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_ttglep_withMVA_LT.root",
      "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_ttgsemilepfromt_withMVA_LT.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v3/plotTree_ttgsemilepfromtbar_withMVA_LT.root"
    };

  for (unsigned i=0; i< oldfilename_list.size();i++)
    RemoveTrainingEvents( oldfilename_list.at(i) , "all_highMx" );

  //for (unsigned i=0; i< oldfilename_list_LT.size();i++)
  //  RemoveTrainingEvents( oldfilename_list_LT.at(i) , "all_lowMx" );

}
