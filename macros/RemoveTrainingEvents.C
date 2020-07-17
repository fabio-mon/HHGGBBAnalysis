

///////////////////////////////////////////////////////////////////////////////////////////////                                                           
// remove training for new mass categories:  highMx --> mtot>=480 ; low_Mx --> 350<mtot<480
///////////////////////////////////////////////////////////////////////////////////////////////                                                           
void RemoveTrainingEvents3(TString oldfilename,TString treename)
{
  cout<<"Reading file "<<oldfilename.Data()<<endl;
  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile(oldfilename.Data());
  TTree *oldtree = (TTree*)oldfile->Get(treename.Data());
  Long64_t nentries = oldtree->GetEntries();
  float event;
  float evWeight;
  float mtot;
  float HHTagger;
  float HHTagger_v16;
  float HHTagger_v17;
  float HHTagger_v18;
  float HHTagger_v19;
  float HHTagger_v20;
  float HHTagger_v7bis;
  //float HHTagger_v7bisreducedtrain;
  float ttHTagger;
  float ttHTagger_v4;
  //float ttHTagger_v4reducedtrain;
  int cut_based_ct;
  float dibjet_leadbtaglevel;
  float dibjet_subleadbtaglevel;

  oldtree->SetBranchAddress("event",&event);
  oldtree->SetBranchAddress("mtot",&mtot);
  oldtree->SetBranchAddress("dibjet_leadbtaglevel",&dibjet_leadbtaglevel);
  oldtree->SetBranchAddress("dibjet_subleadbtaglevel",&dibjet_subleadbtaglevel);
  oldtree->SetBranchAddress("evWeight",&evWeight);
  //oldtree->SetBranchAddress("HHTagger",&HHTagger);
  oldtree->SetBranchAddress("HHTagger_v7bis",&HHTagger_v7bis);
  //  oldtree->SetBranchAddress("HHTagger_v16",&HHTagger_v16);
  //  oldtree->SetBranchAddress("HHTagger_v17",&HHTagger_v17);
  oldtree->SetBranchAddress("HHTagger_v18",&HHTagger_v18);
  //  oldtree->SetBranchAddress("HHTagger_v19",&HHTagger_v19);
  //  oldtree->SetBranchAddress("HHTagger_v20",&HHTagger_v20);
  //oldtree->SetBranchAddress("HHTagger_v7bisreducedtrain",&HHTagger_v7bisreducedtrain);
  //oldtree->SetBranchAddress("ttHTagger",&ttHTagger);
  oldtree->SetBranchAddress("ttHTagger_v4",&ttHTagger_v4);
  //oldtree->SetBranchAddress("ttHTagger_v4reducedtrain",&ttHTagger_v4reducedtrain);
  oldtree->SetBranchAddress("cut_based_ct",&cut_based_ct);

  //Create a new file + a clone of old tree in new file
  TString newfilename(oldfilename);
  if(newfilename.Contains("/"))
    newfilename.Insert(newfilename.Last('/')+1,"subset1_");
  else
    newfilename.Insert(0,"subset1_");
  TFile *newfile = new TFile(newfilename.Data(),"recreate");
  TTree *newtree_lowMx = oldtree->CloneTree(0);
  newtree_lowMx->SetName("all_lowMx");
  newtree_lowMx->SetTitle("all_lowMx");
  TTree *newtree_highMx = oldtree->CloneTree(0);
  newtree_highMx->SetName("all_highMx");
  newtree_highMx->SetTitle("all_highMx");
  for (Long64_t i=0;i<nentries; i++) 
  {
    oldtree->GetEntry(i);

    /*
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //for training events reduced to 1/3
    ///////////////////////////////////////////////////////////////////////////////////////////////
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
    ///////////////////////////////////////////////////////////////////////////////////////////////
    */

    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // requiring 2 btag loose as preselection 
    ///////////////////////////////////////////////////////////////////////////////////////////////    
    evWeight *= 2.;
    if(oldfilename.Contains("plotTree_gg_withMVA.root"))
       evWeight /= 1.1;
    //if(oldfilename.Contains("plotTree_HHggbb_withMVA.root"))
    //   evWeight /= 1.078223;
    if(((int)event)%2==0)  // filtering unwanted entries.
    {
      cut_based_ct=-1;
      if(!(dibjet_subleadbtaglevel>=4 && dibjet_leadbtaglevel>=4))
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
      if(mtot>480)
	newtree_highMx->Fill();
      else
	if(mtot>350 && mtot<=480)
	  newtree_lowMx->Fill();
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /*
    evWeight *= 2.;
    if(oldfilename.Contains("plotTree_gg_withMVA.root"))
       evWeight *= 1.1;
    if(((int)event)%2!=0)  //In this case revert train and test!!!!!!!
    {
      cut_based_ct=-1;
      if(!(dibjet_subleadbtaglevel>=4 && dibjet_leadbtaglevel>=4))
	cut_based_ct=-1;
      else
      {
	if(ttHTagger_v4<-0.2 || HHTagger_v20<0.5) 
	  cut_based_ct=-1;
	else
	  if(HHTagger_v20>0.5 && HHTagger_v20<0.93)
	    cut_based_ct=1;
	  else
	    if(HHTagger_v20>0.93)
	      cut_based_ct=0;
      }
      if(mtot>480)
	newtree_highMx->Fill();
      else
	if(mtot>350 && mtot<=480)
	  newtree_lowMx->Fill();
    }
    */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  }
  newtree_highMx->AutoSave();
  newtree_lowMx->AutoSave();
  delete oldfile;
  delete newfile;
}

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
  float HHTagger_v16;
  float HHTagger_v17;
  float HHTagger_v7bis;
  //float HHTagger_v7bisreducedtrain;
  float ttHTagger;
  float ttHTagger_v4;
  //float ttHTagger_v4reducedtrain;
  int cut_based_ct;
  int dibjet_leadbtaglevel;
  int dibjet_subleadbtaglevel;

  oldtree->SetBranchAddress("event",&event);
  oldtree->SetBranchAddress("dibjet_leadbtaglevel",&dibjet_leadbtaglevel);
  oldtree->SetBranchAddress("dibjet_subleadbtaglevel",&dibjet_subleadbtaglevel);
  oldtree->SetBranchAddress("evWeight",&evWeight);
  //oldtree->SetBranchAddress("HHTagger",&HHTagger);
  oldtree->SetBranchAddress("HHTagger_v7bis",&HHTagger_v7bis);
  oldtree->SetBranchAddress("HHTagger_v16",&HHTagger_v16);
  oldtree->SetBranchAddress("HHTagger_v17",&HHTagger_v17);
  //oldtree->SetBranchAddress("HHTagger_v7bisreducedtrain",&HHTagger_v7bisreducedtrain);
  //oldtree->SetBranchAddress("ttHTagger",&ttHTagger);
  oldtree->SetBranchAddress("ttHTagger_v4",&ttHTagger_v4);
  //oldtree->SetBranchAddress("ttHTagger_v4reducedtrain",&ttHTagger_v4reducedtrain);
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
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //for training events reduced to 1/3
    ///////////////////////////////////////////////////////////////////////////////////////////////
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
    ///////////////////////////////////////////////////////////////////////////////////////////////
    */

    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // requiring 2 btag loose as preselection 
    ///////////////////////////////////////////////////////////////////////////////////////////////    
    evWeight *= 2.;
    if(((int)event)%2==0)  // filtering unwanted entries.
    {
      if(!(dibjet_subleadbtaglevel>=4 && dibjet_leadbtaglevel>=4))
      {
	cut_based_ct=-1;
      	newtree->Fill();
      }
      else
      {
	if(ttHTagger_v4<-0.2 || HHTagger_v17<0.5) 
	  cut_based_ct=-1;
	else
	  if(HHTagger_v17>0.5 && HHTagger_v17<0.93)
	    cut_based_ct=1;
	  else
	    if(HHTagger_v7bis>0.93)
	      cut_based_ct=0;
	newtree->Fill(); // For those entry passing the filter, store the modified version of the data (and the untouched part too).    
      }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    /*
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // default
    ///////////////////////////////////////////////////////////////////////////////////////////////    
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
    ///////////////////////////////////////////////////////////////////////////////////////////////
    */

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
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_HHggbb_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_VH_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_bbH_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_ggH_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_gg_withMVA.root",
    //"/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_gjet_withMVA.root",
    //"/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_qcd_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_qqH_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_ttH_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_tt_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_ttgg_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_ttghad_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_ttglep_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_ttgsemilepfromt_withMVA.root",
    "/afs/cern.ch/user/f/fmonti/public/HHGGBBAnalysis/output/PhotonLoose_mbb_80-200_MTD_v7/plotTree_ttgsemilepfromtbar_withMVA.root"
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
    RemoveTrainingEvents3( oldfilename_list.at(i) , "all_highMx" );

  //for (unsigned i=0; i< oldfilename_list_LT.size();i++)
  //  RemoveTrainingEvents( oldfilename_list_LT.at(i) , "all_lowMx" );

}
