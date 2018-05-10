void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction, const bool& verbosity);

void draw_invmass()
{
   //gStyle->SetOptStat(0);  
   gStyle->SetOptTitle(0);
   TLatex *title = new TLatex();
   title->SetNDC(); 

   TChain *ch_delphes = new TChain("plotTree","plotTree");
   TChain *ch_fullsim_phase2 = new TChain("plotTree","plotTree");
   TChain *ch_fullsim_phase1 = new TChain("bbggSelectionTree","bbggSelectionTree");

   ch_delphes->Add("/eos/user/f/fmonti/HHGGBB/data/DelphesDump4/plotTree/DelphesDump4_v8.root");
   ch_fullsim_phase2->Add("/eos/user/f/fmonti/HHGGBB/data/FullsimPU200_Dump2/plotTree/FullsimPU200_Dump2_v8.root");
   ch_fullsim_phase1->Add("/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/Mar292018_ForUpgrade_ttHBDT/Signal/Hadd/output_GluGluToHHTo2B2G_node_10_13TeV-madgraph_0.root");

   TF1* fitfunc = new TF1("gaus_fitfunc","gaus(0)");
   float* vals = new float[4];
   float max, min;



   TCanvas *c1 = new TCanvas("c_delphes","c_delphes",400,400);
   c1->cd();
   ch_delphes->Draw("dipho_mass>>h_delphes(600,110,140)","fabs(dipho_leadEta)>1.6 && fabs(dipho_subleadEta)>1.6 && dipho_leadPt>40 && dipho_subleadPt>30");
   TH1F* h_delphes = (TH1F*)gDirectory->Get("h_delphes");
   h_delphes->SetTitle("Delphes, PhaseII");
   h_delphes->SetLineColor(1);
   h_delphes->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
   title->DrawLatex(0.1,0.93,"Delphes, PhaseII");
   fitfunc->SetParameters(300,125,2);
   fitfunc->SetRange( h_delphes->GetMean()-h_delphes->GetRMS()*0.8 , h_delphes->GetMean()+h_delphes->GetRMS()*0.8 );
   h_delphes->Fit(fitfunc,"R");
   FindSmallestInterval(vals,h_delphes,0.68,true);
    min = vals[2];
    max = vals[3]; 
   cout<<"\n\nh_delphes"<<endl;
   cout<<"Mean="<<h_delphes->GetMean()<<" +/- "<<h_delphes->GetMeanError()<<endl;
   cout<<"\tRMS="<<h_delphes->GetRMS()<<endl;
   cout<<"\tsigma_fit="<<fitfunc->GetParameter(2)<<endl;
   cout<<"\tSmallestInt="<<(max-min)*0.5<<endl;
   //   c1->Print("Delphes_PhaseII.pdf");

   TCanvas *c2 = new TCanvas("c_fullsim_phase2","c_fullsim_phase2",400,400);
   c2->cd();
   ch_fullsim_phase2->Draw("dipho_mass>>h_fullsim_phase2(600,110,140)","fabs(dipho_leadEta)>1.6 && fabs(dipho_subleadEta)>1.6 && dipho_leadPt>40 && dipho_subleadPt>30");
   TH1F* h_fullsim_phase2 = (TH1F*)gDirectory->Get("h_fullsim_phase2");
   h_fullsim_phase2->SetTitle("Fullsim, PhaseII");
   h_fullsim_phase2->SetLineColor(4);
   h_fullsim_phase2->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
   title->DrawLatex(0.1,0.93,"Fullsim, PhaseII");
   fitfunc->SetRange( h_fullsim_phase2->GetMean()-h_fullsim_phase2->GetRMS()*0.8 , h_fullsim_phase2->GetMean()+h_fullsim_phase2->GetRMS()*0.8 );
   h_fullsim_phase2->Fit(fitfunc,"R");
   FindSmallestInterval(vals,h_fullsim_phase2,0.68,true);
    min = vals[2];
    max = vals[3]; 
   cout<<"\n\nh_fullsim_phase2"<<endl;
   cout<<"Mean="<<h_fullsim_phase2->GetMean()<<" +/- "<<h_fullsim_phase2->GetMeanError()<<endl;
   cout<<"\tRMS="<<h_fullsim_phase2->GetRMS()<<endl;
   cout<<"\tsigma_fit="<<fitfunc->GetParameter(2)<<endl;
   cout<<"\tSmallestInt="<<(max-min)*0.5<<endl;
   //   c2->Print("Fullsim_PhaseII.pdf");

   TCanvas *c3 = new TCanvas("c_fullsim_phase1","c_fullsim_phase1",400,400);
   c3->cd();
   ch_fullsim_phase1->Draw("diphotonCandidate.M()>>h_fullsim_phase1(600,110,140)","fabs(leadingPhoton.Eta())>1.6 && fabs(subleadingPhoton.Eta())>1.6 && leadingPhoton.Pt()>40 && subleadingPhoton.Pt()>30");
   TH1F* h_fullsim_phase1 = (TH1F*)gDirectory->Get("h_fullsim_phase1");
   h_fullsim_phase1->SetTitle("Fullsim, PhaseI");
   h_fullsim_phase1->SetLineColor(3);
   h_fullsim_phase1->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
   title->DrawLatex(0.1,0.93,"Fullsim, PhaseI");
   fitfunc->SetRange( h_fullsim_phase1->GetMean()-h_fullsim_phase1->GetRMS()*0.8 , h_fullsim_phase1->GetMean()+h_fullsim_phase1->GetRMS()*0.8 );
   h_fullsim_phase1->Fit(fitfunc,"R");
   //   c3->Print("Fullsim_PhaseII.pdf");
   FindSmallestInterval(vals,h_fullsim_phase1,0.68,true);
    min = vals[2];
    max = vals[3]; 
   cout<<"\n\nh_fullsim_phase1"<<endl;
   cout<<"Mean="<<h_fullsim_phase1->GetMean()<<" +/- "<<h_fullsim_phase1->GetMeanError()<<endl;
   cout<<"\tRMS="<<h_fullsim_phase1->GetRMS()<<endl;
   cout<<"\tsigma_fit="<<fitfunc->GetParameter(2)<<endl;
   cout<<"\tSmallestInt="<<(max-min)*0.5<<endl;

   TCanvas *c4 = new TCanvas("c_dipho_mass_comparison","c_dipho_mass_comparison");
   c4->cd();
   h_fullsim_phase1->DrawNormalized("E1");
   h_fullsim_phase2->DrawNormalized("E1SAME");
   h_delphes->DrawNormalized("E1SAME");

   

}


void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction, const bool& verbosity)
{
  float integralMax = fraction * histo->Integral();
  
  int N = histo -> GetNbinsX();
  int M1 = 0;
  int M2 = 0;
  for(int bin1 = 0; bin1 < N; ++bin1)
  {
    if( histo->GetBinContent(bin1+1) > 0. && M1 == 0 ) M1 = bin1-1;
    if( histo->GetBinContent(bin1+1) > 0. ) M2 = bin1+2;
  }
  
  std::map<int,float> binCenters;
  std::map<int,float> binContents;
  std::map<int,float> binIntegrals;
  for(int bin1 = M1; bin1 < M2; ++bin1)
  {
    binCenters[bin1] = histo->GetBinCenter(bin1+1);
    binContents[bin1] = histo->GetBinContent(bin1+1);
    
    for(int bin2 = M1; bin2 <= bin1; ++bin2)
      binIntegrals[bin1] += binContents[bin2];
  }
  
  float min = 0.;
  float max = 0.;
  float delta = 999999.;
  for(int bin1 = M1; bin1 < M2; ++bin1)
  {
    for(int bin2 = bin1+1; bin2 < M2; ++bin2)
    {
      if( (binIntegrals[bin2]-binIntegrals[bin1]) < integralMax ) continue;
      
      float tmpMin = histo -> GetBinCenter(bin1+1);
      float tmpMax = histo -> GetBinCenter(bin2+1);
      
      if( (tmpMax-tmpMin) < delta )
      {
        delta = (tmpMax - tmpMin);
        min = tmpMin;
        max = tmpMax;
      }
      
      break;
    }
  }
  
  TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
  for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
  {
    if( smallHisto->GetBinCenter(bin) < min )
      smallHisto -> SetBinContent(bin,0);
    
    if( smallHisto->GetBinCenter(bin) > max )
      smallHisto -> SetBinContent(bin,0);
  }
  smallHisto -> SetFillColor(kYellow);
  
  float mean = smallHisto -> GetMean();
  float meanErr = smallHisto -> GetMeanError();  
  
  ret[0] = mean;
  ret[1] = meanErr;
  ret[2] = min;
  ret[3] = max;
}
