/*** find effective sigma ***/
void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction, const bool& verbosity = false)
{
  float integralMax = fraction * histo->Integral();
  
  // find first and last bin with non-null content
  int N = histo -> GetNbinsX();
  int M1 = 1;
  int M2 = 1;
  for(int bin1 = 1; bin1 <= N; ++bin1)
  {
    if( histo->GetBinContent(bin1) > 0. && M1 == 1 ) M1 = bin1;
    if( histo->GetBinContent(bin1) > 0. )            M2 = bin1;
  }

  std::map<int,float> binCenters;
  std::map<int,float> binContents;
  std::map<int,float> binIntegrals;
  for(int bin1 = M1; bin1 <= M2; ++bin1)
  {
    if( histo->GetBinContent(bin1) == 0 ) continue;

    binCenters[bin1] = histo->GetBinCenter(bin1);
    binContents[bin1] = histo->GetBinContent(bin1);

    for(int bin2 = M1; bin2 <= bin1; ++bin2)
      binIntegrals[bin1] += histo->GetBinContent(bin2);
  }
  
  float min = 0.;
  float max = 0.;
  float delta = 999999.;
  for(std::map<int,float>::const_iterator mapIt1 = binIntegrals.begin(); mapIt1 != binIntegrals.end(); ++mapIt1)
  {
    for(std::map<int,float>::const_iterator mapIt2 = ++binIntegrals.begin(); mapIt2 != binIntegrals.end(); ++mapIt2)
    {
      if( (mapIt2->second-mapIt1->second) < integralMax ) continue;

      float tmpMin = binCenters[mapIt1->first];
      float tmpMax = binCenters[mapIt2->first];

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
  float RMS = smallHisto -> GetRMS();
  float RMSErr = smallHisto -> GetRMSError();

  ret[0] = mean;
  ret[1] = meanErr;
  ret[2] = RMS;
  ret[3] = RMSErr;
  ret[4] = min;
  ret[5] = max;
}



double ComputeSignificance(TH1F* h_sig, TH1F* h_bkg, const int& mode)
{
  double significance = 0.;

  for(int bin = 1; bin <= h_sig->GetNbinsX(); ++bin)
  {
    float S = h_sig->GetBinContent(bin);
    float B = h_bkg->GetBinContent(bin);

    if( B > 0 )
    {
      if( mode == 1 ) significance += S*S / B;
      if( mode == 2 ) significance += 2.*((S+B)*log(1.+S/B)-S);
    }
  }

  return sqrt(significance);
}



float lumin = 3000;

void drawFinalPlot(std::string folder,std::string plotfolder,float effSigma=-1)
{
  std::map<std::string,TTree*> trees;
  
  TFile* inFile;
  
  inFile = TFile::Open(Form("%s/plotTree_HHggbb_withMVA_LT.root",folder.c_str()),"READ");
  trees["LT_sig"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_HHggbb_withMVA.root",folder.c_str()),"READ");
  trees["HT_sig"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_ggH_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_res_ggH"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_ggH_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_res_ggH"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_qqH_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_res_qqH"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_qqH_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_res_qqH"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_VH_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_res_VH"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_VH_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_res_VH"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_ttH_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_res_ttH"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_ttH_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_res_ttH"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_bbH_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_res_bbH"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_bbH_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_res_bbH"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_gg_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_nonres_gg"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_gg_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_nonres_gg"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_tt_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_nonres_tt"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_tt_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_nonres_tt"] = (TTree*)( inFile->Get("all_highMx") );

  inFile = TFile::Open(Form("%s/plotTree_ttgg_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_nonres_ttgg"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_ttgg_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_nonres_ttgg"] = (TTree*)( inFile->Get("all_highMx") );

  inFile = TFile::Open(Form("%s/plotTree_ttghad_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_nonres_ttghad"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_ttghad_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_nonres_ttghad"] = (TTree*)( inFile->Get("all_highMx") );

  inFile = TFile::Open(Form("%s/plotTree_ttglep_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_nonres_ttglep"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_ttglep_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_nonres_ttglep"] = (TTree*)( inFile->Get("all_highMx") );

  inFile = TFile::Open(Form("%s/plotTree_ttgsemilepfromt_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_nonres_ttgsemilepfromt"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_ttgsemilepfromt_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_nonres_ttgsemilepfromt"] = (TTree*)( inFile->Get("all_highMx") );

  inFile = TFile::Open(Form("%s/plotTree_ttgsemilepfromtbar_LT_withMVA.root",folder.c_str()),"READ");
  trees["LT_bkg_nonres_ttgsemilepfromtbar"] = (TTree*)( inFile->Get("all_lowMx") );
  inFile = TFile::Open(Form("%s/plotTree_ttgsemilepfromtbar_withMVA.root",folder.c_str()),"READ");
  trees["HT_bkg_nonres_ttgsemilepfromtbar"] = (TTree*)( inFile->Get("all_highMx") );
  
  // inFile = TFile::Open(Form("%s/plotTree_gjet_LT_withMVA.root",folder.c_str()),"READ");
  // trees["LT_bkg_nonres_gjet"] = (TTree*)( inFile->Get("all_lowMx") );
  // trees["HT_bkg_nonres_gjet"] = (TTree*)( inFile->Get("all_highMx") );
  
  // inFile = TFile::Open(Form("%s/plotTree_qcd.root",folder.c_str()),"READ");
  // trees["LT_bkg_nonres_qcd"] = (TTree*)( inFile->Get("all_lowMx") );
  // trees["HT_bkg_nonres_qcd"] = (TTree*)( inFile->Get("all_highMx") );
  
  
  std::vector<std::string> labels;
  labels.push_back("HT");
  //labels.push_back("LT");
  
  std::vector<std::string> cats;
  //high purity
  //cats.push_back("*2.* (cut_based_ct >= 0 && event%2==0 && mjj<190 && ttHTagger>-0.5 && HHTagger>=0.9)");
  //medium purity
  //cats.push_back("*2.* (cut_based_ct >= 0 && event%2==0 && mjj<190 && ttHTagger>-0.5 && HHTagger>0.5 && HHTagger<0.9)");
  //high purity
  cats.push_back("*2.* (event%2==0 && dibjet_leadbtaglevel>=4 && dibjet_subleadbtaglevel>=4 && abs(dibjet_subleadEta)<=2.4 && abs(dibjet_leadEta)<=2.4 && ttHTagger_v4>-0.3 && HHTagger_v18>0.93 && mjj>90 && mjj<190 && mtot>350 && mtot<480)");
  cats.push_back("*2.* (event%2==0 && dibjet_leadbtaglevel>=4 && dibjet_subleadbtaglevel>=4 && abs(dibjet_subleadEta)<=2.4 && abs(dibjet_leadEta)<=2.4 && ttHTagger_v4>-0.3 && HHTagger_v18>0.65 && HHTagger_v18<0.93 && mjj>90 && mjj<190 && mtot>350 && mtot<480)");
  
  for(auto label : labels)
  {
    int catIt = 0;
    for(auto cat : cats)
    {
      TH1F* h1_bkg_all_mgg      = new TH1F(Form("h1_bkg_all_mgg__%s_%d",label.c_str(),catIt),    "",  40,100.,180.);
      TH1F* h1_bkg_all_fit_mgg  = new TH1F(Form("h1_bkg_all_fit_mgg__%s_%d",label.c_str(),catIt),"",  40,100.,180.);
      TH1F* h1_bkg_nonres_mgg   = new TH1F(Form("h1_bkg_nonres_mgg__%s_%d",label.c_str(),catIt), "",  40,100.,180.);
      TH1F* h1_bkg_res_mgg      = new TH1F(Form("h1_bkg_res_mgg__%s_%d",label.c_str(),catIt),    "",  40,100.,180.);
      TH1F* h1_sig_mgg          = new TH1F(Form("h1_sig_mgg__%s_%d",label.c_str(),catIt),        "",  40,100.,180.);
      TH1F* h1_sig_mgg_fine     = new TH1F(Form("h1_sig_mgg_fine__%s_%d",label.c_str(),catIt),   "", 250,100.,150.);
      
      for(auto tree : trees)
      {
        std::size_t found_label = (tree.first).find(label);
        if( found_label == std::string::npos ) continue;
	cout<<tree.first<<endl;
        std::size_t found_sig = (tree.first).find("sig");
        if( found_sig != std::string::npos )
        {
          (tree.second) -> Draw(Form("mgg >>+ h1_sig_mgg__%s_%d",     label.c_str(),catIt),Form("evWeight * %f %s",lumin,cat.c_str()),"goff");
          (tree.second) -> Draw(Form("mgg >>+ h1_sig_mgg_fine__%s_%d",label.c_str(),catIt),Form("evWeight * %f %s",lumin,cat.c_str()),"goff");
        }
        
        std::size_t found_bkg_res = (tree.first).find("bkg_res");
        if( found_bkg_res != std::string::npos )
        {
          (tree.second) -> Draw(Form("mgg >>+ h1_bkg_res_mgg__%s_%d",label.c_str(),catIt),Form("evWeight * %f %s",lumin,cat.c_str()),"goff");
        }
 
       std::size_t found_bkg_nonres = (tree.first).find("bkg_nonres");
        if( found_bkg_nonres != std::string::npos )
        {
          (tree.second) -> Draw(Form("mgg >>+ h1_bkg_nonres_mgg__%s_%d",label.c_str(),catIt),Form("evWeight * %f %s",lumin,cat.c_str()),"goff");
          (tree.second) -> Draw(Form("mgg >>+ h1_bkg_all_mgg__%s_%d",   label.c_str(),catIt),Form("evWeight * %f %s",lumin,cat.c_str()),"goff");
        }
      }

      float n_bkg_nonres = h1_bkg_nonres_mgg->Integral();
      float n_sig = h1_sig_mgg->Integral();
      float n_bkg_res = h1_bkg_res_mgg->Integral();
      float* vals = new float[6];
      FindSmallestInterval(vals,h1_sig_mgg_fine,0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      
      TF1* func_bkg_res_mgg = new TF1(Form("func_bkg_res_mgg__%s_%d",label.c_str(),catIt),"[0]*1./([2]*sqrt(3.14159))*exp(-1.*(x-[1])*(x-[1])/(2.*[2]*[2]))",100.,150.);
      func_bkg_res_mgg -> SetParameters(n_bkg_res,125.,effSigma);
      func_bkg_res_mgg -> SetLineColor(65);
      func_bkg_res_mgg -> SetLineWidth(2);
      func_bkg_res_mgg -> SetLineStyle(7);
      //func_bkg_res_mgg -> Draw("same");
      
      TF1* func_sig_mgg = new TF1(Form("func_sig_mgg__%s_%d",label.c_str(),catIt),"[0]*1./([2]*sqrt(3.14159))*exp(-1.*(x-[1])*(x-[1])/(2.*[2]*[2]))",100.,150.);
      func_sig_mgg -> SetParameters(n_sig,125.,effSigma);
      func_sig_mgg -> SetLineColor(50);
      func_sig_mgg -> SetLineWidth(2);
      func_sig_mgg -> SetLineStyle(7);
      //func_sig_mgg -> Draw("same");

      if(effSigma>0)
	{
	  for(int bin = 1; bin <= h1_sig_mgg->GetNbinsX(); ++bin)
	    {
	      float low = h1_sig_mgg->GetBinLowEdge(bin);
	      float hig = h1_sig_mgg->GetBinLowEdge(bin) + h1_sig_mgg->GetBinWidth(bin);
	      h1_bkg_res_mgg -> SetBinContent(bin,func_bkg_res_mgg->Integral(low,hig));
	      h1_sig_mgg -> SetBinContent(bin,func_sig_mgg->Integral(low,hig));
	    }
	}

      h1_bkg_all_mgg->Add(h1_bkg_res_mgg);
      


      TCanvas *cc = new TCanvas();
      h1_bkg_nonres_mgg -> SetLineColor(kGreen);
      //h1_bkg_nonres_mgg -> Draw();

      TCanvas* c1 = new TCanvas(Form("%s_cat%d",label.c_str(),catIt),Form("%s_cat%d",label.c_str(),catIt));      
      TH1F* hPad = (TH1F*)( gPad->DrawFrame(100.,0.,150.,h1_bkg_all_mgg->GetMaximum()*2.));
      hPad -> SetTitle(Form(";m_{#gamma#gamma} [GeV];events / %.1f GeV",h1_bkg_all_mgg->GetBinWidth(1)));
      hPad -> Draw();
      h1_bkg_all_mgg -> Draw("P,same");
     
      h1_bkg_res_mgg -> SetLineColor(65);
      h1_bkg_res_mgg -> SetLineWidth(2);
      h1_bkg_res_mgg -> Draw("hist,same");
      
      h1_sig_mgg -> SetLineColor(50);
      h1_sig_mgg -> SetLineWidth(2);
      h1_sig_mgg -> Draw("hist,same");

      TF1* fit_bkg_nonres_mgg = new TF1(Form("fit_bkg_nonres_mgg__%s_%d",label.c_str(),catIt),"[0]*exp(-[1]*x)",105.,145.);
      fit_bkg_nonres_mgg->SetParameters( 7.91464e+02,2.09758e-02);
      h1_bkg_nonres_mgg -> Fit(fit_bkg_nonres_mgg,"LQNRS+");
      fit_bkg_nonres_mgg -> SetLineColor(kBlack);
      fit_bkg_nonres_mgg -> SetLineWidth(2);
      fit_bkg_nonres_mgg -> SetLineStyle(7);
      fit_bkg_nonres_mgg -> Draw("same");
      
      for(int bin = 1; bin <= h1_bkg_all_mgg->GetNbinsX(); ++bin)
      {
        double binCenter = h1_bkg_all_mgg -> GetBinCenter(bin);
        double binContent = h1_bkg_res_mgg -> GetBinContent(bin);
        double fitContent = fit_bkg_nonres_mgg -> Eval(binCenter);
        h1_bkg_all_fit_mgg -> SetBinContent(bin,fitContent+binContent);
      }

      double significance = ComputeSignificance(h1_sig_mgg,h1_bkg_all_fit_mgg,2);
      
      std::cout << "\n\n " << label << " - cat. " << catIt << std::endl;
      std::cout << ">>>>>> bkg_all:      " << std::fixed << std::setprecision(0) << std::setw(6) << h1_bkg_all_mgg->Integral()     << " ev." << std::endl;
      std::cout << ">>>>>> bkg_nonres:   " << std::fixed << std::setprecision(0) << std::setw(6) << fit_bkg_nonres_mgg->Eval(125.) << " ev. / GeV" << std::endl;
      std::cout << ">>>>>> bkg_res:      " << std::fixed << std::setprecision(1) << std::setw(6) << n_bkg_res                      << " ev."       << std::endl;
      std::cout << ">>>>>> sig:          " << std::fixed << std::setprecision(1) << std::setw(6) << n_sig                          << " ev."       << std::endl;
      std::cout << ">>>>>> effSigma:     " << std::fixed << std::setprecision(2) << std::setw(6) << sigma                          << " GeV"       << std::endl;
      std::cout << ">>>>>> significance: " << std::fixed << std::setprecision(3) << std::setw(6) << significance                                   << std::endl;
      gPad -> Update();
      
      ++catIt;
      
      system(Form("mkdir -p %s",plotfolder.c_str()));
      c1 -> Print(Form("%s/c_mgg_%s_%d.png",plotfolder.c_str(),label.c_str(),catIt));
      c1 -> Print(Form("%s/c_mgg_%s_%d.pdf",plotfolder.c_str(),label.c_str(),catIt));
    }
  }
}
