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

void drawFinalPlot(const std::string& folder)
{
  std::map<std::string,TTree*> trees;
  
  TFile* inFile;
  
  inFile = TFile::Open(Form("%s/plotTree_HHggbb.root",folder.c_str()),"READ");
  trees["LT_sig"] = (TTree*)( inFile->Get("all_lowMx") );
  trees["HT_sig"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_ggH.root",folder.c_str()),"READ");
  trees["LT_bkg_res_ggH"] = (TTree*)( inFile->Get("all_lowMx") );
  trees["HT_bkg_res_ggH"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_qqH.root",folder.c_str()),"READ");
  trees["LT_bkg_res_qqH"] = (TTree*)( inFile->Get("all_lowMx") );
  trees["HT_bkg_res_qqH"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_VH.root",folder.c_str()),"READ");
  trees["LT_bkg_res_VH"] = (TTree*)( inFile->Get("all_lowMx") );
  trees["HT_bkg_res_VH"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_ttH.root",folder.c_str()),"READ");
  trees["LT_bkg_res_ttH"] = (TTree*)( inFile->Get("all_lowMx") );
  trees["HT_bkg_res_ttH"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_bbH.root",folder.c_str()),"READ");
  trees["LT_bkg_res_bbH"] = (TTree*)( inFile->Get("all_lowMx") );
  trees["HT_bkg_res_bbH"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_gg.root",folder.c_str()),"READ");
  trees["LT_bkg_nonres_gg"] = (TTree*)( inFile->Get("all_lowMx") );
  trees["HT_bkg_nonres_gg"] = (TTree*)( inFile->Get("all_highMx") );
  
  inFile = TFile::Open(Form("%s/plotTree_ttg.root",folder.c_str()),"READ");
  trees["LT_bkg_nonres_ttg"] = (TTree*)( inFile->Get("all_lowMx") );
  trees["HT_bkg_nonres_ttg"] = (TTree*)( inFile->Get("all_highMx") );
  
  // inFile = TFile::Open(Form("%s/plotTree_gjet.root",folder.c_str()),"READ");
  // trees["LT_bkg_nonres_gjet"] = (TTree*)( inFile->Get("all_lowMx") );
  // trees["HT_bkg_nonres_gjet"] = (TTree*)( inFile->Get("all_highMx") );
  
  // inFile = TFile::Open(Form("%s/plotTree_qcd.root",folder.c_str()),"READ");
  // trees["LT_bkg_nonres_qcd"] = (TTree*)( inFile->Get("all_lowMx") );
  // trees["HT_bkg_nonres_qcd"] = (TTree*)( inFile->Get("all_highMx") );
  
  
  std::vector<std::string> labels;
  labels.push_back("HT");
  labels.push_back("LT");
  
  std::vector<std::string> cats;
  cats.push_back("* (cut_based_ct == 0)");
  cats.push_back("* (cut_based_ct == 1)");
  
  
  for(auto label : labels)
  {
    int catIt = 0;
    for(auto cat : cats)
    {
      TH1F* h1_bkg_all_mgg      = new TH1F(Form("h1_bkg_all_mgg__%s_%d",label.c_str(),catIt),    "",  50,100.,150.);
      TH1F* h1_bkg_all_fit_mgg  = new TH1F(Form("h1_bkg_all_fit_mgg__%s_%d",label.c_str(),catIt),"",  50,100.,150.);
      TH1F* h1_bkg_nonres_mgg   = new TH1F(Form("h1_bkg_nonres_mgg__%s_%d",label.c_str(),catIt), "",  50,100.,150.);
      TH1F* h1_bkg_res_mgg      = new TH1F(Form("h1_bkg_res_mgg__%s_%d",label.c_str(),catIt),    "",  50,100.,150.);
      TH1F* h1_sig_mgg          = new TH1F(Form("h1_sig_mgg__%s_%d",label.c_str(),catIt),        "",  50,100.,150.);
      TH1F* h1_sig_mgg_fine     = new TH1F(Form("h1_sig_mgg_fine__%s_%d",label.c_str(),catIt),   "",2000,100.,150.);
      
      for(auto tree : trees)
      {
        std::size_t found_label = (tree.first).find(label);
        if( found_label == std::string::npos ) continue;
        
        std::size_t found_sig = (tree.first).find("sig");
        if( found_sig != std::string::npos )
        {
          (tree.second) -> Draw(Form("mgg >>+ h1_sig_mgg__%s_%d",     label.c_str(),catIt),Form("evWeight * %f %s",lumin,cat.c_str()),"goff");
          (tree.second) -> Draw(Form("mgg >>+ h1_sig_mgg_fine__%s_%d",label.c_str(),catIt),Form("evWeight * %f %s",lumin,cat.c_str()),"goff");
        }
        
        std::size_t found_bkg_nonres = (tree.first).find("bkg_nonres");
        if( found_bkg_nonres != std::string::npos )
        {
          (tree.second) -> Draw(Form("mgg >>+ h1_bkg_nonres_mgg__%s_%d",label.c_str(),catIt),Form("evWeight * %f %s",lumin,cat.c_str()),"goff");
          (tree.second) -> Draw(Form("mgg >>+ h1_bkg_all_mgg__%s_%d",   label.c_str(),catIt),Form("evWeight * %f %s",lumin,cat.c_str()),"goff");
        }
        
        std::size_t found_bkg_res = (tree.first).find("bkg_res");
        if( found_bkg_res != std::string::npos )
        {
          (tree.second) -> Draw(Form("mgg >>+ h1_bkg_res_mgg__%s_%d",label.c_str(),catIt),Form("evWeight * %f %s",lumin,cat.c_str()),"goff");
          (tree.second) -> Draw(Form("mgg >>+ h1_bkg_all_mgg__%s_%d",label.c_str(),catIt),Form("evWeight * %f %s",lumin,cat.c_str()),"goff");
        }
      }
      
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
      
      TF1* fit_bkg_nonres_mgg = new TF1(Form("fit_bkg_nonres_mgg__%s_%d",label.c_str(),catIt),"expo",100.,150.);
      h1_bkg_nonres_mgg -> Fit(fit_bkg_nonres_mgg,"QNRS+");
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
      std::cout << ">>>>>> bkg_nonres:   " << std::fixed << std::setprecision(0) << std::setw(6) << fit_bkg_nonres_mgg->Eval(125.) << " ev. / GeV" << std::endl;
      std::cout << ">>>>>> bkg_res:      " << std::fixed << std::setprecision(1) << std::setw(6) << h1_bkg_res_mgg->Integral()     << " ev."       << std::endl;
      std::cout << ">>>>>> sig:          " << std::fixed << std::setprecision(1) << std::setw(6) << h1_sig_mgg->Integral()         << " ev."       << std::endl;
      std::cout << ">>>>>> significance: " << std::fixed << std::setprecision(2) << std::setw(6) << significance                                    << std::endl;
      gPad -> Update();
      
      ++catIt;
    }
  }
}
