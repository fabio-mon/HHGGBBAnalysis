void drawROC()
{
  TFile* inFile = TFile::Open("ROC.root","READ");
  
  TGraph* g_ROC_cutBased = (TGraph*)( inFile->Get("g_ROC_cutBased") );
  TGraph* g_ROC_mva = (TGraph*)( inFile->Get("g_ROC_mva") );
  TGraph* g_ROC_mva_new = (TGraph*)( inFile->Get("g_ROC_mva_new") );
  
  TCanvas* c1 = new TCanvas();
  c1 -> cd();
  
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,1.05,1.05) );
  hPad -> SetTitle(Form(";#varepsilon_{signal};1 - #varepsilon_{bkg}"));
  hPad -> SetLineColor(kWhite);
  hPad -> Draw();
  
  g_ROC_mva -> SetLineColor(kBlue);
  g_ROC_mva -> SetLineWidth(2);
  g_ROC_mva -> Draw("C,same");
  
  g_ROC_mva_new -> SetLineColor(kRed);
  g_ROC_mva_new -> SetLineWidth(2);
  g_ROC_mva_new -> Draw("C,same");
  
  g_ROC_cutBased -> Draw("P,same");
  
  TLegend* legend = new TLegend(0.20,0.50-0.04*3,0.50,0.50);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  legend -> AddEntry(g_ROC_mva,"diphoton MVA (2016)","L");
  legend -> AddEntry(g_ROC_mva_new,"diphoton MVA (new)","L");
  legend -> AddEntry(g_ROC_cutBased,"cut-based selection","P");
  legend -> Draw("same");
  
  c1 -> Print("ROC_withCutBased.pdf");
}
