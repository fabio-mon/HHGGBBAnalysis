#ifndef HH_REWEIGHT_COMPONENTS_H
#define HH_REWEIGHT_COMPONENTS_H

#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <utility>
#include <memory>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

class HH_reweight_components{
    
    public:
        HH_reweight_components (std::string comp_file);
        ~HH_reweight_components();
        double functionGF(double kl, double kt, double c2, double cg, double c2g, std::array<double, 15> const &A);
        double get_weight(double kl, double mHH); // neglect cos theta
        double get_weight(double kl, double mHH, double costh){return get_weight(kl, mHH);} // a simple bypass for compatibility - neglect cos theta
        double get_xs_ratio(double kl) {return functionGF(kl, 1, 0, 0, 0, A_integralXS);}
    private:
        std::pair<int,int> get_bin_idx(double x, double y, TH2* histo);
        int get_bin_idx(double x, TH1* histo);

        std::unique_ptr<TFile> f_comp;
        TH1* HH_14TeV_histo_SM;
        TH1* HH_14TeV_histo_tri;
        TH1* HH_14TeV_histo_box;
        TH1* HH_14TeV_histo_interf;

        std::array<double, 15> A_integralXS { {
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
	      } };
};

HH_reweight_components::HH_reweight_components(std::string comp_file)
{
    f_comp   = std::unique_ptr<TFile> (new TFile(comp_file.c_str()));

    // read out the histograms
    HH_14TeV_histo_SM     = (TH1*) f_comp->Get("h_SM");
    HH_14TeV_histo_tri    = (TH1*) f_comp->Get("h_tri");
    HH_14TeV_histo_box    = (TH1*) f_comp->Get("h_box");
    HH_14TeV_histo_interf = (TH1*) f_comp->Get("h_interf");
}

HH_reweight_components::~HH_reweight_components()
{
    // will this doubly delete the ptrs?
    f_comp->Close();
}

double HH_reweight_components::functionGF(double kl, double kt, double c2, double cg, double c2g, std::array<double, 15> const &A)
{
    return ( A[0]*pow(kt,4) + A[1]*pow(c2,2) + (A[2]*pow(kt,2) + A[3]*pow(cg,2))*pow(kl,2) + A[4]*pow(c2g,2) + ( A[5]*c2 + A[6]*kt*kl )*pow(kt,2) + (A[7]*kt*kl + A[8]*cg*kl )*c2 + A[9]*c2*c2g + (A[10]*cg*kl + A[11]*c2g)*pow(kt,2)+ (A[12]*kl*cg + A[13]*c2g )*kt*kl + A[14]*cg*c2g*kl );
}

// double HH_reweight_components::get_weight(double kl, double mHH, double costh)
double HH_reweight_components::get_weight(double kl, double mHH)
{
    // 1) get the yields in the specific bin
    int ibin = get_bin_idx(mHH, HH_14TeV_histo_SM);
    double SM      = HH_14TeV_histo_SM->GetBinContent(ibin);
    double tri     = HH_14TeV_histo_tri->GetBinContent(ibin);
    double box     = HH_14TeV_histo_box->GetBinContent(ibin);
    double interf  = HH_14TeV_histo_interf->GetBinContent(ibin);

    // scale the components according to klambda
    // box unscaled by kl
    tri *= kl*kl;
    interf *= kl;

    double tot_xs = box + tri + interf;
    double r = tot_xs / SM;
    return r;
}

int HH_reweight_components::get_bin_idx(double x, TH1* histo)
{
    int ix = histo->GetXaxis()->FindBin(x);
    if (ix < 1) ix = 1;
    if (ix > histo->GetNbinsX()) ix = histo->GetNbinsX();
    return ix;
}

std::pair<int,int> HH_reweight_components::get_bin_idx(double x, double y, TH2* histo)
{
    int ix = histo->GetXaxis()->FindBin(x);
    int iy = histo->GetYaxis()->FindBin(y);
    if (ix < 1) ix = 1;
    if (iy < 1) iy = 1;
    if (ix > histo->GetNbinsX()) ix = histo->GetNbinsX();
    if (iy > histo->GetNbinsY()) iy = histo->GetNbinsY();
    return std::make_pair(ix,iy);
}

#endif
