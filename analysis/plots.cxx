#include "ROOT/RDataFrame.hxx"
#include "TH3D.h"
#include "TH1.h"
#include "TF1.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "THnSparse.h"  // Needed for THnSparseD
#include "TArrayD.h"
#include <string>
#include <TText.h>
#include <TPaveText.h>  // add at top of file if not already included
#include <TLatex.h>

#include "TLorentzVector.h"

void plot_delta_P(ROOT::RDF::RNode rdf,const std::string& output_folder) {
    TCanvas canvas("c1", "delta_P", 800, 600);
    auto hist = rdf.Histo1D(ROOT::RDF::TH1DModel("delta_P", "delta_P (rec - gen); delta P (GeV); Events", 100, -0.5, 0.5), "delta_p");
    hist->Draw();
    canvas.SaveAs((output_folder + "delta_P.pdf").c_str());
    std::cout << "Saved 1D histogram as delta_P.pdf" << std::endl;
}

void plot_momenta_components(ROOT::RDF::RNode rdf, const std::string& output_folder) { // do not use loops, the graphs are too different for slicing and loopiong will lead to lazy eval
    TCanvas canvas("c2", "momenta_components", 800, 600);
    canvas.Divide(3,2);
    canvas.cd(1);
    auto hist1 = rdf.Histo1D(ROOT::RDF::TH1DModel("px_piplus_gen", "px_piplus_gen; px_piplus_gen (GeV); Events", 100, -2, 2), "px_piplus_gen");
    hist1->Draw();
    canvas.cd(2);
    auto hist2 = rdf.Histo1D(ROOT::RDF::TH1DModel("py_piplus_gen", "py_piplus_gen; py_piplus_gen (GeV); Events", 100, -2, 2), "py_piplus_gen");
    hist2->Draw();
    canvas.cd(3);
    auto hist3 = rdf.Histo1D(ROOT::RDF::TH1DModel("pz_piplus_gen", "pz_piplus_gen; pz_piplus_gen (GeV); Events", 100, 0, 8), "pz_piplus_gen");
    hist3->Draw();
    canvas.cd(4);
    auto hist4 = rdf.Histo1D(ROOT::RDF::TH1DModel("px_piplus_rec", "px_piplus_rec; px_piplus_rec (GeV); Events", 100, -2, 2), "px_piplus_rec");
    hist4->Draw();
    canvas.cd(5);
    auto hist5 = rdf.Histo1D(ROOT::RDF::TH1DModel("py_piplus_rec", "py_piplus_rec; py_piplus_rec (GeV); Events", 100, -2, 2), "py_piplus_rec");
    hist5->Draw();
    canvas.cd(6);
    auto hist6 = rdf.Histo1D(ROOT::RDF::TH1DModel("pz_piplus_rec", "pz_piplus_rec; pz_piplus_rec (GeV); Events", 100, 0, 8), "pz_piplus_rec");
    hist6->Draw();

    canvas.SaveAs((output_folder + "momenta_components.pdf").c_str());
    std::cout << "Saved 1D histogram as momenta_components.pdf" << std::endl;
}

void plot_delta_P_VS_P_rec_FD(ROOT::RDF::RNode rdf, const std::string& output_folder) {
    rdf = rdf.Filter("detector == \"FD\" ");
    //rdf = rdf.Filter("Theta_rec < 27");
    TCanvas canvas("c5", "delta P VS P_rec in FD", 800, 600);
    auto hist2D = rdf.Histo2D(
        ROOT::RDF::TH2DModel("delta_P_VS_P_rec_FD", "delta P vs P_rec in FD;  P_rec (GeV); delta P (GeV)", 200, 0, 6, 200, -0.1, 0.1),
        "p_piplus_rec", "delta_p"
    );
    hist2D->Draw("COLZ");
    canvas.SaveAs((output_folder + "delta_P_VS_P_rec_FD.pdf").c_str());
    std::cout << "Saved 2D histogram as delta_P_VS_P_rec_FD.pdf" << std::endl;
}

void delta_P_VS_P_rec_FD_CD(ROOT::RDF::RNode rdf, const std::string& output_folder) {
    TCanvas canvas("c", "delta_P_VS_P_rec_FD_CD", 800, 600);
    canvas.Divide(1,2);
    canvas.cd(1);
    auto hist2D_1 = rdf.Filter("detector == \"FD\"").Histo2D(
        ROOT::RDF::TH2DModel("delta_P_VS_P_rec_FD", "delta P vs P_rec in FD;  P_rec (GeV); delta P (GeV)", 200, 0, 6, 200, -0.1, 0.1),
        "p_piplus_rec", "delta_p"
    );
    hist2D_1->Draw("COLZ");
    canvas.cd(2);
    auto hist2D_2 = rdf.Filter("detector == \"CD\"").Histo2D(
        ROOT::RDF::TH2DModel("delta_P_VS_P_rec_CD", "delta P vs P_rec in CD;  P_rec (GeV); delta P (GeV)", 200, 0, 2.5, 200, -0.1, 0.1),
        "p_piplus_rec", "delta_p"
    );
    hist2D_2->Draw("COLZ");

    canvas.SaveAs((output_folder + "delta_P_VS_P_rec_FD_CD.pdf").c_str());
    std::cout << "Saved 2D histogram as delta_P_VS_P_rec_FD_CD.pdf" << std::endl;
}

void plot_delta_P_VS_P_rec_FD_Theta_below_above(ROOT::RDF::RNode rdf, const std::string& output_folder){
    auto rdf_above = rdf.Filter("(Theta_rec > 40) && detector == \"FD\" ");
    auto rdf_below = rdf.Filter("(Theta_rec < 25) && detector == \"FD\" ");
    TCanvas canvas("c1", "delta_P", 800, 600);
    canvas.Divide(1,2);
    canvas.cd(1);
    auto hist2D_above = rdf_above.Histo2D(
        ROOT::RDF::TH2DModel("delta_P_VS_P_rec_above_Theta", "delta P vs P_rec for theta > 40 deg;  P_rec (GeV); delta P (GeV)", 100, 0, 6, 100, -0.1, 0.1),
        "p_piplus_rec", "delta_p"
    );
    hist2D_above->Draw("COLZ");
    canvas.cd(2);
    auto hist2D_below = rdf_below.Histo2D(
        ROOT::RDF::TH2DModel("delta_P_VS_P_rec_below_Theta", "delta P vs P_rec for theta < 25 deg;  P_rec (GeV); delta P (GeV)", 100, 0, 6, 100, -0.1, 0.1),
        "p_piplus_rec", "delta_p"
    );
    hist2D_below->Draw("COLZ");
    canvas.SaveAs((output_folder + "delta_P_VS_P_rec_FD_Theta_high_low.pdf").c_str());
    std::cout << "Saved 1D histogram as delta_P_VS_P_rec_FD_Theta_high_low.pdf" << std::endl;
}

void Theta_VS_momentum_FD_CD(ROOT::RDF::RNode rdf, const std::string& output_folder) {
    TCanvas canvas("c8", "Theta VS momentum FD CD", 800, 600);
    canvas.Divide(2,2);
    canvas.cd(1);
    auto hist1 = rdf.Filter("detector == \"FD\"").Histo2D(
        ROOT::RDF::TH2DModel("Theta_gen_VS_P_gen_FD", "Theta_gen VS P_gen in FD; P_gen (GeV); Theta_gen (deg)", 100, 0, 6, 100, 0, 100),
        "p_piplus_gen", "Theta_gen"
    );
    hist1->Draw("COLZ");
    canvas.cd(2);
    auto hist2 = rdf.Filter("detector == \"FD\"").Histo2D(
        ROOT::RDF::TH2DModel("Theta_rec_VS_P_rec_FD", "Theta_rec VS P_rec in FD;  P_rec (GeV); Theta_rec (deg);", 100, 0, 6, 100, 0, 100),
        "p_piplus_rec", "Theta_rec"
    );
    hist2->Draw("COLZ");
    canvas.cd(3);
    auto hist3 = rdf.Filter("detector == \"CD\"").Histo2D(
        ROOT::RDF::TH2DModel("Theta_gen_VS_P_gen_CD", "Theta_gen VS P_gen in CD; P_gen (GeV); Theta_gen (deg); ", 100, 0, 6, 100, 0, 100),
        "p_piplus_gen", "Theta_gen"
    );
    hist3->Draw("COLZ");
    canvas.cd(4);
    auto hist4 = rdf.Filter("detector == \"CD\"").Histo2D(
        ROOT::RDF::TH2DModel("Theta_rec_VS_P_rec_CD", "Theta_rec VS P_rec in CD;  P_rec (GeV); Theta_rec (deg);",  100, 0, 6, 100, 0, 100),
        "p_piplus_rec", "Theta_rec"
    );
    hist4->Draw("COLZ");

    canvas.SaveAs((output_folder + "Theta_VS_momentum_FD_CD.pdf").c_str());
    std::cout << "Saved 2D histogram as Theta_VS_momentum_FD_CD.pdf" << std::endl;
}

void Phi_VS_momentum_FD_CD(ROOT::RDF::RNode rdf, const std::string& output_folder) {
    TCanvas canvas("c8", "Phi VS momentum FD CD", 800, 600);
    canvas.Divide(2,2);
    canvas.cd(1);
    auto hist1 = rdf.Filter("detector == \"FD\"").Histo2D(
        ROOT::RDF::TH2DModel("Phi_gen_VS_P_gen_FD", "Phi_gen VS P_gen in FD;  Phi_gen (deg); P_gen (GeV)", 100, -200, 200, 100, 0, 6),
        "Phi_gen", "p_piplus_gen"
    );
    hist1->Draw("COLZ");
    canvas.cd(2);
    auto hist2 = rdf.Filter("detector == \"FD\"").Histo2D(
        ROOT::RDF::TH2DModel("Phi_rec_VS_P_rec_FD", "Phi_rec VS P_rec in FD;  Phi_rec (deg); P_rec (GeV)", 100, -200, 200, 100, 0, 6),
        "Phi_rec", "p_piplus_rec"
    );
    hist2->Draw("COLZ");
    canvas.cd(3);
    auto hist3 = rdf.Filter("detector == \"CD\"").Histo2D(
        ROOT::RDF::TH2DModel("Phi_gen_VS_P_gen_CD", "Phi_gen VS P_gen in CD; Phi_gen (deg); P_gen (GeV)", 100, -200, 200, 100, 0, 6),
        "Phi_gen", "p_piplus_gen"
    );
    hist3->Draw("COLZ");
    canvas.cd(4);
    auto hist4 = rdf.Filter("detector == \"CD\"").Histo2D(
        ROOT::RDF::TH2DModel("Phi_rec_VS_P_rec_CD", "Phi_rec VS P_rec in CD; Phi_rec (deg); P_rec (GeV)", 100, -200, 200, 100, 0, 6),
        "Phi_rec", "p_piplus_rec"
    );
    hist4->Draw("COLZ");

    canvas.SaveAs((output_folder + "Phi_VS_momentum_FD_CD.pdf").c_str());
    std::cout << "Saved 2D histogram as Phi_VS_momentum_FD_CD.pdf" << std::endl;
}

void Phi_VS_Theta_FD_CD(ROOT::RDF::RNode rdf, const std::string& output_folder) {
    TCanvas canvas("c10", "Phi VS Theta FD CD", 800, 600);
    canvas.Divide(2,2);
    canvas.cd(1);
    auto hist1 = rdf.Filter("detector == \"FD\"").Histo2D(
        ROOT::RDF::TH2DModel("Phi_gen_VS_Theta_gen_FD", "Phi_gen VS Theta_gen in FD;  Phi_gen (deg); Theta_gen (deg)", 100, -200, 200, 100, 0, 100),
        "Phi_gen", "Theta_gen"
    );
    hist1->Draw("COLZ");
    canvas.cd(2);
    auto hist2 = rdf.Filter("detector == \"FD\"").Histo2D(
        ROOT::RDF::TH2DModel("Phi_rec_VS_Theta_rec_FD", "Phi_rec VS Theta_rec in FD;  Phi_rec (deg); Theta_rec (deg)", 100, -200, 200, 100, 0, 100),
        "Phi_rec", "Theta_rec"
    );
    hist2->Draw("COLZ");
    canvas.cd(3);
    auto hist3 = rdf.Filter("detector == \"CD\"").Histo2D(
        ROOT::RDF::TH2DModel("Phi_gen_VS_Theta_gen_CD", "Phi_gen VS Theta_gen in CD; Phi_gen (deg); Theta_gen (deg)", 100, -200, 200, 100, 0, 100),
        "Phi_gen", "Theta_gen"
    );
    hist3->Draw("COLZ");
    canvas.cd(4);
    auto hist4 = rdf.Filter("detector == \"CD\"").Histo2D(
        ROOT::RDF::TH2DModel("Phi_rec_VS_Theta_rec_CD", "Phi_rec VS Theta_rec in CD; Phi_rec (deg); Theta_rec (deg)", 100, -200, 200, 100, 0, 100),
        "Phi_rec", "Theta_rec"
    );
    hist4->Draw("COLZ");

    canvas.SaveAs((output_folder + "Phi_VS_Theta_FD_CD.pdf").c_str());
    std::cout << "Saved 2D histogram as Phi_VS_Theta_FD_CD.pdf" << std::endl;
}

// Unite all FD sectors: build one graph over momentum bins and fit a single curve
// Now: loop over theta bins (like momentum bins) and do this per-theta-bin for pi+
void delta_P_VS_P_rec_FD_unified_1D(ROOT::RDF::RNode rdf,
                                   const std::string& output_folder,
                                   const std::vector<double>& theta_bins,
                                   const std::vector<double>& momentum_bins,
                                   const bool normalized) {
  std::string dp_Or_dpp = normalized ? "dp_norm" : "delta_p";


  
  const size_t num_p_bins = momentum_bins.size() - 1;

  // Safety
  if (theta_bins.size() < 2) {
    std::cerr << "ERROR: theta_bins must have at least 2 edges.\n";
    return;
  }

  // Loop over theta bins
  for (size_t tbin = 0; tbin + 1 < theta_bins.size(); ++tbin) {
    const double th_lo = theta_bins[tbin];
    const double th_hi = theta_bins[tbin + 1];
    const double th_ctr = 0.5 * (th_lo + th_hi);

    // Tag for unique names / files
    TString thetaTag = Form("theta_%g_%g", th_lo, th_hi);

    ROOT::RDF::RNode rdf_filtered =
        rdf.Filter(Form("detector == \"FD\" && Theta_rec >= %f && Theta_rec < %f", th_lo, th_hi));

    // 2D histogram over ALL sectors: X = p_rec, Y = Δp (or Δp/p)
    TString h2name = Form("h2_unified_%s_%s", thetaTag.Data(), dp_Or_dpp.c_str());
    auto h2 = normalized
        ? rdf_filtered.Histo2D(
              ROOT::RDF::TH2DModel(h2name.Data(),
                                   Form("dp/p vs P_{rec} (FD, all sectors, %s);P_{rec} (GeV/c);dp/p",
                                        thetaTag.Data()),
                                   100, 0.0, 6.0, 100, -0.05, 0.05),
              "p_piplus_rec", "dp_norm")
        : rdf_filtered.Histo2D(
              ROOT::RDF::TH2DModel(h2name.Data(),
                                   Form("delta P vs P_{rec} (FD, all sectors, %s);P_{rec} (GeV/c);delta P (GeV/c)",
                                        thetaTag.Data()),
                                   100, 0.0, 6.0, 100, -0.05, 0.05),
              "p_piplus_rec", "delta_p");

    // Slices canvas (show all momentum-bin projections)
    const int nCols = 6;
    const int nRows = 6;
    TCanvas* cSlices = new TCanvas(Form("unified_%s_slices", thetaTag.Data()),
                                   Form("Unified %s slices (FD, all sectors, %s)",
                                        dp_Or_dpp.c_str(), thetaTag.Data()),
                                   1400, 900);
    cSlices->Divide(nCols, nRows);

    // Graph of mean Δp (or Δp/p) vs momentum-bin center
    TGraphErrors* gAll = new TGraphErrors();
    gAll->SetName(Form("gUnified_%s_%s", thetaTag.Data(), dp_Or_dpp.c_str()));
    gAll->SetTitle(Form("FD (all sectors): Mean %s vs Momentum Bin (%s);Momentum Bin Center (GeV/c);Mean %s",
                        dp_Or_dpp.c_str(), thetaTag.Data(), dp_Or_dpp.c_str()));

    for (size_t bin_idx = 0; bin_idx < num_p_bins; ++bin_idx) {
      double p_low = momentum_bins[bin_idx];
      double p_high = momentum_bins[bin_idx + 1];
      double p_center = 0.5 * (p_low + p_high);

      // --- Project Y using X-bin indices (excludes under/overflow) ---
      TAxis* xax = h2->GetXaxis();
      const int nbx = xax->GetNbins();
      const double eps = 1e-9;

      int ix_lo = xax->FindFixBin(p_low  + eps);
      int ix_hi = xax->FindFixBin(p_high - eps);
      ix_lo = std::max(1, std::min(nbx, ix_lo));
      ix_hi = std::max(1, std::min(nbx, ix_hi));
      if (ix_lo > ix_hi) continue;

      TH1* hY = h2->ProjectionY(
          Form("unified_%s_%s_bin%zu", thetaTag.Data(), dp_Or_dpp.c_str(), bin_idx + 1),
          ix_lo, ix_hi
      );

      // Style and draw slice
      cSlices->cd((int)bin_idx + 1);
      hY->SetTitle(Form("Theta %.1f-%.1f deg, P_{rec} %.2f-%.2f GeV/c; %s; Counts",
                        th_lo, th_hi, p_low, p_high, dp_Or_dpp.c_str()));
      hY->Draw();

      // --- Adaptive, momentum-dependent two-step Gaussian fit ---

      int imax = hY->GetMaximumBin();
      double mode = hY->GetBinCenter(imax);

      Double_t probs[3] = {0.16, 0.50, 0.84}, q[3] = {0,0,0};
      hY->GetQuantiles(3, q, probs);
      double sigma68 = 0.5 * (q[2] - q[0]);
      if (!(sigma68 > 0) || !std::isfinite(sigma68)) {
        sigma68 = hY->GetRMS();
        if (!(sigma68 > 0) || !std::isfinite(sigma68))
          sigma68 = 3.0 * hY->GetBinWidth(1);
      }

      if (hY->GetEntries() > 0) {
        double binsPerSigma = sigma68 / hY->GetBinWidth(1);
        if (binsPerSigma < 6.0 && hY->GetNbinsX() >= 80) hY->Rebin(2);
      }

      auto k_init_for_p = [](double p){
        double k = 1.3 + 0.25 * std::max(0.0, std::min((p - 1.5), 4.0));
        return std::min(2.3, k);
      };
      double k1 = k_init_for_p(p_center);

      double xLo = std::max(hY->GetXaxis()->GetXmin(), mode - k1 * sigma68);
      double xHi = std::min(hY->GetXaxis()->GetXmax(), mode + k1 * sigma68);
      if (xLo >= xHi) {
        xLo = mode - 1.5 * sigma68;
        xHi = mode + 1.5 * sigma68;
      }

      TF1* fit_init = new TF1(Form("gaus_init_unified_%s_%zu", thetaTag.Data(), bin_idx), "gaus", xLo, xHi);
      hY->Fit(fit_init, "RQ0");

      double mu    = fit_init->GetParameter(1);
      double sigma = std::abs(fit_init->GetParameter(2));
      if (!(sigma > 0) || !std::isfinite(sigma)) sigma = sigma68;

      double k2 = 1.25;
      double rLo = std::max(hY->GetXaxis()->GetXmin(), mu - k2 * sigma);
      double rHi = std::min(hY->GetXaxis()->GetXmax(), mu + k2 * sigma);
      if (rLo >= rHi) { rLo = mu - 1.2 * sigma; rHi = mu + 1.2 * sigma; }

      TF1* fit_refined = new TF1(Form("gaus_refined_unified_%s_%zu", thetaTag.Data(), bin_idx), "gaus", rLo, rHi);
      fit_refined->SetLineWidth(1); 
      hY->Fit(fit_refined, "RQ");

      double mean     = fit_refined->GetParameter(1);
      double mean_err = 0.0001;

      int ip = gAll->GetN();
      gAll->SetPoint(ip, p_center, mean);
      gAll->SetPointError(ip, 0.0, mean_err);

      gPad->Update();
    }

    cSlices->SaveAs((output_folder + std::string(thetaTag.Data()) +
                     Form("_UNIFIED_slices_%s.pdf", dp_Or_dpp.c_str())).c_str());
    delete cSlices;

    // Summary canvas (single panel)
    TCanvas* cSummary = new TCanvas(Form("unified_%s_summary", thetaTag.Data()),
                                    Form("Unified mean %s vs momentum (FD, all sectors, %s)",
                                         dp_Or_dpp.c_str(), thetaTag.Data()),
                                    1400, 900);
    gAll->SetMarkerStyle(20);
    gAll->SetMarkerSize(1.0);
    gAll->SetMarkerColor(kBlack);
    gAll->SetLineColor(kBlack);
    gAll->Draw("AP");
    gPad->Update();

    // Autoscale axes from points, but ALWAYS include y=0
    double xmin_pts = 1e9, xmax_pts = -1e9, ymin_pts = 1e9, ymax_pts = -1e9;
    for (int k = 0; k < gAll->GetN(); ++k) {
      double xp, yp; gAll->GetPoint(k, xp, yp);
      if (!std::isfinite(xp) || !std::isfinite(yp)) continue;
      xmin_pts = std::min(xmin_pts, xp);
      xmax_pts = std::max(xmax_pts, xp);
      ymin_pts = std::min(ymin_pts, yp);
      ymax_pts = std::max(ymax_pts, yp);
    }
    if (xmin_pts < xmax_pts) {
      double xpad = 0.05 * (xmax_pts - xmin_pts);
      double xmin_auto = xmin_pts - xpad;
      double xmax_auto = xmax_pts + xpad;

      ymin_pts = std::min(ymin_pts, 0.0);
      ymax_pts = std::max(ymax_pts, 0.0);
      double ypad = 0.10 * std::max(1e-6, ymax_pts - ymin_pts);
      double ymin_auto = ymin_pts - ypad;
      double ymax_auto = ymax_pts + ypad;

      gAll->GetXaxis()->SetLimits(xmin_auto, xmax_auto);
      gAll->GetYaxis()->SetRangeUser(ymin_auto, ymax_auto);
      if (TH1* fr = gAll->GetHistogram()) {
        fr->GetXaxis()->SetLimits(xmin_auto, xmax_auto);
        fr->SetMinimum(ymin_auto);
        fr->SetMaximum(ymax_auto);
      }
      gPad->Update();
    }

    gPad->SetGrid();

   // // Fit range (keep your original behavior, generalized)
   // double xmin_fit = 0.25, xmax_fit = 3.0;
   // if (th_hi <= 27.0) { xmin_fit = 0.25; xmax_fit = 2.0; }

   // // --- Fit unified data with f(p) = A/(B + C*sqrt(p) + D*p + E*p^2) ---
   // double pmin = 1e9, pmax = -1e9, pL = 0, yL = 0, pR = 0, yR = 0;
   // for (int k = 0; k < gAll->GetN(); ++k) {
   //   double xp, yp; gAll->GetPoint(k, xp, yp);
   //   if (!std::isfinite(xp) || !std::isfinite(yp)) continue;
   //   if (xp < xmin_fit || xp > xmax_fit) continue;
   //   if (xp < pmin) { pmin = xp; pL = xp; yL = yp; }
   //   if (xp > pmax) { pmax = xp; pR = xp; yR = yp; }
   // }

   // double A0 = (std::isfinite(pR * yR) ? pR * yR : -1e-2);
   // if (!std::isfinite(A0)) A0 = -1e-2;

   // auto seedB = [&](double p, double y) {
   //   return (std::abs(y) > 1e-12) ? (A0 / y - p) : (-p + 0.05);
   // };
   // double B0 = 0.5 * (seedB(pL, yL) + seedB(pR, yR));
   // if (!std::isfinite(B0)) B0 = 0.1;

   // double C0 = 0.0;
   // double D0 = 1.0;
   // double E0 = 0.0;

   // const double eps = 1e-3;
   // double Bmin = eps,  Bmax = 20.0;
   // double Cmin = -5.0, Cmax =  5.0;
   // double Dmin =  0.0, Dmax =  5.0;
   // double Emin =  0.0, Emax =  5.0;

   // TF1* fitFunc = new TF1(Form("fit_unified_ABCDsqE_%s", thetaTag.Data()),
   //     "[0]/([1] + [2]*sqrt(x) + [3]*x + [4]*x*x)", xmin_fit, xmax_fit);
   // fitFunc->SetParNames("A","B","C","D","E");
   // fitFunc->SetParameters(A0, B0, C0, D0, E0);
   // fitFunc->SetParLimits(1, Bmin, Bmax);
   // fitFunc->SetParLimits(2, Cmin, Cmax);
   // fitFunc->SetParLimits(3, Dmin, Dmax);
   // fitFunc->SetParLimits(4, Emin, Emax);

   // fitFunc->FixParameter(4, 0.0);
   // gAll->Fit(fitFunc, "RQ");
   // fitFunc->ReleaseParameter(4);
   // gAll->Fit(fitFunc, "RQ");

   // fitFunc->SetLineColor(kBlue);
   // fitFunc->SetLineStyle(1);
   // //fitFunc->Draw("SAME");

   // double A = fitFunc->GetParameter(0);
   // double B = fitFunc->GetParameter(1);
   // double C = fitFunc->GetParameter(2);
   // double D = fitFunc->GetParameter(3);
   // double E = fitFunc->GetParameter(4);

   // double eA = fitFunc->GetParError(0);
   // double eB = fitFunc->GetParError(1);
   // double eC = fitFunc->GetParError(2);
   // double eD = fitFunc->GetParError(3);
   // double eE = fitFunc->GetParError(4);

   // double chi2 = fitFunc->GetChisquare();
   // int    ndf  = fitFunc->GetNDF();

   // TLatex latex;
   // latex.SetTextFont(42);
   // latex.SetTextSize(0.038);
   // latex.SetNDC();
   // latex.DrawLatex(0.55, 0.36, Form("A = %.3e #pm %.1e", A, eA));
   // latex.DrawLatex(0.55, 0.32, Form("B = %.3e #pm %.1e", B, eB));
   // latex.DrawLatex(0.55, 0.28, Form("C = %.3e #pm %.1e", C, eC));
   // latex.DrawLatex(0.55, 0.24, Form("D = %.3e #pm %.1e", D, eD));
   // latex.DrawLatex(0.55, 0.20, Form("E = %.3e #pm %.1e", E, eE));
   // latex.DrawLatex(0.55, 0.16, Form("#chi^{2}/NDF = %.1f / %d = %.2f", chi2, ndf, chi2 / ndf));

    double xlo = gAll->GetXaxis()->GetXmin();
    double xhi = gAll->GetXaxis()->GetXmax();
    TLine* zeroLine = new TLine(xlo, 0.0, xhi, 0.0);
    zeroLine->SetLineColor(kRed);
    zeroLine->SetLineStyle(2);
    zeroLine->SetLineWidth(2);
    zeroLine->Draw("SAME");

    cSummary->SaveAs((output_folder + std::string(thetaTag.Data()) +
                      "_mean_" + dp_Or_dpp + "_vs_momentum_bin_UNIFIED.pdf").c_str());
    delete cSummary;

    // Keep graph object around only if you need it later; otherwise delete
    delete gAll;
  }
}

