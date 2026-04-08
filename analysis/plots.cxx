#include "ROOT/RDataFrame.hxx"
#include "TH3D.h"
#include "TH2D.h"
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
#include <map>
#include <utility>
#include "TLorentzVector.h"

using ThetaRange = std::pair<double,double>;
using ThetaToPBins = std::map<ThetaRange, std::vector<double>>;

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

void delta_P_VS_P_rec_FD_CD(ROOT::RDF::RNode rdf, const std::string& output_folder, bool no_low_bin_cont = true) {
    rdf = rdf.Filter("Vz_cut == true && Vxy_cut == true");
    TCanvas canvas("c", "delta_P_VS_P_rec_FD_CD", 800, 600);
    canvas.Divide(1,2);
    canvas.cd(1);
    auto hist2D_1 = rdf.Filter("detector == \"FD\"").Histo2D(
        ROOT::RDF::TH2DModel("delta_P_VS_P_rec_FD", "delta P vs P_rec in FD;  P_rec (GeV); delta P (GeV)", 200, 0, 10, 200, -0.1, 0.1),
        "p_piplus_rec", "delta_p"
    );
    if (no_low_bin_cont) {
        hist2D_1->SetMinimum(5); // Set minimum to 0.5 to avoid showing empty bins in log scale
    }
    hist2D_1->Draw("COLZ");
    canvas.cd(2);
    auto hist2D_2 = rdf.Filter("detector == \"CD\"").Histo2D(
        ROOT::RDF::TH2DModel("delta_P_VS_P_rec_CD", "delta P vs P_rec in CD;  P_rec (GeV); delta P (GeV)", 200, 0, 8, 200, -0.1, 0.1),
        "p_piplus_rec", "delta_p"
    );
    if (no_low_bin_cont) {
        hist2D_2->SetMinimum(5); // Set minimum to 0.5 to avoid showing empty bins in log scale
    }
    hist2D_2->Draw("COLZ");

    canvas.SaveAs((output_folder + "delta_P_VS_P_rec_FD_CD.pdf").c_str());
    std::cout << "Saved 2D histogram as delta_P_VS_P_rec_FD_CD.pdf" << std::endl;
}

void plot_delta_P_VS_P_rec_FD_Theta_below_above(ROOT::RDF::RNode rdf, const std::string& output_folder) {
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

void Theta_VS_momentum_FD_CD(ROOT::RDF::RNode rdf, const std::string& output_folder,  bool logz = true, bool no_low_bin_cont = true) {
    TCanvas canvas("c8", "Theta VS momentum FD CD", 800, 600);
    canvas.Divide(2,2);
    canvas.cd(1);
    rdf = rdf.Filter("Vz_cut == true && Vxy_cut == true");
    auto hist1 = rdf.Filter("detector == \"FD\"").Histo2D(
        ROOT::RDF::TH2DModel("Theta_gen_VS_P_gen_FD", "Theta_gen VS P_gen in FD; P_gen (GeV); Theta_gen (deg)", 200, 0, 10, 200, 0, 140),
        "p_piplus_gen", "Theta_gen"
    );
    if (logz) gPad->SetLogz();
    if (no_low_bin_cont) {
        hist1->SetMinimum(5); // Set minimum to 0.5 to avoid showing empty bins in log scale
    }
    hist1->Draw("COLZ");
    canvas.cd(2);
    auto hist2 = rdf.Filter("detector == \"FD\"").Histo2D(
        ROOT::RDF::TH2DModel("Theta_rec_VS_P_rec_FD", "Theta_rec VS P_rec in FD;  P_rec (GeV); Theta_rec (deg);", 200, 0, 10, 200, 0, 140),
        "p_piplus_rec", "Theta_rec"
    );
    if (logz) gPad->SetLogz();
    if (no_low_bin_cont) {
        hist2->SetMinimum(5); // Set minimum to 0.5 to avoid showing empty bins in log scale
    }
    hist2->Draw("COLZ");
    canvas.cd(3);
    auto hist3 = rdf.Filter("detector == \"CD\"").Histo2D(
        ROOT::RDF::TH2DModel("Theta_gen_VS_P_gen_CD", "Theta_gen VS P_gen in CD; P_gen (GeV); Theta_gen (deg); ", 200, 0, 10, 200, 0, 180),
        "p_piplus_gen", "Theta_gen"
    );
    if (logz) gPad->SetLogz();
    if (no_low_bin_cont) {
        hist3->SetMinimum(5); // Set minimum to 0.5 to avoid showing empty bins in log scale
    }
    hist3->Draw("COLZ");
    canvas.cd(4);
    auto hist4 = rdf.Filter("detector == \"CD\"").Histo2D(
        ROOT::RDF::TH2DModel("Theta_rec_VS_P_rec_CD", "Theta_rec VS P_rec in CD;  P_rec (GeV); Theta_rec (deg);",  200, 0, 10, 200, 0, 180),
        "p_piplus_rec", "Theta_rec"
    );
    if (logz) gPad->SetLogz();
    if (no_low_bin_cont) {
        hist4->SetMinimum(5); // Set minimum to 0.5 to avoid showing empty bins in log scale
    }
    hist4->Draw("COLZ");

    canvas.SaveAs((output_folder + "Theta_VS_momentum_FD_CD_Vz_Vxy_cut_setmin5.pdf").c_str());
    std::cout << "Saved 2D histogram as Theta_VS_momentum_FD_CD.pdf" << std::endl;
}

void Theta_VS_momentum_FD_gen_cuts(ROOT::RDF::RNode rdf, const std::string& output_folder, double theta_gen_min, double theta_gen_max, double p_gen_min, double p_gen_max, bool logz = true) {
    TCanvas canvas("c8", "Theta VS momentum FD", 800, 600);
    rdf = rdf.Filter("pz_piplus_rec > 0");
    canvas.Divide(1,2);
    canvas.cd(1);
    auto hist1 = rdf.Filter("detector == \"FD\"").Filter("Theta_gen >"+std::to_string(theta_gen_min)).Filter("Theta_gen < "+std::to_string(theta_gen_max)).Filter("p_piplus_gen > "+std::to_string(p_gen_min)).Filter("p_piplus_gen < "+std::to_string(p_gen_max))
    .Histo2D(
        ROOT::RDF::TH2DModel("Theta_gen_VS_P_gen_FD", "Theta_gen VS P_gen in FD; P_gen (GeV); Theta_gen (deg)", 200, 0, 6.5, 200, theta_gen_min*0.9, theta_gen_max*1.1),
        "p_piplus_gen", "Theta_gen"
    );
    if (logz) gPad->SetLogz();
    hist1->Draw("COLZ");
    canvas.cd(2);
    auto hist2 = rdf.Filter("detector == \"FD\"").Filter("Theta_gen >"+std::to_string(theta_gen_min)).Filter("Theta_gen < "+std::to_string(theta_gen_max)).Filter("p_piplus_gen > "+std::to_string(p_gen_min)).Filter("p_piplus_gen < "+std::to_string(p_gen_max))
    .Histo2D(
        ROOT::RDF::TH2DModel("Theta_rec_VS_P_rec_FD", "Theta_rec VS P_rec in FD;  P_rec (GeV); Theta_rec (deg);", 200, 0, 6.5, 200, 0, 150),
        "p_piplus_rec", "Theta_rec"
    );
    if (logz) gPad->SetLogz();
    hist2->Draw("COLZ");

    canvas.SaveAs((output_folder + "Theta_VS_momentum_FD_theta_gen_" + std::to_string(theta_gen_min) + "-" + std::to_string(theta_gen_max) + "_p_gen_" + std::to_string(p_gen_min) + "-" + std::to_string(p_gen_max) + "_cut.pdf").c_str());
    std::cout << "Saved 2D histogram as Theta_VS_momentum_FD_theta_gen_" + std::to_string(theta_gen_min) + "_" + std::to_string(theta_gen_max) + "_p_gen_" + std::to_string(p_gen_min) + "_" + std::to_string(p_gen_max) + "_cut.pdf" << std::endl;

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
void delta_P_VS_P_rec_unified_1D(ROOT::RDF::RNode rdf,
                                   const std::string& output_folder,
                                   const ThetaToPBins& theta_to_momentum_bins,
                                   const bool normalized,
                                   const std::string detector,
                                   double phi_low = -180.0, double phi_high = 180.0){

  std::string dp_Or_dpp = normalized ? "dp_norm" : "delta_p";

// Loop over theta bins provided by the dictionary
for (const auto& kv : theta_to_momentum_bins) {
    const double th_lo = kv.first.first;
    const double th_hi = kv.first.second;
    const auto& momentum_bins = kv.second;

    if (momentum_bins.size() < 2) continue; // need at least 2 edges
    const size_t num_p_bins = momentum_bins.size() - 1;

    TString thetaTag = Form("theta_%g_%g", th_lo, th_hi);

    ROOT::RDF::RNode rdf_filtered =rdf.Filter(Form("detector == \"%s\" && Phi_rec >= %f && Phi_rec <= %f  && Theta_rec >= %f && Theta_rec < %f", detector.c_str(), phi_low, phi_high, th_lo, th_hi))
                                    .Filter("Vz_cut == true") // Apply Vz cut to select good events
                                    //.Filter("DC_fiducial_cut_piplus == true") // Apply DC fiducial cut for pi+
                                    .Filter("Vxy_cut == true"); // Apply Vxy cut to select good events

    std::cout<<"Detector = "<<detector.c_str()<<std::endl;



    // 2D histogram over ALL sectors: X = p_rec, Y = Δp (or Δp/p)
    TString h2name = Form("h2_unified_%s_%s", thetaTag.Data(), dp_Or_dpp.c_str());
    const int ny = 500; // default
    auto h2 = normalized
        ? rdf_filtered.Histo2D(
              ROOT::RDF::TH2DModel(h2name.Data(),
                                   Form("dp/p vs P_{rec} ( all sectors, %s);P_{rec} (GeV/c);dp/p",
                                        thetaTag.Data()),
                                   100, 0.0, 9.0, ny, -0.5, 0.3),
              "p_piplus_rec", "dp_norm")
        : rdf_filtered.Histo2D(
              ROOT::RDF::TH2DModel(h2name.Data(),
                                   Form("delta P vs P_{rec} ( all sectors, %s);P_{rec} (GeV/c);delta P (GeV/c)",
                                        thetaTag.Data()),
                                   100, 0.0, 9.0, ny, -0.5, 0.3),
              "p_piplus_rec", "delta_p");

    // Slices canvas (show all momentum-bin projections)
    const int nCols = 7;
    const int nRows = 7;
    TCanvas* cSlices = new TCanvas(Form("unified_%s_slices", thetaTag.Data()),
                                   Form("Unified %s slices (all sectors, %s)",
                                        dp_Or_dpp.c_str(), thetaTag.Data()),
                                   1400, 900);
    cSlices->Divide(nCols, nRows);

    // Graph of mean Δp (or Δp/p) vs momentum-bin center
    TGraphErrors* gAll = new TGraphErrors();
    gAll->SetName(Form("gUnified_%s_%s", thetaTag.Data(), dp_Or_dpp.c_str()));
    gAll->SetTitle(Form(" (all sectors): Mean %s vs Momentum Bin (%s);Momentum Bin Center (GeV/c);Mean %s",
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
      double xmin = -0.15, xmax = 0.15;     // default
      if (th_lo >= 25  && p_high >= 7.4) { xmin = -0.15; xmax = 0.30; }  
      else if (th_lo >= 37  && p_high >= 1.6) { xmin = -0.3; xmax = 0.15; }
      else if (th_lo >= 40 && p_high >=1.0) {xmin = -0.3; xmax = 0.15;}

      


                // --- Two-step Gaussian fit ---
      cSlices->cd((int)bin_idx + 1);
      gPad->Clear();
      hY->GetXaxis()->SetRangeUser(xmin, xmax);
      hY->Draw();

      // 1st fit: exploratory over the whole displayed region
      TF1* fit_explore = new TF1(
          Form("gaus_explore_unified_%s_%zu", thetaTag.Data(), bin_idx),
          "gaus", xmin, xmax
      );
      fit_explore->SetParameters(hY->GetMaximum(), hY->GetMean(), hY->GetRMS());
      fit_explore->SetLineWidth(1);
      hY->Fit(fit_explore, "RQ0");  // 0 = do not draw exploratory fit

      double mean1  = fit_explore->GetParameter(1);
      double sigma1 = std::abs(fit_explore->GetParameter(2));

      // fallback in case exploratory fit is garbage
      if (!(sigma1 > 0) || !std::isfinite(sigma1)) sigma1 = hY->GetRMS();
      if (!(sigma1 > 0) || !std::isfinite(sigma1)) sigma1 = 2.0 * hY->GetBinWidth(1);
      if (!std::isfinite(mean1)) mean1 = hY->GetMean();

      // 2nd fit: actual displayed fit within mean ± 2 sigma
      double fit_lo = std::max(xmin, mean1 - 1.0 * sigma1);
      double fit_hi = std::min(xmax, mean1 + 1.0 * sigma1);

      if (fit_lo >= fit_hi) {
          fit_lo = xmin;
          fit_hi = xmax;
      }

      TF1* fit_final = new TF1(
          Form("gaus_final_unified_%s_%zu", thetaTag.Data(), bin_idx),
          "gaus", fit_lo, fit_hi
      );
      fit_final->SetParameters(hY->GetMaximum(), mean1, sigma1);
      fit_final->SetLineWidth(1);

      hY->Fit(fit_final, "RQ");   // this one is shown

      double mean     = fit_final->GetParameter(1);
      double mean_err = fit_final->GetParError(1);
      if (!(mean_err > 0) || !std::isfinite(mean_err)) mean_err = 0.0001;
      //double mean_err = 0.0001;


      int ip = gAll->GetN();
      gAll->SetPoint(ip, p_center, mean);
      gAll->SetPointError(ip, 0.0, mean_err);

      // Re-apply zoom just in case Fit/Update messes with it (cheap and safe)
      hY->GetXaxis()->SetRangeUser(xmin, xmax);
      gPad->Update();
    }

    cSlices->SaveAs((output_folder + detector + "/" + detector + "_" + std::string(thetaTag.Data()) +
                     Form("_UNIFIED_slices_%s_phi_%.0f-%.0f_Vxyz_cut.pdf", dp_Or_dpp.c_str(), phi_low, phi_high)).c_str());
    delete cSlices;

    // Summary canvas (single panel)
    TCanvas* cSummary = new TCanvas(Form("unified_%s_summary", thetaTag.Data()),
                                    Form("Unified mean %s vs momentum (all sectors, %s)",
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

    cSummary->SaveAs((output_folder + detector + "/" + detector + "_" + std::string(thetaTag.Data()) +
                      "_mean_" + dp_Or_dpp + "_vs_momentum_bin_UNIFIED_phi_" + std::to_string(int(phi_low)) + "-" + std::to_string(int(phi_high)) + "_Vxyz_cut.pdf").c_str());
    delete cSummary;

    // Keep graph object around only if you need it later; otherwise delete
    delete gAll;
  }
}

void plot_X_vs_Y_piplus_FD(ROOT::RDF::RNode rdf,
                           const std::string& output_folder,
                           bool logz = true){
auto rdf_with_cuts = rdf.Filter("detector == \"FD\"").Filter("DC_fiducial_cut_piplus == true");
auto rdf_wo_cuts = rdf.Filter("detector == \"FD\"");

TCanvas canvas("c_X_vs_Y_piplus_FD", "X vs Y for pi+ in FD", 900, 700);

canvas.Divide(1,2);
canvas.cd(1);
auto hist1 = rdf_with_cuts.Histo2D(
    ROOT::RDF::TH2DModel("X_vs_Y_piplus_FD_cuts", "X vs Y for pi+ in FD with DC fiducial cuts; X (cm); Y (cm)", 200, -200, 200, 200, -200, 200),
    "x1_piplus", "y1_piplus"
);
if (logz) gPad->SetLogz();
hist1->Draw("COLZ");
canvas.cd(2);
auto hist2 = rdf_wo_cuts.Histo2D(
    ROOT::RDF::TH2DModel("X_vs_Y_piplus_FD_no_cuts", "X vs Y for pi+ in FD without DC fiducial cuts; X (cm); Y (cm)", 200, -200, 200, 200, -200, 200),
    "x1_piplus", "y1_piplus"
);
if (logz) gPad->SetLogz();
hist2->Draw("COLZ");

canvas.SaveAs((output_folder + "X_vs_Y_piplus_FD_comparison.pdf").c_str());
                           }

/////////////////////////////////////PHASE SPACE HISTOGRAMS/////////////////////////////////////
// 1) For each momentum bin: 2D histogram Theta VS Phi (reconstructed)
void Theta_VS_Phi_per_P_bin(ROOT::RDF::RNode rdf,
                           const std::string& output_folder,
                           const std::string& p_col     = "p_piplus_rec",
                           const std::string& theta_col = "Theta_rec",
                           const std::string& phi_col   = "Phi_rec",
                           int nP = 9, double pMin = 0.5, double pMax = 5.0,
                           int nPhi = 200, double phiMin = -180.0, double phiMax = 180.0,
                           int nTheta = 200, double thMin = 0.0, double thMax = 50.0)
{
     rdf = rdf.Filter("detector==\"FD\"");

    // X=phi, Y=theta, Z=p
    auto h3 = rdf.Histo3D(
        ROOT::RDF::TH3DModel("h3_phi_theta_p__ThetaVsPhiPerP",
                             "Phase space;#phi [deg];#theta [deg];p [GeV]",
                             nPhi, phiMin, phiMax,
                             nTheta, thMin, thMax,
                             nP, pMin, pMax),
        phi_col, theta_col, p_col
    );

    TCanvas c("c_ThetaVsPhiPerP", "Theta vs Phi per P bin", 900, 700);

    const std::string pdf = output_folder + "Theta_VS_Phi_per_P_bin.pdf";
    c.SaveAs((pdf + "[").c_str());

    TAxis* zax = h3->GetZaxis();
    for (int iz = 1; iz <= zax->GetNbins(); ++iz) {
        const double plo = zax->GetBinLowEdge(iz);
        const double phi = zax->GetBinUpEdge(iz);

        h3->GetZaxis()->SetRange(iz, iz);

        TH1*  proj = h3->Project3D("yx"); // x=phi, y=theta
        TH2D* h2   = (TH2D*)proj->Clone(Form("h2_ThetaVsPhi_pbin_%d", iz));
        delete proj;

        h2->SetDirectory(nullptr);
        h2->SetTitle(Form("#theta vs #phi (p #in [%.3f, %.3f] GeV);#phi [deg];#theta [deg]", plo, phi));
        h2->Draw("COLZ");
        c.SaveAs(pdf.c_str());
        delete h2;

        h3->GetZaxis()->SetRange(0, 0);
    }

    c.SaveAs((pdf + "]").c_str());
}


// 2) For each theta bin: 2D histogram P VS Phi (reconstructed)
void P_VS_Phi_per_Theta_bin(ROOT::RDF::RNode rdf,
                           const std::string& output_folder,
                           const std::string& p_col     = "p_piplus_rec",
                           const std::string& theta_col = "Theta_rec",
                           const std::string& phi_col   = "Phi_rec",
                           int nTheta = 40, double thMin = 5, double thMax = 45,
                           int nPhi = 200, double phiMin = -180.0, double phiMax = 180.0,
                           int nP = 200, double pMin = 0.0, double pMax = 5.0)
{
     rdf = rdf.Filter("detector==\"FD\"");

    // X=phi, Y=theta, Z=p
    auto h3 = rdf.Histo3D(
        ROOT::RDF::TH3DModel("h3_phi_theta_p__PvsPhiPerTheta",
                             "Phase space;#phi [deg];#theta [deg];p [GeV]",
                             nPhi, phiMin, phiMax,
                             nTheta, thMin, thMax,
                             nP, pMin, pMax),
        phi_col, theta_col, p_col
    );

    TCanvas c("c_PvsPhiPerTheta", "P vs Phi per Theta bin", 900, 700);

    const std::string pdf = output_folder + "P_VS_Phi_per_Theta_bin.pdf";
    c.SaveAs((pdf + "[").c_str());

    TAxis* yax = h3->GetYaxis();
    for (int iy = 1; iy <= yax->GetNbins(); ++iy) {
        const double tlo = yax->GetBinLowEdge(iy);
        const double thi = yax->GetBinUpEdge(iy);

        h3->GetYaxis()->SetRange(iy, iy);

        TH1*  proj = h3->Project3D("zx"); 
        TH2D* h2   = (TH2D*)proj->Clone(Form("h2_PvsPhi_thetabin_%d", iy));
        delete proj;

        h2->SetDirectory(nullptr);
        h2->SetTitle(Form("p vs #phi (#theta #in [%.3f, %.3f] deg);#phi [deg];p [GeV]", tlo, thi));
        h2->Draw("COLZ");
        c.SaveAs(pdf.c_str());
        delete h2;

        h3->GetYaxis()->SetRange(0, 0);
    }

    c.SaveAs((pdf + "]").c_str());
}


// 3) For each phi bin: 2D histogram P VS Theta (reconstructed)
void P_VS_Theta_per_Phi_bin(ROOT::RDF::RNode rdf,
                           const std::string& output_folder,
                           const std::string& p_col     = "p_piplus_rec",
                           const std::string& theta_col = "Theta_rec",
                           const std::string& phi_col   = "Phi_rec",
                           
                           int nPhi = 72, double phiMin = -180.0, double phiMax = 180.0,
                           int nTheta = 200, double thMin = 0.0, double thMax = 60.0,
                           int nP = 200, double pMin = 0.0, double pMax = 5.0)
{
    
    rdf = rdf.Filter("detector==\"FD\"");

    // X=phi, Y=theta, Z=p
    auto h3 = rdf.Histo3D(
        ROOT::RDF::TH3DModel("h3_phi_theta_p__PvsThetaPerPhi",
                             "Phase space;#phi [deg];#theta [deg];p [GeV]",
                             nPhi, phiMin, phiMax,
                             nTheta, thMin, thMax,
                             nP, pMin, pMax),
        phi_col, theta_col, p_col
    );

    TCanvas c("c_PvsThetaPerPhi", "P vs Theta per Phi bin", 900, 700);

    const std::string pdf = output_folder + "P_VS_Theta_per_Phi_bin.pdf";
    c.SaveAs((pdf + "[").c_str());

    TAxis* xax = h3->GetXaxis();
    for (int ix = 1; ix <= xax->GetNbins(); ++ix) {
        const double plo = xax->GetBinLowEdge(ix);
        const double phi = xax->GetBinUpEdge(ix);

        h3->GetXaxis()->SetRange(ix, ix);

        TH1*  proj = h3->Project3D("yz"); // x=p, y = theta
        TH2D* h2   = (TH2D*)proj->Clone(Form("h2_PvsTheta_phibin_%d", ix));
        delete proj;

        h2->SetDirectory(nullptr);
        h2->SetTitle(Form("p vs #theta (#phi #in [%.3f, %.3f] deg);p [GeV];#theta [deg]", plo, phi));
        h2->Draw("COLZ");
        c.SaveAs(pdf.c_str());
        delete h2;

        h3->GetXaxis()->SetRange(0, 0);
    }

    c.SaveAs((pdf + "]").c_str());
}

//////////////////////////////////////////////////////////// Delta P sliced////////////////////////////////////////////////////////////////
// 1) 2D delta p VS P_rec for different theta slices
void deltaP_VS_Prec_per_Theta_bin(ROOT::RDF::RNode rdf,
                                 const std::string& output_folder,
                                 const std::string& p_col     = "p_piplus_rec",
                                 const std::string& dp_col    = "delta_p",
                                 const std::string& theta_col = "Theta_rec",
                                 int nTheta = 40, double thMin = 5.0, double thMax = 45.0,
                                 int nP = 200, double pMin = 0.0, double pMax = 5.0,
                                 int nDP = 200, double dpMin = -0.1, double dpMax = 0.1)
{
    rdf = rdf.Filter("detector==\"FD\"");

    // X=theta, Y=p_rec, Z=delta_p
    auto h3 = rdf.Histo3D(
        ROOT::RDF::TH3DModel("h3_theta_p_dp__dPvsP_perTheta",
                             "Phase space;#theta [deg];p_{rec} [GeV];#Deltap [GeV]",
                             nTheta, thMin, thMax,
                             nP, pMin, pMax,
                             nDP, dpMin, dpMax),
        theta_col, p_col, dp_col
    );

    TCanvas c("c_dPvsP_perTheta", "delta p vs p_rec per theta bin", 900, 700);

    const std::string pdf = output_folder + "deltaP_VS_Prec_per_Theta_bin.pdf";
    c.SaveAs((pdf + "[").c_str());

    TAxis* xax = h3->GetXaxis();
    for (int ix = 1; ix <= xax->GetNbins(); ++ix) {
        const double tlo = xax->GetBinLowEdge(ix);
        const double thi = xax->GetBinUpEdge(ix);

        h3->GetXaxis()->SetRange(ix, ix);

        // Project to (p_rec, delta_p): y,z in TH3 -> "yz"
        TH1*  proj = h3->Project3D("zy"); // x=p_rec, y=delta_p
        TH2D* h2   = (TH2D*)proj->Clone(Form("h2_dPvsP_theta_%d", ix));
        delete proj;

        h2->SetDirectory(nullptr);
        h2->SetTitle(Form("#Deltap vs p_{rec} (#theta #in [%.2f, %.2f] deg); p_{rec} [GeV];#Deltap [GeV];", tlo, thi));
        h2->Draw("COLZ");
        c.SaveAs(pdf.c_str());
        delete h2;

        h3->GetXaxis()->SetRange(0, 0);
    }

    c.SaveAs((pdf + "]").c_str());
}


// 2) 2D delta p VS P_rec for different phi slices
void deltaP_VS_Prec_per_Phi_bin(ROOT::RDF::RNode rdf,
                               const std::string& output_folder,
                               const std::string& p_col   = "p_piplus_rec",
                               const std::string& dp_col  = "delta_p",
                               const std::string& phi_col = "Phi_rec",
                               int nPhi = 36, double phiMin = -180.0, double phiMax = 180.0,
                               int nP = 200, double pMin = 0.0, double pMax = 5.0,
                               int nDP = 200, double dpMin = -0.1, double dpMax = 0.1)
{
    rdf = rdf.Filter("detector==\"FD\"");

    // X=phi, Y=p_rec, Z=delta_p
    auto h3 = rdf.Histo3D(
        ROOT::RDF::TH3DModel("h3_phi_p_dp__dPvsP_perPhi",
                             "Phase space;#phi [deg];p_{rec} [GeV];#Deltap [GeV]",
                             nPhi, phiMin, phiMax,
                             nP, pMin, pMax,
                             nDP, dpMin, dpMax),
        phi_col, p_col, dp_col
    );

    TCanvas c("c_dPvsP_perPhi", "delta p vs p_rec per phi bin", 900, 700);

    const std::string pdf = output_folder + "deltaP_VS_Prec_per_Phi_bin.pdf";
    c.SaveAs((pdf + "[").c_str());

    TAxis* xax = h3->GetXaxis();
    for (int ix = 1; ix <= xax->GetNbins(); ++ix) {
        const double plo = xax->GetBinLowEdge(ix);
        const double phi = xax->GetBinUpEdge(ix);

        h3->GetXaxis()->SetRange(ix, ix);

        // Project to (p_rec, delta_p): y,z -> "yz"
        TH1*  proj = h3->Project3D("zy"); // x=p_rec, y=delta_p
        TH2D* h2   = (TH2D*)proj->Clone(Form("h2_dPvsP_phi_%d", ix));
        delete proj;

        h2->SetDirectory(nullptr);
        h2->SetTitle(Form("#Deltap vs p_{rec} (#phi #in [%.1f, %.1f] deg); p_{rec} [GeV];#Deltap [GeV];", plo, phi));
        h2->Draw("COLZ");
        c.SaveAs(pdf.c_str());
        delete h2;

        h3->GetXaxis()->SetRange(0, 0);
    }

    c.SaveAs((pdf + "]").c_str());
}

void plot_theta_vs_vz_pion(ROOT::RDF::RNode rdf, const std::string& output_folder, bool logz = true) {

    auto df_FD = rdf.Filter("detector == \"FD\"");
    auto df_CD = rdf.Filter("detector == \"CD\"");

    // Book actions first
    auto h_FD = df_FD.Histo2D(
        {"h_theta_vs_vz_piplus_FD",
         "FD #pi^{+} #theta vs vertex z;v_{z} (cm);#theta_{#pi^{+}} (deg)",
         200, -40.0, 40.0, 200, 0.0, 150.0},
        "vz_piplus", "Theta_rec"
    );

    auto h_CD = df_CD.Histo2D(
        {"h_theta_vs_vz_piplus_CD",
         "CD #pi^{+} #theta vs vertex z;v_{z} (cm);#theta_{#pi^{+}} (deg)",
         200, -20.0, 20.0, 200, 0.0, 180.0},
        "vz_piplus", "Theta_rec"
    );

    // Trigger one event loop for both histograms
    h_FD.GetValue();

    TCanvas c("c_theta_vs_vz_piplus", "c_theta_vs_vz_piplus", 900, 1000);
    c.Divide(1, 2);

    c.cd(1);
    if (logz) gPad->SetLogz();
    h_FD->Draw("COLZ");

    c.cd(2);
    if (logz) gPad->SetLogz();
    h_CD->Draw("COLZ");

    c.SaveAs((output_folder + "/theta_vs_vz_piplus.pdf").c_str());
}

void plot_Vx_VS_Vy_piplus(ROOT::RDF::RNode rdf, const std::string& output_folder, bool logz = false) {

    auto df_FD = rdf.Filter("detector == \"FD\"");
    auto df_CD = rdf.Filter("detector == \"CD\"");

    // Book actions first
    auto h_FD = df_FD.Histo2D(
        {"h_vx_vs_vy_piplus_FD",
         "FD #pi^{+} vertex x vs y;v_{x} (cm);v_{y} (cm)",
         200, -3, 3, 200, -3, 3},
        "vx_piplus", "vy_piplus"
    );

    auto h_CD = df_CD.Histo2D(
        {"h_vx_vs_vy_piplus_CD",
         "CD #pi^{+} vertex x vs y;v_{x} (cm);v_{y} (cm)",
         200, -0.1, -0.04, 200, -0.1, -0.2},
        "vx_piplus", "vy_piplus"
    );

    // Trigger one event loop for both histograms
    h_FD.GetValue();

    TCanvas c("c_vx_vs_vy_piplus", "c_vx_vs_vy_piplus", 900, 1000);
    c.Divide(1, 2);

    c.cd(1);
    if (logz) gPad->SetLogz();
    //h_FD->SetMarkerStyle(20);
    //h_FD->SetMarkerSize(0.4);
    //h_FD->SetMarkerColor(kBlack);
    h_FD->Draw("SCAT");

    c.cd(2);
    if (logz) gPad->SetLogz();
    //h_CD->GetXaxis()->SetRangeUser(-0.2, 0.2);
    //h_CD->GetYaxis()->SetRangeUser(-0.26, -0.06);
    //h_CD->SetMarkerStyle(20);
    //h_CD->SetMarkerSize(0.4);
    //h_CD->SetMarkerColor(kBlack);
    h_CD->Draw("SCAT");

    std::string log_suffix = logz ? "_logz" : "";

    c.SaveAs((output_folder + "/vx_vs_vy_piplus" + log_suffix + "_SCAT.pdf").c_str());
}