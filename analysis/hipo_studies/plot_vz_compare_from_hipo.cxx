// plot_vz_compare_from_hipo.cxx
//
// Run:
// clas12root -l -b -q 'plot_vz_compare_from_hipo.cxx+("high_loss_dp_negative.hipo","low_loss_dp_peak.hipo","vz_compare.png")'

#include <iostream>
#include <fstream>
#include <string>

#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

#include "clas12reader.h"
#include "HipoChain.h"

bool file_exists(const char* path)
{
    std::ifstream f(path);
    return f.good();
}

int fill_vz_hist_from_hipo(const char* input_hipo,
                           TH1D* h_vz,
                           bool require_exactly_one_pip = true)
{
    if (!file_exists(input_hipo)) {
        std::cerr << "ERROR: Cannot open file: " << input_hipo << "\n";
        return -1;
    }

    clas12root::HipoChain chain;
    chain.Add(input_hipo);
    chain.db()->turnOffQADB();

    auto config_c12 = chain.GetC12Reader();
    auto& c12 = chain.C12ref();

    int n_events_scanned = 0;
    int n_events_filled  = 0;
    int n_bad_pip_count  = 0;

    while (chain.Next()) {
        n_events_scanned++;

        auto pips = c12->getByID(211);

        if (require_exactly_one_pip && pips.size() != 1) {
            n_bad_pip_count++;
            continue;
        }

        if (pips.size() < 1) {
            n_bad_pip_count++;
            continue;
        }

        // If exactly one pi+ is required, this is that pi+.
        // If not, this uses the first reconstructed pi+.
        auto& pip = pips[0];

        const double vz = pip->par()->getVz();

        h_vz->Fill(vz);
        n_events_filled++;
    }

    std::cout << "\nFile: " << input_hipo << "\n";
    std::cout << "Events scanned:        " << n_events_scanned << "\n";
    std::cout << "Events filled:         " << n_events_filled  << "\n";
    std::cout << "Skipped bad pi count:  " << n_bad_pip_count  << "\n";

    return n_events_filled;
}

void plot_vz_compare_from_hipo(const char* high_loss_hipo,
                               const char* low_loss_hipo,
                               const char* output_png = "vz_compare_high_low_loss.png",
                               int nBins = 160,
                               double vzMin = -20.0,
                               double vzMax = 10.0,
                               bool normalize = false)
{
    std::cout << "High-loss HIPO: " << high_loss_hipo << "\n";
    std::cout << "Low-loss HIPO:  " << low_loss_hipo  << "\n";
    std::cout << "Output PNG:     " << output_png     << "\n";

    TH1D* h_high = new TH1D("h_vz_high",
                            "Reconstructed #pi^{+} vertex z;v_{z}^{rec}(#pi^{+}) [cm];Counts",
                            nBins, vzMin, vzMax);

    TH1D* h_low = new TH1D("h_vz_low",
                           "Reconstructed #pi^{+} vertex z;v_{z}^{rec}(#pi^{+}) [cm];Counts",
                           nBins, vzMin, vzMax);

    const int n_high = fill_vz_hist_from_hipo(high_loss_hipo, h_high, true);
    const int n_low  = fill_vz_hist_from_hipo(low_loss_hipo,  h_low,  true);

    if (n_high <= 0 && n_low <= 0) {
        std::cerr << "ERROR: Both histograms are empty. Nothing to plot.\n";
        return;
    }

    if (normalize) {
        if (h_high->Integral() > 0) h_high->Scale(1.0 / h_high->Integral());
        if (h_low->Integral()  > 0) h_low->Scale(1.0 / h_low->Integral());

        h_high->GetYaxis()->SetTitle("Normalized counts");
        h_low->GetYaxis()->SetTitle("Normalized counts");
    }

    gStyle->SetOptStat(1110);

TCanvas* c = new TCanvas("c", "Vz high vs low loss", 1400, 650);
c->Divide(2, 1);

// ---------------- left panel: high loss ----------------
c->cd(1);
gPad->SetTopMargin(0.12);
gPad->SetRightMargin(0.05);
gPad->SetLeftMargin(0.13);
gPad->SetBottomMargin(0.13);

h_high->SetLineColor(kRed + 1);
h_high->SetLineWidth(2);
h_high->SetTitle(Form("High-loss sample;v_{z}^{rec}(#pi^{+}) [cm];%s",
                      normalize ? "Normalized counts" : "Counts"));
h_high->Draw("hist");

TLatex lat_high;
lat_high.SetNDC();
lat_high.SetTextSize(0.045);
lat_high.DrawLatex(0.16, 0.92, Form("High loss, N = %d", n_high));

// ---------------- right panel: low loss ----------------
c->cd(2);
gPad->SetTopMargin(0.12);
gPad->SetRightMargin(0.05);
gPad->SetLeftMargin(0.13);
gPad->SetBottomMargin(0.13);

h_low->SetLineColor(kBlue + 1);
h_low->SetLineWidth(2);
h_low->SetTitle(Form("Low-loss sample;v_{z}^{rec}(#pi^{+}) [cm];%s",
                     normalize ? "Normalized counts" : "Counts"));
h_low->Draw("hist");

TLatex lat_low;
lat_low.SetNDC();
lat_low.SetTextSize(0.045);
lat_low.DrawLatex(0.16, 0.92, Form("Low loss, N = %d", n_low));

c->SaveAs(output_png);

    std::cout << "\nSaved plot: " << output_png << "\n";
}