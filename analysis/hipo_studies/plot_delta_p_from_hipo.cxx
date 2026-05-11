// clas12root -l -b -q 'plot_delta_p_from_hipo.cxx+("ced_FD_high_loss.hipo","dp_FD_high_loss.png")'
#include <iostream>
#include <cmath>
#include <limits>

#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"

#include "clas12reader.h"
#include "HipoChain.h"

void plot_delta_p_from_hipo(const char* input_hipo,
                            const char* output_png = "delta_p_check.png",
                            int nBins = 200,
                            double dpMin = -0.30,
                            double dpMax = 0.30)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";
    std::cout << "Output PNG: " << output_png << "\n";

    clas12root::HipoChain chain;
    chain.Add(input_hipo);
    chain.db()->turnOffQADB();

    auto config_c12 = chain.GetC12Reader();
    auto& c12 = chain.C12ref();

    TH1D* h_dp = new TH1D("h_dp",
                          "#Delta p sanity check;#Delta p = p_{rec} - p_{gen} (GeV);Counts",
                          nBins, dpMin, dpMax);

    int n_events_scanned = 0;
    int n_events_filled  = 0;
    int n_bad_rec_pip    = 0;
    int n_bad_gen_pip    = 0;

    while (chain.Next()) {
        n_events_scanned++;

        auto pips = c12->getByID(211);

        // For your skimmed diagnostic HIPO files, this should usually be exactly 1.
        if (pips.size() != 1) {
            n_bad_rec_pip++;
            continue;
        }

        auto& pip = pips[0];

        const double p_rec = pip->getP();

        auto mc = c12->mcparts();

        double p_gen = std::numeric_limits<double>::quiet_NaN();
        int n_gen_pip = 0;

        for (int j = 0; j < mc->getRows(); j++) {
            if (mc->getPid(j) != 211) continue;

            const double px = mc->getPx(j);
            const double py = mc->getPy(j);
            const double pz = mc->getPz(j);

            p_gen = std::sqrt(px * px + py * py + pz * pz);
            n_gen_pip++;
        }

        if (n_gen_pip != 1 || !std::isfinite(p_gen)) {
            n_bad_gen_pip++;
            continue;
        }

        const double delta_p = p_rec - p_gen;

        h_dp->Fill(delta_p);
        n_events_filled++;
    }

    std::cout << "\nDone reading.\n"
              << "Events scanned:             " << n_events_scanned << "\n"
              << "Events filled:              " << n_events_filled  << "\n"
              << "Skipped: REC pi+ != 1:      " << n_bad_rec_pip    << "\n"
              << "Skipped: GEN pi+ != 1:      " << n_bad_gen_pip    << "\n";

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c", "delta p check", 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    h_dp->SetLineWidth(2);
    h_dp->Draw("hist");

    TLine* zero = new TLine(0.0, 0.0, 0.0, h_dp->GetMaximum() * 1.05);
    zero->SetLineStyle(2);
    zero->SetLineWidth(2);
    zero->Draw("same");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.15, 0.92, Form("Events filled: %d", n_events_filled));

    c->SaveAs(output_png);

    std::cout << "Saved plot: " << output_png << "\n";
}