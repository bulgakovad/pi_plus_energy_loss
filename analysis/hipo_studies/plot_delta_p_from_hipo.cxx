// plot_delta_p_from_hipo.cxx
//
// Example:
// clas12root -l -b -q 'plot_delta_p_from_hipo.cxx+("FD_posTrack_noVz_allDp.hipo","dp_posTrack.png")'

#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>

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

    std::ifstream test(input_hipo);
    if (!test.good()) {
        std::cerr << "ERROR: Cannot open input HIPO file:\n"
                  << "  " << input_hipo << "\n";
        return;
    }

    clas12root::HipoChain chain;
    chain.Add(input_hipo);
    chain.db()->turnOffQADB();

    auto config_c12 = chain.GetC12Reader();
    auto& c12 = chain.C12ref();

    TH1D* h_dp = new TH1D(
        "h_dp",
        "#Delta p for positive REC tracks;#Delta p = p_{rec}(positive track) - p_{gen}(#pi^{+}) [GeV];Counts",
        nBins, dpMin, dpMax
    );

    int n_events_scanned = 0;
    int n_events_with_good_gen_pip = 0;
    int n_positive_tracks_found = 0;
    int n_hist_fills = 0;

    int n_bad_gen_pip = 0;
    int n_events_no_positive_track = 0;

    while (chain.Next()) {
        n_events_scanned++;

        // ------------------------------------------------------------
        // Find generated pi+
        // ------------------------------------------------------------
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

        n_events_with_good_gen_pip++;

        // ------------------------------------------------------------
        // Loop over all reconstructed positive tracks
        // ------------------------------------------------------------
        auto rec_particles = c12->getDetParticles();

        int n_positive_this_event = 0;

        for (auto& part : rec_particles) {
            const int charge = part->par()->getCharge();

            if (charge <= 0) continue;

            const int pid_rec = part->par()->getPid();
            const int status  = part->par()->getStatus();

            const double p_rec = part->getP();
            const double delta_p = p_rec - p_gen;

            h_dp->Fill(delta_p);

            n_positive_tracks_found++;
            n_positive_this_event++;
            n_hist_fills++;

            if (n_hist_fills <= 20) {
                std::cout << "FILL"
                          << "  run=" << c12->runconfig()->getRun()
                          << "  event=" << c12->runconfig()->getEvent()
                          << "  rec_pid=" << pid_rec
                          << "  charge=" << charge
                          << "  status=" << status
                          << "  p_rec=" << p_rec
                          << "  p_gen_pi+=" << p_gen
                          << "  dp=" << delta_p
                          << "\n";
            }
        }

        if (n_positive_this_event == 0) {
            n_events_no_positive_track++;
        }
    }

    std::cout << "\nDone reading.\n"
              << "Events scanned:                  " << n_events_scanned << "\n"
              << "Events with good GEN pi+:         " << n_events_with_good_gen_pip << "\n"
              << "Positive REC tracks found:        " << n_positive_tracks_found << "\n"
              << "Histogram fills:                  " << n_hist_fills << "\n"
              << "Skipped: GEN pi+ != 1:            " << n_bad_gen_pip << "\n"
              << "Events with no positive REC track:" << n_events_no_positive_track << "\n";

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c", "delta p positive tracks", 900, 700);
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
    lat.DrawLatex(0.15, 0.92, Form("Positive-track fills: %d", n_hist_fills));

    c->SaveAs(output_png);

    std::cout << "Saved plot: " << output_png << "\n";
}