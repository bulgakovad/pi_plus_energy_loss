// plot_delta_p_rec_minus_gen.cxx
//
// Purpose:
//   Read one HIPO file containing:
//     - MC::Lund generated particles
//     - REC::Particle reconstructed particles
//
//   For each event:
//     1. Find generated pi+ from MC::Lund, pid = 211
//     2. Compute p_gen = sqrt(px^2 + py^2 + pz^2)
//     3. Check REC::Particle exists
//     4. Require trigger electron: pid = 11 and status < 0
//     5. Require at least one reconstructed pi+: pid = 211
//     6. If multiple REC pi+, choose the one whose p_rec is closest to p_gen
//     7. Fill delta p = p_rec - p_gen
//
// Run:
//   clas12root -l -b -q 'plot_delta_p_rec_minus_gen.cxx+("input.hipo")'
//
// Optional:
//   clas12root -l -b -q 'plot_delta_p_rec_minus_gen.cxx+("input.hipo","delta_p.png","delta_p.root")'

#include <iostream>
#include <cmath>
#include <string>
#include <limits>

#include "hipo4/reader.h"
#include "hipo4/dictionary.h"
#include "hipo4/event.h"
#include "hipo4/bank.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"


static double mag3(double px, double py, double pz)
{
    return std::sqrt(px * px + py * py + pz * pz);
}


void plot_delta_p_rec_minus_gen(const char* input_hipo,
                                const char* output_png  = "delta_p_rec_minus_gen.png",
                                const char* output_root = "delta_p_rec_minus_gen.root",
                                int n_bins = 200,
                                double dp_min = -0.5,
                                double dp_max = 0.5)
{
    std::cout << "\nInput HIPO file: " << input_hipo << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank mc_lund(factory.getSchema("MC::Lund"));
    hipo::bank rec_part(factory.getSchema("REC::Particle"));

    TH1D* h_dp = new TH1D("h_delta_p",
                          "#Delta p = p_{rec} - p_{gen};#Delta p [GeV];Counts",
                          n_bins, dp_min, dp_max);

    Long64_t n_events_total = 0;
    Long64_t n_events_with_gen_pip = 0;
    Long64_t n_events_no_rec_bank = 0;
    Long64_t n_events_with_trigger_e = 0;
    Long64_t n_events_with_rec_pip = 0;
    Long64_t n_events_filled = 0;

    while (reader.next(event)) {
        n_events_total++;

        event.getStructure(mc_lund);
        event.getStructure(rec_part);

        // ------------------------------------------------------------
        // 1. Find generated pi+ in MC::Lund
        // ------------------------------------------------------------
        bool found_gen_pip = false;
        double p_gen = std::numeric_limits<double>::quiet_NaN();

        for (int i = 0; i < mc_lund.getRows(); i++) {
            int pid = mc_lund.getInt("pid", i);

            if (pid != 211)
                continue;

            double px = mc_lund.getFloat("px", i);
            double py = mc_lund.getFloat("py", i);
            double pz = mc_lund.getFloat("pz", i);

            p_gen = mag3(px, py, pz);
            found_gen_pip = true;

            // Your generated sample should have exactly one pi+,
            // so once we find it, stop.
            break;
        }

        if (!found_gen_pip)
            continue;

        n_events_with_gen_pip++;

        // ------------------------------------------------------------
        // 2. Require REC::Particle bank
        // ------------------------------------------------------------
        if (rec_part.getRows() <= 0) {
            n_events_no_rec_bank++;
            continue;
        }

        // ------------------------------------------------------------
        // 3. Scan REC::Particle
        //    Need:
        //      - trigger electron: pid = 11 and status < 0
        //      - at least one reconstructed pi+: pid = 211
        // ------------------------------------------------------------
        bool has_trigger_e = false;
        bool has_rec_pip = false;

        double best_p_rec = std::numeric_limits<double>::quiet_NaN();
        double best_abs_diff = std::numeric_limits<double>::max();

        int n_rec_pip = 0;

        for (int i = 0; i < rec_part.getRows(); i++) {
            int pid = rec_part.getInt("pid", i);
            int status = rec_part.getShort("status", i);

            double px = rec_part.getFloat("px", i);
            double py = rec_part.getFloat("py", i);
            double pz = rec_part.getFloat("pz", i);

            double p_rec = mag3(px, py, pz);

            if (pid == 11 && status < 0) {
                has_trigger_e = true;
            }

            if (pid == 211) {
                has_rec_pip = true;
                n_rec_pip++;

                double abs_diff = std::abs(p_rec - p_gen);

                if (abs_diff < best_abs_diff) {
                    best_abs_diff = abs_diff;
                    best_p_rec = p_rec;
                }
            }
        }

        if (!has_trigger_e)
            continue;

        n_events_with_trigger_e++;

        if (!has_rec_pip)
            continue;

        n_events_with_rec_pip++;

        // ------------------------------------------------------------
        // 4. Fill delta p
        // ------------------------------------------------------------
        double delta_p = best_p_rec - p_gen;
        h_dp->Fill(delta_p);

        n_events_filled++;
    }

    // ------------------------------------------------------------
    // Print summary
    // ------------------------------------------------------------
    std::cout << "\n========== Summary ==========\n";
    std::cout << "Total events scanned:              " << n_events_total << "\n";
    std::cout << "Events with gen pi+:               " << n_events_with_gen_pip << "\n";
    std::cout << "Events with missing/empty REC::Particle bank:   " << n_events_no_rec_bank << " missing/empty\n";
    std::cout << "Events with trigger electron:       " << n_events_with_trigger_e << "\n";
    std::cout << "Events with >=1 REC pi+:            " << n_events_with_rec_pip << "\n";
    std::cout << "Events filled in delta-p histogram: " << n_events_filled << "\n";

    if (n_events_filled > 0) {
        std::cout << "\nHistogram stats:\n";
        std::cout << "Mean delta p:  " << h_dp->GetMean() << " GeV\n";
        std::cout << "RMS delta p:   " << h_dp->GetRMS()  << " GeV\n";
    }

    // ------------------------------------------------------------
    // Save ROOT output
    // ------------------------------------------------------------
    TFile* fout = new TFile(output_root, "RECREATE");
    h_dp->Write();
    fout->Close();

    // ------------------------------------------------------------
    // Draw PNG/PDF-style plot
    // ------------------------------------------------------------
    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c", "delta p", 900, 700);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.12);

    h_dp->SetLineWidth(2);
    h_dp->Draw("hist");

    TLine* line_zero = new TLine(0.0, 0.0, 0.0, h_dp->GetMaximum() * 1.05);
    line_zero->SetLineStyle(2);
    line_zero->SetLineWidth(2);
    line_zero->Draw("same");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.16, 0.86, Form("Filled events: %lld", n_events_filled));


    c->SaveAs(output_png);

    std::cout << "\nSaved histogram ROOT file: " << output_root << "\n";
    std::cout << "Saved plot:                " << output_png << "\n\n";
}