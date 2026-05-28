// plot_step2p1_phi_vs_deltaP_from_hipo.cxx
//
// Purpose:
//   Take a Step-2.1 HIPO file from skim_three_step_positive_track_study_singlepass.cxx
//   and plot:
//
//      phi_rec(selected pi+) VS delta_p
//
//   where selected pi+ is reconstructed using the SAME logic as in the 5-step skim:
//
//      1. require exactly one MC::Lund pi+
//      2. among requested-detector positive tracks, choose the one closest to p_Lund(pi+)
//      3. require selected track inside theta/p bin
//      4. require selected track pid == 211
//      5. fill phi_rec(selected pi+) vs delta_p
//
// Output:
//   PNG only
//
// Example:
//
// clas12root -l -b -q 'plot_step2p1_phi_vs_deltaP_from_hipo.cxx+("FD_theta38_39_p1p0_1p2_step2p1_selectedTrackIsPip.hipo","step2p1_phi_vs_deltaP_38_39.png","FD",38,39,1.0,1.2,200,-180,180,200,-0.5,0.5,-1)'

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"

#include "hipo4/reader.h"
#include "hipo4/event.h"
#include "hipo4/bank.h"
#include "hipo4/dictionary.h"

struct TrackInfo {
    int index = -1;
    int pid = 0;
    int charge = 0;
    int status = 0;

    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    double p = 0.0;
    double theta_deg = 0.0;
    double phi_deg = 0.0;
};

static bool get_schema_safe(hipo::dictionary& factory,
                            const char* bank_name,
                            hipo::schema& schema)
{
    try {
        schema = factory.getSchema(bank_name);
        return true;
    }
    catch (...) {
        return false;
    }
}

static bool detector_status_pass(int status_raw, const std::string& det)
{
    const int status = std::abs(status_raw);

    if (det == "FD") return status < 4000;
    if (det == "CD") return status >= 4000;

    return false;
}

static double momentum(double px, double py, double pz)
{
    return std::sqrt(px * px + py * py + pz * pz);
}

static double theta_deg_from_p(double px, double py, double pz)
{
    const double p = momentum(px, py, pz);

    if (p <= 0.0 || !std::isfinite(p)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double c = pz / p;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    return std::acos(c) * TMath::RadToDeg();
}

static double phi_deg_from_px_py(double px, double py)
{
    double phi = std::atan2(py, px) * TMath::RadToDeg();

    if (phi >= 180.0) phi -= 360.0;
    if (phi < -180.0) phi += 360.0;

    return phi;
}

static bool track_in_theta_p_bin(const TrackInfo& trk,
                                 double theta_min_deg,
                                 double theta_max_deg,
                                 double p_min,
                                 double p_max)
{
    if (!std::isfinite(trk.p) || !std::isfinite(trk.theta_deg)) {
        return false;
    }

    if (trk.theta_deg < theta_min_deg || trk.theta_deg >= theta_max_deg) {
        return false;
    }

    if (trk.p < p_min || trk.p >= p_max) {
        return false;
    }

    return true;
}

static bool get_mclund_pip_momentum(hipo::bank& mcLund,
                                    double& p_lund,
                                    int& n_lund_pip)
{
    p_lund = std::numeric_limits<double>::quiet_NaN();
    n_lund_pip = 0;

    for (int i = 0; i < mcLund.getRows(); i++) {
        const int pid = mcLund.getInt("pid", i);

        if (pid != 211) continue;

        const double px = mcLund.getFloat("px", i);
        const double py = mcLund.getFloat("py", i);
        const double pz = mcLund.getFloat("pz", i);

        p_lund = momentum(px, py, pz);
        n_lund_pip++;
    }

    return (n_lund_pip == 1 && std::isfinite(p_lund));
}

void plot_step2p1_phi_vs_deltaP_from_hipo(const char* input_hipo,
                                          const std::string& output_png = "step2p1_phi_vs_deltaP.png",
                                          const char* detector = "FD",
                                          double theta_min_deg = 38.0,
                                          double theta_max_deg = 39.0,
                                          double p_min = 1.0,
                                          double p_max = 1.2,
                                          int nPhi = 200,
                                          double phiMin = -180.0,
                                          double phiMax = 180.0,
                                          int nDP = 200,
                                          double dpMin = -0.5,
                                          double dpMax = 0.5,
                                          int max_events = -1)
{
    std::string det = detector;

    if (det != "FD" && det != "CD") {
        std::cerr << "ERROR: detector must be \"FD\" or \"CD\", got: "
                  << det << "\n";
        return;
    }

    std::cout << "Input HIPO:      " << input_hipo << "\n";
    std::cout << "Output PNG:      " << output_png << "\n";
    std::cout << "Detector:        " << det << "\n";
    std::cout << "Theta bin:       [" << theta_min_deg << ", " << theta_max_deg << ") deg\n";
    std::cout << "Momentum bin:    [" << p_min << ", " << p_max << ") GeV\n";
    std::cout << "Plot:            phi_rec(selected pi+) vs delta_p\n";
    std::cout << "Selection logic: Step-1 chosen positive track closest to MC::Lund pi+ momentum\n";

    hipo::reader reader;

    try {
        reader.open(input_hipo);
    }
    catch (...) {
        std::cerr << "ERROR: could not open input HIPO:\n"
                  << "  " << input_hipo << "\n";
        return;
    }

    if (!reader.is_open()) {
        std::cerr << "ERROR: reader is not open for file:\n"
                  << "  " << input_hipo << "\n";
        return;
    }

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::schema rec_schema;
    hipo::schema mclund_schema;

    if (!get_schema_safe(factory, "REC::Particle", rec_schema)) {
        std::cerr << "ERROR: REC::Particle bank not found.\n";
        return;
    }

    if (!get_schema_safe(factory, "MC::Lund", mclund_schema)) {
        std::cerr << "ERROR: MC::Lund bank not found.\n";
        return;
    }

    hipo::event event;
    hipo::bank recPart(rec_schema);
    hipo::bank mcLund(mclund_schema);

    TH2D* h_phi_dp = new TH2D(
        "h_step2p1_phi_vs_deltaP",
        "Step 2.1: selected #pi^{+};#phi_{rec}(selected #pi^{+}) [deg];#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV]",
        nPhi, phiMin, phiMax,
        nDP, dpMin, dpMax
    );

    long long n_events_scanned = 0;
    long long n_events_filled = 0;

    long long n_fail_lund_pip = 0;
    long long n_fail_no_detector_positive = 0;
    long long n_fail_selected_outside_bin = 0;
    long long n_fail_selected_not_pip = 0;
    long long n_fail_bad_kinematics = 0;

    while (reader.next()) {
        reader.read(event);

        n_events_scanned++;

        if (max_events > 0 && n_events_scanned > max_events) {
            break;
        }

        event.getStructure(recPart);
        event.getStructure(mcLund);

        double p_lund = std::numeric_limits<double>::quiet_NaN();
        int n_lund_pip = 0;

        if (!get_mclund_pip_momentum(mcLund, p_lund, n_lund_pip)) {
            n_fail_lund_pip++;
            continue;
        }

        std::vector<TrackInfo> detector_positive_tracks;

        for (int i = 0; i < recPart.getRows(); i++) {
            const int pid = recPart.getInt("pid", i);
            const int charge = recPart.getInt("charge", i);
            const int status = recPart.getInt("status", i);

            if (charge <= 0) continue;
            if (!detector_status_pass(status, det)) continue;

            const double px = recPart.getFloat("px", i);
            const double py = recPart.getFloat("py", i);
            const double pz = recPart.getFloat("pz", i);

            TrackInfo trk;
            trk.index = i;
            trk.pid = pid;
            trk.charge = charge;
            trk.status = status;
            trk.px = px;
            trk.py = py;
            trk.pz = pz;
            trk.p = momentum(px, py, pz);
            trk.theta_deg = theta_deg_from_p(px, py, pz);
            trk.phi_deg = phi_deg_from_px_py(px, py);

            if (!std::isfinite(trk.p) ||
                !std::isfinite(trk.theta_deg) ||
                !std::isfinite(trk.phi_deg)) {
                continue;
            }

            detector_positive_tracks.push_back(trk);
        }

        if (detector_positive_tracks.empty()) {
            n_fail_no_detector_positive++;
            continue;
        }

        TrackInfo selected_step1;
        double best_abs_dp = std::numeric_limits<double>::infinity();

        for (const auto& trk : detector_positive_tracks) {
            const double abs_dp = std::fabs(trk.p - p_lund);

            if (abs_dp < best_abs_dp) {
                best_abs_dp = abs_dp;
                selected_step1 = trk;
            }
        }

        if (!track_in_theta_p_bin(selected_step1,
                                  theta_min_deg,
                                  theta_max_deg,
                                  p_min,
                                  p_max)) {
            n_fail_selected_outside_bin++;
            continue;
        }

        if (selected_step1.pid != 211) {
            n_fail_selected_not_pip++;
            continue;
        }

        const double delta_p = selected_step1.p - p_lund;
        const double phi = selected_step1.phi_deg;

        if (!std::isfinite(delta_p) || !std::isfinite(phi)) {
            n_fail_bad_kinematics++;
            continue;
        }

        h_phi_dp->Fill(phi, delta_p);
        n_events_filled++;
    }

    std::cout << "\nDone reading.\n";
    std::cout << "Events scanned:                         " << n_events_scanned << "\n";
    std::cout << "Events filled:                          " << n_events_filled << "\n";
    std::cout << "Histogram entries:                      " << h_phi_dp->GetEntries() << "\n";
    std::cout << "Fail exactly 1 MC::Lund pi+:             " << n_fail_lund_pip << "\n";
    std::cout << "Fail no requested-detector positive:     " << n_fail_no_detector_positive << "\n";
    std::cout << "Fail selected outside theta/p bin:       " << n_fail_selected_outside_bin << "\n";
    std::cout << "Fail selected track not pi+:             " << n_fail_selected_not_pip << "\n";
    std::cout << "Fail bad kinematics:                    " << n_fail_bad_kinematics << "\n";

    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("c_step2p1_phi_vs_deltaP",
                             "Step 2.1 phi vs deltaP",
                             1200, 900);

    c->SetTopMargin(0.10);
    c->SetRightMargin(0.16);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    h_phi_dp->Draw("COLZ");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.15, 0.93, Form("Entries = %.0f", h_phi_dp->GetEntries()));

    std::string outname = output_png;

    if (outname.size() < 4 || outname.substr(outname.size() - 4) != ".png") {
        outname += ".png";
    }

    c->SaveAs(outname.c_str());

    std::cout << "\nWrote:\n";
    std::cout << "  " << outname << "\n";

    delete c;
    delete h_phi_dp;
}