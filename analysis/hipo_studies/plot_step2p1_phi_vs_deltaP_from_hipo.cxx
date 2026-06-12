// plot_step2p1_phi_vs_deltaP_from_hipo.cxx
//
// Purpose:
//   Take a .dat file containing paths to HIPO files and reproduce the skim-code
//   selection logic up to Step 2.1, then plot selected pi+ delta_p information.
//
// Selection logic reproduced from skim_three_step_positive_track_study_singlepass.cxx:
//
// Base:
//   - exactly 1 reconstructed electron total: pid == 11
//   - that one electron must be trigger electron: status < 0
//   - no negative hadrons: no REC::Particle with charge < 0 except pid == 11
//   - exactly 1 generated pi+ from MC::Lund
//
// Step 1:
//   - among requested-detector positive tracks, choose the one closest to p_Lund(pi+)
//   - require selected track inside theta/p bin
//
// Step 2.1:
//   - require selected Step-1 track pid == 211
//
// Plot modes:
//
//   plot_mode = "2D":
//      make one PNG:
//         phi(selected pi+) VS delta_p
//
//   plot_mode = "1D":
//      make one multipage PDF:
//         one delta_p histogram per phi bin
//         phi bins: [-180, 180) deg, bin size 2.5 deg
//         total pages: 144
//
// Phi definitions:
//
//   phi_definition = "REC":
//      phi = atan2(py_rec, px_rec)
//
//   phi_definition = "DC":
//      phi = atan2(y_DC, x_DC)
//      x_DC,y_DC are taken from REC::Traj for selected_step1.index,
//      detector = dc_detector_id, layer = dc_layer
//
// Examples:
//
// 2D, REC phi:
//
// clas12root -l -b -q 'plot_step2p1_phi_vs_deltaP_from_hipo.cxx+("good_hipo.dat","step2p1_phi_vs_deltaP_REC_38_39.png","FD",38,39,1.0,1.2,"REC",6,6,"2D",200,-180,180,200,-0.5,0.5,-1)'
//
// 2D, DC phi, exact DC layer 6:
//
// clas12root -l -b -q 'plot_step2p1_phi_vs_deltaP_from_hipo.cxx+("good_hipo.dat","step2p1_phi_vs_deltaP_DC_layer6_38_39.png","FD",38,39,1.0,1.2,"DC",6,6,"2D",200,-180,180,200,-0.5,0.5,-1)'
//
// 2D, DC phi, first available DC layer:
//
// clas12root -l -b -q 'plot_step2p1_phi_vs_deltaP_from_hipo.cxx+("good_hipo.dat","step2p1_phi_vs_deltaP_DC_firstLayer_38_39.png","FD",38,39,1.0,1.2,"DC",6,-1,"2D",200,-180,180,200,-0.5,0.5,-1)'
//
// 1D phi-bin histograms, DC phi:
//
// clas12root -l -b -q 'plot_step2p1_phi_vs_deltaP_from_hipo.cxx+("good_hipo.dat","step2p1_deltaP_per_phiBin_DC_38_39.pdf","FD",38,39,1.0,1.2,"DC",6,6,"1D",200,-180,180,200,-0.5,0.5,-1)'

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLine.h"

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

static bool ends_with(const std::string& s, const std::string& suffix)
{
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static bool file_exists(const std::string& path)
{
    std::ifstream f(path);
    return f.good();
}

static std::vector<std::string> get_hipo_files_from_dat(const std::string& dat_file)
{
    std::vector<std::string> hipo_files;

    std::ifstream fin(dat_file);

    if (!fin.is_open()) {
        std::cerr << "ERROR: cannot open .dat file: " << dat_file << "\n";
        return hipo_files;
    }

    std::string line;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (!ends_with(line, ".hipo")) continue;

        hipo_files.push_back(line);
    }

    return hipo_files;
}

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

static bool get_phi_dc_from_rec_traj(hipo::bank& recTraj,
                                     int pindex,
                                     int dc_detector_id,
                                     int dc_layer,
                                     double& phi_dc_deg)
{
    phi_dc_deg = std::numeric_limits<double>::quiet_NaN();

    bool found = false;
    int best_layer = 999999;

    for (int i = 0; i < recTraj.getRows(); i++) {
        const int row_pindex = recTraj.getInt("pindex", i);
        if (row_pindex != pindex) continue;

        const int detector = recTraj.getInt("detector", i);
        if (detector != dc_detector_id) continue;

        const int layer = recTraj.getInt("layer", i);

        // If dc_layer > 0, require exact layer.
        // If dc_layer <= 0, use the smallest available DC layer for this pindex.
        if (dc_layer > 0) {
            if (layer != dc_layer) continue;
        }
        else {
            if (layer >= best_layer) continue;
        }

        const double x = recTraj.getFloat("x", i);
        const double y = recTraj.getFloat("y", i);

        phi_dc_deg = std::atan2(y, x) * TMath::RadToDeg();

        if (phi_dc_deg >= 180.0) phi_dc_deg -= 360.0;
        if (phi_dc_deg < -180.0) phi_dc_deg += 360.0;

        best_layer = layer;
        found = std::isfinite(phi_dc_deg);
    }

    return found;
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

static std::string force_extension(const std::string& name,
                                   const std::string& ext)
{
    if (name.size() >= ext.size() &&
        name.substr(name.size() - ext.size()) == ext) {
        return name;
    }

    return name + ext;
}

static int get_phi_bin_1d(double phi_deg,
                          double phi_min,
                          double phi_bin_width,
                          int n_phi_bins)
{
    if (!std::isfinite(phi_deg)) return -1;

    if (phi_deg >= 180.0) phi_deg -= 360.0;
    if (phi_deg < -180.0) phi_deg += 360.0;

    const int bin = static_cast<int>(std::floor((phi_deg - phi_min) / phi_bin_width));

    if (bin < 0 || bin >= n_phi_bins) return -1;

    return bin;
}

void plot_step2p1_phi_vs_deltaP_from_hipo(const char* input_dat_file,
                                          const std::string& output_file = "step2p1_phi_vs_deltaP.png",
                                          const char* detector = "FD",
                                          double theta_min_deg = 38.0,
                                          double theta_max_deg = 39.0,
                                          double p_min = 1.0,
                                          double p_max = 1.2,
                                          const char* phi_definition = "REC",
                                          int dc_detector_id = 6,
                                          int dc_layer = 6,
                                          const char* plot_mode = "2D",
                                          int nPhi = 200,
                                          double phiMin = -180.0,
                                          double phiMax = 180.0,
                                          int nDP = 200,
                                          double dpMin = -0.5,
                                          double dpMax = 0.5,
                                          int max_events = -1)
{
    std::string det = detector;
    std::string phi_def = phi_definition;
    std::string mode = plot_mode;

    if (det != "FD" && det != "CD") {
        std::cerr << "ERROR: detector must be \"FD\" or \"CD\", got: "
                  << det << "\n";
        return;
    }

    if (phi_def != "REC" && phi_def != "DC") {
        std::cerr << "ERROR: phi_definition must be \"REC\" or \"DC\", got: "
                  << phi_def << "\n";
        return;
    }

    if (mode != "2D" && mode != "1D") {
        std::cerr << "ERROR: plot_mode must be \"2D\" or \"1D\", got: "
                  << mode << "\n";
        return;
    }

    std::cout << "Input .dat file: " << input_dat_file << "\n";
    std::cout << "Output file:     " << output_file << "\n";
    std::cout << "Detector:        " << det << "\n";
    std::cout << "Theta bin:       [" << theta_min_deg << ", " << theta_max_deg << ") deg\n";
    std::cout << "Momentum bin:    [" << p_min << ", " << p_max << ") GeV\n";
    std::cout << "Phi definition:  " << phi_def << "\n";
    std::cout << "Plot mode:       " << mode << "\n";
    std::cout << "Max events:      " << max_events << " (-1 or 0 means all)\n";

    if (phi_def == "DC") {
        std::cout << "DC detector id:  " << dc_detector_id << "\n";
        std::cout << "DC layer:        " << dc_layer
                  << " (-1 means smallest available DC layer)\n";
    }

    if (mode == "1D") {
        std::cout << "1D phi bins:     [-180, 180) deg, bin width 2.5 deg, N=144\n";
    }

    std::cout << "Selection logic: skim-code base + Step 1 + Step 2.1\n";

    auto hipo_files_raw = get_hipo_files_from_dat(input_dat_file);

    if (hipo_files_raw.empty()) {
        std::cerr << "ERROR: no .hipo files found in " << input_dat_file << "\n";
        return;
    }

    std::vector<std::string> hipo_files;

    for (const auto& f : hipo_files_raw) {
        if (!file_exists(f)) {
            std::cerr << "WARNING: skipping unreadable file:\n"
                      << "  " << f << "\n";
            continue;
        }

        hipo_files.push_back(f);
    }

    if (hipo_files.empty()) {
        std::cerr << "ERROR: no readable .hipo files left.\n";
        return;
    }

    std::cout << "Readable HIPO files: " << hipo_files.size() << "\n";

    TH2D* h_phi_dp = nullptr;

    const double phi1D_min = -180.0;
    const double phi1D_width = 2.5;
    const int n_phi_1d_bins = 144;

    std::vector<TH1D*> h_dp_phi_bins;

    if (mode == "2D") {
        h_phi_dp = new TH2D(
            "h_step2p1_phi_vs_deltaP",
            Form("Step 2.1: selected #pi^{+}, #phi definition = %s;#phi(selected #pi^{+}) [deg];#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV]",
                 phi_def.c_str()),
            nPhi, phiMin, phiMax,
            nDP, dpMin, dpMax
        );
    }
    else if (mode == "1D") {
        h_dp_phi_bins.reserve(n_phi_1d_bins);

        for (int ibin = 0; ibin < n_phi_1d_bins; ibin++) {
            const double lo = phi1D_min + ibin * phi1D_width;
            const double hi = lo + phi1D_width;

            TH1D* h = new TH1D(
                Form("h_dp_phi_bin_%03d", ibin),
                Form("Step 2.1: %.1f #leq #phi < %.1f deg, #phi definition = %s;#Delta p = p_{rec} - p_{Lund} [GeV];Counts",
                     lo, hi, phi_def.c_str()),
                nDP, dpMin, dpMax
            );

            h_dp_phi_bins.push_back(h);
        }
    }

    long long n_events_scanned = 0;
    long long n_events_filled = 0;

    long long n_files_opened = 0;
    long long n_files_skipped_open = 0;
    long long n_files_skipped_missing_rec = 0;
    long long n_files_skipped_missing_lund = 0;
    long long n_files_skipped_missing_traj = 0;

    // Skim-code-matching cut-flow counters
    long long n_fail_base_electron = 0;
    long long n_fail_negative_hadrons = 0;
    long long n_fail_lund_pip = 0;
    long long n_fail_no_detector_positive = 0;
    long long n_fail_selected_outside_bin = 0;
    long long n_fail_selected_not_pip = 0;

    long long n_base_one_trigger_electron = 0;
    long long n_base_no_negative_hadrons = 0;
    long long n_base_mclund_pip = 0;
    long long n_step1_has_detector_positive = 0;
    long long n_step1_pass = 0;
    long long n_step21_pass = 0;

    // Plotting-only failure counters
    long long n_fail_bad_kinematics = 0;
    long long n_fail_phi_bin = 0;

    for (size_t ifile = 0; ifile < hipo_files.size(); ifile++) {
        const std::string& input_hipo = hipo_files[ifile];

        std::cout << "\nOpening file " << (ifile + 1)
                  << " / " << hipo_files.size()
                  << ": " << input_hipo << "\n";

        hipo::reader reader;

        try {
            reader.open(input_hipo.c_str());
        }
        catch (...) {
            std::cerr << "WARNING: failed to open, skipping:\n"
                      << "  " << input_hipo << "\n";
            n_files_skipped_open++;
            continue;
        }

        if (!reader.is_open()) {
            std::cerr << "WARNING: reader is not open, skipping:\n"
                      << "  " << input_hipo << "\n";
            n_files_skipped_open++;
            continue;
        }

        n_files_opened++;

        hipo::dictionary factory;
        reader.readDictionary(factory);

        hipo::schema rec_schema;
        hipo::schema mclund_schema;
        hipo::schema rec_traj_schema;

        if (!get_schema_safe(factory, "REC::Particle", rec_schema)) {
            std::cerr << "WARNING: REC::Particle bank not found, skipping file:\n"
                      << "  " << input_hipo << "\n";
            n_files_skipped_missing_rec++;
            continue;
        }

        if (!get_schema_safe(factory, "MC::Lund", mclund_schema)) {
            std::cerr << "WARNING: MC::Lund bank not found, skipping file:\n"
                      << "  " << input_hipo << "\n";
            n_files_skipped_missing_lund++;
            continue;
        }

        const bool have_rec_traj =
            get_schema_safe(factory, "REC::Traj", rec_traj_schema);

        if (phi_def == "DC" && !have_rec_traj) {
            std::cerr << "WARNING: phi_definition=\"DC\" requested, but REC::Traj missing. Skipping file:\n"
                      << "  " << input_hipo << "\n";
            n_files_skipped_missing_traj++;
            continue;
        }

        hipo::event event;
        hipo::bank recPart(rec_schema);
        hipo::bank mcLund(mclund_schema);
        hipo::bank recTraj(rec_traj_schema);

        while (reader.next()) {
            reader.read(event);

            n_events_scanned++;

            if (max_events > 0 && n_events_scanned > max_events) {
                break;
            }

            event.getStructure(recPart);
            event.getStructure(mcLund);

            if (have_rec_traj) {
                event.getStructure(recTraj);
            }

            int n_electron_total = 0;
            int n_trigger_electron = 0;
            int n_negative_hadrons_total = 0;

            std::vector<TrackInfo> detector_positive_tracks;

            for (int i = 0; i < recPart.getRows(); i++) {
                const int pid = recPart.getInt("pid", i);
                const int charge = recPart.getInt("charge", i);
                const int status = recPart.getInt("status", i);

                const double px = recPart.getFloat("px", i);
                const double py = recPart.getFloat("py", i);
                const double pz = recPart.getFloat("pz", i);

                if (pid == 11) {
                    n_electron_total++;

                    if (status < 0) {
                        n_trigger_electron++;
                    }
                }

                if (charge < 0 && pid != 11) {
                    n_negative_hadrons_total++;
                }

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

                if (charge > 0 && detector_status_pass(status, det)) {
                    if (!std::isfinite(trk.p) ||
                        !std::isfinite(trk.theta_deg) ||
                        !std::isfinite(trk.phi_deg)) {
                        continue;
                    }

                    detector_positive_tracks.push_back(trk);
                }
            }

            // --------------------------------------------------------
            // Base selection: same as skim code.
            // exactly one REC electron total, and it is trigger electron.
            // --------------------------------------------------------
            if (!(n_electron_total == 1 && n_trigger_electron == 1)) {
                n_fail_base_electron++;
                continue;
            }

            n_base_one_trigger_electron++;

            // --------------------------------------------------------
            // Base selection: no negative hadrons.
            // --------------------------------------------------------
            if (n_negative_hadrons_total != 0) {
                n_fail_negative_hadrons++;
                continue;
            }

            n_base_no_negative_hadrons++;

            // --------------------------------------------------------
            // Base selection: exactly one MC::Lund pi+.
            // --------------------------------------------------------
            double p_lund = std::numeric_limits<double>::quiet_NaN();
            int n_lund_pip = 0;

            if (!get_mclund_pip_momentum(mcLund, p_lund, n_lund_pip)) {
                n_fail_lund_pip++;
                continue;
            }

            n_base_mclund_pip++;

            // --------------------------------------------------------
            // Step 1: choose requested-detector positive track closest
            // to p_Lund(pi+), then require theta/p bin.
            // --------------------------------------------------------
            if (detector_positive_tracks.empty()) {
                n_fail_no_detector_positive++;
                continue;
            }

            n_step1_has_detector_positive++;

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

            n_step1_pass++;

            // --------------------------------------------------------
            // Step 2.1: selected Step-1 track is pi+.
            // --------------------------------------------------------
            if (selected_step1.pid != 211) {
                n_fail_selected_not_pip++;
                continue;
            }

            n_step21_pass++;

            // --------------------------------------------------------
            // Plotting part.
            // For REC phi, this should not add extra rejection.
            // For DC phi, this can reject events if no matching REC::Traj hit exists.
            // --------------------------------------------------------
            const double delta_p = selected_step1.p - p_lund;

            double phi = std::numeric_limits<double>::quiet_NaN();

            if (phi_def == "REC") {
                phi = selected_step1.phi_deg;
            }
            else if (phi_def == "DC") {
                const bool got_phi_dc =
                    get_phi_dc_from_rec_traj(recTraj,
                                              selected_step1.index,
                                              dc_detector_id,
                                              dc_layer,
                                              phi);

                if (!got_phi_dc) {
                    n_fail_bad_kinematics++;
                    continue;
                }
            }

            if (!std::isfinite(delta_p) || !std::isfinite(phi)) {
                n_fail_bad_kinematics++;
                continue;
            }

            if (mode == "2D") {
                h_phi_dp->Fill(phi, delta_p);
            }
            else if (mode == "1D") {
                const int phi_bin =
                    get_phi_bin_1d(phi, phi1D_min, phi1D_width, n_phi_1d_bins);

                if (phi_bin < 0) {
                    n_fail_phi_bin++;
                    continue;
                }

                h_dp_phi_bins[phi_bin]->Fill(delta_p);
            }

            n_events_filled++;
        }

        if (max_events > 0 && n_events_scanned > max_events) {
            std::cout << "\nReached max_events = " << max_events << "\n";
            break;
        }
    }

    std::cout << "\nDone reading.\n";
    std::cout << "Files in .dat:                          " << hipo_files.size() << "\n";
    std::cout << "Files opened:                           " << n_files_opened << "\n";
    std::cout << "Files skipped open/read error:           " << n_files_skipped_open << "\n";
    std::cout << "Files skipped missing REC::Particle:     " << n_files_skipped_missing_rec << "\n";
    std::cout << "Files skipped missing MC::Lund:          " << n_files_skipped_missing_lund << "\n";
    std::cout << "Files skipped missing REC::Traj:         " << n_files_skipped_missing_traj << "\n";

    std::cout << "\nSkim-code-matching cut flow:\n";
    std::cout << "Events scanned:                         " << n_events_scanned << "\n";
    std::cout << "Base: exactly 1 trigger electron only:  " << n_base_one_trigger_electron << "\n";
    std::cout << "Base: no negative hadrons:              " << n_base_no_negative_hadrons << "\n";
    std::cout << "Base: exactly 1 MC::Lund pi+:           " << n_base_mclund_pip << "\n";
    std::cout << "Step 1: has detector-positive track:    " << n_step1_has_detector_positive << "\n";
    std::cout << "Step 1 pass:                            " << n_step1_pass << "\n";
    std::cout << "Step 2.1 pass before plotting:          " << n_step21_pass << "\n";

    std::cout << "\nPlotting result:\n";
    std::cout << "Events filled into histogram:           " << n_events_filled << "\n";

    if (mode == "2D") {
        std::cout << "Histogram entries:                      " << h_phi_dp->GetEntries() << "\n";
    }
    else if (mode == "1D") {
        double total_entries_1d = 0.0;

        for (int i = 0; i < n_phi_1d_bins; i++) {
            total_entries_1d += h_dp_phi_bins[i]->GetEntries();
        }

        std::cout << "Total 1D histogram entries:             " << total_entries_1d << "\n";
    }

    std::cout << "\nFailure counters:\n";
    std::cout << "Fail base electron requirement:          " << n_fail_base_electron << "\n";
    std::cout << "Fail no-negative-hadrons requirement:    " << n_fail_negative_hadrons << "\n";
    std::cout << "Fail exactly 1 MC::Lund pi+:             " << n_fail_lund_pip << "\n";
    std::cout << "Fail no requested-detector positive:     " << n_fail_no_detector_positive << "\n";
    std::cout << "Fail selected outside theta/p bin:       " << n_fail_selected_outside_bin << "\n";
    std::cout << "Fail selected track not pi+:             " << n_fail_selected_not_pip << "\n";
    std::cout << "Fail bad/DC phi/kinematics:              " << n_fail_bad_kinematics << "\n";
    std::cout << "Fail phi-bin assignment:                 " << n_fail_phi_bin << "\n";

    gStyle->SetOptStat(0);

    if (mode == "2D") {
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
        lat.DrawLatex(0.15, 0.88, Form("Step 2.1 pass = %lld", n_step21_pass));
        lat.DrawLatex(0.15, 0.83, Form("#phi definition = %s", phi_def.c_str()));

        std::string outname = force_extension(output_file, ".png");

        c->SaveAs(outname.c_str());

        std::cout << "\nWrote:\n";
        std::cout << "  " << outname << "\n";

        delete c;
        delete h_phi_dp;
    }
    else if (mode == "1D") {
        gStyle->SetOptStat(1110);

        std::string outname = force_extension(output_file, ".pdf");

        TCanvas* c = new TCanvas("c_deltaP_per_phi_bin",
                                 "Delta p per phi bin",
                                 900, 700);

        c->SetTopMargin(0.10);
        c->SetRightMargin(0.05);
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.12);

        c->Print((outname + "[").c_str());

        for (int ibin = 0; ibin < n_phi_1d_bins; ibin++) {
            const double lo = phi1D_min + ibin * phi1D_width;
            const double hi = lo + phi1D_width;

            c->Clear();
            c->SetTopMargin(0.10);
            c->SetRightMargin(0.05);
            c->SetLeftMargin(0.12);
            c->SetBottomMargin(0.12);

            TH1D* h = h_dp_phi_bins[ibin];

            h->SetLineWidth(2);
            h->SetTitle(Form("Step 2.1: %.1f #leq #phi < %.1f deg;#Delta p = p_{rec} - p_{Lund} [GeV];Counts",
                             lo, hi));

            h->Draw("hist");

            TLine* zero = new TLine(0.0, 0.0, 0.0, h->GetMaximum() * 1.05);
            zero->SetLineStyle(2);
            zero->SetLineWidth(2);
            zero->Draw("same");

            TLatex lat;
            lat.SetNDC();
            lat.SetTextSize(0.035);
            lat.DrawLatex(0.15, 0.92,
                          Form("#phi bin %d / %d: %.1f #leq #phi < %.1f deg",
                               ibin + 1, n_phi_1d_bins, lo, hi));
            lat.DrawLatex(0.15, 0.86,
                          Form("Entries = %.0f", h->GetEntries()));
            lat.DrawLatex(0.15, 0.80,
                          Form("#phi definition = %s", phi_def.c_str()));
            lat.DrawLatex(0.15, 0.74,
                          Form("Step 2.1 pass total = %lld", n_step21_pass));

            c->Print(outname.c_str());

            delete zero;
        }

        c->Print((outname + "]").c_str());

        std::cout << "\nWrote multipage PDF:\n";
        std::cout << "  " << outname << "\n";
        std::cout << "Pages: " << n_phi_1d_bins << "\n";

        delete c;

        for (auto* h : h_dp_phi_bins) {
            delete h;
        }
    }
}