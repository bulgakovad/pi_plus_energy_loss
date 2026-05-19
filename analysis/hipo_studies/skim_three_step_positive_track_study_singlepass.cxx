// skim_three_step_positive_track_study_singlepass.cxx
//
// Single-pass version.
//
// Reads each input event ONCE and can write the same event to several output files.
//
// Common event requirement:
//   - exactly 1 trigger electron: pid = 11 and status < 0
//   - no other negative reconstructed tracks
//   - exactly 1 generated pi+ from MC::Lund
//   - NO vertex cuts
//
// Track bin:
//   - FD/CD selected by abs(REC::Particle.status)
//   - theta_rec in [theta_min, theta_max)
//   - p_rec in [p_min, p_max)
//
// Step 1:
//   - event must have at least one positive reconstructed track
//   - ALL positive reconstructed tracks must be labeled as requested detector
//       FD: abs(status) < 4000
//       CD: abs(status) >= 4000
//   - ALL positive reconstructed tracks must be inside the theta/p bin
//   - if several positive tracks exist, choose the one closest to MC::Lund pi+ momentum
//   - delta_p is computed using this selected track
//
// Step 2, ADVISOR-EXACT:
//   - from Step 1, keep events where the Step-1 selected track is identified as pi+
//   - multiple positive tracks are still allowed
//   - delta_p is computed using the same Step-1 selected track
//
// Step 2b:
//   - exactly 1 positive reconstructed track total
//   - PID can be anything
//   - that one positive track must pass detector/theta/p bin
//
// Step 3, ADVISOR-EXACT:
//   - from Step 2, keep events with no other positive tracks
//   - therefore exactly 1 positive track total, and it is the Step-1/Step-2 selected pi+
//
// Important:
//   This script uses MC::Lund ONLY for generated truth.
//   MC::Particle is intentionally not used and not used as fallback.
//
// Test:
//
// clas12root -l -b -q 'skim_three_step_positive_track_study_singlepass.cxx+("good_hipo.dat","FD_step1_test.hipo","FD_step2_test.hipo","FD_step2b_test.hipo","FD_step3_test.hipo","FD_step1_test.png","FD_step2_test.png","FD_step2b_test.png","FD_step3_test.png","FD",38,39,1.0,1.2,-0.5,0.5,200,10,"FD_selection_stats_test.txt")'
//
// Full:
//
// clas12root -l -b -q 'skim_three_step_positive_track_study_singlepass.cxx+("good_hipo.dat","FD_step1_allPositiveTracksInBin.hipo","FD_step2_selectedTrackIsPip.hipo","FD_step2b_onlyOnePositive.hipo","FD_step3_selectedPipOnlyPositive.hipo","FD_step1_dp.png","FD_step2_dp.png","FD_step2b_dp.png","FD_step3_dp.png","FD",38,39,1.0,1.2,-0.5,0.5,200,-1,"FD_selection_stats.txt")'

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>

#include "TMath.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLatex.h"

#include "hipo4/reader.h"
#include "hipo4/writer.h"
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
};

static bool ends_with(const std::string& s, const std::string& suffix)
{
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static std::vector<std::string> get_hipo_files_from_dat(const std::string& dat_file)
{
    std::vector<std::string> hipo_files;

    std::ifstream fin(dat_file);
    if (!fin.is_open()) {
        std::cerr << "ERROR: Cannot open .dat file: " << dat_file << std::endl;
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

static bool file_exists(const std::string& path)
{
    std::ifstream f(path);
    return f.good();
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

static bool track_in_bin(const TrackInfo& trk,
                         const std::string& det,
                         double theta_min_deg,
                         double theta_max_deg,
                         double p_min,
                         double p_max)
{
    if (!detector_status_pass(trk.status, det)) {
        return false;
    }

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

static void save_hist_png(TH1D* h,
                          const char* output_png,
                          const char* title,
                          long long n_written)
{
    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas(Form("c_%s", h->GetName()), title, 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    h->SetLineWidth(2);
    h->SetTitle(title);
    h->Draw("hist");

    TLine* zero = new TLine(0.0, 0.0, 0.0, h->GetMaximum() * 1.05);
    zero->SetLineStyle(2);
    zero->SetLineWidth(2);
    zero->Draw("same");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.15, 0.92, Form("Written events: %lld", n_written));

    c->SaveAs(output_png);

    delete c;
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

static bool reached_cap(long long n_written, int max_events_to_write)
{
    return (max_events_to_write > 0 && n_written >= max_events_to_write);
}

static double pct(long long numerator, long long denominator)
{
    if (denominator <= 0) return 0.0;
    return 100.0 * static_cast<double>(numerator) / static_cast<double>(denominator);
}

static void print_stat_line(std::ostream& os,
                            const std::string& label,
                            long long count,
                            long long previous,
                            long long total)
{
    os << std::left << std::setw(62) << label
       << std::right << std::setw(14) << count
       << "   "
       << std::fixed << std::setprecision(2)
       << std::setw(8) << pct(count, previous) << "% of previous"
       << "   "
       << std::fixed << std::setprecision(2)
       << std::setw(8) << pct(count, total) << "% of scanned"
       << "\n";
}

void skim_three_step_positive_track_study_singlepass(const char* input_dat_file,
                                                     const char* output_step1_hipo,
                                                     const char* output_step2_hipo,
                                                     const char* output_step2b_hipo,
                                                     const char* output_step3_hipo,
                                                     const char* output_step1_png,
                                                     const char* output_step2_png,
                                                     const char* output_step2b_png,
                                                     const char* output_step3_png,
                                                     const char* detector = "FD",
                                                     double theta_min_deg = 38.0,
                                                     double theta_max_deg = 39.0,
                                                     double p_min = 1.0,
                                                     double p_max = 1.2,
                                                     double dp_min = -0.50,
                                                     double dp_max = 0.50,
                                                     int n_dp_bins = 200,
                                                     int max_events_to_write = 1000,
                                                     const char* output_stats_txt = "selection_stats.txt")
{
    std::cout << "Macro started: SINGLE-PASS VERSION\n";
    std::cout << "Input .dat file:    " << input_dat_file << "\n";
    std::cout << "Detector:           " << detector << "\n";
    std::cout << "Theta bin:          [" << theta_min_deg << ", " << theta_max_deg << ") deg\n";
    std::cout << "Momentum bin:       [" << p_min << ", " << p_max << ") GeV\n";
    std::cout << "Delta-p plot range: [" << dp_min << ", " << dp_max << "] GeV\n";
    std::cout << "Max events/output:  " << max_events_to_write << "\n";
    std::cout << "Generated truth:    MC::Lund ONLY\n";
    std::cout << "Stats output:       " << output_stats_txt << "\n";

    std::string det = detector;

    if (det != "FD" && det != "CD") {
        std::cerr << "ERROR: detector must be \"FD\" or \"CD\", got: "
                  << det << "\n";
        return;
    }

    auto hipo_files_raw = get_hipo_files_from_dat(input_dat_file);

    if (hipo_files_raw.empty()) {
        std::cerr << "ERROR: No .hipo files found in .dat file: "
                  << input_dat_file << "\n";
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
        std::cerr << "ERROR: No readable .hipo files left.\n";
        return;
    }

    std::cout << "\nFound " << hipo_files.size() << " readable .hipo files\n";

    for (size_t i = 0; i < hipo_files.size(); i++) {
        if (i < 10) {
            std::cout << "  [" << i << "] " << hipo_files[i] << "\n";
        }
    }

    if (hipo_files.size() > 10) {
        std::cout << "  ...\n";
    }

    TH1D* h_step1 = new TH1D(
        "h_step1",
        "Step 1: all positive tracks in detector/theta/p bin;#Delta p = p_{rec}(selected positive) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_dp_bins, dp_min, dp_max
    );

    TH1D* h_step2 = new TH1D(
        "h_step2",
        "Step 2: Step-1 selected track is #pi^{+};#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_dp_bins, dp_min, dp_max
    );

    TH1D* h_step2b = new TH1D(
        "h_step2b",
        "Step 2b: exactly one positive track;#Delta p = p_{rec}(only positive) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_dp_bins, dp_min, dp_max
    );

    TH1D* h_step3 = new TH1D(
        "h_step3",
        "Step 3: Step-2 #pi^{+} is the only positive track;#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_dp_bins, dp_min, dp_max
    );

    hipo::writer writer1;
    hipo::writer writer2;
    hipo::writer writer2b;
    hipo::writer writer3;

    bool writers_opened = false;

    long long scanned_events = 0;

    // Common cut-flow.
    long long n_common_trigger_pass = 0;
    long long n_common_no_other_negative_pass = 0;
    long long n_common_lund_pip_pass = 0;

    // Step-1 internal cut-flow after common cuts.
    long long n_step1_has_positive = 0;
    long long n_step1_all_positive_detector = 0;
    long long n_step1_all_positive_theta_p = 0;

    // Step pass counts before output cap.
    long long n_step1_pass = 0;
    long long n_step2_pass = 0;
    long long n_step2b_pass = 0;
    long long n_step3_pass = 0;

    // Written counts after output cap.
    long long n_step1_written = 0;
    long long n_step2_written = 0;
    long long n_step2b_written = 0;
    long long n_step3_written = 0;

    // Fail counters.
    long long n_fail_trigger = 0;
    long long n_fail_other_negative = 0;
    long long n_fail_lund_pip = 0;
    long long n_fail_missing_mclund = 0;

    long long n_step1_fail_no_positive = 0;
    long long n_step1_fail_positive_not_detector = 0;
    long long n_step1_fail_positive_not_in_theta_p_bin = 0;

    long long n_step2_fail_selected_track_not_pip = 0;
    long long n_step3_fail_extra_positive_tracks = 0;

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
            continue;
        }

        hipo::dictionary factory;
        reader.readDictionary(factory);

        if (!writers_opened) {
            writer1.getDictionary()  = factory;
            writer2.getDictionary()  = factory;
            writer2b.getDictionary() = factory;
            writer3.getDictionary()  = factory;

            writer1.open(output_step1_hipo);
            writer2.open(output_step2_hipo);
            writer2b.open(output_step2b_hipo);
            writer3.open(output_step3_hipo);

            writers_opened = true;

            std::cout << "Opened output files:\n"
                      << "  " << output_step1_hipo << "\n"
                      << "  " << output_step2_hipo << "\n"
                      << "  " << output_step2b_hipo << "\n"
                      << "  " << output_step3_hipo << "\n";
        }

        hipo::schema rec_schema;
        hipo::schema mclund_schema;

        try {
            rec_schema = factory.getSchema("REC::Particle");
        }
        catch (...) {
            std::cerr << "WARNING: REC::Particle missing, skipping file:\n"
                      << "  " << input_hipo << "\n";
            continue;
        }

        bool have_mclund_schema = false;

        try {
            mclund_schema = factory.getSchema("MC::Lund");
            have_mclund_schema = true;
        }
        catch (...) {
            have_mclund_schema = false;
        }

        if (!have_mclund_schema) {
            n_fail_missing_mclund++;
            std::cerr << "WARNING: MC::Lund missing, skipping file:\n"
                      << "  " << input_hipo << "\n";
            continue;
        }

        hipo::event event;
        hipo::bank recPart(rec_schema);
        hipo::bank mcLund(mclund_schema);

        while (reader.next()) {
            reader.read(event);
            scanned_events++;

            if (scanned_events % 50000 == 0) {
                std::cout << "Scanned " << scanned_events
                          << " events | written:"
                          << " step1=" << n_step1_written
                          << " step2=" << n_step2_written
                          << " step2b=" << n_step2b_written
                          << " step3=" << n_step3_written
                          << "\n";
            }

            event.getStructure(recPart);
            event.getStructure(mcLund);

            int n_trigger_e = 0;
            int n_other_negative = 0;
            int n_positive_tracks_total = 0;

            std::vector<TrackInfo> positive_tracks;

            for (int i = 0; i < recPart.getRows(); i++) {
                const int pid = recPart.getInt("pid", i);
                const int charge = recPart.getInt("charge", i);
                const int status = recPart.getInt("status", i);

                const double px = recPart.getFloat("px", i);
                const double py = recPart.getFloat("py", i);
                const double pz = recPart.getFloat("pz", i);

                const bool is_trigger_e = (pid == 11 && status < 0);

                if (is_trigger_e) {
                    n_trigger_e++;
                }
                else if (charge < 0) {
                    n_other_negative++;
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

                if (charge > 0) {
                    n_positive_tracks_total++;
                    positive_tracks.push_back(trk);
                }
            }

            if (n_trigger_e != 1) {
                n_fail_trigger++;
                continue;
            }

            n_common_trigger_pass++;

            if (n_other_negative != 0) {
                n_fail_other_negative++;
                continue;
            }

            n_common_no_other_negative_pass++;

            double p_lund = std::numeric_limits<double>::quiet_NaN();
            int n_lund_pip = 0;

            if (!get_mclund_pip_momentum(mcLund, p_lund, n_lund_pip)) {
                n_fail_lund_pip++;
                continue;
            }

            n_common_lund_pip_pass++;

            bool pass_step1 = false;
            bool pass_step2 = false;
            bool pass_step2b = false;
            bool pass_step3 = false;

            TrackInfo trk_step1;
            TrackInfo trk_step2;
            TrackInfo trk_step2b;
            TrackInfo trk_step3;

            // --------------------------------------------------------
            // Step 1 cut-flow:
            //   after common cuts, require:
            //   1) at least one positive track
            //   2) all positive tracks are requested detector
            //   3) all positive tracks are inside theta/p bin
            //   4) choose closest positive track to p_Lund(pi+)
            // --------------------------------------------------------
            bool step1_all_positive_tracks_in_requested_detector = true;
            bool step1_all_positive_tracks_in_theta_p_bin = true;

            if (positive_tracks.empty()) {
                step1_all_positive_tracks_in_requested_detector = false;
                step1_all_positive_tracks_in_theta_p_bin = false;
                n_step1_fail_no_positive++;
            }
            else {
                n_step1_has_positive++;
            }

            for (const auto& pos : positive_tracks) {
                if (!detector_status_pass(pos.status, det)) {
                    step1_all_positive_tracks_in_requested_detector = false;
                }

                if (!track_in_bin(pos, det, theta_min_deg, theta_max_deg, p_min, p_max)) {
                    step1_all_positive_tracks_in_theta_p_bin = false;
                }
            }

            if (!positive_tracks.empty()) {
                if (!step1_all_positive_tracks_in_requested_detector) {
                    n_step1_fail_positive_not_detector++;
                }
                else {
                    n_step1_all_positive_detector++;

                    if (!step1_all_positive_tracks_in_theta_p_bin) {
                        n_step1_fail_positive_not_in_theta_p_bin++;
                    }
                    else {
                        n_step1_all_positive_theta_p++;

                        double best_abs_dp = std::numeric_limits<double>::infinity();

                        for (const auto& pos : positive_tracks) {
                            const double abs_dp = std::abs(pos.p - p_lund);

                            if (abs_dp < best_abs_dp) {
                                best_abs_dp = abs_dp;
                                trk_step1 = pos;
                                pass_step1 = true;
                            }
                        }
                    }
                }
            }

            // --------------------------------------------------------
            // Step 2, ADVISOR-EXACT:
            //
            // From Step 1, keep events where the Step-1 selected track
            // is identified as a pion.
            // --------------------------------------------------------
            if (pass_step1) {
                if (trk_step1.pid == 211) {
                    pass_step2 = true;
                    trk_step2 = trk_step1;
                }
                else {
                    n_step2_fail_selected_track_not_pip++;
                }
            }

            // --------------------------------------------------------
            // Step 2b:
            //
            // From Step 1's loose context, require exactly one positive
            // track. That one positive track must be in detector/theta/p bin.
            // --------------------------------------------------------
            if (n_positive_tracks_total == 1 && positive_tracks.size() == 1) {
                const auto& only_positive = positive_tracks[0];

                if (track_in_bin(only_positive, det, theta_min_deg, theta_max_deg, p_min, p_max)) {
                    pass_step2b = true;
                    trk_step2b = only_positive;
                }
            }

            // --------------------------------------------------------
            // Step 3, ADVISOR-EXACT:
            //
            // From Step 2, require no other positive tracks.
            // --------------------------------------------------------
            if (pass_step2) {
                if (n_positive_tracks_total == 1) {
                    pass_step3 = true;
                    trk_step3 = trk_step2;
                }
                else {
                    n_step3_fail_extra_positive_tracks++;
                }
            }

            if (pass_step1) n_step1_pass++;
            if (pass_step2) n_step2_pass++;
            if (pass_step2b) n_step2b_pass++;
            if (pass_step3) n_step3_pass++;

            if (pass_step1 && !reached_cap(n_step1_written, max_events_to_write)) {
                const double dp = trk_step1.p - p_lund;

                h_step1->Fill(dp);
                writer1.addEvent(event);
                n_step1_written++;

                if (n_step1_written <= 20) {
                    std::cout << "STEP1 KEEP"
                              << " rec_index=" << trk_step1.index
                              << " pid=" << trk_step1.pid
                              << " charge=" << trk_step1.charge
                              << " status=" << trk_step1.status
                              << " theta=" << trk_step1.theta_deg
                              << " p_rec=" << trk_step1.p
                              << " p_lund=" << p_lund
                              << " dp=" << dp
                              << " n_pos=" << n_positive_tracks_total
                              << "\n";
                }
            }

            if (pass_step2 && !reached_cap(n_step2_written, max_events_to_write)) {
                const double dp = trk_step2.p - p_lund;

                h_step2->Fill(dp);
                writer2.addEvent(event);
                n_step2_written++;

                if (n_step2_written <= 20) {
                    std::cout << "STEP2 KEEP"
                              << " rec_index=" << trk_step2.index
                              << " pid=" << trk_step2.pid
                              << " charge=" << trk_step2.charge
                              << " status=" << trk_step2.status
                              << " theta=" << trk_step2.theta_deg
                              << " p_rec=" << trk_step2.p
                              << " p_lund=" << p_lund
                              << " dp=" << dp
                              << " n_pos=" << n_positive_tracks_total
                              << "\n";
                }
            }

            if (pass_step2b && !reached_cap(n_step2b_written, max_events_to_write)) {
                const double dp = trk_step2b.p - p_lund;

                h_step2b->Fill(dp);
                writer2b.addEvent(event);
                n_step2b_written++;

                if (n_step2b_written <= 20) {
                    std::cout << "STEP2B KEEP"
                              << " rec_index=" << trk_step2b.index
                              << " pid=" << trk_step2b.pid
                              << " charge=" << trk_step2b.charge
                              << " status=" << trk_step2b.status
                              << " theta=" << trk_step2b.theta_deg
                              << " p_rec=" << trk_step2b.p
                              << " p_lund=" << p_lund
                              << " dp=" << dp
                              << " n_pos=" << n_positive_tracks_total
                              << "\n";
                }
            }

            if (pass_step3 && !reached_cap(n_step3_written, max_events_to_write)) {
                const double dp = trk_step3.p - p_lund;

                h_step3->Fill(dp);
                writer3.addEvent(event);
                n_step3_written++;

                if (n_step3_written <= 20) {
                    std::cout << "STEP3 KEEP"
                              << " rec_index=" << trk_step3.index
                              << " pid=" << trk_step3.pid
                              << " charge=" << trk_step3.charge
                              << " status=" << trk_step3.status
                              << " theta=" << trk_step3.theta_deg
                              << " p_rec=" << trk_step3.p
                              << " p_lund=" << p_lund
                              << " dp=" << dp
                              << " n_pos=" << n_positive_tracks_total
                              << "\n";
                }
            }

            if (max_events_to_write > 0 &&
                n_step1_written  >= max_events_to_write &&
                n_step2_written  >= max_events_to_write &&
                n_step2b_written >= max_events_to_write &&
                n_step3_written  >= max_events_to_write) {

                std::cout << "\nReached max_events_to_write for all four outputs.\n";
                break;
            }
        }

        if (max_events_to_write > 0 &&
            n_step1_written  >= max_events_to_write &&
            n_step2_written  >= max_events_to_write &&
            n_step2b_written >= max_events_to_write &&
            n_step3_written  >= max_events_to_write) {

            break;
        }
    }

    if (writers_opened) {
        writer1.close();
        writer2.close();
        writer2b.close();
        writer3.close();
    }
    else {
        std::cerr << "ERROR: writers were never opened. No output HIPO files written.\n";
    }

    save_hist_png(
        h_step1,
        output_step1_png,
        "Step 1: all positive tracks in detector/theta/p bin;#Delta p = p_{rec}(selected positive) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_step1_written
    );

    save_hist_png(
        h_step2,
        output_step2_png,
        "Step 2: Step-1 selected track is #pi^{+};#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_step2_written
    );

    save_hist_png(
        h_step2b,
        output_step2b_png,
        "Step 2b: exactly one positive track;#Delta p = p_{rec}(only positive) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_step2b_written
    );

    save_hist_png(
        h_step3,
        output_step3_png,
        "Step 3: Step-2 #pi^{+} is the only positive track;#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_step3_written
    );

    // ------------------------------------------------------------
    // Write stats file.
    // ------------------------------------------------------------
    std::ofstream stats(output_stats_txt);

    if (!stats.is_open()) {
        std::cerr << "ERROR: Could not open stats output file: "
                  << output_stats_txt << "\n";
    }
    else {
        stats << "Selection statistics\n";
        stats << "====================\n\n";

        stats << "Input .dat file:       " << input_dat_file << "\n";
        stats << "Detector:              " << det << "\n";
        stats << "Theta bin:             [" << theta_min_deg << ", " << theta_max_deg << ") deg\n";
        stats << "Momentum bin:          [" << p_min << ", " << p_max << ") GeV\n";
        stats << "Delta-p plot range:    [" << dp_min << ", " << dp_max << "] GeV\n";
        stats << "Generated truth:       MC::Lund only\n";
        stats << "Max events/output:     " << max_events_to_write << "\n\n";

        stats << "Input files\n";
        stats << "-----------\n";
        stats << "Readable HIPO files:   " << hipo_files.size() << "\n";
        stats << "Files skipped missing MC::Lund: " << n_fail_missing_mclund << "\n\n";

        stats << "Common cut flow\n";
        stats << "---------------\n";
        print_stat_line(stats, "Scanned events:", scanned_events, scanned_events, scanned_events);
        print_stat_line(stats, "Events with exactly 1 trigger electron:", n_common_trigger_pass, scanned_events, scanned_events);
        print_stat_line(stats, "Events with no other negative tracks:", n_common_no_other_negative_pass, n_common_trigger_pass, scanned_events);
        print_stat_line(stats, "Events with exactly 1 MC::Lund pi+:", n_common_lund_pip_pass, n_common_no_other_negative_pass, scanned_events);

        stats << "\nStep 1 cut flow after common cuts\n";
        stats << "---------------------------------\n";
        print_stat_line(stats, "Events with at least 1 positive track:", n_step1_has_positive, n_common_lund_pip_pass, scanned_events);
        print_stat_line(stats, ("Events where all positive tracks are " + det + ":"), n_step1_all_positive_detector, n_step1_has_positive, scanned_events);
        print_stat_line(stats, "Events where all positive tracks are in theta/p bin:", n_step1_all_positive_theta_p, n_step1_all_positive_detector, scanned_events);
        print_stat_line(stats, "Step 1 pass:", n_step1_pass, n_step1_all_positive_theta_p, scanned_events);

        stats << "\nAdvisor-requested steps\n";
        stats << "-----------------------\n";
        print_stat_line(stats, "Step 1: choose closest positive track to MC::Lund pi+:", n_step1_pass, n_common_lund_pip_pass, scanned_events);
        print_stat_line(stats, "Step 2: Step-1 selected track is pi+:", n_step2_pass, n_step1_pass, scanned_events);
        print_stat_line(stats, "Step 2b: exactly one positive track:", n_step2b_pass, n_step1_pass, scanned_events);
        print_stat_line(stats, "Step 3: Step-2 pi+ is only positive track:", n_step3_pass, n_step2_pass, scanned_events);

        stats << "\nWritten counts after output cap\n";
        stats << "-------------------------------\n";
        stats << "Step 1 written:   " << n_step1_written  << "   " << output_step1_hipo  << "   " << output_step1_png  << "\n";
        stats << "Step 2 written:   " << n_step2_written  << "   " << output_step2_hipo  << "   " << output_step2_png  << "\n";
        stats << "Step 2b written:  " << n_step2b_written << "   " << output_step2b_hipo << "   " << output_step2b_png << "\n";
        stats << "Step 3 written:   " << n_step3_written  << "   " << output_step3_hipo  << "   " << output_step3_png  << "\n";

        stats << "\nFailure counters\n";
        stats << "----------------\n";
        stats << "Fail exactly 1 trigger electron:                  " << n_fail_trigger << "\n";
        stats << "Fail no-other-negative-tracks requirement:         " << n_fail_other_negative << "\n";
        stats << "Fail exactly 1 MC::Lund pi+ check:                 " << n_fail_lund_pip << "\n";
        stats << "Step 1 fail: no positive tracks:                   " << n_step1_fail_no_positive << "\n";
        stats << "Step 1 fail: at least one positive not " << det << ":              " << n_step1_fail_positive_not_detector << "\n";
        stats << "Step 1 fail: all positives " << det << ", but one outside theta/p: " << n_step1_fail_positive_not_in_theta_p_bin << "\n";
        stats << "Step 2 fail: selected Step-1 track is not pi+:      " << n_step2_fail_selected_track_not_pip << "\n";
        stats << "Step 3 fail: Step-2 event has extra positives:      " << n_step3_fail_extra_positive_tracks << "\n";

        stats << "\nHistogram statistics\n";
        stats << "--------------------\n";
        stats << "Step 1 entries:   " << h_step1->GetEntries()
              << "   mean=" << h_step1->GetMean()
              << "   stddev=" << h_step1->GetStdDev() << "\n";
        stats << "Step 2 entries:   " << h_step2->GetEntries()
              << "   mean=" << h_step2->GetMean()
              << "   stddev=" << h_step2->GetStdDev() << "\n";
        stats << "Step 2b entries:  " << h_step2b->GetEntries()
              << "   mean=" << h_step2b->GetMean()
              << "   stddev=" << h_step2b->GetStdDev() << "\n";
        stats << "Step 3 entries:   " << h_step3->GetEntries()
              << "   mean=" << h_step3->GetMean()
              << "   stddev=" << h_step3->GetStdDev() << "\n";

        stats.close();
    }

    std::cout << "\n============================================================\n";
    std::cout << "DONE SINGLE-PASS STUDY\n";
    std::cout << "Scanned events:                       " << scanned_events << "\n";
    std::cout << "Common: exactly 1 trigger electron:   " << n_common_trigger_pass << "\n";
    std::cout << "Common: no other negative tracks:     " << n_common_no_other_negative_pass << "\n";
    std::cout << "Common: exactly 1 MC::Lund pi+:       " << n_common_lund_pip_pass << "\n";

    std::cout << "\nStep 1 cut flow:\n";
    std::cout << "  At least 1 positive track:          " << n_step1_has_positive << "\n";
    std::cout << "  All positive tracks are " << det << ":          " << n_step1_all_positive_detector << "\n";
    std::cout << "  All positives in theta/p bin:       " << n_step1_all_positive_theta_p << "\n";

    std::cout << "\nPass counts before output cap:\n";
    std::cout << "  Step 1 pass:                        " << n_step1_pass << "\n";
    std::cout << "  Step 2 pass:                        " << n_step2_pass << "\n";
    std::cout << "  Step 2b pass:                       " << n_step2b_pass << "\n";
    std::cout << "  Step 3 pass:                        " << n_step3_pass << "\n";

    std::cout << "\nStep 1 diagnostic counts:\n";
    std::cout << "  Step 1 fail: no positive tracks:                         "
              << n_step1_fail_no_positive << "\n";
    std::cout << "  Step 1 fail: at least one positive not labeled as "
              << det << ":              "
              << n_step1_fail_positive_not_detector << "\n";
    std::cout << "  Step 1 fail: positives all labeled "
              << det << ", but at least one outside theta/p bin: "
              << n_step1_fail_positive_not_in_theta_p_bin << "\n";

    std::cout << "\nStep 2 diagnostic counts:\n";
    std::cout << "  Step 2 fail: Step-1 selected track is not pi+:             "
              << n_step2_fail_selected_track_not_pip << "\n";

    std::cout << "\nStep 3 diagnostic counts:\n";
    std::cout << "  Step 3 fail: Step-2 event has extra positive tracks:       "
              << n_step3_fail_extra_positive_tracks << "\n";

    std::cout << "\nWritten counts:\n";
    std::cout << "  Step 1 written:                     " << n_step1_written
              << "  " << output_step1_hipo << "  " << output_step1_png << "\n";
    std::cout << "  Step 2 written:                     " << n_step2_written
              << "  " << output_step2_hipo << "  " << output_step2_png << "\n";
    std::cout << "  Step 2b written:                    " << n_step2b_written
              << "  " << output_step2b_hipo << "  " << output_step2b_png << "\n";
    std::cout << "  Step 3 written:                     " << n_step3_written
              << "  " << output_step3_hipo << "  " << output_step3_png << "\n";

    std::cout << "\nFailures:\n";
    std::cout << "  fail trigger electron count:        " << n_fail_trigger << "\n";
    std::cout << "  fail other negative tracks:         " << n_fail_other_negative << "\n";
    std::cout << "  fail MC::Lund pi+ check:            " << n_fail_lund_pip << "\n";
    std::cout << "  files skipped missing MC::Lund:     " << n_fail_missing_mclund << "\n";

    std::cout << "\nStats file written: " << output_stats_txt << "\n";

    delete h_step1;
    delete h_step2;
    delete h_step2b;
    delete h_step3;
}