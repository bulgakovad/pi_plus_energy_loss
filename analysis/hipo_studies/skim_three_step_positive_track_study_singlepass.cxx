// skim_three_step_positive_track_study_singlepass.cxx
//
// Single-pass version for updated advisor request.
//
// Base selection:
//   - exactly 1 reconstructed electron total
//   - that one electron must be trigger electron: pid = 11 and status < 0
//   - no other REC electrons with any status
//   - NO negative hadrons allowed
//       i.e. no REC tracks with charge < 0 except the trigger electron
//   - exactly 1 generated pi+ from MC::Lund
//   - NO vertex cuts
//
// Step 1:
//   A. choose events with at least one positive track in requested detector, usually FD
//   B. among requested-detector positive tracks, choose the one closest to MC::Lund pi+ momentum
//   C. require this selected positive track to be inside theta/p bin
//   D. make delta p histogram and save HIPO
//
// Step 2.1 after Step 1:
//   - selected Step-1 track is identified as pi+
//   - other positive tracks are allowed anywhere
//
// Step 2.2 after Step 1:
//   - exactly one positive track total in REC::Particle
//   - since this is after Step 1, that one track is the selected detector track in bin
//
// Step 3.1 after Step 2.1:
//   - selected pion track has no other positive tracks
//
// Step 3.2 after Step 2.2:
//   - the only positive track from Step 2.2 is identified as pion
//
// With the no-negative-hadrons veto added back, Step 3.1/3.2 should be closer
// to the old Step 3 logic.
//
// Generated truth:
//   - MC::Lund ONLY
//   - MC::Particle is not used
//
// Example test:
// clas12root -l -b -q 'skim_three_step_positive_track_study_singlepass.cxx+("good_hipo.dat","FD_theta38_39_p1p0_1p2_test","FD",38,39,1.0,1.2,-0.5,0.1,200,1000)'
//
// Example full:
// clas12root -l -b -q 'skim_three_step_positive_track_study_singlepass.cxx+("good_hipo.dat","FD_theta38_39_p1p0_1p2","FD",38,39,1.0,1.2,-0.5,0.5,200,-1)'

#include <iostream>
#include <fstream>
#include <sstream>
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
        std::cerr << "ERROR: Cannot open .dat file: " << dat_file << "\n";
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
    os << std::left << std::setw(72) << label
       << std::right << std::setw(14) << count
       << "   "
       << std::fixed << std::setprecision(2)
       << std::setw(8) << pct(count, previous) << "% of previous"
       << "   "
       << std::fixed << std::setprecision(2)
       << std::setw(8) << pct(count, total) << "% of scanned"
       << "\n";
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

void skim_three_step_positive_track_study_singlepass(const char* input_dat_file,
                                                     const char* output_prefix = "FD_theta38_39_p1p0_1p2",
                                                     const char* detector = "FD",
                                                     double theta_min_deg = 38.0,
                                                     double theta_max_deg = 39.0,
                                                     double p_min = 1.0,
                                                     double p_max = 1.2,
                                                     double dp_min = -0.50,
                                                     double dp_max = 0.50,
                                                     int n_dp_bins = 200,
                                                     int max_events_to_write = 1000)
{
    std::string prefix = output_prefix;
    std::string det = detector;

    if (det != "FD" && det != "CD") {
        std::cerr << "ERROR: detector must be \"FD\" or \"CD\", got: "
                  << det << "\n";
        return;
    }

    const std::string out_step1_hipo = prefix + "_step1.hipo";
    const std::string out_21_hipo    = prefix + "_step2p1_selectedTrackIsPip.hipo";
    const std::string out_22_hipo    = prefix + "_step2p2_onlyOnePositive.hipo";
    const std::string out_31_hipo    = prefix + "_step3p1_selectedPipOnlyPositive.hipo";
    const std::string out_32_hipo    = prefix + "_step3p2_onlyPositiveIsPip.hipo";

    const std::string out_step1_png = prefix + "_step1_dp.png";
    const std::string out_21_png    = prefix + "_step2p1_dp.png";
    const std::string out_22_png    = prefix + "_step2p2_dp.png";
    const std::string out_31_png    = prefix + "_step3p1_dp.png";
    const std::string out_32_png    = prefix + "_step3p2_dp.png";

    const std::string out_stats_txt = prefix + "_selection_stats.txt";

    std::cout << "Macro started: SINGLE-PASS UPDATED 5-OUTPUT VERSION\n";
    std::cout << "Input .dat file:    " << input_dat_file << "\n";
    std::cout << "Output prefix:      " << prefix << "\n";
    std::cout << "Detector for Step1 selected track: " << det << "\n";
    std::cout << "Theta bin:          [" << theta_min_deg << ", " << theta_max_deg << ") deg\n";
    std::cout << "Momentum bin:       [" << p_min << ", " << p_max << ") GeV\n";
    std::cout << "Delta-p plot range: [" << dp_min << ", " << dp_max << "] GeV\n";
    std::cout << "Max events/output:  " << max_events_to_write << "\n";
    std::cout << "Base selection:     exactly 1 trigger electron and NO negative hadrons\n";
    std::cout << "Generated truth:    MC::Lund ONLY\n\n";

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

    std::cout << "Found " << hipo_files.size() << " readable .hipo files\n";

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
        "Step 1: selected FD positive in #theta/p bin;#Delta p = p_{rec}(selected positive) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_dp_bins, dp_min, dp_max
    );

    TH1D* h_21 = new TH1D(
        "h_21",
        "Step 2.1: Step-1 selected track is #pi^{+};#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_dp_bins, dp_min, dp_max
    );

    TH1D* h_22 = new TH1D(
        "h_22",
        "Step 2.2: exactly one positive track after Step 1;#Delta p = p_{rec}(only positive) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_dp_bins, dp_min, dp_max
    );

    TH1D* h_31 = new TH1D(
        "h_31",
        "Step 3.1: Step-2.1 #pi^{+} is the only positive track;#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_dp_bins, dp_min, dp_max
    );

    TH1D* h_32 = new TH1D(
        "h_32",
        "Step 3.2: Step-2.2 only positive track is #pi^{+};#Delta p = p_{rec}(only positive #pi^{+}) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_dp_bins, dp_min, dp_max
    );

    hipo::writer writer1;
    hipo::writer writer21;
    hipo::writer writer22;
    hipo::writer writer31;
    hipo::writer writer32;

    bool writers_opened = false;

    long long scanned_events = 0;

    // Base cut-flow.
    long long n_base_one_electron_total = 0;
    long long n_base_one_trigger_electron = 0;
    long long n_base_no_negative_hadrons = 0;
    long long n_base_mclund_pip = 0;

    // Step 1 cut-flow.
    long long n_step1_has_fd_positive = 0;
    long long n_step1_selected_in_theta_p = 0;

    // Step pass counts before output cap.
    long long n_step1_pass = 0;
    long long n_21_pass = 0;
    long long n_22_pass = 0;
    long long n_31_pass = 0;
    long long n_32_pass = 0;

    // Written counts after output cap.
    long long n_step1_written = 0;
    long long n_21_written = 0;
    long long n_22_written = 0;
    long long n_31_written = 0;
    long long n_32_written = 0;

    // Fail counters.
    long long n_fail_base_electron = 0;
    long long n_fail_negative_hadrons = 0;
    long long n_fail_lund_pip = 0;
    long long n_fail_missing_mclund = 0;
    long long n_step1_fail_no_fd_positive = 0;
    long long n_step1_fail_selected_not_in_bin = 0;
    long long n_21_fail_selected_not_pip = 0;
    long long n_22_fail_not_exactly_one_positive = 0;
    long long n_31_fail_extra_positive = 0;
    long long n_32_fail_only_positive_not_pip = 0;

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
            writer21.getDictionary() = factory;
            writer22.getDictionary() = factory;
            writer31.getDictionary() = factory;
            writer32.getDictionary() = factory;

            writer1.open(out_step1_hipo.c_str());
            writer21.open(out_21_hipo.c_str());
            writer22.open(out_22_hipo.c_str());
            writer31.open(out_31_hipo.c_str());
            writer32.open(out_32_hipo.c_str());

            writers_opened = true;

            std::cout << "Opened output files:\n"
                      << "  " << out_step1_hipo << "\n"
                      << "  " << out_21_hipo << "\n"
                      << "  " << out_22_hipo << "\n"
                      << "  " << out_31_hipo << "\n"
                      << "  " << out_32_hipo << "\n";
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
                          << " | written:"
                          << " step1=" << n_step1_written
                          << " 2.1=" << n_21_written
                          << " 2.2=" << n_22_written
                          << " 3.1=" << n_31_written
                          << " 3.2=" << n_32_written
                          << "\n";
            }

            event.getStructure(recPart);
            event.getStructure(mcLund);

            int n_electron_total = 0;
            int n_trigger_electron = 0;
            int n_positive_tracks_total = 0;
            int n_negative_hadrons_total = 0;

            std::vector<TrackInfo> positive_tracks;
            std::vector<TrackInfo> fd_positive_tracks;

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

                // Added back to match the old "only trigger electron is negative" logic:
                // reject negative non-electron tracks, e.g. pi-, K-, anti-proton.
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

                if (charge > 0) {
                    n_positive_tracks_total++;
                    positive_tracks.push_back(trk);

                    if (detector_status_pass(status, det)) {
                        fd_positive_tracks.push_back(trk);
                    }
                }
            }

            // --------------------------------------------------------
            // Base selection:
            // exactly one REC electron total, and it is trigger electron.
            // --------------------------------------------------------
            if (!(n_electron_total == 1 && n_trigger_electron == 1)) {
                n_fail_base_electron++;
                continue;
            }

            n_base_one_electron_total++;
            n_base_one_trigger_electron++;

            // --------------------------------------------------------
            // Added base selection:
            // no negative hadrons allowed.
            // The trigger electron is the only negative reconstructed track.
            // --------------------------------------------------------
            if (n_negative_hadrons_total != 0) {
                n_fail_negative_hadrons++;
                continue;
            }

            n_base_no_negative_hadrons++;

            double p_lund = std::numeric_limits<double>::quiet_NaN();
            int n_lund_pip = 0;

            if (!get_mclund_pip_momentum(mcLund, p_lund, n_lund_pip)) {
                n_fail_lund_pip++;
                continue;
            }

            n_base_mclund_pip++;

            bool pass_step1 = false;
            bool pass_21 = false;
            bool pass_22 = false;
            bool pass_31 = false;
            bool pass_32 = false;

            TrackInfo selected_step1;

            // --------------------------------------------------------
            // Step 1:
            // choose requested-detector positive track closest to p_Lund(pi+),
            // then require selected track inside theta/p bin.
            // --------------------------------------------------------
            if (fd_positive_tracks.empty()) {
                n_step1_fail_no_fd_positive++;
            }
            else {
                n_step1_has_fd_positive++;

                double best_abs_dp = std::numeric_limits<double>::infinity();

                for (const auto& trk : fd_positive_tracks) {
                    const double abs_dp = std::abs(trk.p - p_lund);

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
                    n_step1_fail_selected_not_in_bin++;
                }
                else {
                    pass_step1 = true;
                    n_step1_selected_in_theta_p++;
                }
            }

            // --------------------------------------------------------
            // Step 2.1:
            // from Step 1, selected track is pion.
            // --------------------------------------------------------
            if (pass_step1) {
                if (selected_step1.pid == 211) {
                    pass_21 = true;
                }
                else {
                    n_21_fail_selected_not_pip++;
                }
            }

            // --------------------------------------------------------
            // Step 2.2:
            // from Step 1, exactly one positive track total.
            // --------------------------------------------------------
            if (pass_step1) {
                if (n_positive_tracks_total == 1) {
                    pass_22 = true;
                }
                else {
                    n_22_fail_not_exactly_one_positive++;
                }
            }

            // --------------------------------------------------------
            // Step 3.1:
            // from Step 2.1, no other positive tracks.
            // --------------------------------------------------------
            if (pass_21) {
                if (n_positive_tracks_total == 1) {
                    pass_31 = true;
                }
                else {
                    n_31_fail_extra_positive++;
                }
            }

            // --------------------------------------------------------
            // Step 3.2:
            // from Step 2.2, only positive track is identified as pion.
            // --------------------------------------------------------
            if (pass_22) {
                if (selected_step1.pid == 211) {
                    pass_32 = true;
                }
                else {
                    n_32_fail_only_positive_not_pip++;
                }
            }

            if (pass_step1) n_step1_pass++;
            if (pass_21)    n_21_pass++;
            if (pass_22)    n_22_pass++;
            if (pass_31)    n_31_pass++;
            if (pass_32)    n_32_pass++;

            if (pass_step1 && !reached_cap(n_step1_written, max_events_to_write)) {
                const double dp = selected_step1.p - p_lund;

                h_step1->Fill(dp);
                writer1.addEvent(event);
                n_step1_written++;

                if (n_step1_written <= 20) {
                    std::cout << "STEP1 KEEP"
                              << " idx=" << selected_step1.index
                              << " pid=" << selected_step1.pid
                              << " charge=" << selected_step1.charge
                              << " status=" << selected_step1.status
                              << " theta=" << selected_step1.theta_deg
                              << " p_rec=" << selected_step1.p
                              << " p_lund=" << p_lund
                              << " dp=" << dp
                              << " n_pos_total=" << n_positive_tracks_total
                              << " n_fd_pos=" << fd_positive_tracks.size()
                              << " n_neg_hadrons=" << n_negative_hadrons_total
                              << "\n";
                }
            }

            if (pass_21 && !reached_cap(n_21_written, max_events_to_write)) {
                const double dp = selected_step1.p - p_lund;

                h_21->Fill(dp);
                writer21.addEvent(event);
                n_21_written++;
            }

            if (pass_22 && !reached_cap(n_22_written, max_events_to_write)) {
                const double dp = selected_step1.p - p_lund;

                h_22->Fill(dp);
                writer22.addEvent(event);
                n_22_written++;
            }

            if (pass_31 && !reached_cap(n_31_written, max_events_to_write)) {
                const double dp = selected_step1.p - p_lund;

                h_31->Fill(dp);
                writer31.addEvent(event);
                n_31_written++;
            }

            if (pass_32 && !reached_cap(n_32_written, max_events_to_write)) {
                const double dp = selected_step1.p - p_lund;

                h_32->Fill(dp);
                writer32.addEvent(event);
                n_32_written++;
            }

            if (max_events_to_write > 0 &&
                n_step1_written >= max_events_to_write &&
                n_21_written    >= max_events_to_write &&
                n_22_written    >= max_events_to_write &&
                n_31_written    >= max_events_to_write &&
                n_32_written    >= max_events_to_write) {

                std::cout << "\nReached max_events_to_write for all five outputs.\n";
                break;
            }
        }

        if (max_events_to_write > 0 &&
            n_step1_written >= max_events_to_write &&
            n_21_written    >= max_events_to_write &&
            n_22_written    >= max_events_to_write &&
            n_31_written    >= max_events_to_write &&
            n_32_written    >= max_events_to_write) {

            break;
        }
    }

    if (writers_opened) {
        writer1.close();
        writer21.close();
        writer22.close();
        writer31.close();
        writer32.close();
    }
    else {
        std::cerr << "ERROR: writers were never opened. No HIPO files written.\n";
    }

    save_hist_png(
        h_step1,
        out_step1_png.c_str(),
        "Step 1: selected FD positive in #theta/p bin;#Delta p = p_{rec}(selected positive) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_step1_written
    );

    save_hist_png(
        h_21,
        out_21_png.c_str(),
        "Step 2.1: Step-1 selected track is #pi^{+};#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_21_written
    );

    save_hist_png(
        h_22,
        out_22_png.c_str(),
        "Step 2.2: exactly one positive track after Step 1;#Delta p = p_{rec}(only positive) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_22_written
    );

    save_hist_png(
        h_31,
        out_31_png.c_str(),
        "Step 3.1: Step-2.1 #pi^{+} is the only positive track;#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_31_written
    );

    save_hist_png(
        h_32,
        out_32_png.c_str(),
        "Step 3.2: Step-2.2 only positive track is #pi^{+};#Delta p = p_{rec}(only positive #pi^{+}) - p_{Lund}(#pi^{+}) [GeV];Counts",
        n_32_written
    );

    std::ofstream stats(out_stats_txt);

    if (!stats.is_open()) {
        std::cerr << "ERROR: could not open stats output file: "
                  << out_stats_txt << "\n";
    }
    else {
        stats << "Updated three-step / five-output positive-track study\n";
        stats << "Negative-hadron veto enabled\n";
        stats << "====================================================\n\n";

        stats << "Input .dat file:       " << input_dat_file << "\n";
        stats << "Output prefix:         " << prefix << "\n";
        stats << "Detector used for Step-1 candidate choice: " << det << "\n";
        stats << "Theta bin:             [" << theta_min_deg << ", " << theta_max_deg << ") deg\n";
        stats << "Momentum bin:          [" << p_min << ", " << p_max << ") GeV\n";
        stats << "Delta-p plot range:    [" << dp_min << ", " << dp_max << "] GeV\n";
        stats << "Base selection:        exactly 1 trigger electron, no other electrons, no negative hadrons\n";
        stats << "Generated truth:       MC::Lund only\n";
        stats << "Max events/output:     " << max_events_to_write << "\n";
        stats << "Readable HIPO files:   " << hipo_files.size() << "\n";
        stats << "Files skipped missing MC::Lund: " << n_fail_missing_mclund << "\n\n";

        stats << "Base selection cut flow\n";
        stats << "-----------------------\n";
        print_stat_line(stats, "Scanned events:", scanned_events, scanned_events, scanned_events);
        print_stat_line(stats, "Exactly 1 REC electron total and it is trigger electron:", n_base_one_trigger_electron, scanned_events, scanned_events);
        print_stat_line(stats, "No negative hadrons:", n_base_no_negative_hadrons, n_base_one_trigger_electron, scanned_events);
        print_stat_line(stats, "Exactly 1 MC::Lund pi+:", n_base_mclund_pip, n_base_no_negative_hadrons, scanned_events);

        stats << "\nStep 1 cut flow after base selection\n";
        stats << "------------------------------------\n";
        print_stat_line(stats, "Events with at least one requested-detector positive track:", n_step1_has_fd_positive, n_base_mclund_pip, scanned_events);
        print_stat_line(stats, "Selected requested-detector positive track inside theta/p bin:", n_step1_selected_in_theta_p, n_step1_has_fd_positive, scanned_events);
        print_stat_line(stats, "Step 1 pass:", n_step1_pass, n_step1_selected_in_theta_p, scanned_events);

        stats << "\nAdvisor-requested outputs\n";
        stats << "-------------------------\n";
        print_stat_line(stats, "Step 1: selected detector-positive in theta/p bin:", n_step1_pass, n_base_mclund_pip, scanned_events);
        print_stat_line(stats, "Step 2.1: Step-1 selected track is pi+:", n_21_pass, n_step1_pass, scanned_events);
        print_stat_line(stats, "Step 2.2: exactly one positive track after Step 1:", n_22_pass, n_step1_pass, scanned_events);
        print_stat_line(stats, "Step 3.1: Step-2.1 pion is only positive track:", n_31_pass, n_21_pass, scanned_events);
        print_stat_line(stats, "Step 3.2: Step-2.2 only positive track is pion:", n_32_pass, n_22_pass, scanned_events);

        stats << "\nConsistency check\n";
        stats << "-----------------\n";
        stats << "Step 3.1 pass: " << n_31_pass << "\n";
        stats << "Step 3.2 pass: " << n_32_pass << "\n";
        stats << "Difference Step3.1 - Step3.2: " << (n_31_pass - n_32_pass) << "\n";

        stats << "\nWritten counts after output cap\n";
        stats << "-------------------------------\n";
        stats << "Step 1 written:    " << n_step1_written << "   " << out_step1_hipo << "   " << out_step1_png << "\n";
        stats << "Step 2.1 written:  " << n_21_written    << "   " << out_21_hipo    << "   " << out_21_png << "\n";
        stats << "Step 2.2 written:  " << n_22_written    << "   " << out_22_hipo    << "   " << out_22_png << "\n";
        stats << "Step 3.1 written:  " << n_31_written    << "   " << out_31_hipo    << "   " << out_31_png << "\n";
        stats << "Step 3.2 written:  " << n_32_written    << "   " << out_32_hipo    << "   " << out_32_png << "\n";

        stats << "\nFailure counters\n";
        stats << "----------------\n";
        stats << "Fail base electron requirement:                    " << n_fail_base_electron << "\n";
        stats << "Fail no-negative-hadrons requirement:              " << n_fail_negative_hadrons << "\n";
        stats << "Fail exactly 1 MC::Lund pi+ check:                 " << n_fail_lund_pip << "\n";
        stats << "Step 1 fail: no requested-detector positive track: " << n_step1_fail_no_fd_positive << "\n";
        stats << "Step 1 fail: selected positive outside theta/p:    " << n_step1_fail_selected_not_in_bin << "\n";
        stats << "Step 2.1 fail: selected Step-1 track is not pi+:   " << n_21_fail_selected_not_pip << "\n";
        stats << "Step 2.2 fail: not exactly one positive track:     " << n_22_fail_not_exactly_one_positive << "\n";
        stats << "Step 3.1 fail: Step-2.1 has extra positives:       " << n_31_fail_extra_positive << "\n";
        stats << "Step 3.2 fail: only positive is not pi+:           " << n_32_fail_only_positive_not_pip << "\n";

        stats << "\nHistogram statistics\n";
        stats << "--------------------\n";
        stats << "Step 1 entries:    " << h_step1->GetEntries()
              << "   mean=" << h_step1->GetMean()
              << "   stddev=" << h_step1->GetStdDev() << "\n";
        stats << "Step 2.1 entries:  " << h_21->GetEntries()
              << "   mean=" << h_21->GetMean()
              << "   stddev=" << h_21->GetStdDev() << "\n";
        stats << "Step 2.2 entries:  " << h_22->GetEntries()
              << "   mean=" << h_22->GetMean()
              << "   stddev=" << h_22->GetStdDev() << "\n";
        stats << "Step 3.1 entries:  " << h_31->GetEntries()
              << "   mean=" << h_31->GetMean()
              << "   stddev=" << h_31->GetStdDev() << "\n";
        stats << "Step 3.2 entries:  " << h_32->GetEntries()
              << "   mean=" << h_32->GetMean()
              << "   stddev=" << h_32->GetStdDev() << "\n";

        stats.close();
    }

    std::cout << "\n============================================================\n";
    std::cout << "DONE UPDATED SINGLE-PASS STUDY\n";
    std::cout << "Scanned events:                         " << scanned_events << "\n";
    std::cout << "Base: exactly 1 trigger electron only:  " << n_base_one_trigger_electron << "\n";
    std::cout << "Base: no negative hadrons:              " << n_base_no_negative_hadrons << "\n";
    std::cout << "Base: exactly 1 MC::Lund pi+:           " << n_base_mclund_pip << "\n";
    std::cout << "Step 1 pass:                            " << n_step1_pass << "\n";
    std::cout << "Step 2.1 pass:                          " << n_21_pass << "\n";
    std::cout << "Step 2.2 pass:                          " << n_22_pass << "\n";
    std::cout << "Step 3.1 pass:                          " << n_31_pass << "\n";
    std::cout << "Step 3.2 pass:                          " << n_32_pass << "\n";
    std::cout << "Step 3.1 - Step 3.2:                    " << (n_31_pass - n_32_pass) << "\n";

    std::cout << "\nFailures:\n";
    std::cout << "  fail base electron requirement:       " << n_fail_base_electron << "\n";
    std::cout << "  fail no-negative-hadrons requirement: " << n_fail_negative_hadrons << "\n";
    std::cout << "  fail MC::Lund pi+ check:              " << n_fail_lund_pip << "\n";

    std::cout << "\nOutput files:\n";
    std::cout << "  " << out_step1_hipo << "   " << out_step1_png << "\n";
    std::cout << "  " << out_21_hipo    << "   " << out_21_png << "\n";
    std::cout << "  " << out_22_hipo    << "   " << out_22_png << "\n";
    std::cout << "  " << out_31_hipo    << "   " << out_31_png << "\n";
    std::cout << "  " << out_32_hipo    << "   " << out_32_png << "\n";
    std::cout << "  " << out_stats_txt  << "\n";

    delete h_step1;
    delete h_21;
    delete h_22;
    delete h_31;
    delete h_32;
}