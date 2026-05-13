// skim_posTrack_and_pip_for_ced.cxx
//
// Makes TWO output HIPO files:
//
// 1) Positive-track topology sample:
//    - exactly 1 trigger electron: pid = 11 and status < 0
//    - no other reconstructed electrons
//    - exactly 1 positively charged reconstructed track: charge > 0
//      This can be pi+, proton, K+, mis-ID, etc.
//    - optional gammas allowed
//    - nothing else allowed
//
// 2) REC pi+ topology sample:
//    - exactly 1 trigger electron: pid = 11 and status < 0
//    - no other reconstructed electrons
//    - exactly 1 reconstructed pi+: pid = 211
//    - optional gammas allowed
//    - nothing else allowed
//
// Common cuts:
//    - detector = FD or CD, based on abs(status)
//    - theta_rec bin
//    - p_rec bin
//    - exactly 1 generated pi+
//    - NO vertex cuts
//    - NO delta_p cut
//
// Example:
//
// clas12root -l -b -q 'skim_posTrack_and_pip_for_ced.cxx+(
//     "/home/bulgakov/projects/momcor/pi_plus_energy_loss/utils/hipo2root/pi_plus_toy_FD_CD_full.dat",
//     "FD_posTrack_noVz_allDp.hipo",
//     "FD_pip_noVz_allDp.hipo",
//     "FD",
//     38,39,
//     1.0,1.2,
//     -1
// )'

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

#include "TMath.h"

#include "clas12reader.h"
#include "HipoChainWriter.h"

using namespace clas12;

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

static bool detector_status_pass(int status_raw, const std::string& det)
{
    const int status = std::abs(status_raw);

    if (det == "FD") {
        return status < 4000;
    }

    if (det == "CD") {
        return status >= 4000;
    }

    return false;
}

static int process_one_selection(const std::vector<std::string>& hipo_files,
                                 const char* output_hipo,
                                 const std::string& det,
                                 bool require_rec_pip,
                                 double theta_min_deg,
                                 double theta_max_deg,
                                 double p_min,
                                 double p_max,
                                 int max_events_to_write)
{
    std::cout << "\n============================================================\n";

    if (require_rec_pip) {
        std::cout << "PASS 2: writing REC pi+ topology sample\n";
    }
    else {
        std::cout << "PASS 1: writing positive-track topology sample\n";
    }

    std::cout << "Output file: " << output_hipo << "\n";

    if (require_rec_pip) {
        std::cout << "Topology:\n"
                  << "  exactly 1 trigger electron: pid = 11 and status < 0\n"
                  << "  no other reconstructed electrons allowed\n"
                  << "  exactly 1 reconstructed pi+: pid = 211\n"
                  << "  optional gammas allowed\n"
                  << "  no other REC::Particle entries allowed\n";
    }
    else {
        std::cout << "Topology:\n"
                  << "  exactly 1 trigger electron: pid = 11 and status < 0\n"
                  << "  no other reconstructed electrons allowed\n"
                  << "  exactly 1 positively charged reconstructed track: charge > 0\n"
                  << "  optional gammas allowed\n"
                  << "  no other REC::Particle entries allowed\n";
    }

    std::cout << "Cuts:\n"
              << "  no vertex cuts applied\n"
              << "  detector: " << det << "\n"
              << "  theta_rec in [" << theta_min_deg << ", " << theta_max_deg << ") deg\n"
              << "  p_rec     in [" << p_min << ", " << p_max << ") GeV\n"
              << "  generated pi+ check: exactly 1 generated pi+\n"
              << "  delta_p cut: NONE\n"
              << "  max events = " << max_events_to_write
              << "  (-1 or 0 means write all)\n";

    clas12root::HipoChainWriter chain(output_hipo);

    for (const auto& f : hipo_files) {
        chain.Add(f.c_str());
    }

    chain.GetWriter().writeSpecialBanks(true);
    chain.db()->turnOffQADB();

    auto config_c12 = chain.GetC12Reader();
    auto& c12 = chain.C12ref();

    long long scanned_events = 0;
    long long written_events = 0;

    long long n_one_trigger_e = 0;
    long long n_no_other_e = 0;
    long long n_topology_pass = 0;

    long long n_selected_candidates = 0;
    long long n_detector_pass = 0;
    long long n_theta_pass = 0;
    long long n_p_pass = 0;
    long long n_gen_pip_pass = 0;
    long long n_dp_computed = 0;

    long long n_fail_trigger_e = 0;
    long long n_fail_other_e = 0;
    long long n_fail_topology = 0;
    long long n_fail_detector = 0;
    long long n_fail_theta = 0;
    long long n_fail_p = 0;
    long long n_fail_gen_pip = 0;

    while (chain.Next()) {
        scanned_events++;

        if (scanned_events % 50000 == 0) {
            std::cout << "Scanned " << scanned_events
                      << " events, written " << written_events << "\n";
        }

        auto rec_particles = c12->getDetParticles();
        const int n_rec_total = rec_particles.size();

        int n_trigger_e = 0;
        int n_other_e = 0;
        int n_gamma = 0;

        int n_positive_tracks = 0;
        int n_rec_pip = 0;
        int n_forbidden = 0;

        clas12::region_part_ptr selected_track = nullptr;

        for (auto& part : rec_particles) {
            const int pid = part->par()->getPid();
            const int status_part = part->par()->getStatus();
            const int charge = part->par()->getCharge();

            const bool is_trigger_e = (pid == 11 && status_part < 0);
            const bool is_other_e   = (pid == 11 && status_part >= 0);
            const bool is_gamma     = (pid == 22);
            const bool is_pip       = (pid == 211);
            const bool is_positive_track = (charge > 0);

            if (is_trigger_e) {
                n_trigger_e++;
                continue;
            }

            if (is_other_e) {
                n_other_e++;
                continue;
            }

            if (is_gamma) {
                n_gamma++;
                continue;
            }

            if (require_rec_pip) {
                if (is_pip) {
                    n_rec_pip++;
                    selected_track = part;
                    continue;
                }

                n_forbidden++;
                continue;
            }
            else {
                if (is_positive_track) {
                    n_positive_tracks++;
                    selected_track = part;
                    continue;
                }

                n_forbidden++;
                continue;
            }
        }

        if (n_trigger_e != 1) {
            n_fail_trigger_e++;
            continue;
        }
        n_one_trigger_e++;

        if (n_other_e != 0) {
            n_fail_other_e++;
            continue;
        }
        n_no_other_e++;

        if (require_rec_pip) {
            if (n_rec_pip != 1 || n_forbidden != 0 || selected_track == nullptr) {
                n_fail_topology++;
                continue;
            }
        }
        else {
            if (n_positive_tracks != 1 || n_forbidden != 0 || selected_track == nullptr) {
                n_fail_topology++;
                continue;
            }
        }

        n_topology_pass++;
        n_selected_candidates++;

        const int pid_rec = selected_track->par()->getPid();
        const int charge_rec = selected_track->par()->getCharge();
        const int status_raw = selected_track->par()->getStatus();

        if (!detector_status_pass(status_raw, det)) {
            n_fail_detector++;
            continue;
        }
        n_detector_pass++;

        const double p_rec = selected_track->getP();
        const double theta_rec_deg = selected_track->getTheta() * TMath::RadToDeg();

        const double vx = selected_track->par()->getVx();
        const double vy = selected_track->par()->getVy();
        const double vz = selected_track->par()->getVz();

        if (theta_rec_deg < theta_min_deg || theta_rec_deg >= theta_max_deg) {
            n_fail_theta++;
            continue;
        }
        n_theta_pass++;

        if (p_rec < p_min || p_rec >= p_max) {
            n_fail_p++;
            continue;
        }
        n_p_pass++;

        // ------------------------------------------------------------
        // Generated pi+ check.
        // Keep this inline to avoid clas12reader type-name issues.
        // ------------------------------------------------------------
        double p_gen = std::numeric_limits<double>::quiet_NaN();
        int n_gen_pip = 0;

        auto mc = c12->mcparts();

        for (int j = 0; j < mc->getRows(); j++) {
            if (mc->getPid(j) != 211) continue;

            const double px_gen = mc->getPx(j);
            const double py_gen = mc->getPy(j);
            const double pz_gen = mc->getPz(j);

            p_gen = std::sqrt(px_gen * px_gen +
                              py_gen * py_gen +
                              pz_gen * pz_gen);

            n_gen_pip++;
        }

        if (n_gen_pip != 1 || !std::isfinite(p_gen)) {
            n_fail_gen_pip++;
            continue;
        }

        n_gen_pip_pass++;

        const double delta_p = p_rec - p_gen;
        n_dp_computed++;

        std::cout << "KEEP"
                  << "  sample=" << (require_rec_pip ? "REC_pip" : "positive_track")
                  << "  run=" << c12->runconfig()->getRun()
                  << "  event=" << c12->runconfig()->getEvent()
                  << "  det=" << det
                  << "  rec_pid=" << pid_rec
                  << "  rec_charge=" << charge_rec
                  << "  status=" << status_raw
                  << "  theta=" << theta_rec_deg
                  << "  p_rec=" << p_rec
                  << "  p_gen_pi+=" << p_gen
                  << "  dp=" << delta_p
                  << "  vx=" << vx
                  << "  vy=" << vy
                  << "  vz=" << vz
                  << "  n_rec_total=" << n_rec_total
                  << "  n_gamma=" << n_gamma
                  << "\n";

        chain.WriteEvent();
        written_events++;

        if (max_events_to_write > 0 && written_events >= max_events_to_write) {
            std::cout << "\nReached max_events_to_write = "
                      << max_events_to_write << "\n";
            break;
        }
    }

    std::cout << "\nSummary for output: " << output_hipo << "\n"
              << "Scanned events:                         " << scanned_events << "\n"
              << "Events with exactly 1 trigger electron: " << n_one_trigger_e << "\n"
              << "Events with no other electrons:         " << n_no_other_e << "\n"
              << "Events passing REC topology:            " << n_topology_pass << "\n"
              << "Selected-track candidates:              " << n_selected_candidates << "\n"
              << "After detector-status cut:              " << n_detector_pass << "\n"
              << "After theta cut:                        " << n_theta_pass << "\n"
              << "After p cut:                            " << n_p_pass << "\n"
              << "After generated pi+ check:              " << n_gen_pip_pass << "\n"
              << "After delta_p computed:                 " << n_dp_computed << "\n"
              << "Written events:                         " << written_events << "\n"
              << "\nFailures:\n"
              << "  fail trigger electron count:          " << n_fail_trigger_e << "\n"
              << "  fail other-electron veto:             " << n_fail_other_e << "\n"
              << "  fail REC topology:                    " << n_fail_topology << "\n"
              << "  fail detector-status cut:             " << n_fail_detector << "\n"
              << "  fail theta cut:                       " << n_fail_theta << "\n"
              << "  fail p cut:                           " << n_fail_p << "\n"
              << "  fail generated pi+ check:             " << n_fail_gen_pip << "\n";

    return written_events;
}

void skim_posTrack_and_pip_for_ced(const char* input_dat_file,
                                   const char* output_pos_track_hipo,
                                   const char* output_pip_hipo,
                                   const char* detector = "FD",
                                   double theta_min_deg = 38.0,
                                   double theta_max_deg = 39.0,
                                   double p_min = 1.0,
                                   double p_max = 1.2,
                                   int max_events_to_write = 1000)
{
    std::cout << "Macro started\n";
    std::cout << "Input .dat file:              " << input_dat_file << "\n";
    std::cout << "Output positive-track HIPO:   " << output_pos_track_hipo << "\n";
    std::cout << "Output REC pi+ HIPO:          " << output_pip_hipo << "\n";

    std::string det = detector;

    if (det != "FD" && det != "CD") {
        std::cerr << "ERROR: detector must be \"FD\" or \"CD\", got: "
                  << det << std::endl;
        return;
    }

    auto hipo_files = get_hipo_files_from_dat(input_dat_file);

    if (hipo_files.empty()) {
        std::cerr << "ERROR: No .hipo files found in .dat file: "
                  << input_dat_file << std::endl;
        return;
    }

    std::cout << "\nFound " << hipo_files.size() << " .hipo files\n";

    for (size_t i = 0; i < hipo_files.size(); i++) {
        if (i < 10) {
            std::cout << "  [" << i << "] " << hipo_files[i] << "\n";
        }
    }

    if (hipo_files.size() > 10) {
        std::cout << "  ...\n";
    }

    std::cout << "\nGlobal selection parameters:\n"
              << "  detector:        " << det << "\n"
              << "  theta_rec:       [" << theta_min_deg << ", " << theta_max_deg << ") deg\n"
              << "  p_rec:           [" << p_min << ", " << p_max << ") GeV\n"
              << "  vertex cuts:     NONE\n"
              << "  delta_p cut:     NONE\n"
              << "  max per output:  " << max_events_to_write << "\n";

    const int n_pos = process_one_selection(hipo_files,
                                            output_pos_track_hipo,
                                            det,
                                            false,
                                            theta_min_deg,
                                            theta_max_deg,
                                            p_min,
                                            p_max,
                                            max_events_to_write);

    const int n_pip = process_one_selection(hipo_files,
                                            output_pip_hipo,
                                            det,
                                            true,
                                            theta_min_deg,
                                            theta_max_deg,
                                            p_min,
                                            p_max,
                                            max_events_to_write);

    std::cout << "\n============================================================\n"
              << "DONE BOTH SKIMS\n"
              << "Positive-track output: " << output_pos_track_hipo
              << "  written events = " << n_pos << "\n"
              << "REC pi+ output:        " << output_pip_hipo
              << "  written events = " << n_pip << "\n";
}