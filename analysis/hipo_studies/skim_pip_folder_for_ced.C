#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <algorithm>

#include "TMath.h"
#include "TVector3.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"

#include "clas12reader.h"
#include "HipoChainWriter.h"

#include <fstream>

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

        // Remove Windows carriage return if present
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        // Optional safety: only accept .hipo paths
        if (!ends_with(line, ".hipo")) continue;

        hipo_files.push_back(line);
    }

    return hipo_files;
}

void skim_pip_folder_for_ced(const char* input_dat_file,
                             const char* output_hipo,
                             const char* detector,
                             double theta_min_deg = 38.0,
                             double theta_max_deg = 39.0,
                             double p_min = 1.8,
                             double p_max = 2.0,
                             double dp_min = -0.05,
                             double dp_max = 0.05,
                             int max_events_to_write = 500
                             )
{
    std::cout << "Macro started\n";
    std::cout << "Input .dat file: " << input_dat_file << "\n";
    std::cout << "Output file:  " << output_hipo << "\n\n";

    std::string det = detector;

    if (det != "FD" && det != "CD") {
        std::cerr << "ERROR: detector must be \"FD\" or \"CD\", got: "
                  << det << std::endl;
        return;
    }

    std::cout << "Detector:        " << det << "\n\n";

    auto hipo_files = get_hipo_files_from_dat(input_dat_file);

    if (hipo_files.empty()) {
        std::cerr << "ERROR: No .hipo files found in .dat file: "<< input_dat_file << std::endl;
        return;
    }

    std::cout << "Found " << hipo_files.size() << " .hipo files\n";

    for (size_t i = 0; i < hipo_files.size(); i++) {
        if (i < 10) {
            std::cout << "  [" << i << "] " << hipo_files[i] << "\n";
        }
    }

    if (hipo_files.size() > 10) {
        std::cout << "  ...\n";
    }

    std::cout << "\nSelection:\n"
          << "  exactly 1 reconstructed electron with status < 0\n"
          << "  exactly 1 reconstructed pi+\n"
          << "  pion theta_rec in [" << theta_min_deg << ", " << theta_max_deg << ") deg\n"
          << "  pion p_rec     in [" << p_min << ", " << p_max << ") GeV\n"
          << "  pion Vz        in (-10, 2) cm\n"
          << "  pion |Vx| < 2 cm, |Vy| < 2 cm\n"
          << "  delta_p        in [" << dp_min << ", " << dp_max << ") GeV\n"
          << "  detector: " << det << "\n";
            

    clas12root::HipoChainWriter chain(output_hipo);

    for (const auto& f : hipo_files) {
        chain.Add(f.c_str());
    }

    chain.GetWriter().writeSpecialBanks(true);
    chain.db()->turnOffQADB();

    auto config_c12 = chain.GetC12Reader();
    auto& c12 = chain.C12ref();

    int scanned_events = 0;
    int written_events = 0;

    long long n_events_with_one_trigger_e = 0;
    long long n_events_with_one_pip = 0;

    long long n_rec_pip_candidates = 0;
    long long n_detector_status_candidates = 0;
    long long n_theta_pass = 0;
    long long n_p_pass = 0;
    long long n_vz_pass = 0;
    long long n_vxy_pass = 0;
    long long n_mc_match_pass = 0;
    long long n_dp_pass = 0;

   while (chain.Next()) {
    scanned_events++;

    // ------------------------------------------------------------
    // Require one reconstructed electron with negative status and maybe also some other electrons
    // ------------------------------------------------------------
    auto electrons = c12->getByID(11);

    int n_negative_status_e = 0;

    for (auto& ele : electrons) {
        if (ele->par()->getStatus() < 0) {
            n_negative_status_e++;
        }
    }

    if (n_negative_status_e != 1) continue;

    n_events_with_one_trigger_e++;

    // ------------------------------------------------------------
    // Require exactly one reconstructed pi+
    // ------------------------------------------------------------
    auto pips = c12->getByID(211);

    if (pips.size() != 1) continue;

    n_events_with_one_pip++;

    for (auto& pip : pips) {
        n_rec_pip_candidates++;

           const int status_raw = pip->par()->getStatus();
           const int status = std::abs(status_raw);

            bool detector_pass = false;

            if (det == "FD") {
                detector_pass = (status < 4000);
            }
            else if (det == "CD") {
                detector_pass = (status >= 4000);
            }

            if (!detector_pass) continue;
            n_detector_status_candidates++;

           const double p_rec = pip->getP();
            const double theta_rec_deg = pip->getTheta() * TMath::RadToDeg();

            const double vx_pip = pip->par()->getVx();
            const double vy_pip = pip->par()->getVy();
            const double vz_pip = pip->par()->getVz();

            if (theta_rec_deg < theta_min_deg || theta_rec_deg >= theta_max_deg) continue;
            n_theta_pass++;

            if (p_rec < p_min || p_rec >= p_max) continue;
            n_p_pass++;

            // Same as: vz_piplus > -10 && vz_piplus < 2
            if (!(vz_pip > -10.0 && vz_pip < 2.0)) continue;
            n_vz_pass++;

            // Same as: abs(vx_piplus) < 2.0 && abs(vy_piplus) < 2.0
            if (!(std::abs(vx_pip) < 2.0 && std::abs(vy_pip) < 2.0)) continue;
            n_vxy_pass++;

            auto mc = c12->mcparts();

            double p_gen = std::numeric_limits<double>::quiet_NaN();
            int n_gen_pip = 0;

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

            // Your files should have exactly one generated pi+.
            // If not, skip because something is not what we think it is.
            if (n_gen_pip != 1) continue;
            if (!std::isfinite(p_gen)) continue;

            n_mc_match_pass++;

            const double delta_p = p_rec - p_gen;

            if (delta_p < dp_min || delta_p >= dp_max) continue;
            n_dp_pass++;

            std::cout << "KEEP"
                      << "  run=" << c12->runconfig()->getRun()
                      << "  event=" << c12->runconfig()->getEvent()
                      << "  p_rec=" << std::fixed << std::setprecision(4) << p_rec
                      << "  p_gen=" << p_gen
                      << "  dp=" << delta_p
                      << "  theta=" << theta_rec_deg
                      << "  status=" << status_raw
                      << "  n_gen_pip=" << n_gen_pip
                      << "\n";

            chain.WriteEvent();
            written_events++;

            break; // write event once, even if multiple candidates pass
        }

        if (max_events_to_write > 0 && written_events >= max_events_to_write) {
            std::cout << "\nReached max_events_to_write = "
                      << max_events_to_write << "\n";
            break;
        }

        if (scanned_events % 50000 == 0) {
            std::cout << "Scanned " << scanned_events
                      << " events, written " << written_events << "\n";
        }
    }

  std::cout << "\nDone.\n"
          << "Scanned events:                         " << scanned_events << "\n"
          << "Events with exactly 1 trigger electron: " << n_events_with_one_trigger_e << "\n"
          << "Events with exactly 1 REC pi+:          " << n_events_with_one_pip << "\n"
          << "REC pi+ candidates:                     " << n_rec_pip_candidates << "\n"
          << det << "-status pi+ candidates:               " << n_detector_status_candidates << "\n"
          << "After theta cut:                        " << n_theta_pass << "\n"
          << "After p cut:                            " << n_p_pass << "\n"
          << "After pion Vz cut:                      " << n_vz_pass << "\n"
          << "After pion Vxy cut:                     " << n_vxy_pass << "\n"
          << "After generated pi+ check:              " << n_mc_match_pass << "\n"
          << "After delta_p cut:                      " << n_dp_pass << "\n"
          << "Written events:                         " << written_events << "\n"
          << "Output file:                            " << output_hipo << "\n";
}