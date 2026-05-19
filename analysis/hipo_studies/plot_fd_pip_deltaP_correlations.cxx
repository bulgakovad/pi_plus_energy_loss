// plot_fd_pip_deltaP_correlations.cxx
//
// Selection:
//   - exactly 1 trigger electron: pid = 11 and status < 0
//   - exactly 1 REC pi+: pid = 211
//   - REC pi+ is FD-like: abs(status) < 4000
//   - REC pi+ is inside theta/p bin
//   - exactly 1 GEN pi+: pid = 211
//   - anything else in REC::Particle is allowed
//
// Plots:
//   deltaP = p_rec - p_gen
//   deltaTheta = theta_rec - theta_gen
//   deltaVx = vx_rec - vx_gen
//   deltaVy = vy_rec - vy_gen
//   deltaVz = vz_rec - vz_gen
//
// Output PNGs:
//   <prefix>_deltaP_vs_deltaTheta.png
//   <prefix>_deltaP_vs_deltaVx.png
//   <prefix>_deltaP_vs_deltaVy.png
//   <prefix>_deltaP_vs_deltaVz.png
//
// Example:
//
// clas12root -l -b -q 'plot_fd_pip_deltaP_correlations.cxx+(
//   "good_hipo.dat",
//   "FD_pip_theta38_39_p1p0_1p2",
//   38,39,
//   1.0,1.2,
//   -0.5,0.5,
//   -5,5,
//   -2,2,
//   -2,2,
//   -10,10,
//   200,
//   200,
//   -1
// )'

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

#include "TMath.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLatex.h"

#include "hipo4/reader.h"
#include "hipo4/event.h"
#include "hipo4/bank.h"
#include "hipo4/dictionary.h"

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

static bool is_fd_status(int status_raw)
{
    return std::abs(status_raw) < 4000;
}

static void save_2d_png(TH2D* h,
                        const char* output_png,
                        const char* title,
                        long long n_filled,
                        bool logz = true)
{
    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas(Form("c_%s", h->GetName()), title, 950, 750);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.10);

    if (logz) {
        c->SetLogz();
    }

    h->SetTitle(title);
    h->Draw("COLZ");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.032);
    lat.DrawLatex(0.15, 0.93, Form("Filled events: %lld", n_filled));

    c->SaveAs(output_png);

    delete c;
}

void plot_fd_pip_deltaP_correlations(const char* input_dat_file,
                                     const char* output_prefix = "FD_pip_deltaP_corr",
                                     double theta_min_deg = 38.0,
                                     double theta_max_deg = 39.0,
                                     double p_min = 1.0,
                                     double p_max = 1.2,
                                     double dp_min = -0.50,
                                     double dp_max = 0.50,
                                     double dtheta_min = -5.0,
                                     double dtheta_max = 5.0,
                                     double dvx_min = -2.0,
                                     double dvx_max = 2.0,
                                     double dvy_min = -2.0,
                                     double dvy_max = 2.0,
                                     double dvz_min = -10.0,
                                     double dvz_max = 10.0,
                                     int n_x_bins = 200,
                                     int n_dp_bins = 200,
                                     int max_events_to_fill = -1)
{
    std::cout << "Macro started\n";
    std::cout << "Input .dat file: " << input_dat_file << "\n";
    std::cout << "Output prefix:   " << output_prefix << "\n";
    std::cout << "Selection:\n"
              << "  exactly 1 trigger electron: pid=11, status<0\n"
              << "  exactly 1 REC pi+: pid=211\n"
              << "  REC pi+ FD-like: abs(status)<4000\n"
              << "  REC pi+ theta in [" << theta_min_deg << ", " << theta_max_deg << ") deg\n"
              << "  REC pi+ p     in [" << p_min << ", " << p_max << ") GeV\n"
              << "  exactly 1 GEN pi+: pid=211\n"
              << "  everything else allowed\n"
              << "  max fills = " << max_events_to_fill << " (-1 or 0 means all)\n\n";

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

    TH2D* h_dp_vs_dtheta = new TH2D(
        "h_dp_vs_dtheta",
        "#Delta p vs #Delta#theta;#Delta#theta = #theta_{rec} - #theta_{gen} [deg];#Delta p = p_{rec} - p_{gen} [GeV]",
        n_x_bins, dtheta_min, dtheta_max,
        n_dp_bins, dp_min, dp_max
    );

    TH2D* h_dp_vs_dvx = new TH2D(
        "h_dp_vs_dvx",
        "#Delta p vs #Delta v_{x};#Delta v_{x} = v_{x}^{rec} - v_{x}^{gen} [cm];#Delta p = p_{rec} - p_{gen} [GeV]",
        n_x_bins, dvx_min, dvx_max,
        n_dp_bins, dp_min, dp_max
    );

    TH2D* h_dp_vs_dvy = new TH2D(
        "h_dp_vs_dvy",
        "#Delta p vs #Delta v_{y};#Delta v_{y} = v_{y}^{rec} - v_{y}^{gen} [cm];#Delta p = p_{rec} - p_{gen} [GeV]",
        n_x_bins, dvy_min, dvy_max,
        n_dp_bins, dp_min, dp_max
    );

    TH2D* h_dp_vs_dvz = new TH2D(
        "h_dp_vs_dvz",
        "#Delta p vs #Delta v_{z};#Delta v_{z} = v_{z}^{rec} - v_{z}^{gen} [cm];#Delta p = p_{rec} - p_{gen} [GeV]",
        n_x_bins, dvz_min, dvz_max,
        n_dp_bins, dp_min, dp_max
    );

    long long scanned_events = 0;
    long long filled_events = 0;

    long long n_fail_trigger = 0;
    long long n_fail_rec_pip = 0;
    long long n_fail_fd = 0;
    long long n_fail_theta_p = 0;
    long long n_fail_gen_pip = 0;

    long long n_missing_runconfig_bank = 0;

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

        hipo::schema rec_particle_schema;
        hipo::schema mc_schema;
        hipo::schema runconfig_schema;

        const bool have_rec_particle =
            get_schema_safe(factory, "REC::Particle", rec_particle_schema);

        if (!have_rec_particle) {
            std::cerr << "WARNING: REC::Particle missing, skipping file:\n"
                      << "  " << input_hipo << "\n";
            continue;
        }

        bool have_mc = false;

        if (get_schema_safe(factory, "MC::Lund", mc_schema)) {
            have_mc = true;
        }
        else if (get_schema_safe(factory, "MC::Particle", mc_schema)) {
            have_mc = true;
        }

        if (!have_mc) {
            std::cerr << "WARNING: neither MC::Lund nor MC::Particle exists, skipping file:\n"
                      << "  " << input_hipo << "\n";
            continue;
        }

        const bool have_runconfig =
            get_schema_safe(factory, "RUN::config", runconfig_schema);

        if (!have_runconfig) {
            n_missing_runconfig_bank++;
        }

        hipo::event event;

        hipo::bank rec_particle(rec_particle_schema);
        hipo::bank mc_bank(mc_schema);
        hipo::bank runconfig(runconfig_schema);

        while (reader.next()) {
            reader.read(event);
            scanned_events++;

            if (scanned_events % 50000 == 0) {
                std::cout << "Scanned " << scanned_events
                          << ", filled " << filled_events << "\n";
            }

            event.getStructure(rec_particle);
            event.getStructure(mc_bank);
            if (have_runconfig) event.getStructure(runconfig);

            int n_trigger_e = 0;
            int n_rec_pip = 0;
            int pip_rec_row = -1;

            for (int i = 0; i < rec_particle.getRows(); i++) {
                const int pid = rec_particle.getInt("pid", i);
                const int status = rec_particle.getInt("status", i);

                if (pid == 11 && status < 0) {
                    n_trigger_e++;
                }

                if (pid == 211) {
                    n_rec_pip++;
                    pip_rec_row = i;
                }
            }

            if (n_trigger_e != 1) {
                n_fail_trigger++;
                continue;
            }

            if (n_rec_pip != 1 || pip_rec_row < 0) {
                n_fail_rec_pip++;
                continue;
            }

            const int pip_status = rec_particle.getInt("status", pip_rec_row);

            if (!is_fd_status(pip_status)) {
                n_fail_fd++;
                continue;
            }

            const double pip_rec_px = rec_particle.getFloat("px", pip_rec_row);
            const double pip_rec_py = rec_particle.getFloat("py", pip_rec_row);
            const double pip_rec_pz = rec_particle.getFloat("pz", pip_rec_row);

            const double p_rec = momentum(pip_rec_px, pip_rec_py, pip_rec_pz);
            const double theta_rec_deg =
                theta_deg_from_p(pip_rec_px, pip_rec_py, pip_rec_pz);

            if (theta_rec_deg < theta_min_deg || theta_rec_deg >= theta_max_deg ||
                p_rec < p_min || p_rec >= p_max) {
                n_fail_theta_p++;
                continue;
            }

            const double vx_rec = rec_particle.getFloat("vx", pip_rec_row);
            const double vy_rec = rec_particle.getFloat("vy", pip_rec_row);
            const double vz_rec = rec_particle.getFloat("vz", pip_rec_row);

            int n_gen_pip = 0;
            int pip_gen_row = -1;

            for (int i = 0; i < mc_bank.getRows(); i++) {
                const int pid = mc_bank.getInt("pid", i);

                if (pid == 211) {
                    n_gen_pip++;
                    pip_gen_row = i;
                }
            }

            if (n_gen_pip != 1 || pip_gen_row < 0) {
                n_fail_gen_pip++;
                continue;
            }

            const double pip_gen_px = mc_bank.getFloat("px", pip_gen_row);
            const double pip_gen_py = mc_bank.getFloat("py", pip_gen_row);
            const double pip_gen_pz = mc_bank.getFloat("pz", pip_gen_row);

            const double p_gen = momentum(pip_gen_px, pip_gen_py, pip_gen_pz);
            const double theta_gen_deg =
                theta_deg_from_p(pip_gen_px, pip_gen_py, pip_gen_pz);

            const double vx_gen = mc_bank.getFloat("vx", pip_gen_row);
            const double vy_gen = mc_bank.getFloat("vy", pip_gen_row);
            const double vz_gen = mc_bank.getFloat("vz", pip_gen_row);

            if (!std::isfinite(p_rec) ||
                !std::isfinite(p_gen) ||
                !std::isfinite(theta_rec_deg) ||
                !std::isfinite(theta_gen_deg)) {
                continue;
            }

            const double delta_p = p_rec - p_gen;
            const double delta_theta = theta_rec_deg - theta_gen_deg;
            const double delta_vx = vx_rec - vx_gen;
            const double delta_vy = vy_rec - vy_gen;
            const double delta_vz = vz_rec - vz_gen;

            h_dp_vs_dtheta->Fill(delta_theta, delta_p);
            h_dp_vs_dvx->Fill(delta_vx, delta_p);
            h_dp_vs_dvy->Fill(delta_vy, delta_p);
            h_dp_vs_dvz->Fill(delta_vz, delta_p);

            filled_events++;

            if (filled_events <= 20) {
                int run = -1;
                int evnum = -1;

                if (have_runconfig && runconfig.getRows() > 0) {
                    run = runconfig.getInt("run", 0);
                    evnum = runconfig.getInt("event", 0);
                }

                std::cout << "FILL"
                          << " run=" << run
                          << " event=" << evnum
                          << " p_rec=" << p_rec
                          << " p_gen=" << p_gen
                          << " dp=" << delta_p
                          << " theta_rec=" << theta_rec_deg
                          << " theta_gen=" << theta_gen_deg
                          << " dtheta=" << delta_theta
                          << " dvx=" << delta_vx
                          << " dvy=" << delta_vy
                          << " dvz=" << delta_vz
                          << "\n";
            }

            if (max_events_to_fill > 0 &&
                filled_events >= max_events_to_fill) {
                std::cout << "\nReached max_events_to_fill = "
                          << max_events_to_fill << "\n";
                break;
            }
        }

        if (max_events_to_fill > 0 &&
            filled_events >= max_events_to_fill) {
            break;
        }
    }

    const std::string prefix = output_prefix;

    const std::string png_dtheta =
        prefix + "_deltaP_vs_deltaTheta.png";

    const std::string png_dvx =
        prefix + "_deltaP_vs_deltaVx.png";

    const std::string png_dvy =
        prefix + "_deltaP_vs_deltaVy.png";

    const std::string png_dvz =
        prefix + "_deltaP_vs_deltaVz.png";

    save_2d_png(
        h_dp_vs_dtheta,
        png_dtheta.c_str(),
        "#Delta p vs #Delta#theta;#Delta#theta = #theta_{rec} - #theta_{gen} [deg];#Delta p = p_{rec} - p_{gen} [GeV]",
        filled_events,
        true
    );

    save_2d_png(
        h_dp_vs_dvx,
        png_dvx.c_str(),
        "#Delta p vs #Delta v_{x};#Delta v_{x} = v_{x}^{rec} - v_{x}^{gen} [cm];#Delta p = p_{rec} - p_{gen} [GeV]",
        filled_events,
        true
    );

    save_2d_png(
        h_dp_vs_dvy,
        png_dvy.c_str(),
        "#Delta p vs #Delta v_{y};#Delta v_{y} = v_{y}^{rec} - v_{y}^{gen} [cm];#Delta p = p_{rec} - p_{gen} [GeV]",
        filled_events,
        true
    );

    save_2d_png(
        h_dp_vs_dvz,
        png_dvz.c_str(),
        "#Delta p vs #Delta v_{z};#Delta v_{z} = v_{z}^{rec} - v_{z}^{gen} [cm];#Delta p = p_{rec} - p_{gen} [GeV]",
        filled_events,
        true
    );

    std::cout << "\nDone.\n"
              << "Scanned events:                         " << scanned_events << "\n"
              << "Filled events:                          " << filled_events << "\n"
              << "Fail exactly 1 trigger electron:        " << n_fail_trigger << "\n"
              << "Fail exactly 1 REC pi+:                 " << n_fail_rec_pip << "\n"
              << "Fail FD status for REC pi+:             " << n_fail_fd << "\n"
              << "Fail theta/p bin for REC pi+:           " << n_fail_theta_p << "\n"
              << "Fail exactly 1 GEN pi+:                 " << n_fail_gen_pip << "\n"
              << "Files missing RUN::config:              " << n_missing_runconfig_bank << "\n"
              << "Saved:                                  " << png_dtheta << "\n"
              << "Saved:                                  " << png_dvx << "\n"
              << "Saved:                                  " << png_dvy << "\n"
              << "Saved:                                  " << png_dvz << "\n";

    delete h_dp_vs_dtheta;
    delete h_dp_vs_dvx;
    delete h_dp_vs_dvy;
    delete h_dp_vs_dvz;
}