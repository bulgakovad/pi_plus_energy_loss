// export_fd_pip_bin_to_csv.cxx
//
// Purpose:
//   Loop over HIPO files listed in a .dat file and export one CSV row per event.
//
// Event selection:
//   - exactly 1 trigger electron: pid = 11 and status < 0
//   - exactly 1 reconstructed pi+: pid = 211
//   - reconstructed pi+ is FD-like: abs(status) < 4000
//   - reconstructed pi+ is inside theta/p bin
//   - exactly 1 generated pi+: pid = 211
//   - everything else is allowed
//
// Output CSV columns:
//   run,event,
//   rec pi+ px,py,pz,
//   gen pi+ px,py,pz,
//   theta_gen,theta_rec,
//   rec pi+ vx,vy,vz,
//   gen pi+ vx,vy,vz,
//   REC::Track detector list for pi+ pindex,
//   REC::Traj detector list for pi+ pindex,
//   REC::Scintillator detector:component list for pi+ pindex,
//   full REC::Particle pid list
//
// Example:
//
// clas12root -l -b -q 'export_fd_pip_bin_to_csv.cxx+("good_hipo.dat","fd_pip_theta38_39_p1p0_1p2.csv", 38,39, 1.0, 1.2, -1 )'

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

#include "TMath.h"

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

static std::string csv_quote(const std::string& s)
{
    std::string out = "\"";

    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out += c;
    }

    out += "\"";
    return out;
}

static std::string join_int_list(const std::vector<int>& values)
{
    std::ostringstream ss;

    for (size_t i = 0; i < values.size(); i++) {
        if (i > 0) ss << ",";
        ss << values[i];
    }

    return ss.str();
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

// For REC::Track and REC::Traj:
// only store the detector list matched by pindex.
static std::string get_detector_list(hipo::bank& bank,
                                     int pindex)
{
    std::ostringstream ss;
    bool first = true;

    for (int i = 0; i < bank.getRows(); i++) {
        const int row_pindex = bank.getInt("pindex", i);

        if (row_pindex != pindex) continue;

        const int detector = bank.getInt("detector", i);

        if (!first) ss << ",";
        ss << detector;
        first = false;
    }

    return ss.str();
}

// For REC::Scintillator:
// store detector:component matched by pindex.
static std::string get_scint_detector_component_list(hipo::bank& bank,
                                                     int pindex)
{
    std::ostringstream ss;
    bool first = true;

    for (int i = 0; i < bank.getRows(); i++) {
        const int row_pindex = bank.getInt("pindex", i);

        if (row_pindex != pindex) continue;

        const int detector = bank.getInt("detector", i);
        const int component = bank.getShort("component", i);

        if (!first) ss << ",";
        ss << detector << ":" << component;
        first = false;
    }

    return ss.str();
}

void export_fd_pip_bin_to_csv(const char* input_dat_file,
                              const char* output_csv,
                              double theta_min_deg = 38.0,
                              double theta_max_deg = 39.0,
                              double p_min = 1.0,
                              double p_max = 1.2,
                              int max_events_to_write = -1)
{
    std::cout << "Macro started\n";
    std::cout << "Input .dat file: " << input_dat_file << "\n";
    std::cout << "Output CSV:      " << output_csv << "\n";
    std::cout << "Selection:\n"
              << "  exactly 1 trigger electron: pid=11, status<0\n"
              << "  exactly 1 REC pi+: pid=211\n"
              << "  REC pi+ FD-like: abs(status)<4000\n"
              << "  REC pi+ theta in [" << theta_min_deg << ", " << theta_max_deg << ") deg\n"
              << "  REC pi+ p     in [" << p_min << ", " << p_max << ") GeV\n"
              << "  exactly 1 GEN pi+: pid=211\n"
              << "  everything else allowed\n"
              << "  max rows = " << max_events_to_write << " (-1 or 0 means all)\n";

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

    std::ofstream fout(output_csv);

    if (!fout.is_open()) {
        std::cerr << "ERROR: cannot open output CSV: " << output_csv << "\n";
        return;
    }

    fout << "run,event,"
         << "pip_rec_px,pip_rec_py,pip_rec_pz,"
         << "pip_lund_px,pip_lund_py,pip_lund_pz,"
         << "theta_lund_deg,theta_rec_deg,"
         << "pip_rec_vx,pip_rec_vy,pip_rec_vz,"
         << "pip_lund_vx,pip_lund_vy,pip_lund_vz,"
         << "rec_track_detector_list,"
         << "rec_traj_detector_list,"
         << "rec_scintillator_detector_component,"
         << "rec_particle_pid_list\n";

    long long scanned_events = 0;
    long long written_rows = 0;

    long long n_fail_trigger = 0;
    long long n_fail_rec_pip = 0;
    long long n_fail_fd = 0;
    long long n_fail_theta_p = 0;
    long long n_fail_gen_pip = 0;

    long long n_missing_rec_track_bank = 0;
    long long n_missing_rec_traj_bank = 0;
    long long n_missing_rec_scint_bank = 0;
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
        hipo::schema rec_track_schema;
        hipo::schema rec_traj_schema;
        hipo::schema rec_scint_schema;

        const bool have_rec_particle =
            get_schema_safe(factory, "REC::Particle", rec_particle_schema);

        if (!have_rec_particle) {
            std::cerr << "WARNING: REC::Particle missing, skipping file:\n"
                      << "  " << input_hipo << "\n";
            continue;
        }

        const bool have_mc_lund =
        get_schema_safe(factory, "MC::Lund", mc_schema);

        if (!have_mc_lund) {
            std::cerr << "WARNING: MC::Lund missing, skipping file:\n"
                      << "  " << input_hipo << "\n";
            continue;
        }

        const bool have_runconfig =
            get_schema_safe(factory, "RUN::config", runconfig_schema);

        const bool have_rec_track =
            get_schema_safe(factory, "REC::Track", rec_track_schema);

        const bool have_rec_traj =
            get_schema_safe(factory, "REC::Traj", rec_traj_schema);

        const bool have_rec_scint =
            get_schema_safe(factory, "REC::Scintillator", rec_scint_schema);

        if (!have_runconfig) n_missing_runconfig_bank++;
        if (!have_rec_track) n_missing_rec_track_bank++;
        if (!have_rec_traj) n_missing_rec_traj_bank++;
        if (!have_rec_scint) n_missing_rec_scint_bank++;

        hipo::event event;

        hipo::bank rec_particle(rec_particle_schema);
        hipo::bank mc_bank(mc_schema);

        hipo::bank runconfig(runconfig_schema);
        hipo::bank rec_track(rec_track_schema);
        hipo::bank rec_traj(rec_traj_schema);
        hipo::bank rec_scint(rec_scint_schema);

        while (reader.next()) {
            reader.read(event);
            scanned_events++;

            if (scanned_events % 50000 == 0) {
                std::cout << "Scanned " << scanned_events
                          << ", written rows " << written_rows << "\n";
            }

            event.getStructure(rec_particle);
            event.getStructure(mc_bank);

            if (have_runconfig) event.getStructure(runconfig);
            if (have_rec_track) event.getStructure(rec_track);
            if (have_rec_traj) event.getStructure(rec_traj);
            if (have_rec_scint) event.getStructure(rec_scint);

            int run = -1;
            int evnum = -1;

            if (have_runconfig && runconfig.getRows() > 0) {
                run = runconfig.getInt("run", 0);
                evnum = runconfig.getInt("event", 0);
            }

            int n_trigger_e = 0;
            int n_rec_pip = 0;
            int pip_rec_row = -1;

            std::vector<int> pid_list;

            for (int i = 0; i < rec_particle.getRows(); i++) {
                const int pid = rec_particle.getInt("pid", i);
                const int status = rec_particle.getInt("status", i);

                pid_list.push_back(pid);

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

            const double pip_rec_p = momentum(pip_rec_px, pip_rec_py, pip_rec_pz);
            const double theta_rec_deg =
                theta_deg_from_p(pip_rec_px, pip_rec_py, pip_rec_pz);

            if (theta_rec_deg < theta_min_deg || theta_rec_deg >= theta_max_deg ||
                pip_rec_p < p_min || pip_rec_p >= p_max) {
                n_fail_theta_p++;
                continue;
            }

            const double pip_rec_vx = rec_particle.getFloat("vx", pip_rec_row);
            const double pip_rec_vy = rec_particle.getFloat("vy", pip_rec_row);
            const double pip_rec_vz = rec_particle.getFloat("vz", pip_rec_row);

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

            const double theta_gen_deg =
                theta_deg_from_p(pip_gen_px, pip_gen_py, pip_gen_pz);

            const double pip_gen_vx = mc_bank.getFloat("vx", pip_gen_row);
            const double pip_gen_vy = mc_bank.getFloat("vy", pip_gen_row);
            const double pip_gen_vz = mc_bank.getFloat("vz", pip_gen_row);

            std::string track_detectors = "";
            std::string traj_detectors = "";
            std::string scint_dc = "";

            if (have_rec_track) {
                track_detectors = get_detector_list(rec_track, pip_rec_row);
            }

            if (have_rec_traj) {
                traj_detectors = get_detector_list(rec_traj, pip_rec_row);
            }

            if (have_rec_scint) {
                scint_dc = get_scint_detector_component_list(rec_scint, pip_rec_row);
            }

            const std::string pid_string = join_int_list(pid_list);

            fout << run << ","
                 << evnum << ","
                 << pip_rec_px << ","
                 << pip_rec_py << ","
                 << pip_rec_pz << ","
                 << pip_gen_px << ","
                 << pip_gen_py << ","
                 << pip_gen_pz << ","
                 << theta_gen_deg << ","
                 << theta_rec_deg << ","
                 << pip_rec_vx << ","
                 << pip_rec_vy << ","
                 << pip_rec_vz << ","
                 << pip_gen_vx << ","
                 << pip_gen_vy << ","
                 << pip_gen_vz << ","
                 << csv_quote(track_detectors) << ","
                 << csv_quote(traj_detectors) << ","
                 << csv_quote(scint_dc) << ","
                 << csv_quote(pid_string) << "\n";

            written_rows++;

            if (max_events_to_write > 0 &&
                written_rows >= max_events_to_write) {
                std::cout << "\nReached max_events_to_write = "
                          << max_events_to_write << "\n";
                break;
            }
        }

        if (max_events_to_write > 0 &&
            written_rows >= max_events_to_write) {
            break;
        }
    }

    fout.close();

    std::cout << "\nDone.\n"
              << "Scanned events:                         " << scanned_events << "\n"
              << "Written CSV rows:                       " << written_rows << "\n"
              << "Fail exactly 1 trigger electron:        " << n_fail_trigger << "\n"
              << "Fail exactly 1 REC pi+:                 " << n_fail_rec_pip << "\n"
              << "Fail FD status for REC pi+:             " << n_fail_fd << "\n"
              << "Fail theta/p bin for REC pi+:           " << n_fail_theta_p << "\n"
              << "Fail exactly 1 GEN pi+:                 " << n_fail_gen_pip << "\n"
              << "Files missing RUN::config:              " << n_missing_runconfig_bank << "\n"
              << "Files missing REC::Track:               " << n_missing_rec_track_bank << "\n"
              << "Files missing REC::Traj:                " << n_missing_rec_traj_bank << "\n"
              << "Files missing REC::Scintillator:        " << n_missing_rec_scint_bank << "\n"
              << "Output CSV:                             " << output_csv << "\n";
}