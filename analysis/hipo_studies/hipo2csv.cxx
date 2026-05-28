// hipo2csv.cxx
//
// Purpose:
//   Loop over HIPO files listed in a .dat file and export one CSV row per event.
//
// No selection:
//   - every event is written
//   - no electron requirement
//   - no pion requirement
//   - no FD/CD requirement
//   - no theta/p bin cut
//   - no MC::Lund pi+ requirement
//
// Output:
//   run,event,
//   REC::Particle lists:
//      pid, charge, status,
//      px, py, pz,
//      vx, vy, vz,
//      theta_deg, p
//   MC::Lund lists:
//      pid,
//      px, py, pz,
//      vx, vy, vz,
//      theta_deg, p
//   REC::Track list:
//      pindex:detector
//   REC::Traj list:
//      pindex:detector:layer
//   REC::Scintillator list:
//      pindex:detector:component
//
// Example:
//
// clas12root -l -b -q 'hipo2csv.cxx+("good_hipo.dat","all_events.csv",-1)'
//
// Test first 10 events:
//
// clas12root -l -b -q 'hipo2csv.cxx+("good_hipo.dat","test_10_events.csv",10)'

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

static std::string join_double_list(const std::vector<double>& values)
{
    std::ostringstream ss;
    ss.setf(std::ios::fixed);
    ss.precision(8);

    for (size_t i = 0; i < values.size(); i++) {
        if (i > 0) ss << ",";
        ss << values[i];
    }

    return ss.str();
}

// REC::Track:
// store all rows as pindex:detector
static std::string get_track_list(hipo::bank& bank)
{
    std::ostringstream ss;
    bool first = true;

    for (int i = 0; i < bank.getRows(); i++) {
        const int pindex = bank.getInt("pindex", i);
        const int detector = bank.getInt("detector", i);

        if (!first) ss << ",";
        ss << pindex << ":" << detector;
        first = false;
    }

    return ss.str();
}

// REC::Traj:
// store all rows as pindex:detector:layer
static std::string get_traj_list(hipo::bank& bank)
{
    std::ostringstream ss;
    bool first = true;

    for (int i = 0; i < bank.getRows(); i++) {
        const int pindex = bank.getInt("pindex", i);
        const int detector = bank.getInt("detector", i);
        const int layer = bank.getInt("layer", i);

        if (!first) ss << ",";
        ss << pindex << ":" << detector << ":" << layer;
        first = false;
    }

    return ss.str();
}

// REC::Scintillator:
// store all rows as pindex:detector:component
static std::string get_scintillator_list(hipo::bank& bank)
{
    std::ostringstream ss;
    bool first = true;

    for (int i = 0; i < bank.getRows(); i++) {
        const int pindex = bank.getInt("pindex", i);
        const int detector = bank.getInt("detector", i);
        const int component = bank.getShort("component", i);

        if (!first) ss << ",";
        ss << pindex << ":" << detector << ":" << component;
        first = false;
    }

    return ss.str();
}

void hipo2csv(const char* input_dat_file,
              const char* output_csv,
              int max_events_to_write = -1)
{
    std::cout << "Macro started\n";
    std::cout << "Input .dat file: " << input_dat_file << "\n";
    std::cout << "Output CSV:      " << output_csv << "\n";
    std::cout << "Selection:       none, writing all events\n";
    std::cout << "Max rows:        " << max_events_to_write
              << " (-1 or 0 means all)\n";

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
         << "rec_particle_pid_list,"
         << "rec_particle_charge_list,"
         << "rec_particle_status_list,"
         << "rec_particle_px_list,"
         << "rec_particle_py_list,"
         << "rec_particle_pz_list,"
         << "rec_particle_p_list,"
         << "rec_particle_theta_deg_list,"
         << "rec_particle_vx_list,"
         << "rec_particle_vy_list,"
         << "rec_particle_vz_list,"
         << "mc_lund_pid_list,"
         << "mc_lund_px_list,"
         << "mc_lund_py_list,"
         << "mc_lund_pz_list,"
         << "mc_lund_p_list,"
         << "mc_lund_theta_deg_list,"
         << "mc_lund_vx_list,"
         << "mc_lund_vy_list,"
         << "mc_lund_vz_list,"
         << "rec_track_pindex_detector_list,"
         << "rec_traj_pindex_detector_layer_list,"
         << "rec_scintillator_pindex_detector_component_list\n";

    long long scanned_events = 0;
    long long written_rows = 0;

    long long n_missing_rec_particle_bank_files = 0;
    long long n_missing_mc_lund_bank_files = 0;
    long long n_missing_runconfig_bank_files = 0;
    long long n_missing_rec_track_bank_files = 0;
    long long n_missing_rec_traj_bank_files = 0;
    long long n_missing_rec_scint_bank_files = 0;

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
        hipo::schema mc_lund_schema;
        hipo::schema runconfig_schema;
        hipo::schema rec_track_schema;
        hipo::schema rec_traj_schema;
        hipo::schema rec_scint_schema;

        const bool have_rec_particle =
            get_schema_safe(factory, "REC::Particle", rec_particle_schema);

        const bool have_mc_lund =
            get_schema_safe(factory, "MC::Lund", mc_lund_schema);

        const bool have_runconfig =
            get_schema_safe(factory, "RUN::config", runconfig_schema);

        const bool have_rec_track =
            get_schema_safe(factory, "REC::Track", rec_track_schema);

        const bool have_rec_traj =
            get_schema_safe(factory, "REC::Traj", rec_traj_schema);

        const bool have_rec_scint =
            get_schema_safe(factory, "REC::Scintillator", rec_scint_schema);

        if (!have_rec_particle) n_missing_rec_particle_bank_files++;
        if (!have_mc_lund) n_missing_mc_lund_bank_files++;
        if (!have_runconfig) n_missing_runconfig_bank_files++;
        if (!have_rec_track) n_missing_rec_track_bank_files++;
        if (!have_rec_traj) n_missing_rec_traj_bank_files++;
        if (!have_rec_scint) n_missing_rec_scint_bank_files++;

        hipo::event event;

        hipo::bank rec_particle(rec_particle_schema);
        hipo::bank mc_lund(mc_lund_schema);
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

            if (have_rec_particle) event.getStructure(rec_particle);
            if (have_mc_lund) event.getStructure(mc_lund);
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

            std::vector<int> rec_pid_list;
            std::vector<int> rec_charge_list;
            std::vector<int> rec_status_list;

            std::vector<double> rec_px_list;
            std::vector<double> rec_py_list;
            std::vector<double> rec_pz_list;
            std::vector<double> rec_p_list;
            std::vector<double> rec_theta_list;

            std::vector<double> rec_vx_list;
            std::vector<double> rec_vy_list;
            std::vector<double> rec_vz_list;

            if (have_rec_particle) {
                for (int i = 0; i < rec_particle.getRows(); i++) {
                    const int pid = rec_particle.getInt("pid", i);
                    const int charge = rec_particle.getInt("charge", i);
                    const int status = rec_particle.getInt("status", i);

                    const double px = rec_particle.getFloat("px", i);
                    const double py = rec_particle.getFloat("py", i);
                    const double pz = rec_particle.getFloat("pz", i);

                    const double p = momentum(px, py, pz);
                    const double theta = theta_deg_from_p(px, py, pz);

                    const double vx = rec_particle.getFloat("vx", i);
                    const double vy = rec_particle.getFloat("vy", i);
                    const double vz = rec_particle.getFloat("vz", i);

                    rec_pid_list.push_back(pid);
                    rec_charge_list.push_back(charge);
                    rec_status_list.push_back(status);

                    rec_px_list.push_back(px);
                    rec_py_list.push_back(py);
                    rec_pz_list.push_back(pz);
                    rec_p_list.push_back(p);
                    rec_theta_list.push_back(theta);

                    rec_vx_list.push_back(vx);
                    rec_vy_list.push_back(vy);
                    rec_vz_list.push_back(vz);
                }
            }

            std::vector<int> mc_pid_list;

            std::vector<double> mc_px_list;
            std::vector<double> mc_py_list;
            std::vector<double> mc_pz_list;
            std::vector<double> mc_p_list;
            std::vector<double> mc_theta_list;

            std::vector<double> mc_vx_list;
            std::vector<double> mc_vy_list;
            std::vector<double> mc_vz_list;

            if (have_mc_lund) {
                for (int i = 0; i < mc_lund.getRows(); i++) {
                    const int pid = mc_lund.getInt("pid", i);

                    const double px = mc_lund.getFloat("px", i);
                    const double py = mc_lund.getFloat("py", i);
                    const double pz = mc_lund.getFloat("pz", i);

                    const double p = momentum(px, py, pz);
                    const double theta = theta_deg_from_p(px, py, pz);

                    const double vx = mc_lund.getFloat("vx", i);
                    const double vy = mc_lund.getFloat("vy", i);
                    const double vz = mc_lund.getFloat("vz", i);

                    mc_pid_list.push_back(pid);

                    mc_px_list.push_back(px);
                    mc_py_list.push_back(py);
                    mc_pz_list.push_back(pz);
                    mc_p_list.push_back(p);
                    mc_theta_list.push_back(theta);

                    mc_vx_list.push_back(vx);
                    mc_vy_list.push_back(vy);
                    mc_vz_list.push_back(vz);
                }
            }

            std::string track_list = "";
            std::string traj_list = "";
            std::string scint_list = "";

            if (have_rec_track) {
                track_list = get_track_list(rec_track);
            }

            if (have_rec_traj) {
                traj_list = get_traj_list(rec_traj);
            }

            if (have_rec_scint) {
                scint_list = get_scintillator_list(rec_scint);
            }

            fout << run << ","
                 << evnum << ","

                 << csv_quote(join_int_list(rec_pid_list)) << ","
                 << csv_quote(join_int_list(rec_charge_list)) << ","
                 << csv_quote(join_int_list(rec_status_list)) << ","

                 << csv_quote(join_double_list(rec_px_list)) << ","
                 << csv_quote(join_double_list(rec_py_list)) << ","
                 << csv_quote(join_double_list(rec_pz_list)) << ","
                 << csv_quote(join_double_list(rec_p_list)) << ","
                 << csv_quote(join_double_list(rec_theta_list)) << ","

                 << csv_quote(join_double_list(rec_vx_list)) << ","
                 << csv_quote(join_double_list(rec_vy_list)) << ","
                 << csv_quote(join_double_list(rec_vz_list)) << ","

                 << csv_quote(join_int_list(mc_pid_list)) << ","

                 << csv_quote(join_double_list(mc_px_list)) << ","
                 << csv_quote(join_double_list(mc_py_list)) << ","
                 << csv_quote(join_double_list(mc_pz_list)) << ","
                 << csv_quote(join_double_list(mc_p_list)) << ","
                 << csv_quote(join_double_list(mc_theta_list)) << ","

                 << csv_quote(join_double_list(mc_vx_list)) << ","
                 << csv_quote(join_double_list(mc_vy_list)) << ","
                 << csv_quote(join_double_list(mc_vz_list)) << ","

                 << csv_quote(track_list) << ","
                 << csv_quote(traj_list) << ","
                 << csv_quote(scint_list) << "\n";

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
              << "Files missing REC::Particle:            " << n_missing_rec_particle_bank_files << "\n"
              << "Files missing MC::Lund:                 " << n_missing_mc_lund_bank_files << "\n"
              << "Files missing RUN::config:              " << n_missing_runconfig_bank_files << "\n"
              << "Files missing REC::Track:               " << n_missing_rec_track_bank_files << "\n"
              << "Files missing REC::Traj:                " << n_missing_rec_traj_bank_files << "\n"
              << "Files missing REC::Scintillator:        " << n_missing_rec_scint_bank_files << "\n"
              << "Output CSV:                             " << output_csv << "\n";
}