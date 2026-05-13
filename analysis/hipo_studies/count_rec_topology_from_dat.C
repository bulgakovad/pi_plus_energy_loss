#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>

#include "hipo4/reader.h"
#include "hipo4/dictionary.h"
#include "hipo4/event.h"
#include "hipo4/bank.h"

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

static std::string pid_name(int pid)
{
    if (pid == 11) return "e-";
    if (pid == -11) return "e+";
    if (pid == 22) return "gamma";
    if (pid == 211) return "pi+";
    if (pid == -211) return "pi-";
    if (pid == 111) return "pi0";
    if (pid == 2212) return "proton";
    if (pid == -2212) return "anti-proton";
    if (pid == 2112) return "neutron";
    if (pid == -2112) return "anti-neutron";
    if (pid == 321) return "K+";
    if (pid == -321) return "K-";
    if (pid == 311) return "K0";
    if (pid == 310) return "K0S";
    if (pid == 130) return "K0L";
    if (pid == 45) return "deuteron";
    return "other/unknown";
}

void count_rec_topology_from_dat(const char* input_dat_file,
                                 const char* output_txt = "rec_topology_counts.txt",
                                 long long max_events_to_scan = -1)
{
    std::cout << "Input .dat file: " << input_dat_file << "\n";
    std::cout << "Output report:   " << output_txt << "\n";

    auto hipo_files = get_hipo_files_from_dat(input_dat_file);

    if (hipo_files.empty()) {
        std::cerr << "ERROR: no .hipo files found in .dat file\n";
        return;
    }

    std::cout << "Found " << hipo_files.size() << " HIPO files\n";

    std::ofstream fout(output_txt);
    if (!fout.is_open()) {
        std::cerr << "ERROR: cannot open output report file: " << output_txt << "\n";
        return;
    }

    auto print_both = [&](const std::string& s) {
        std::cout << s;
        fout << s;
    };

    long long n_files_opened = 0;
    long long n_files_failed = 0;

    long long n_events = 0;
    long long n_events_no_rec_particles = 0;

    long long n_events_zero_pip = 0;
    long long n_events_one_pip = 0;
    long long n_events_two_or_more_pip = 0;

    long long n_events_zero_charged_pions = 0;
    long long n_events_one_charged_pion = 0;
    long long n_events_two_or_more_charged_pions = 0;

    long long n_events_exactly_one_pip_and_nothing_else = 0;

    long long n_events_exactly_one_trigger_e = 0;
    long long n_events_trigger_e_and_pip_only = 0;
    long long n_events_trigger_e_pip_plus_extra = 0;

    long long n_events_with_proton = 0;
    long long n_events_with_kaon = 0;
    long long n_events_with_gamma = 0;
    long long n_events_with_unexpected_pid = 0;

    long long n_events_with_extra_pid_not_e_or_pip = 0;

    std::map<int, long long> pid_row_counts;
    std::map<int, long long> pid_event_counts;
    std::map<int, long long> pip_mult_counts;
    std::map<int, long long> charged_pion_mult_counts;
    std::map<int, long long> rec_mult_counts;

    for (size_t ifile = 0; ifile < hipo_files.size(); ifile++) {
        const std::string& file = hipo_files[ifile];

        std::cout << "Opening [" << ifile + 1 << "/" << hipo_files.size()
                  << "]: " << file << "\n";

        hipo::reader reader;

        try {
            reader.open(file.c_str());
        }
        catch (...) {
            std::cerr << "WARNING: failed to open file, skipping: " << file << "\n";
            n_files_failed++;
            continue;
        }

        n_files_opened++;

        hipo::dictionary factory;
        reader.readDictionary(factory);

        hipo::event event;
        hipo::bank recParticle(factory.getSchema("REC::Particle"));

        while (reader.next()) {
            reader.read(event);
            event.getStructure(recParticle);

            n_events++;

            const int n_rec = recParticle.getRows();
            rec_mult_counts[n_rec]++;

            if (n_rec == 0) {
                n_events_no_rec_particles++;
            }

            int n_e_minus = 0;
            int n_trigger_e = 0;
            int n_pip = 0;
            int n_pim = 0;
            int n_proton = 0;
            int n_kaon = 0;
            int n_gamma = 0;
            int n_extra_not_e_or_pip = 0;
            int n_unexpected = 0;

            std::map<int, int> event_pid_seen;

            for (int i = 0; i < n_rec; i++) {
                const int pid = recParticle.getInt("pid", i);
                const int status = recParticle.getShort("status", i);

                pid_row_counts[pid]++;
                event_pid_seen[pid]++;

                if (pid == 11) {
                    n_e_minus++;
                    if (status < 0) n_trigger_e++;
                }

                if (pid == 211) n_pip++;
                if (pid == -211) n_pim++;

                if (pid == 2212) n_proton++;

                if (pid == 321 || pid == -321 || pid == 311 || pid == 310 || pid == 130) {
                    n_kaon++;
                }

                if (pid == 22) n_gamma++;

                // "Extra relative to clean e- + pi+ topology"
                if (!(pid == 11 || pid == 211)) {
                    n_extra_not_e_or_pip++;
                }

                // Broad suspicious/unexpected bucket for your toy e- + pi+ study.
                // This is not proof of mis-ID, just reconstructed PID content
                // beyond the expected electron and pi+.
                if (!(pid == 11 || pid == 211)) {
                    n_unexpected++;
                }
            }

            for (const auto& kv : event_pid_seen) {
                pid_event_counts[kv.first]++;
            }

            pip_mult_counts[n_pip]++;

            const int n_charged_pions = n_pip + n_pim;
            charged_pion_mult_counts[n_charged_pions]++;

            if (n_pip == 0) n_events_zero_pip++;
            else if (n_pip == 1) n_events_one_pip++;
            else n_events_two_or_more_pip++;

            if (n_charged_pions == 0) n_events_zero_charged_pions++;
            else if (n_charged_pions == 1) n_events_one_charged_pion++;
            else n_events_two_or_more_charged_pions++;

            if (n_rec == 1 && n_pip == 1) {
                n_events_exactly_one_pip_and_nothing_else++;
            }

            if (n_trigger_e == 1) {
                n_events_exactly_one_trigger_e++;
            }

            if (n_rec == 2 && n_trigger_e == 1 && n_pip == 1) {
                n_events_trigger_e_and_pip_only++;
            }

            if (n_rec > 2 && n_trigger_e == 1 && n_pip == 1) {
                n_events_trigger_e_pip_plus_extra++;
            }

            if (n_proton > 0) n_events_with_proton++;
            if (n_kaon > 0) n_events_with_kaon++;
            if (n_gamma > 0) n_events_with_gamma++;
            if (n_unexpected > 0) n_events_with_unexpected_pid++;
            if (n_extra_not_e_or_pip > 0) n_events_with_extra_pid_not_e_or_pip++;

            if (max_events_to_scan > 0 && n_events >= max_events_to_scan) {
                break;
            }
        }

        if (max_events_to_scan > 0 && n_events >= max_events_to_scan) {
            break;
        }
    }

    print_both("\n==================== SUMMARY ====================\n");

    {
        std::ostringstream ss;
        ss << "Files listed:                              " << hipo_files.size() << "\n"
           << "Files opened:                              " << n_files_opened << "\n"
           << "Files failed:                              " << n_files_failed << "\n"
           << "Events scanned:                            " << n_events << "\n"
           << "Events with no REC::Particle rows:         " << n_events_no_rec_particles << "\n";
        print_both(ss.str());
    }

    print_both("\n--- pi+ multiplicity ---\n");
    {
        std::ostringstream ss;
        ss << "Events with 0 rec pi+:                     " << n_events_zero_pip << "\n"
           << "Events with exactly 1 rec pi+:             " << n_events_one_pip << "\n"
           << "Events with 2+ rec pi+:                    " << n_events_two_or_more_pip << "\n"
           << "Events with exactly 1 rec pi+ and nothing else: "
           << n_events_exactly_one_pip_and_nothing_else << "\n";
        print_both(ss.str());
    }

    print_both("\n--- charged pion multiplicity, pi+ + pi- ---\n");
    {
        std::ostringstream ss;
        ss << "Events with 0 charged pions:               " << n_events_zero_charged_pions << "\n"
           << "Events with exactly 1 charged pion:        " << n_events_one_charged_pion << "\n"
           << "Events with 2+ charged pions:              " << n_events_two_or_more_charged_pions << "\n";
        print_both(ss.str());
    }

    print_both("\n--- trigger-electron + pi+ topology ---\n");
    {
        std::ostringstream ss;
        ss << "Events with exactly 1 trigger e- status<0: " << n_events_exactly_one_trigger_e << "\n"
           << "Events with trigger e- + pi+ only:         " << n_events_trigger_e_and_pip_only << "\n"
           << "Events with trigger e- + pi+ + extras:     " << n_events_trigger_e_pip_plus_extra << "\n";
        print_both(ss.str());
    }

    print_both("\n--- common extra reconstructed PIDs ---\n");
    {
        std::ostringstream ss;
        ss << "Events with proton pid=2212:               " << n_events_with_proton << "\n"
           << "Events with kaon pid=321/-321/311/310/130: " << n_events_with_kaon << "\n"
           << "Events with gamma pid=22:                  " << n_events_with_gamma << "\n"
           << "Events with any PID not e- or pi+:         " << n_events_with_extra_pid_not_e_or_pip << "\n"
           << "Events with unexpected reconstructed PID:  " << n_events_with_unexpected_pid << "\n";
        print_both(ss.str());
    }

    print_both("\n--- rec pi+ multiplicity table ---\n");
    for (const auto& kv : pip_mult_counts) {
        std::ostringstream ss;
        ss << "n_pi+ = " << std::setw(3) << kv.first
           << "   events = " << kv.second << "\n";
        print_both(ss.str());
    }

    print_both("\n--- charged pion multiplicity table ---\n");
    for (const auto& kv : charged_pion_mult_counts) {
        std::ostringstream ss;
        ss << "n_pi+ + n_pi- = " << std::setw(3) << kv.first
           << "   events = " << kv.second << "\n";
        print_both(ss.str());
    }

    print_both("\n--- REC::Particle multiplicity table ---\n");
    for (const auto& kv : rec_mult_counts) {
        std::ostringstream ss;
        ss << "n_REC_particles = " << std::setw(3) << kv.first
           << "   events = " << kv.second << "\n";
        print_both(ss.str());
    }

    print_both("\n--- REC::Particle PID row counts ---\n");
    for (const auto& kv : pid_row_counts) {
        std::ostringstream ss;
        ss << "pid = " << std::setw(8) << kv.first
           << "   rows = " << std::setw(10) << kv.second
           << "   " << pid_name(kv.first) << "\n";
        print_both(ss.str());
    }

    print_both("\n--- REC::Particle PID event counts ---\n");
    for (const auto& kv : pid_event_counts) {
        std::ostringstream ss;
        ss << "pid = " << std::setw(8) << kv.first
           << "   events with PID = " << std::setw(10) << kv.second
           << "   " << pid_name(kv.first) << "\n";
        print_both(ss.str());
    }

    print_both("\nNOTE:\n");
    print_both("  The script counts reconstructed REC::Particle PID content.\n");
    print_both("  'Unexpected' or 'extra' PID does not prove truth-level misidentification by itself.\n");
    print_both("  It means the reconstruction produced particles beyond the expected e- + pi+ topology.\n");

    fout.close();

    std::cout << "\nSaved report: " << output_txt << "\n";
}