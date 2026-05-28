// count_piplus_multiplicity_in_hipo.cxx
//
// Purpose:
//   Read one HIPO file and count REC::Particle pi+ and positive-track multiplicities.
//
// Positive-track definition used here:
//   REC::Particle.charge > 0
//
// Additional chi2pid diagnostic:
//   chi2pid is read from REC::Particle, not REC::Track.
//   Counts events with multiple positive tracks and valid chi2pid != 9999.
//
// Examples:
//
// clas12root -l -b -q 'count_piplus_multiplicity_in_hipo.cxx+("FD_theta38_39_p1p0_1p2_step2p1_selectedTrackIsPip.hipo")'
// clas12root -l -b -q 'count_piplus_multiplicity_in_hipo.cxx+("FD_theta20_25_p1p0_1p2_step2p1_selectedTrackIsPip.hipo")'

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include "hipo4/reader.h"
#include "hipo4/event.h"
#include "hipo4/bank.h"
#include "hipo4/dictionary.h"

static void print_count_line(const char* label,
                             long long count,
                             long long total_events)
{
    const double frac = (total_events > 0)
        ? 100.0 * static_cast<double>(count) / static_cast<double>(total_events)
        : 0.0;

    std::cout << std::left << std::setw(70) << label
              << std::right << std::setw(12) << count
              << "   "
              << std::fixed << std::setprecision(3)
              << std::setw(9) << frac << " %\n";
}

void count_piplus_multiplicity_in_hipo(const char* input_hipo)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";
    std::cout << "Positive-track definition: REC::Particle.charge > 0\n";
    std::cout << "Valid PID chi2 definition: REC::Particle.chi2pid != 9999\n";

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

    hipo::schema rec_particle_schema;

    try {
        rec_particle_schema = factory.getSchema("REC::Particle");
    }
    catch (...) {
        std::cerr << "ERROR: REC::Particle bank not found in file:\n"
                  << "  " << input_hipo << "\n";
        return;
    }

    hipo::event event;
    hipo::bank rec_particle(rec_particle_schema);

    long long n_events_scanned = 0;

    long long n_events_zero_pip = 0;
    long long n_events_one_pip = 0;
    long long n_events_one_pip_only_positive = 0;
    long long n_events_one_pip_extra_positive = 0;
    long long n_events_gt1_pip = 0;

    long long n_total_rec_pip_seen = 0;
    long long n_total_positive_rows_seen = 0;

    // New counters
    long long n_events_multiple_positive = 0;

    // This is probably the one you want:
    // event has >1 positive REC::Particle rows and at least one positive row has chi2pid != 9999.
    long long n_events_multiple_positive_any_valid_chi2pid = 0;

    // Stricter diagnostic:
    // event has >1 positive REC::Particle rows and all positive rows have chi2pid != 9999.
    long long n_events_multiple_positive_all_valid_chi2pid = 0;

    while (reader.next()) {
        reader.read(event);
        event.getStructure(rec_particle);

        n_events_scanned++;

        int n_pip_this_event = 0;
        int n_positive_this_event = 0;
        int n_positive_with_valid_chi2pid_this_event = 0;

        for (int i = 0; i < rec_particle.getRows(); i++) {
            const int pid = rec_particle.getInt("pid", i);
            const int charge = rec_particle.getInt("charge", i);

            // chi2pid is in REC::Particle.
            // In CLAS12 cooked files, 9999 is the usual "not assigned / invalid" sentinel.
            const float chi2pid = rec_particle.getFloat("chi2pid", i);

            if (pid == 211) {
                n_pip_this_event++;
            }

            if (charge > 0) {
                n_positive_this_event++;

                if (std::fabs(chi2pid - 9999.0f) > 1.0e-6f) {
                    n_positive_with_valid_chi2pid_this_event++;
                }
            }
        }

        n_total_rec_pip_seen += n_pip_this_event;
        n_total_positive_rows_seen += n_positive_this_event;

        if (n_positive_this_event > 1) {
            n_events_multiple_positive++;

            if (n_positive_with_valid_chi2pid_this_event > 0) {
                n_events_multiple_positive_any_valid_chi2pid++;
            }

            if (n_positive_with_valid_chi2pid_this_event == n_positive_this_event) {
                n_events_multiple_positive_all_valid_chi2pid++;
            }
        }

        if (n_pip_this_event == 0) {
            n_events_zero_pip++;
        }
        else if (n_pip_this_event == 1) {
            n_events_one_pip++;

            if (n_positive_this_event == 1) {
                n_events_one_pip_only_positive++;
            }
            else if (n_positive_this_event > 1) {
                n_events_one_pip_extra_positive++;
            }
            else {
                // This should not happen for normal REC bookkeeping:
                // pid=211 should have positive charge.
                std::cout << "WARNING: event has exactly 1 pid=211 but "
                          << n_positive_this_event
                          << " positive REC::Particle rows. scanned_event="
                          << n_events_scanned << "\n";
            }
        }
        else {
            n_events_gt1_pip++;
        }
    }

    std::cout << "\nEvents scanned:                              "
              << n_events_scanned << "\n";

    std::cout << "\nREC pi+ / positive-track multiplicity per event:\n";


    print_count_line("Events with 1 REC pi+:", n_events_one_pip, n_events_scanned);
    print_count_line("  Out of those: events with 1 REC pi+ and no other positive tracks:",
                     n_events_one_pip_only_positive, n_events_scanned);
    print_count_line("  Out of those: events with 1 REC pi+ but extra positive tracks:",
                     n_events_one_pip_extra_positive, n_events_scanned);
    print_count_line("Events with >1 REC pi+:", n_events_gt1_pip, n_events_scanned);

    std::cout << "\nPositive-track chi2pid diagnostic:\n";

    print_count_line("Events with multiple positive tracks:",
                     n_events_multiple_positive,
                     n_events_scanned);

    print_count_line("Events with multiple positive tracks and ANY positive chi2pid != 9999:",
                     n_events_multiple_positive_any_valid_chi2pid,
                     n_events_scanned);

    print_count_line("Events with multiple positive tracks and ALL positive chi2pid != 9999:",
                     n_events_multiple_positive_all_valid_chi2pid,
                     n_events_scanned);

std::cout << "\nRatios inside multiple-positive-track sample:\n";

const double ratio_any_inside_multi_pos =
    (n_events_multiple_positive > 0)
    ? 100.0 * static_cast<double>(n_events_multiple_positive_any_valid_chi2pid) /
      static_cast<double>(n_events_multiple_positive)
    : 0.0;

std::cout << "Ratio ANY = 100 * "
          << "(events with multiple positive tracks and ANY positive chi2pid != 9999)"
          << " / (events with multiple positive tracks)\n";

std::cout << "Ratio ANY: "
          << n_events_multiple_positive_any_valid_chi2pid
          << " / "
          << n_events_multiple_positive
          << " = "
          << std::fixed << std::setprecision(3)
          << ratio_any_inside_multi_pos
          << " %\n";

const double ratio_all_inside_multi_pos =
    (n_events_multiple_positive > 0)
    ? 100.0 * static_cast<double>(n_events_multiple_positive_all_valid_chi2pid) /
      static_cast<double>(n_events_multiple_positive)
    : 0.0;

std::cout << "\nRatio ALL = 100 * "
          << "(events with multiple positive tracks and ALL positive chi2pid != 9999)"
          << " / (events with multiple positive tracks)\n";

std::cout << "Ratio ALL: "
          << n_events_multiple_positive_all_valid_chi2pid
          << " / "
          << n_events_multiple_positive
          << " = "
          << std::fixed << std::setprecision(3)
          << ratio_all_inside_multi_pos
          << " %\n";

    std::cout << "\nConditional fractions inside exactly-1-REC-pi+ sample:\n";

    if (n_events_one_pip > 0) {
        const double frac_only_positive =
            100.0 * static_cast<double>(n_events_one_pip_only_positive) /
            static_cast<double>(n_events_one_pip);

        const double frac_extra_positive =
            100.0 * static_cast<double>(n_events_one_pip_extra_positive) /
            static_cast<double>(n_events_one_pip);

        std::cout << "Among events with 1 REC pi+:\n";

        std::cout << "  no other positive tracks:                 "
                  << n_events_one_pip_only_positive
                  << "   "
                  << std::fixed << std::setprecision(3)
                  << frac_only_positive << " %\n";

        std::cout << "  extra positive tracks:                    "
                  << n_events_one_pip_extra_positive
                  << "   "
                  << std::fixed << std::setprecision(3)
                  << frac_extra_positive << " %\n";
    }
    else {
        std::cout << "No events with exactly 1 REC pi+.\n";
    }

    const long long sanity_sum =
        n_events_zero_pip +
        n_events_one_pip +
        n_events_gt1_pip;

    std::cout << "\nSanity check:\n";
    std::cout << "0 pi+ + 1 pi+ + >1 pi+ events:              "
              << sanity_sum << "\n";
    std::cout << "Events scanned:                             "
              << n_events_scanned << "\n";

    if (sanity_sum != n_events_scanned) {
        std::cout << "WARNING: pi+ multiplicity-bin sum does not match scanned events.\n";
    }

    if (n_events_one_pip_only_positive + n_events_one_pip_extra_positive != n_events_one_pip) {
        std::cout << "WARNING: exactly-1-pi+ subdivision does not match exactly-1-pi+ count.\n";
    }

    std::cout << "\nTotal row counters:\n";
    std::cout << "Total REC pi+ rows seen:                    "
              << n_total_rec_pip_seen << "\n";
    std::cout << "Total positive REC::Particle rows seen:     "
              << n_total_positive_rows_seen << "\n";
}