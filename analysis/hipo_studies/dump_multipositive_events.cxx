// dump_multipositive_events.cxx
//
// Purpose:
//   Read one Step-2.1 HIPO file.
//   Find first N events with multiple positive REC::Particle tracks.
//   Print full REC::Particle bank and REC::Track bank for those events.
//
// Positive-track definition:
//   REC::Particle.charge > 0
//
// Matching:
//   REC::Track.pindex points to REC::Particle row index.
//
// Example:
//
// clas12root -l -b -q 'dump_multipositive_events.cxx+("FD_theta20_25_p1p0_1p2_step2p1_selectedTrackIsPip.hipo",10,"multipositive_events_dump_20-25_step_2p1.txt")'
// clas12root -l -b -q 'dump_multipositive_events.cxx+("FD_theta38_39_p1p0_1p2_step2p1_selectedTrackIsPip.hipo",10,"multipositive_events_dump_38-39_step_2p1.txt")'


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "hipo4/reader.h"
#include "hipo4/event.h"
#include "hipo4/bank.h"
#include "hipo4/dictionary.h"

static bool has_entry(const hipo::schema& sch, const char* name)
{
    return sch.getEntryOrder(name) >= 0;
}

static double momentum(double px, double py, double pz)
{
    return std::sqrt(px * px + py * py + pz * pz);
}

static double theta_deg(double px, double py, double pz)
{
    const double p = momentum(px, py, pz);

    if (p <= 0.0 || !std::isfinite(p)) {
        return -999.0;
    }

    double c = pz / p;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    return std::acos(c) * 180.0 / M_PI;
}

static void print_rec_particle_bank(std::ostream& out,
                                    hipo::bank& rec_particle,
                                    const hipo::schema& rec_particle_schema)
{
    const bool has_beta    = has_entry(rec_particle_schema, "beta");
    const bool has_chi2pid = has_entry(rec_particle_schema, "chi2pid");

    out << "\nREC::Particle bank\n";
    out << "---------------------------------------------------------------------------------------------\n";

    out << std::setw(4)  << "row"
        << std::setw(8)  << "pid"
        << std::setw(8)  << "charge"
        << std::setw(10) << "status"
        << std::setw(12) << "p"
        << std::setw(12) << "theta"
        << std::setw(12) << "px"
        << std::setw(12) << "py"
        << std::setw(12) << "pz"
        << std::setw(10) << "vx"
        << std::setw(10) << "vy"
        << std::setw(10) << "vz";

    if (has_beta) {
        out << std::setw(10) << "beta";
    }

    if (has_chi2pid) {
        out << std::setw(12) << "chi2pid";
    }

    out << "   note\n";

    for (int i = 0; i < rec_particle.getRows(); i++) {
        const int pid    = rec_particle.getInt("pid", i);
        const int charge = rec_particle.getInt("charge", i);
        const int status = rec_particle.getInt("status", i);

        const double px = rec_particle.getFloat("px", i);
        const double py = rec_particle.getFloat("py", i);
        const double pz = rec_particle.getFloat("pz", i);

        const double p  = momentum(px, py, pz);
        const double th = theta_deg(px, py, pz);

        const double vx = rec_particle.getFloat("vx", i);
        const double vy = rec_particle.getFloat("vy", i);
        const double vz = rec_particle.getFloat("vz", i);

        out << std::setw(4)  << i
            << std::setw(8)  << pid
            << std::setw(8)  << charge
            << std::setw(10) << status
            << std::setw(12) << std::fixed << std::setprecision(5) << p
            << std::setw(12) << std::fixed << std::setprecision(3) << th
            << std::setw(12) << std::fixed << std::setprecision(5) << px
            << std::setw(12) << std::fixed << std::setprecision(5) << py
            << std::setw(12) << std::fixed << std::setprecision(5) << pz
            << std::setw(10) << std::fixed << std::setprecision(3) << vx
            << std::setw(10) << std::fixed << std::setprecision(3) << vy
            << std::setw(10) << std::fixed << std::setprecision(3) << vz;

        if (has_beta) {
            out << std::setw(10)
                << std::fixed << std::setprecision(4)
                << rec_particle.getFloat("beta", i);
        }

        if (has_chi2pid) {
            out << std::setw(12)
                << std::fixed << std::setprecision(3)
                << rec_particle.getFloat("chi2pid", i);
        }

        if (charge > 0) {
            out << "   POSITIVE";
        }

        if (pid == 211) {
            out << " pi+";
        }

        if (pid == 11 && status < 0) {
            out << " trigger_e";
        }

        out << "\n";
    }
}

static void print_rec_track_bank(std::ostream& out,
                                 hipo::bank& rec_track,
                                 hipo::bank& rec_particle,
                                 const hipo::schema& rec_track_schema)
{
    const bool has_detector = has_entry(rec_track_schema, "detector");
    const bool has_sector   = has_entry(rec_track_schema, "sector");
    const bool has_q        = has_entry(rec_track_schema, "q");
    const bool has_chi2     = has_entry(rec_track_schema, "chi2");
    const bool has_ndf      = has_entry(rec_track_schema, "NDF");
    const bool has_status   = has_entry(rec_track_schema, "status");

    out << "\nREC::Track bank\n";
    out << "---------------------------------------------------------------------------------------------\n";

    out << std::setw(4)  << "row"
        << std::setw(10) << "pindex";

    if (has_detector) out << std::setw(10) << "detector";
    if (has_sector)   out << std::setw(10) << "sector";
    if (has_q)        out << std::setw(8)  << "q";
    if (has_chi2)     out << std::setw(12) << "chi2";
    if (has_ndf)      out << std::setw(8)  << "NDF";
    if (has_status)   out << std::setw(10) << "status";

    out << "   matched REC::Particle\n";

    for (int it = 0; it < rec_track.getRows(); it++) {
        const int pindex = rec_track.getInt("pindex", it);

        out << std::setw(4) << it
            << std::setw(10) << pindex;

        if (has_detector) {
            out << std::setw(10)
                << rec_track.getInt("detector", it);
        }

        if (has_sector) {
            out << std::setw(10)
                << rec_track.getInt("sector", it);
        }

        if (has_q) {
            out << std::setw(8)
                << static_cast<int>(rec_track.getByte("q", it));
        }

        if (has_chi2) {
            out << std::setw(12)
                << std::fixed << std::setprecision(3)
                << rec_track.getFloat("chi2", it);
        }

        if (has_ndf) {
            out << std::setw(8)
                << rec_track.getShort("NDF", it);
        }

        if (has_status) {
            out << std::setw(10)
                << rec_track.getShort("status", it);
        }

        if (pindex >= 0 && pindex < rec_particle.getRows()) {
            const int pid    = rec_particle.getInt("pid", pindex);
            const int charge = rec_particle.getInt("charge", pindex);
            const int status = rec_particle.getInt("status", pindex);

            out << "   row=" << pindex
                << " pid=" << pid
                << " charge=" << charge
                << " status=" << status;

            if (charge > 0) {
                out << " POSITIVE";
            }

            if (pid == 211) {
                out << " pi+";
            }

            if (pid == 11 && status < 0) {
                out << " trigger_e";
            }
        }
        else {
            out << "   invalid pindex";
        }

        out << "\n";
    }
}

void dump_multipositive_events(const char* input_hipo,
                               int max_events_to_print = 10,
                               const char* output_txt = "multipositive_events_dump.txt")
{
    std::cout << "Input HIPO:  " << input_hipo << "\n";
    std::cout << "Output text: " << output_txt << "\n";
    std::cout << "Looking for first " << max_events_to_print
              << " events with >1 positive REC::Particle tracks.\n";

    std::ofstream fout(output_txt);

    if (!fout.is_open()) {
        std::cerr << "ERROR: cannot open output text file: "
                  << output_txt << "\n";
        return;
    }

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
        std::cerr << "ERROR: reader is not open.\n";
        return;
    }

    hipo::dictionary factory;
    reader.readDictionary(factory);

    if (!factory.hasSchema("REC::Particle")) {
        std::cerr << "ERROR: REC::Particle bank missing.\n";
        return;
    }

    if (!factory.hasSchema("REC::Track")) {
        std::cerr << "ERROR: REC::Track bank missing.\n";
        return;
    }

    hipo::schema rec_particle_schema = factory.getSchema("REC::Particle");
    hipo::schema rec_track_schema    = factory.getSchema("REC::Track");

    const bool rec_particle_ok =
        has_entry(rec_particle_schema, "pid") &&
        has_entry(rec_particle_schema, "charge") &&
        has_entry(rec_particle_schema, "status") &&
        has_entry(rec_particle_schema, "px") &&
        has_entry(rec_particle_schema, "py") &&
        has_entry(rec_particle_schema, "pz") &&
        has_entry(rec_particle_schema, "vx") &&
        has_entry(rec_particle_schema, "vy") &&
        has_entry(rec_particle_schema, "vz");

    if (!rec_particle_ok) {
        std::cerr << "ERROR: REC::Particle missing required entries.\n";
        return;
    }

    if (!has_entry(rec_track_schema, "pindex")) {
        std::cerr << "ERROR: REC::Track missing pindex entry.\n";
        return;
    }

    hipo::event event;
    hipo::bank rec_particle(rec_particle_schema);
    hipo::bank rec_track(rec_track_schema);

    long long scanned_events = 0;
    long long multipositive_events = 0;
    long long printed_events = 0;

    while (reader.next()) {
        reader.read(event);
        scanned_events++;

        event.getStructure(rec_particle);
        event.getStructure(rec_track);

        int n_positive = 0;
        int n_pip = 0;

        for (int i = 0; i < rec_particle.getRows(); i++) {
            const int pid = rec_particle.getInt("pid", i);
            const int charge = rec_particle.getInt("charge", i);

            if (charge > 0) {
                n_positive++;
            }

            if (pid == 211) {
                n_pip++;
            }
        }

        if (n_positive <= 1) {
            continue;
        }

        multipositive_events++;

        if (printed_events >= max_events_to_print) {
            continue;
        }

        printed_events++;

        fout << "\n\n";
        fout << "================================================================================\n";
        fout << "MULTI-POSITIVE EVENT " << printed_events << "\n";
        fout << "scanned_event_index = " << scanned_events << "\n";
        fout << "n_positive_RECParticle_charge_gt_0 = " << n_positive << "\n";
        fout << "n_REC_pid_211 = " << n_pip << "\n";
        fout << "================================================================================\n";

        print_rec_particle_bank(fout, rec_particle, rec_particle_schema);
        print_rec_track_bank(fout, rec_track, rec_particle, rec_track_schema);
    }

    fout << "\n\n";
    fout << "================================================================================\n";
    fout << "SUMMARY\n";
    fout << "================================================================================\n";
    fout << "Input HIPO:                         " << input_hipo << "\n";
    fout << "Events scanned:                     " << scanned_events << "\n";
    fout << "Events with >1 positive REC tracks: " << multipositive_events << "\n";
    fout << "Events printed:                     " << printed_events << "\n";

    fout.close();

    std::cout << "\nDone.\n";
    std::cout << "Events scanned:                     " << scanned_events << "\n";
    std::cout << "Events with >1 positive REC tracks: " << multipositive_events << "\n";
    std::cout << "Events printed:                     " << printed_events << "\n";
    std::cout << "Output text:                        " << output_txt << "\n";
}