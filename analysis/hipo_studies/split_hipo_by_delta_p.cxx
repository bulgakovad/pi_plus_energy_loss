// split_hipo_by_delta_p.cxx

#include <iostream>
#include <cmath>
#include <fstream>
#include <limits>

#include "hipo4/reader.h"
#include "hipo4/writer.h"
#include "hipo4/event.h"
#include "hipo4/bank.h"
#include "hipo4/dictionary.h"

void add_schema_if_exists(hipo::dictionary& out_dict,
                          hipo::dictionary& in_dict,
                          const char* schema_name)
{
    try {
        out_dict.addSchema(in_dict.getSchema(schema_name));
        std::cout << "Added schema: " << schema_name << "\n";
    }
    catch (...) {
        std::cout << "WARNING: schema not found, skipping: " << schema_name << "\n";
    }
}

void split_hipo_by_delta_p(const char* input_hipo,
                           const char* output_low  = "low_loss_dp_peak.hipo",
                           const char* output_high = "high_loss_dp_negative.hipo",
                           double lowMin  = -0.02,
                           double lowMax  =  0.02,
                           double highMax = -0.04)
{
    std::cout << "Input HIPO:       " << input_hipo  << "\n";
    std::cout << "Low-loss output:  " << output_low  << "\n";
    std::cout << "High-loss output: " << output_high << "\n";
    std::cout << "Low-loss cut:     " << lowMin << " < dp < " << lowMax << "\n";
    std::cout << "High-loss cut:    dp < " << highMax << "\n";

    std::ifstream test(input_hipo);
    if (!test.good()) {
        std::cerr << "\nERROR: Cannot open input file:\n"
                  << "  " << input_hipo << "\n";
        return;
    }

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    std::cout << "\nInput dictionary loaded.\n";

    hipo::schema rec_schema;
    hipo::schema mc_schema;

    try {
        rec_schema = factory.getSchema("REC::Particle");
    }
    catch (...) {
        std::cerr << "ERROR: REC::Particle schema does not exist in this file.\n";
        return;
    }

    try {
        mc_schema = factory.getSchema("MC::Lund");
    }
    catch (...) {
        std::cerr << "ERROR: MC::Lund schema does not exist in this file.\n";
        std::cerr << "Available schemas:\n";
        factory.show();
        return;
    }

    hipo::writer writer_low;
    hipo::writer writer_high;

    // Minimum schemas needed for CED-ish/event viewing and this analysis.
    // Add more here if your CED output needs additional banks.
    add_schema_if_exists(writer_low.getDictionary(), factory, "RUN::config");
    add_schema_if_exists(writer_low.getDictionary(), factory, "REC::Particle");
    add_schema_if_exists(writer_low.getDictionary(), factory, "REC::Track");
    add_schema_if_exists(writer_low.getDictionary(), factory, "REC::Traj");
    add_schema_if_exists(writer_low.getDictionary(), factory, "REC::Scintillator");
    add_schema_if_exists(writer_low.getDictionary(), factory, "REC::Calorimeter");
    add_schema_if_exists(writer_low.getDictionary(), factory, "MC::Lund");

    add_schema_if_exists(writer_high.getDictionary(), factory, "RUN::config");
    add_schema_if_exists(writer_high.getDictionary(), factory, "REC::Particle");
    add_schema_if_exists(writer_high.getDictionary(), factory, "REC::Track");
    add_schema_if_exists(writer_high.getDictionary(), factory, "REC::Traj");
    add_schema_if_exists(writer_high.getDictionary(), factory, "REC::Scintillator");
    add_schema_if_exists(writer_high.getDictionary(), factory, "REC::Calorimeter");
    add_schema_if_exists(writer_high.getDictionary(), factory, "MC::Lund");

    writer_low.open(output_low);
    writer_high.open(output_high);

    hipo::event event;

    hipo::bank recPart(rec_schema);
    hipo::bank mcLund(mc_schema);

    long n_events_total = 0;
    long n_low = 0;
    long n_high = 0;
    long n_bad_rec_pip = 0;
    long n_bad_gen_pip = 0;
    long n_neither = 0;

    while (reader.next()) {
        reader.read(event);
        n_events_total++;

        event.getStructure(recPart);
        event.getStructure(mcLund);

        int n_rec_pip = 0;
        double p_rec = std::numeric_limits<double>::quiet_NaN();

        for (int i = 0; i < recPart.getRows(); i++) {
            const int pid = recPart.getInt("pid", i);

            if (pid != 211) continue;

            const double px = recPart.getFloat("px", i);
            const double py = recPart.getFloat("py", i);
            const double pz = recPart.getFloat("pz", i);

            p_rec = std::sqrt(px * px + py * py + pz * pz);
            n_rec_pip++;
        }

        if (n_rec_pip != 1 || !std::isfinite(p_rec)) {
            n_bad_rec_pip++;
            continue;
        }

        int n_gen_pip = 0;
        double p_gen = std::numeric_limits<double>::quiet_NaN();

        for (int i = 0; i < mcLund.getRows(); i++) {
            const int pid = mcLund.getInt("pid", i);

            if (pid != 211) continue;

            const double px = mcLund.getFloat("px", i);
            const double py = mcLund.getFloat("py", i);
            const double pz = mcLund.getFloat("pz", i);

            p_gen = std::sqrt(px * px + py * py + pz * pz);
            n_gen_pip++;
        }

        if (n_gen_pip != 1 || !std::isfinite(p_gen)) {
            n_bad_gen_pip++;
            continue;
        }

        const double dp = p_rec - p_gen;

        if (dp > lowMin && dp < lowMax) {
            writer_low.addEvent(event);
            n_low++;
        }
        else if (dp < highMax) {
            writer_high.addEvent(event);
            n_high++;
        }
        else {
            n_neither++;
        }
    }

    writer_low.close();
    writer_high.close();

    std::cout << "\nDone.\n";
    std::cout << "Total events scanned:        " << n_events_total << "\n";
    std::cout << "Written to low-loss file:    " << n_low << "\n";
    std::cout << "Written to high-loss file:   " << n_high << "\n";
    std::cout << "Skipped: REC pi+ != 1:       " << n_bad_rec_pip << "\n";
    std::cout << "Skipped: GEN/LUND pi+ != 1:  " << n_bad_gen_pip << "\n";
    std::cout << "Skipped: neither dp region:  " << n_neither << "\n";

    std::cout << "\nOutput files:\n";
    std::cout << "  " << output_low  << "\n";
    std::cout << "  " << output_high << "\n";
}