// skim_five_step_bin_quality_study_singlepass.cxx
//
// Purpose:
//   Scan HIPO files listed in one .dat file exactly once.
//   Apply the same base selection and five-step logic as
//   skim_three_step_positive_track_study_singlepass.cxx.
//
//   Route each accepted Step-1 selected track into exactly one reconstructed
//   theta/p bin.
//
// Migration modes:
//
//   allow_bin_migration = true
//       The selected reconstructed track determines the theta/p bin.
//       The generated MC::Lund pi+ may originate from any theta/p bin.
//
//   allow_bin_migration = false
//       Require the generated MC::Lund pi+ to fall into the same theta/p
//       bin as the selected reconstructed track.
//
// Automatic output folders:
//
//   FD, migration allowed:
//       bin_quality_studies_FD_YES_migration
//
//   FD, migration rejected:
//       bin_quality_studies_FD_NO_migration
//
//   CD, migration allowed:
//       bin_quality_studies_CD_YES_migration
//
//   CD, migration rejected:
//       bin_quality_studies_CD_NO_migration
//
// For every kinematic bin, write:
//
//   step1_dp.png
//   step2p1_dp.png
//   step2p2_dp.png
//   step3p1_dp.png
//   step3p2_dp.png
//   step1_positive_track_multiplicity.png
//   bin_stats.txt
//
// No HIPO output files are written.
//
// Base selection:
//   - exactly 1 reconstructed electron total: pid == 11
//   - that electron must be trigger electron: status < 0
//   - no negative hadrons:
//       no REC::Particle rows with charge < 0 except pid == 11
//   - exactly 1 generated pi+ from MC::Lund
//
// Step 1:
//   - require at least one positive track in requested detector
//   - among requested-detector positive tracks, choose the one closest
//     to MC::Lund pi+ momentum
//   - assign selected reconstructed track into exactly one theta/p grid bin
//   - if allow_bin_migration == false, require generated pi+ to be in
//     the same theta/p bin
//
// Step 2.1:
//   - selected Step-1 track is identified as pi+: pid == 211
//
// Step 2.2:
//   - exactly one positive REC::Particle track total
//
// Step 3.1:
//   - Step-2.1 selected pion is the only positive track
//
// Step 3.2:
//   - Step-2.2 only positive track is identified as pi+
//
// Examples:
//
// FD, allow migration:
//
// clas12root -l -b -q 'skim_five_step_bin_quality_study_singlepass.cxx+("good_hipo.dat",true,"FD",30,42,1.0,0.4,4.8,0.2,-0.5,0.5,200,-1)'
//
// FD, reject migration:
//
// clas12root -l -b -q 'skim_five_step_bin_quality_study_singlepass.cxx+("good_hipo.dat",false,"FD",30,42,1.0,0.4,4.8,0.2,-0.5,0.5,200,-1)'
//
// CD, allow migration:
//
// clas12root -l -b -q 'skim_five_step_bin_quality_study_singlepass.cxx+("good_hipo.dat",true,"CD",30,42,1.0,0.4,4.8,0.2,-0.5,0.5,200,-1)'
//
// CD, reject migration:
//
// clas12root -l -b -q 'skim_five_step_bin_quality_study_singlepass.cxx+("good_hipo.dat",false,"CD",30,42,1.0,0.4,4.8,0.2,-0.5,0.5,200,-1)'

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>
#include <algorithm>

#include "TMath.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TAxis.h"

#include "hipo4/reader.h"
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


struct BinData {
    int i_theta = -1;
    int i_p = -1;

    double theta_min = 0.0;
    double theta_max = 0.0;
    double p_min = 0.0;
    double p_max = 0.0;

    std::string folder;

    TH1D* h_step1 = nullptr;
    TH1D* h_21 = nullptr;
    TH1D* h_22 = nullptr;
    TH1D* h_31 = nullptr;
    TH1D* h_32 = nullptr;

    // Positive-track multiplicity for events passing Step 1.
    TH1D* h_step1_positive_multiplicity = nullptr;

    // Five-step counts.
    long long n_step1 = 0;
    long long n_21 = 0;
    long long n_22 = 0;
    long long n_31 = 0;
    long long n_32 = 0;

    // Number of events reconstructed into this bin but rejected because
    // the generated pi+ belongs to a different theta/p bin.
    // This remains zero when allow_bin_migration == true.
    long long n_rejected_bin_migration = 0;

    // Step-2.1 REC pi+ / positive-track multiplicity diagnostic.
    long long n_21_zero_rec_pip = 0;
    long long n_21_one_rec_pip = 0;
    long long n_21_one_rec_pip_only_positive = 0;
    long long n_21_one_rec_pip_extra_positive = 0;
    long long n_21_gt1_rec_pip = 0;

    long long n_21_total_rec_pip_rows = 0;
    long long n_21_total_positive_rows = 0;
};


static bool ends_with(const std::string& s,
                      const std::string& suffix)
{
    return s.size() >= suffix.size() &&
           s.compare(
               s.size() - suffix.size(),
               suffix.size(),
               suffix
           ) == 0;
}


static bool file_exists(const std::string& path)
{
    std::ifstream f(path);
    return f.good();
}


static std::vector<std::string>
get_hipo_files_from_dat(const std::string& dat_file)
{
    std::vector<std::string> hipo_files;

    std::ifstream fin(dat_file);

    if (!fin.is_open()) {
        std::cerr << "ERROR: cannot open .dat file: "
                  << dat_file << "\n";

        return hipo_files;
    }

    std::string line;

    while (std::getline(fin, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (line.empty()) continue;
        if (line[0] == '#') continue;
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


static bool detector_status_pass(int status_raw,
                                 const std::string& det)
{
    const int status = std::abs(status_raw);

    if (det == "FD") return status < 4000;
    if (det == "CD") return status >= 4000;

    return false;
}


static double momentum(double px,
                       double py,
                       double pz)
{
    return std::sqrt(
        px * px +
        py * py +
        pz * pz
    );
}


static double theta_deg_from_p(double px,
                               double py,
                               double pz)
{
    const double p =
        momentum(px, py, pz);

    if (p <= 0.0 || !std::isfinite(p)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double c = pz / p;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    return std::acos(c) *
           TMath::RadToDeg();
}


// Extract both generated momentum and generated theta.
// Exactly one MC::Lund pi+ is required.
static bool get_mclund_pip_kinematics(
    hipo::bank& mc_lund,
    double& p_lund,
    double& theta_lund_deg,
    int& n_lund_pip)
{
    p_lund =
        std::numeric_limits<double>::quiet_NaN();

    theta_lund_deg =
        std::numeric_limits<double>::quiet_NaN();

    n_lund_pip = 0;

    for (int i = 0; i < mc_lund.getRows(); i++) {
        const int pid =
            mc_lund.getInt("pid", i);

        if (pid != 211) continue;

        const double px =
            mc_lund.getFloat("px", i);

        const double py =
            mc_lund.getFloat("py", i);

        const double pz =
            mc_lund.getFloat("pz", i);

        p_lund =
            momentum(px, py, pz);

        theta_lund_deg =
            theta_deg_from_p(px, py, pz);

        n_lund_pip++;
    }

    return n_lund_pip == 1 &&
           std::isfinite(p_lund) &&
           std::isfinite(theta_lund_deg);
}


static double pct(long long numerator,
                  long long denominator)
{
    if (denominator <= 0) {
        return 0.0;
    }

    return 100.0 *
           static_cast<double>(numerator) /
           static_cast<double>(denominator);
}


static std::string format_number(double x,
                                 int precision = 6)
{
    std::ostringstream ss;

    ss << std::fixed
       << std::setprecision(precision)
       << x;

    std::string out =
        ss.str();

    while (!out.empty() &&
           out.back() == '0') {

        out.pop_back();
    }

    if (!out.empty() &&
        out.back() == '.') {

        out.pop_back();
    }

    if (out == "-0") {
        out = "0";
    }

    return out;
}


static std::string join_path(const std::string& a,
                             const std::string& b)
{
    if (a.empty()) {
        return b;
    }

    if (a.back() == '/') {
        return a + b;
    }

    return a + "/" + b;
}


static bool create_folder(const std::string& folder)
{
    const int rc =
        gSystem->mkdir(
            folder.c_str(),
            true
        );

    if (rc != 0 &&
        gSystem->AccessPathName(folder.c_str())) {

        std::cerr << "ERROR: could not create folder:\n"
                  << "  " << folder << "\n";

        return false;
    }

    return true;
}


static int get_grid_bin(double value,
                        double grid_min,
                        double grid_max,
                        double step,
                        int n_bins)
{
    if (!std::isfinite(value)) {
        return -1;
    }

    if (value < grid_min ||
        value >= grid_max) {

        return -1;
    }

    const int bin =
        static_cast<int>(
            std::floor(
                (value - grid_min) /
                step
            )
        );

    if (bin < 0 ||
        bin >= n_bins) {

        return -1;
    }

    return bin;
}


static std::string make_plot_title(
    const std::string& step_title,
    const BinData& bin,
    const std::string& detector,
    bool allow_bin_migration,
    const std::string& x_axis,
    const std::string& y_axis)
{
    const std::string migration_label =
        allow_bin_migration
        ? "migration allowed"
        : "no migration";

    std::ostringstream ss;

    ss << step_title
       << ", "
       << detector
       << ", #theta_{rec} #in ["
       << format_number(bin.theta_min)
       << ", "
       << format_number(bin.theta_max)
       << ") deg, p_{rec} #in ["
       << format_number(bin.p_min)
       << ", "
       << format_number(bin.p_max)
       << ") GeV, "
       << migration_label
       << ";"
       << x_axis
       << ";"
       << y_axis;

    return ss.str();
}


static void save_dp_hist_png(
    TH1D* h,
    const std::string& output_png,
    const std::string& title,
    long long n_entries)
{
    gStyle->SetOptStat(0);

    TCanvas* c =
        new TCanvas(
            Form("c_%s", h->GetName()),
            title.c_str(),
            1100,
            750
        );

    c->SetTopMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    h->SetLineWidth(2);
    h->SetTitle(title.c_str());
    h->Draw("hist");

    const double ymax =
        std::max(
            1.0,
            h->GetMaximum() * 1.05
        );

    TLine* zero =
        new TLine(
            0.0,
            0.0,
            0.0,
            ymax
        );

    zero->SetLineStyle(2);
    zero->SetLineWidth(2);
    zero->Draw("same");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.033);

    lat.DrawLatex(
        0.15,
        0.88,
        Form("Events: %lld", n_entries)
    );

    c->SaveAs(output_png.c_str());

    delete zero;
    delete c;
}


static void save_multiplicity_hist_png(
    TH1D* h,
    const std::string& output_png,
    const std::string& title,
    long long n_entries)
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.05);

    TCanvas* c =
        new TCanvas(
            Form("c_%s", h->GetName()),
            title.c_str(),
            1100,
            750
        );

    c->SetTopMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    h->SetLineWidth(2);
    h->SetTitle(title.c_str());

    h->GetXaxis()->SetNdivisions(4, false);
    h->GetXaxis()->SetDecimals(false);
    h->GetXaxis()->CenterLabels(false);

    h->GetXaxis()->SetBinLabel(1, "1");
    h->GetXaxis()->SetBinLabel(2, "2");
    h->GetXaxis()->SetBinLabel(3, "3");
    h->GetXaxis()->SetBinLabel(4, "4");

    h->SetMinimum(0.0);
    h->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.033);

    lat.DrawLatex(
        0.15,
        0.88,
        Form("Step-1 events: %lld", n_entries)
    );

    c->SaveAs(output_png.c_str());

    delete c;
}


static void update_step21_multiplicity(
    BinData& bin,
    int n_rec_pip_this_event,
    int n_positive_this_event)
{
    bin.n_21_total_rec_pip_rows +=
        n_rec_pip_this_event;

    bin.n_21_total_positive_rows +=
        n_positive_this_event;

    if (n_rec_pip_this_event == 0) {
        bin.n_21_zero_rec_pip++;
    }
    else if (n_rec_pip_this_event == 1) {
        bin.n_21_one_rec_pip++;

        if (n_positive_this_event == 1) {
            bin.n_21_one_rec_pip_only_positive++;
        }
        else if (n_positive_this_event > 1) {
            bin.n_21_one_rec_pip_extra_positive++;
        }
    }
    else {
        bin.n_21_gt1_rec_pip++;
    }
}


static void write_bin_stats(
    const BinData& bin,
    const std::string& input_dat_file,
    const std::string& detector,
    bool allow_bin_migration,
    long long global_events_scanned,
    long long global_base_electron,
    long long global_no_negative_hadrons,
    long long global_one_lund_pip,
    long long global_rejected_bin_migration)
{
    const std::string stats_path =
        join_path(
            bin.folder,
            "bin_stats.txt"
        );

    std::ofstream out(stats_path);

    if (!out.is_open()) {
        std::cerr << "ERROR: cannot write stats file:\n"
                  << "  " << stats_path << "\n";

        return;
    }

    out << "Five-step bin-quality study\n";
    out << "===========================\n\n";

    out << "Input .dat file: " << input_dat_file << "\n";
    out << "Detector:        " << detector << "\n";

    out << "Allow bin migration: "
        << (allow_bin_migration ? "true" : "false")
        << "\n";

    if (allow_bin_migration) {
        out << "Migration mode: reconstructed bin only; "
            << "generated pi+ may come from another bin\n";
    }
    else {
        out << "Migration mode: require generated and reconstructed pi+ "
            << "to be in the same theta/p bin\n";
    }

    out << "Theta bin:       ["
        << format_number(bin.theta_min)
        << ", "
        << format_number(bin.theta_max)
        << ") deg\n";

    out << "Momentum bin:    ["
        << format_number(bin.p_min)
        << ", "
        << format_number(bin.p_max)
        << ") GeV\n\n";

    out << "Global base-selection counters\n";
    out << "------------------------------\n";

    out << "Events scanned:                                "
        << global_events_scanned << "\n";

    out << "Exactly 1 REC electron total and trigger:      "
        << global_base_electron << "\n";

    out << "No negative hadrons:                           "
        << global_no_negative_hadrons << "\n";

    out << "Exactly 1 MC::Lund pi+:                        "
        << global_one_lund_pip << "\n";

    out << "Rejected reconstructed/generated bin mismatch:"
        << std::setw(12)
        << global_rejected_bin_migration
        << "\n\n";

    out << "Five-step counts for this theta/p bin\n";
    out << "------------------------------------\n";

    out << std::left
        << std::setw(60)
        << "Rejected migration into this reconstructed bin:"
        << std::right
        << std::setw(12)
        << bin.n_rejected_bin_migration
        << "\n";

    out << std::left
        << std::setw(60)
        << "Step 1: selected detector-positive track in bin:"
        << std::right
        << std::setw(12)
        << bin.n_step1
        << "\n";

    out << std::left
        << std::setw(60)
        << "Step 2.1: selected Step-1 track is pi+:"
        << std::right
        << std::setw(12)
        << bin.n_21
        << "   "
        << std::fixed
        << std::setprecision(3)
        << pct(bin.n_21, bin.n_step1)
        << " % of Step 1\n";

    out << std::left
        << std::setw(60)
        << "Step 2.2: exactly one positive track total:"
        << std::right
        << std::setw(12)
        << bin.n_22
        << "   "
        << std::fixed
        << std::setprecision(3)
        << pct(bin.n_22, bin.n_step1)
        << " % of Step 1\n";

    out << std::left
        << std::setw(60)
        << "Step 3.1: Step-2.1 pion is only positive track:"
        << std::right
        << std::setw(12)
        << bin.n_31
        << "   "
        << std::fixed
        << std::setprecision(3)
        << pct(bin.n_31, bin.n_21)
        << " % of Step 2.1\n";

    out << std::left
        << std::setw(60)
        << "Step 3.2: Step-2.2 only positive track is pi+:"
        << std::right
        << std::setw(12)
        << bin.n_32
        << "   "
        << std::fixed
        << std::setprecision(3)
        << pct(bin.n_32, bin.n_22)
        << " % of Step 2.2\n";

    out << "\nStep-1 positive-track multiplicity\n";
    out << "----------------------------------\n";

    out << "Step-1 events: "
        << bin.n_step1
        << "\n";

    for (int i = 1;
         i <= bin.h_step1_positive_multiplicity->GetNbinsX();
         i++) {

        const long long count =
            static_cast<long long>(
                std::llround(bin.h_step1_positive_multiplicity->GetBinContent(i)
                )
            );

        if (count == 0) {
            continue;
        }

        const int multiplicity =
            static_cast<int>(
                std::llround(
                    bin.h_step1_positive_multiplicity
                       ->GetXaxis()
                       ->GetBinCenter(i)
                )
            );

        out << "Positive tracks = "
            << std::setw(2)
            << multiplicity
            << ": "
            << std::setw(12)
            << count
            << "   "
            << std::fixed
            << std::setprecision(3)
            << pct(count, bin.n_step1)
            << " % of Step 1\n";
    }

    const int multiplicity_overflow_bin =
        bin.h_step1_positive_multiplicity
           ->GetNbinsX() + 1;

    const long long multiplicity_overflow =
        static_cast<long long>(
            std::llround(
                bin.h_step1_positive_multiplicity
                   ->GetBinContent(
                       multiplicity_overflow_bin
                   )
            )
        );

    if (multiplicity_overflow > 0) {
        out << "Positive tracks > 4"
            << ": "
            << multiplicity_overflow
            << "   "
            << std::fixed
            << std::setprecision(3)
            << pct(
                   multiplicity_overflow,
                   bin.n_step1
               )
            << " % of Step 1\n";
    }

    out << "\nStep-2.1 REC pi+ / positive-track multiplicity per event\n";
    out << "-------------------------------------------------------\n";

    out << std::left
        << std::setw(72)
        << "Step-2.1 events:"
        << std::right
        << std::setw(12)
        << bin.n_21
        << "\n";

    out << std::left
        << std::setw(72)
        << "Events with 0 REC pi+:"
        << std::right
        << std::setw(12)
        << bin.n_21_zero_rec_pip
        << "   "
        << std::fixed
        << std::setprecision(3)
        << pct(
               bin.n_21_zero_rec_pip,
               bin.n_21
           )
        << " %\n";

    out << std::left
        << std::setw(72)
        << "Events with 1 REC pi+:"
        << std::right
        << std::setw(12)
        << bin.n_21_one_rec_pip
        << "   "
        << std::fixed
        << std::setprecision(3)
        << pct(
               bin.n_21_one_rec_pip,
               bin.n_21
           )
        << " %\n";

    out << std::left
        << std::setw(72)
        << "  Out of those: 1 REC pi+ and no other positive tracks:"
        << std::right
        << std::setw(12)
        << bin.n_21_one_rec_pip_only_positive
        << "   "
        << std::fixed
        << std::setprecision(3)
        << pct(
               bin.n_21_one_rec_pip_only_positive,
               bin.n_21
           )
        << " %\n";

    out << std::left
        << std::setw(72)
        << "  Out of those: 1 REC pi+ but extra positive tracks:"
        << std::right
        << std::setw(12)
        << bin.n_21_one_rec_pip_extra_positive
        << "   "
        << std::fixed
        << std::setprecision(3)
        << pct(
               bin.n_21_one_rec_pip_extra_positive,
               bin.n_21
           )
        << " %\n";

    out << std::left
        << std::setw(72)
        << "Events with >1 REC pi+:"
        << std::right
        << std::setw(12)
        << bin.n_21_gt1_rec_pip
        << "   "
        << std::fixed
        << std::setprecision(3)
        << pct(
               bin.n_21_gt1_rec_pip,
               bin.n_21
           )
        << " %\n";

    out << "\nConditional fractions inside exactly-1-REC-pi+ sample\n";
    out << "----------------------------------------------------\n";

    out << "No other positive tracks: "
        << bin.n_21_one_rec_pip_only_positive
        << " / "
        << bin.n_21_one_rec_pip
        << " = "
        << std::fixed
        << std::setprecision(3)
        << pct(
               bin.n_21_one_rec_pip_only_positive,
               bin.n_21_one_rec_pip
           )
        << " %\n";

    out << "Extra positive tracks:    "
        << bin.n_21_one_rec_pip_extra_positive
        << " / "
        << bin.n_21_one_rec_pip
        << " = "
        << std::fixed
        << std::setprecision(3)
        << pct(
               bin.n_21_one_rec_pip_extra_positive,
               bin.n_21_one_rec_pip
           )
        << " %\n";

    out << "\nHistogram statistics\n";
    out << "--------------------\n";

    out << "Step 1 delta-p entries:   "
        << bin.h_step1->GetEntries()
        << "   mean="
        << bin.h_step1->GetMean()
        << "   stddev="
        << bin.h_step1->GetStdDev()
        << "\n";

    out << "Step 2.1 delta-p entries: "
        << bin.h_21->GetEntries()
        << "   mean="
        << bin.h_21->GetMean()
        << "   stddev="
        << bin.h_21->GetStdDev()
        << "\n";

    out << "Step 2.2 delta-p entries: "
        << bin.h_22->GetEntries()
        << "   mean="
        << bin.h_22->GetMean()
        << "   stddev="
        << bin.h_22->GetStdDev()
        << "\n";

    out << "Step 3.1 delta-p entries: "
        << bin.h_31->GetEntries()
        << "   mean="
        << bin.h_31->GetMean()
        << "   stddev="
        << bin.h_31->GetStdDev()
        << "\n";

    out << "Step 3.2 delta-p entries: "
        << bin.h_32->GetEntries()
        << "   mean="
        << bin.h_32->GetMean()
        << "   stddev="
        << bin.h_32->GetStdDev()
        << "\n";

    out << "Step 1 positive-multiplicity entries: "
        << bin.h_step1_positive_multiplicity->GetEntries()
        << "   mean="
        << bin.h_step1_positive_multiplicity->GetMean()
        << "   stddev="
        << bin.h_step1_positive_multiplicity->GetStdDev()
        << "\n";

    out.close();
}


void skim_five_step_bin_quality_study_singlepass(
    const char* input_dat_file,
    bool allow_bin_migration = true,
    const char* detector = "FD",
    double theta_grid_min = 30.0,
    double theta_grid_max = 42.0,
    double theta_step = 1.0,
    double p_grid_min = 0.4,
    double p_grid_max = 4.8,
    double p_step = 0.2,
    double dp_min = -0.5,
    double dp_max = 0.5,
    int n_dp_bins = 200,
    long long max_events_to_scan = -1)
{
    const std::string det =
        detector;

    if (det != "FD" &&
        det != "CD") {

        std::cerr << "ERROR: detector must be \"FD\" or \"CD\", got: "
                  << det << "\n";

        return;
    }

    const std::string output_root =
        std::string("bin_quality_studies_") +
        det +
        (
            allow_bin_migration
            ? "_YES_migration"
            : "_NO_migration"
        );

    if (!(theta_grid_max > theta_grid_min) ||
        !(p_grid_max > p_grid_min) ||
        !(theta_step > 0.0) ||
        !(p_step > 0.0)) {

        std::cerr << "ERROR: invalid theta/p grid.\n";
        return;
    }

    const int n_theta_bins =
        static_cast<int>(
            std::llround(
                (theta_grid_max -
                 theta_grid_min) /
                theta_step
            )
        );

    const int n_p_bins =
        static_cast<int>(
            std::llround(
                (p_grid_max -
                 p_grid_min) /
                p_step
            )
        );

    const double theta_span_rebuilt =
        theta_grid_min +
        n_theta_bins * theta_step;

    const double p_span_rebuilt =
        p_grid_min +
        n_p_bins * p_step;

    if (std::fabs(
            theta_span_rebuilt -
            theta_grid_max
        ) > 1.0e-9 ||
        std::fabs(
            p_span_rebuilt -
            p_grid_max
        ) > 1.0e-9) {

        std::cerr
            << "ERROR: grid range is not divisible by grid step.\n";

        return;
    }

    std::cout
        << "Macro started: FIVE-STEP MULTI-BIN QUALITY STUDY\n";

    std::cout << "Input .dat file:       "
              << input_dat_file << "\n";

    std::cout << "Allow bin migration:   "
              << (
                     allow_bin_migration
                     ? "true"
                     : "false"
                 )
              << "\n";

    std::cout << "Output root folder:    "
              << output_root << "\n";

    std::cout << "Detector:              "
              << det << "\n";

    std::cout << "Theta grid:            ["
              << theta_grid_min << ", "
              << theta_grid_max
              << ") deg, step "
              << theta_step << "\n";

    std::cout << "Momentum grid:         ["
              << p_grid_min << ", "
              << p_grid_max
              << ") GeV, step "
              << p_step << "\n";

    std::cout << "Theta bins:            "
              << n_theta_bins << "\n";

    std::cout << "Momentum bins:         "
              << n_p_bins << "\n";

    std::cout << "Total theta/p bins:    "
              << n_theta_bins * n_p_bins
              << "\n";

    std::cout << "Delta-p histogram:     "
              << n_dp_bins
              << " bins in ["
              << dp_min << ", "
              << dp_max << "] GeV\n";

    std::cout << "Max events to scan:    "
              << max_events_to_scan
              << " (-1 or 0 means all)\n\n";

    if (allow_bin_migration) {
        std::cout
            << "Migration behavior: generated pi+ may come from any bin.\n\n";
    }
    else {
        std::cout
            << "Migration behavior: generated pi+ must be in the same "
            << "theta/p bin as the selected reconstructed track.\n\n";
    }

    auto hipo_files_raw =
        get_hipo_files_from_dat(
            input_dat_file
        );

    if (hipo_files_raw.empty()) {
        std::cerr
            << "ERROR: no .hipo files found in "
            << input_dat_file << "\n";

        return;
    }

    std::vector<std::string> hipo_files;

    for (const auto& f : hipo_files_raw) {
        if (!file_exists(f)) {
            std::cerr
                << "WARNING: skipping unreadable file:\n"
                << "  " << f << "\n";

            continue;
        }

        hipo_files.push_back(f);
    }

    if (hipo_files.empty()) {
        std::cerr
            << "ERROR: no readable .hipo files left.\n";

        return;
    }

    if (!create_folder(output_root)) {
        return;
    }

    std::cout << "Readable HIPO files: "
              << hipo_files.size()
              << "\n";

    std::vector<BinData> bins;

    bins.reserve(
        n_theta_bins *
        n_p_bins
    );

    for (int i_theta = 0;
         i_theta < n_theta_bins;
         i_theta++) {

        const double theta_min =
            theta_grid_min +
            i_theta * theta_step;

        const double theta_max =
            theta_min +
            theta_step;

        for (int i_p = 0;
             i_p < n_p_bins;
             i_p++) {

            const double p_min =
                p_grid_min +
                i_p * p_step;

            const double p_max =
                p_min +
                p_step;

            BinData bin;

            bin.i_theta =
                i_theta;

            bin.i_p =
                i_p;

            bin.theta_min =
                theta_min;

            bin.theta_max =
                theta_max;

            bin.p_min =
                p_min;

            bin.p_max =
                p_max;

            const std::string subfolder_name =
                "theta_" +
                format_number(theta_min) +
                "-" +
                format_number(theta_max) +
                "_P_" +
                format_number(p_min) +
                "-" +
                format_number(p_max);

            bin.folder =
                join_path(
                    output_root,
                    subfolder_name
                );

            if (!create_folder(bin.folder)) {
                return;
            }

            const std::string unique =
                "t" +
                std::to_string(i_theta) +
                "_p" +
                std::to_string(i_p);

            bin.h_step1 =
                new TH1D(
                    (
                        "h_step1_" +
                        unique
                    ).c_str(),
                    "",
                    n_dp_bins,
                    dp_min,
                    dp_max
                );

            bin.h_21 =
                new TH1D(
                    (
                        "h_21_" +
                        unique
                    ).c_str(),
                    "",
                    n_dp_bins,
                    dp_min,
                    dp_max
                );

            bin.h_22 =
                new TH1D(
                    (
                        "h_22_" +
                        unique
                    ).c_str(),
                    "",
                    n_dp_bins,
                    dp_min,
                    dp_max
                );

            bin.h_31 =
                new TH1D(
                    (
                        "h_31_" +
                        unique
                    ).c_str(),
                    "",
                    n_dp_bins,
                    dp_min,
                    dp_max
                );

            bin.h_32 =
                new TH1D(
                    (
                        "h_32_" +
                        unique
                    ).c_str(),
                    "",
                    n_dp_bins,
                    dp_min,
                    dp_max
                );

            // Integer-centered bins:
            // 1, 2, 3, ..., 20 positive tracks.
            bin.h_step1_positive_multiplicity = new TH1D(("h_step1_positive_multiplicity_" + unique).c_str(),"",4,0.5,4.5);

            bins.push_back(bin);
        }
    }

    long long n_events_scanned = 0;

    long long n_base_one_trigger_electron = 0;
    long long n_base_no_negative_hadrons = 0;
    long long n_base_one_lund_pip = 0;

    long long n_fail_base_electron = 0;
    long long n_fail_negative_hadrons = 0;
    long long n_fail_lund_pip = 0;
    long long n_fail_no_detector_positive = 0;
    long long n_fail_selected_outside_grid = 0;
    long long n_fail_bin_migration = 0;

    long long n_files_opened = 0;
    long long n_files_skipped_open = 0;
    long long n_files_skipped_missing_rec = 0;
    long long n_files_skipped_missing_lund = 0;

    bool reached_max_events =
        false;

    for (size_t ifile = 0;
         ifile < hipo_files.size() &&
         !reached_max_events;
         ifile++) {

        const std::string& input_hipo =
            hipo_files[ifile];

        std::cout << "\nOpening file "
                  << (ifile + 1)
                  << " / "
                  << hipo_files.size()
                  << ": "
                  << input_hipo
                  << "\n";

        hipo::reader reader;

        try {
            reader.open(
                input_hipo.c_str()
            );
        }
        catch (...) {
            std::cerr
                << "WARNING: failed to open, skipping:\n"
                << "  "
                << input_hipo
                << "\n";

            n_files_skipped_open++;
            continue;
        }

        if (!reader.is_open()) {
            std::cerr
                << "WARNING: reader is not open, skipping:\n"
                << "  "
                << input_hipo
                << "\n";

            n_files_skipped_open++;
            continue;
        }

        n_files_opened++;

        hipo::dictionary factory;

        reader.readDictionary(
            factory
        );

        hipo::schema rec_schema;
        hipo::schema mclund_schema;

        if (!get_schema_safe(
                factory,
                "REC::Particle",
                rec_schema)) {

            std::cerr
                << "WARNING: REC::Particle missing, skipping:\n"
                << "  "
                << input_hipo
                << "\n";

            n_files_skipped_missing_rec++;
            continue;
        }

        if (!get_schema_safe(
                factory,
                "MC::Lund",
                mclund_schema)) {

            std::cerr
                << "WARNING: MC::Lund missing, skipping:\n"
                << "  "
                << input_hipo
                << "\n";

            n_files_skipped_missing_lund++;
            continue;
        }

        hipo::event event;

        hipo::bank rec_part(
            rec_schema
        );

        hipo::bank mc_lund(
            mclund_schema
        );

        while (reader.next()) {
            if (max_events_to_scan > 0 &&
                n_events_scanned >= max_events_to_scan) {

                reached_max_events =
                    true;

                break;
            }

            reader.read(event);

            n_events_scanned++;

            if (n_events_scanned % 50000 == 0) {
                std::cout
                    << "Scanned "
                    << n_events_scanned
                    << " | base electron="
                    << n_base_one_trigger_electron
                    << " | no negative hadrons="
                    << n_base_no_negative_hadrons
                    << " | one MC::Lund pi+="
                    << n_base_one_lund_pip
                    << " | migration rejected="
                    << n_fail_bin_migration
                    << "\n";
            }

            event.getStructure(
                rec_part
            );

            event.getStructure(
                mc_lund
            );

            int n_electron_total = 0;
            int n_trigger_electron = 0;
            int n_negative_hadrons_total = 0;
            int n_positive_tracks_total = 0;
            int n_rec_pip_total = 0;

            std::vector<TrackInfo>
                detector_positive_tracks;

            for (int i = 0;
                 i < rec_part.getRows();
                 i++) {

                const int pid =
                    rec_part.getInt(
                        "pid",
                        i
                    );

                const int charge =
                    rec_part.getInt(
                        "charge",
                        i
                    );

                const int status =
                    rec_part.getInt(
                        "status",
                        i
                    );

                const double px =
                    rec_part.getFloat(
                        "px",
                        i
                    );

                const double py =
                    rec_part.getFloat(
                        "py",
                        i
                    );

                const double pz =
                    rec_part.getFloat(
                        "pz",
                        i
                    );

                if (pid == 11) {
                    n_electron_total++;

                    if (status < 0) {
                        n_trigger_electron++;
                    }
                }

                if (charge < 0 &&
                    pid != 11) {

                    n_negative_hadrons_total++;
                }

                if (pid == 211) {
                    n_rec_pip_total++;
                }

                TrackInfo trk;

                trk.index =
                    i;

                trk.pid =
                    pid;

                trk.charge =
                    charge;

                trk.status =
                    status;

                trk.px =
                    px;

                trk.py =
                    py;

                trk.pz =
                    pz;

                trk.p =
                    momentum(
                        px,
                        py,
                        pz
                    );

                trk.theta_deg =
                    theta_deg_from_p(
                        px,
                        py,
                        pz
                    );

                if (charge > 0) {
                    n_positive_tracks_total++;

                    if (detector_status_pass(
                            status,
                            det)) {

                        detector_positive_tracks
                            .push_back(trk);
                    }
                }
            }

            if (!(n_electron_total == 1 &&
                  n_trigger_electron == 1)) {

                n_fail_base_electron++;
                continue;
            }

            n_base_one_trigger_electron++;

            if (n_negative_hadrons_total != 0) {
                n_fail_negative_hadrons++;
                continue;
            }

            n_base_no_negative_hadrons++;

            double p_lund =
                std::numeric_limits<double>::quiet_NaN();

            double theta_lund_deg =
                std::numeric_limits<double>::quiet_NaN();

            int n_lund_pip = 0;

            if (!get_mclund_pip_kinematics(
                    mc_lund,
                    p_lund,
                    theta_lund_deg,
                    n_lund_pip)) {

                n_fail_lund_pip++;
                continue;
            }

            n_base_one_lund_pip++;

            if (detector_positive_tracks.empty()) {
                n_fail_no_detector_positive++;
                continue;
            }

            TrackInfo selected_step1;

            double best_abs_dp =
                std::numeric_limits<double>::infinity();

            for (const auto& trk :
                 detector_positive_tracks) {

                const double abs_dp =
                    std::fabs(
                        trk.p -
                        p_lund
                    );

                if (abs_dp < best_abs_dp) {
                    best_abs_dp =
                        abs_dp;

                    selected_step1 =
                        trk;
                }
            }

            const int i_theta_rec =
                get_grid_bin(
                    selected_step1.theta_deg,
                    theta_grid_min,
                    theta_grid_max,
                    theta_step,
                    n_theta_bins
                );

            const int i_p_rec =
                get_grid_bin(
                    selected_step1.p,
                    p_grid_min,
                    p_grid_max,
                    p_step,
                    n_p_bins
                );

            if (i_theta_rec < 0 ||
                i_p_rec < 0) {

                n_fail_selected_outside_grid++;
                continue;
            }

            const int flat_bin_rec =
                i_theta_rec *
                n_p_bins +
                i_p_rec;

            BinData& bin =
                bins.at(
                    flat_bin_rec
                );

            if (!allow_bin_migration) {
                const int i_theta_gen =
                    get_grid_bin(
                        theta_lund_deg,
                        theta_grid_min,
                        theta_grid_max,
                        theta_step,
                        n_theta_bins
                    );

                const int i_p_gen =
                    get_grid_bin(
                        p_lund,
                        p_grid_min,
                        p_grid_max,
                        p_step,
                        n_p_bins
                    );

                if (i_theta_gen != i_theta_rec ||
                    i_p_gen != i_p_rec) {

                    n_fail_bin_migration++;
                    bin.n_rejected_bin_migration++;

                    continue;
                }
            }

            const double dp =
                selected_step1.p -
                p_lund;

            // Step 1.
            bin.n_step1++;

            bin.h_step1->Fill(
                dp
            );

            // Positive-track multiplicity for the Step-1 sample.
            // This counts all REC::Particle rows with charge > 0,
            // regardless of their FD/CD status.
            bin.h_step1_positive_multiplicity
               ->Fill(
                   n_positive_tracks_total
               );

            // Step 2.1.
            const bool pass_21 =
                selected_step1.pid == 211;

            if (pass_21) {
                bin.n_21++;

                bin.h_21->Fill(
                    dp
                );

                update_step21_multiplicity(
                    bin,
                    n_rec_pip_total,
                    n_positive_tracks_total
                );
            }

            // Step 2.2.
            const bool pass_22 =
                n_positive_tracks_total == 1;

            if (pass_22) {
                bin.n_22++;

                bin.h_22->Fill(
                    dp
                );
            }

            // Step 3.1.
            const bool pass_31 =
                pass_21 &&
                n_positive_tracks_total == 1;

            if (pass_31) {
                bin.n_31++;

                bin.h_31->Fill(
                    dp
                );
            }

            // Step 3.2.
            const bool pass_32 =
                pass_22 &&
                selected_step1.pid == 211;

            if (pass_32) {
                bin.n_32++;

                bin.h_32->Fill(
                    dp
                );
            }
        }
    }

    std::cout << "\nFinished event loop.\n";

    std::cout << "Files listed in .dat:                    "
              << hipo_files.size()
              << "\n";

    std::cout << "Files opened:                            "
              << n_files_opened
              << "\n";

    std::cout << "Files skipped open error:                "
              << n_files_skipped_open
              << "\n";

    std::cout << "Files skipped missing REC::Particle:     "
              << n_files_skipped_missing_rec
              << "\n";

    std::cout << "Files skipped missing MC::Lund:          "
              << n_files_skipped_missing_lund
              << "\n";

    std::cout << "Events scanned:                          "
              << n_events_scanned
              << "\n";

    std::cout << "Base exactly 1 trigger electron only:    "
              << n_base_one_trigger_electron
              << "\n";

    std::cout << "Base no negative hadrons:                "
              << n_base_no_negative_hadrons
              << "\n";

    std::cout << "Base exactly 1 MC::Lund pi+:             "
              << n_base_one_lund_pip
              << "\n";

    std::cout << "Fail base electron requirement:          "
              << n_fail_base_electron
              << "\n";

    std::cout << "Fail negative-hadron requirement:        "
              << n_fail_negative_hadrons
              << "\n";

    std::cout << "Fail exactly-one-MC-pi+ requirement:     "
              << n_fail_lund_pip
              << "\n";

    std::cout << "Fail no detector-positive track:         "
              << n_fail_no_detector_positive
              << "\n";

    std::cout << "Selected track outside theta/p grid:     "
              << n_fail_selected_outside_grid
              << "\n";

    std::cout << "Rejected reconstructed/gen bin mismatch: "
              << n_fail_bin_migration
              << "\n";

    std::cout << "\nWriting PNG and TXT outputs...\n";

    for (auto& bin : bins) {
        const std::string step1_title =
            make_plot_title(
                "Step 1: selected detector-positive track",
                bin,
                det,
                allow_bin_migration,
                "#Delta p = p_{rec}(selected positive) - p_{Lund}(#pi^{+}) [GeV]",
                "Counts"
            );

        const std::string step21_title =
            make_plot_title(
                "Step 2.1: selected Step-1 track is #pi^{+}",
                bin,
                det,
                allow_bin_migration,
                "#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV]",
                "Counts"
            );

        const std::string step22_title =
            make_plot_title(
                "Step 2.2: exactly one positive track total",
                bin,
                det,
                allow_bin_migration,
                "#Delta p = p_{rec}(selected positive) - p_{Lund}(#pi^{+}) [GeV]",
                "Counts"
            );

        const std::string step31_title =
            make_plot_title(
                "Step 3.1: Step-2.1 pion is the only positive track",
                bin,
                det,
                allow_bin_migration,
                "#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV]",
                "Counts"
            );

        const std::string step32_title =
            make_plot_title(
                "Step 3.2: Step-2.2 only positive track is #pi^{+}",
                bin,
                det,
                allow_bin_migration,
                "#Delta p = p_{rec}(selected #pi^{+}) - p_{Lund}(#pi^{+}) [GeV]",
                "Counts"
            );

        const std::string multiplicity_title =
            make_plot_title(
                "Step 1: positive-track multiplicity",
                bin,
                det,
                allow_bin_migration,
                "Number of positive REC::Particle tracks",
                "Counts"
            );

        save_dp_hist_png(
            bin.h_step1,
            join_path(
                bin.folder,
                "step1_dp.png"
            ),
            step1_title,
            bin.n_step1
        );

        save_dp_hist_png(
            bin.h_21,
            join_path(
                bin.folder,
                "step2p1_dp.png"
            ),
            step21_title,
            bin.n_21
        );

        save_dp_hist_png(
            bin.h_22,
            join_path(
                bin.folder,
                "step2p2_dp.png"
            ),
            step22_title,
            bin.n_22
        );

        save_dp_hist_png(
            bin.h_31,
            join_path(
                bin.folder,
                "step3p1_dp.png"
            ),
            step31_title,
            bin.n_31
        );

        save_dp_hist_png(
            bin.h_32,
            join_path(
                bin.folder,
                "step3p2_dp.png"
            ),
            step32_title,
            bin.n_32
        );

        save_multiplicity_hist_png(
            bin.h_step1_positive_multiplicity,
            join_path(
                bin.folder,
                "step1_positive_track_multiplicity.png"
            ),
            multiplicity_title,
            bin.n_step1
        );

        write_bin_stats(
            bin,
            input_dat_file,
            det,
            allow_bin_migration,
            n_events_scanned,
            n_base_one_trigger_electron,
            n_base_no_negative_hadrons,
            n_base_one_lund_pip,
            n_fail_bin_migration
        );
    }

    std::cout << "\nDone.\n";

    std::cout << "Output folder:\n";
    std::cout << "  "
              << output_root
              << "\n";

    std::cout << "Subfolders written: "
              << bins.size()
              << "\n";

    for (auto& bin : bins) {
        delete bin.h_step1;
        delete bin.h_21;
        delete bin.h_22;
        delete bin.h_31;
        delete bin.h_32;
        delete bin.h_step1_positive_multiplicity;
    }
}