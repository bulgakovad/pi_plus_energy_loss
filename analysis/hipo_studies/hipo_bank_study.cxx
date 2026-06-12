#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <cmath>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TLine.h"

#include "hipo4/reader.h"
#include "hipo4/dictionary.h"
#include "hipo4/event.h"
#include "hipo4/bank.h"

static std::string base_name_no_ext(const std::string& path)
{
    std::string name = path;

    const size_t slash = name.find_last_of("/");
    if (slash != std::string::npos) {
        name = name.substr(slash + 1);
    }

    const size_t dot = name.find_last_of(".");
    if (dot != std::string::npos) {
        name = name.substr(0, dot);
    }

    return name;
}

static const char* track_detector_label(int det)
{
    // For REC::Track, the important expected values are usually:
    // 5 = CVT-like central tracking
    // 6 = DC-like forward tracking
    if (det == 5) return "CVT/CD-like";
    if (det == 6) return "DC/FD-like";
    return "";
}

void plot_RECTrack_detector_pip(const char* input_hipo,
                                const char* output_png = "",
                                int require_exactly_one_rec_pip = 1)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECTrack_detector_pip_" + base_name_no_ext(input_hipo) + ".png";
    }
    else {
        std::string s = outname;

        // If user gave only a filename, put it into plots/.
        // If user gave a path like "some_folder/file.png", keep it as-is.
        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recTrack(factory.getSchema("REC::Track"));

   TH1D* h_det = new TH1D(
    "h_det",
    "REC::Track detector for #pi^{+}-associated rows;REC::Track.detector;Counts",
    11, -0.5, 10.5
    );

    // Label every integer detector-code bin: 0, 1, 2, ..., 20
    for (int det = 0; det <= 10; det++) {
        int bin = h_det->GetXaxis()->FindBin(det);
        h_det->GetXaxis()->SetBinLabel(bin, Form("%d", det));
    }

    h_det->GetXaxis()->LabelsOption("h");   // horizontal labels
    h_det->GetXaxis()->SetLabelSize(0.035);

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_track_for_pip = 0;

    long long n_track_rows_total = 0;
    long long n_track_rows_for_pip = 0;

    std::map<int, long long> det_counts;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recTrack);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i; // this row index is what REC::Track.pindex points to
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_track_for_this_pip = false;

        const int n_tracks = recTrack.getRows();
        n_track_rows_total += n_tracks;

        for (int j = 0; j < n_tracks; j++) {
            const int pidx = recTrack.getShort("pindex", j);

            // Critical association:
            // REC::Track row belongs to the pion iff REC::Track.pindex == pion REC::Particle row.
            if (pidx != pip_pindex) continue;

            const int det = recTrack.getByte("detector", j);

            h_det->Fill(det);
            det_counts[det]++;

            n_track_rows_for_pip++;
            found_track_for_this_pip = true;
        }

        if (!found_track_for_this_pip) {
            n_events_no_track_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                         " << n_events << "\n"
              << "  Events with exactly one REC pi+:      " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:           " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:     " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                          " << n_events_used << "\n"
              << "  Total REC::Track rows:                " << n_track_rows_total << "\n"
              << "  Pion-associated REC::Track rows:      " << n_track_rows_for_pip << "\n"
              << "  Events with no pion REC::Track row:   " << n_events_no_track_for_pip << "\n";

    std::cout << "\nREC::Track.detector counts for pion-associated rows:\n";
    for (const auto& kv : det_counts) {
        std::cout << "  detector = " << std::setw(3) << kv.first
                  << "  count = " << std::setw(8) << kv.second
                  << "  " << track_detector_label(kv.first)
                  << "\n";
    }

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c", "REC::Track detector for pion", 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.14);

    h_det->SetLineWidth(2);
    h_det->SetMinimum(0.0);
    h_det->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(0.16, 0.82, Form("File: %s", base_name_no_ext(input_hipo).c_str()));
    lat.DrawLatex(0.16, 0.77, Form("Events used: %d", n_events_used));
    lat.DrawLatex(0.16, 0.72, Form("#pi^{+}-associated REC::Track rows: %lld", n_track_rows_for_pip));

    lat.DrawLatex(0.68, 0.68, "5 = CVT / CD-like tracking");
    lat.DrawLatex(0.68, 0.63, "6 = DC / FD-like tracking");
    lat.SetTextSize(0.032);

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}

void plot_RECTrack_chi2NDF_pip(const char* input_hipo,
                               const char* output_png = "",
                               int require_exactly_one_rec_pip = 1,
                               int detector_filter = -1,
                               int nBins = 100,
                               double xMin = 0.0,
                               double xMax = 40.0,
                               bool logy = false)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECTrack_chi2NDF_pip_" + base_name_no_ext(input_hipo);

        if (detector_filter > 0) {
            outname += Form("_det%d", detector_filter);
        }

        outname += ".png";
    }
    else {
        std::string s = outname;

        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recTrack(factory.getSchema("REC::Track"));

    std::string title = "REC::Track #chi^{2}/NDF for #pi^{+}-associated rows";

    if (detector_filter > 0) {
        title += Form(", detector = %d", detector_filter);
    }

    title += ";REC::Track #chi^{2}/NDF;Counts";

    TH1D* h_chi2ndf = new TH1D(
        "h_track_chi2ndf",
        title.c_str(),
        nBins, xMin, xMax
    );

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_track_for_pip = 0;

    long long n_track_rows_total = 0;
    long long n_track_rows_for_pip = 0;
    long long n_track_rows_after_detector_filter = 0;
    long long n_chi2ndf_filled = 0;
    long long n_tracks_bad_ndf = 0;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recTrack);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_track_for_this_pip = false;

        const int n_tracks = recTrack.getRows();
        n_track_rows_total += n_tracks;

        for (int j = 0; j < n_tracks; j++) {
            const int pidx = recTrack.getShort("pindex", j);

            // REC::Track row belongs to the reconstructed pi+ iff
            // REC::Track.pindex equals the REC::Particle row index of the pi+.
            if (pidx != pip_pindex) continue;

            found_track_for_this_pip = true;
            n_track_rows_for_pip++;

            const int det = recTrack.getByte("detector", j);

            if (detector_filter > 0 && det != detector_filter) continue;
            n_track_rows_after_detector_filter++;

            const double chi2 = recTrack.getFloat("chi2", j);
            const int ndf = recTrack.getShort("NDF", j);

            if (ndf <= 0) {
                n_tracks_bad_ndf++;
                continue;
            }

            const double chi2_ndf = chi2 / static_cast<double>(ndf);

            h_chi2ndf->Fill(chi2_ndf);
            n_chi2ndf_filled++;
        }

        if (!found_track_for_this_pip) {
            n_events_no_track_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                         " << n_events << "\n"
              << "  Events with exactly one REC pi+:     " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:          " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:    " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                         " << n_events_used << "\n"
              << "  Total REC::Track rows:               " << n_track_rows_total << "\n"
              << "  Pion-associated REC::Track rows:     " << n_track_rows_for_pip << "\n"
              << "  Rows after detector filter:          " << n_track_rows_after_detector_filter << "\n"
              << "  Filled chi2/NDF entries:             " << n_chi2ndf_filled << "\n"
              << "  Tracks skipped because NDF <= 0:     " << n_tracks_bad_ndf << "\n"
              << "  Events with no pion REC::Track row:  " << n_events_no_track_for_pip << "\n";

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas(
        "c_track_chi2ndf",
        "REC::Track chi2/NDF for pion",
        900, 700
    );

    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    if (logy) {
        c->SetLogy();
    }

    h_chi2ndf->SetLineWidth(2);
    h_chi2ndf->SetMinimum(0.0);
    h_chi2ndf->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(
        0.16, 0.82,
        Form("File: %s", base_name_no_ext(input_hipo).c_str())
    );

    lat.DrawLatex(
        0.16, 0.77,
        Form("Events used: %d", n_events_used)
    );

    lat.DrawLatex(
        0.16, 0.72,
        Form("#pi^{+}-associated REC::Track rows: %lld", n_track_rows_for_pip)
    );

    lat.DrawLatex(
        0.16, 0.67,
        Form("Filled #chi^{2}/NDF entries: %lld", n_chi2ndf_filled)
    );

    lat.DrawLatex(
        0.16, 0.62,
        Form("Mean #chi^{2}/NDF: %.3f", h_chi2ndf->GetMean())
    );

    if (detector_filter > 0) {
        lat.DrawLatex(
            0.16, 0.57,
            Form("Detector filter: %d", detector_filter)
        );
    }
    else {
        lat.DrawLatex(0.16, 0.57, "Detector filter: none");
    }

    lat.DrawLatex(0.58, 0.77, "REC::Track detector IDs:");
    lat.DrawLatex(0.58, 0.72, "5 = CVT / CD-like");
    lat.DrawLatex(0.58, 0.67, "6 = DC / FD-like");

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}
void plot_RECTraj_detector_pip(const char* input_hipo,
                               const char* output_png = "",
                               int require_exactly_one_rec_pip = 1)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECTraj_detector_pip_" + base_name_no_ext(input_hipo) + ".png";
    }
    else {
        std::string s = outname;
        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recTraj(factory.getSchema("REC::Traj"));

    TH1D* h_det = new TH1D(
        "h_traj_det",
        "REC::Traj detector for #pi^{+}-associated rows;REC::Traj.detector;Counts",
        31, -0.5, 30.5
    );

    // Label all integer detector-code bins
    for (int det = 0; det <= 30; det++) {
        int bin = h_det->GetXaxis()->FindBin(det);
        h_det->GetXaxis()->SetBinLabel(bin, Form("%d", det));
    }

    h_det->GetXaxis()->LabelsOption("h");
    h_det->GetXaxis()->SetLabelSize(0.030);

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_traj_for_pip = 0;

    long long n_traj_rows_total = 0;
    long long n_traj_rows_for_pip = 0;

    std::map<int, long long> det_counts;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recTraj);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_traj_for_this_pip = false;

        const int n_traj = recTraj.getRows();
        n_traj_rows_total += n_traj;

        for (int j = 0; j < n_traj; j++) {
            const int pidx = recTraj.getShort("pindex", j);

            // Critical association:
            // REC::Traj row belongs to the pion iff REC::Traj.pindex == pion REC::Particle row.
            if (pidx != pip_pindex) continue;

            const int det = recTraj.getByte("detector", j);

            h_det->Fill(det);
            det_counts[det]++;

            n_traj_rows_for_pip++;
            found_traj_for_this_pip = true;
        }

        if (!found_traj_for_this_pip) {
            n_events_no_traj_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                         " << n_events << "\n"
              << "  Events with exactly one REC pi+:      " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:           " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:     " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                          " << n_events_used << "\n"
              << "  Total REC::Traj rows:                 " << n_traj_rows_total << "\n"
              << "  Pion-associated REC::Traj rows:       " << n_traj_rows_for_pip << "\n"
              << "  Events with no pion REC::Traj row:    " << n_events_no_traj_for_pip << "\n";

    std::cout << "\nREC::Traj.detector counts for pion-associated rows:\n";
    for (const auto& kv : det_counts) {
        std::cout << "  detector = " << std::setw(3) << kv.first
                  << "  count = " << std::setw(8) << kv.second
                  << "\n";
    }

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c_traj_det", "REC::Traj detector for pion", 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.14);

    h_det->SetLineWidth(2);
    h_det->SetMinimum(0.0);
    h_det->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(0.16, 0.82, Form("File: %s", base_name_no_ext(input_hipo).c_str()));
    lat.DrawLatex(0.16, 0.77, Form("Events used: %d", n_events_used));
    lat.DrawLatex(0.16, 0.72, Form("#pi^{+}-associated REC::Traj rows: %lld", n_traj_rows_for_pip));

    lat.DrawLatex(0.55, 0.72, "Rows selected by REC::Traj.pindex");
    lat.DrawLatex(0.55, 0.67, "Multiple traj rows per #pi^{+} are expected");

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}

void plot_RECTraj_layer_pip(const char* input_hipo,
                            const char* output_png = "",
                            int require_exactly_one_rec_pip = 1,
                            int detector_filter = -1,
                            int nBins = 5,
                            double xMin = -0.5,
                            double xMax = 4.5)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECTraj_layer_pip_" + base_name_no_ext(input_hipo);

        if (detector_filter > 0) {
            outname += Form("_det%d", detector_filter);
        }

        outname += ".png";
    }
    else {
        std::string s = outname;
        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recTraj(factory.getSchema("REC::Traj"));

    std::string title = "REC::Traj layer for #pi^{+}-associated rows";

    if (detector_filter > 0) {
        title += Form(", detector = %d", detector_filter);
    }

    title += ";REC::Traj.layer;Counts";

    TH1D* h_layer = new TH1D(
        "h_traj_layer",
        title.c_str(),
        nBins, xMin, xMax
    );

    // Label integer layer bins
    for (int layer = int(xMin + 0.5); layer <= int(xMax - 0.5); layer++) {
        if (layer < 0) continue;
        int bin = h_layer->GetXaxis()->FindBin(layer);
        h_layer->GetXaxis()->SetBinLabel(bin, Form("%d", layer));
    }

    h_layer->GetXaxis()->LabelsOption("h");
    h_layer->GetXaxis()->SetLabelSize(0.025);

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_traj_for_pip = 0;

    long long n_traj_rows_total = 0;
    long long n_traj_rows_for_pip = 0;
    long long n_traj_rows_after_detector_filter = 0;

    std::map<int, long long> layer_counts;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recTraj);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_traj_for_this_pip = false;

        const int n_traj = recTraj.getRows();
        n_traj_rows_total += n_traj;

        for (int j = 0; j < n_traj; j++) {
            const int pidx = recTraj.getShort("pindex", j);

            // Critical association:
            // REC::Traj row belongs to pi+ iff REC::Traj.pindex == pion REC::Particle row.
            if (pidx != pip_pindex) continue;

            found_traj_for_this_pip = true;
            n_traj_rows_for_pip++;

            const int det = recTraj.getByte("detector", j);

            if (detector_filter > 0 && det != detector_filter) continue;

            const int layer = recTraj.getByte("layer", j);

            h_layer->Fill(layer);
            layer_counts[layer]++;

            n_traj_rows_after_detector_filter++;
        }

        if (!found_traj_for_this_pip) {
            n_events_no_traj_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                         " << n_events << "\n"
              << "  Events with exactly one REC pi+:      " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:           " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:     " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                          " << n_events_used << "\n"
              << "  Total REC::Traj rows:                 " << n_traj_rows_total << "\n"
              << "  Pion-associated REC::Traj rows:       " << n_traj_rows_for_pip << "\n"
              << "  Rows after detector filter:           " << n_traj_rows_after_detector_filter << "\n"
              << "  Events with no pion REC::Traj row:    " << n_events_no_traj_for_pip << "\n";

    std::cout << "\nREC::Traj.layer counts for pion-associated rows";
    if (detector_filter > 0) {
        std::cout << " with detector = " << detector_filter;
    }
    std::cout << ":\n";

    for (const auto& kv : layer_counts) {
        std::cout << "  layer = " << std::setw(3) << kv.first
                  << "  count = " << std::setw(8) << kv.second
                  << "\n";
    }

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c_traj_layer", "REC::Traj layer for pion", 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.14);

    h_layer->SetLineWidth(2);
    h_layer->SetMinimum(0.0);
    h_layer->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(0.16, 0.82, Form("File: %s", base_name_no_ext(input_hipo).c_str()));
    lat.DrawLatex(0.16, 0.77, Form("Events used: %d", n_events_used));
    lat.DrawLatex(0.16, 0.72, Form("#pi^{+}-associated REC::Traj rows: %lld", n_traj_rows_for_pip));

    if (detector_filter > 0) {
        lat.DrawLatex(0.16, 0.67, Form("Detector filter: %d", detector_filter));
    }
    else {
        lat.DrawLatex(0.16, 0.67, "Detector filter: none");
    }

    lat.DrawLatex(0.52, 0.72, "Rows selected by REC::Traj.pindex");
    lat.DrawLatex(0.52, 0.67, "Multiple traj rows per #pi^{+} are expected");

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}
void plot_RECTraj_xy_pip(const char* input_hipo,
                         const char* output_png = "",
                         int require_exactly_one_rec_pip = 1,
                         int detector_filter = -1,
                         int nBinsX = 300,
                         double xMin = -500.0,
                         double xMax = 500.0,
                         int nBinsY = 300,
                         double yMin = -500.0,
                         double yMax = 500.0,
                         bool logz = true)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECTraj_xy_pip_" + base_name_no_ext(input_hipo);

        if (detector_filter > 0) {
            outname += Form("_det%d", detector_filter);
        }

        outname += ".png";
    }
    else {
        std::string s = outname;
        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recTraj(factory.getSchema("REC::Traj"));

    std::string title = "REC::Traj x vs y for #pi^{+}-associated rows";

    if (detector_filter > 0) {
        title += Form(", detector = %d", detector_filter);
    }

    title += ";REC::Traj.x (cm);REC::Traj.y (cm)";

    TH2D* h_xy = new TH2D(
        "h_traj_xy",
        title.c_str(),
        nBinsX, xMin, xMax,
        nBinsY, yMin, yMax
    );

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_traj_for_pip = 0;

    long long n_traj_rows_total = 0;
    long long n_traj_rows_for_pip = 0;
    long long n_traj_rows_after_detector_filter = 0;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recTraj);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_traj_for_this_pip = false;

        const int n_traj = recTraj.getRows();
        n_traj_rows_total += n_traj;

        for (int j = 0; j < n_traj; j++) {
            const int pidx = recTraj.getShort("pindex", j);

            // Critical association:
            // REC::Traj row belongs to pi+ iff REC::Traj.pindex == pion REC::Particle row.
            if (pidx != pip_pindex) continue;

            found_traj_for_this_pip = true;
            n_traj_rows_for_pip++;

            const int det = recTraj.getByte("detector", j);

            if (detector_filter > 0 && det != detector_filter) continue;

            const double x = recTraj.getFloat("x", j);
            const double y = recTraj.getFloat("y", j);

            h_xy->Fill(x, y);

            n_traj_rows_after_detector_filter++;
        }

        if (!found_traj_for_this_pip) {
            n_events_no_traj_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                         " << n_events << "\n"
              << "  Events with exactly one REC pi+:      " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:           " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:     " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                          " << n_events_used << "\n"
              << "  Total REC::Traj rows:                 " << n_traj_rows_total << "\n"
              << "  Pion-associated REC::Traj rows:       " << n_traj_rows_for_pip << "\n"
              << "  Rows after detector filter:           " << n_traj_rows_after_detector_filter << "\n"
              << "  Events with no pion REC::Traj row:    " << n_events_no_traj_for_pip << "\n";

    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("c_traj_xy", "REC::Traj x vs y for pion", 900, 800);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    if (logz) {
        c->SetLogz();
    }

    h_xy->Draw("COLZ");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(0.14, 0.93, Form("File: %s", base_name_no_ext(input_hipo).c_str()));
    lat.DrawLatex(0.14, 0.88, Form("Events used: %d", n_events_used));
    lat.DrawLatex(0.14, 0.83, Form("#pi^{+}-associated REC::Traj rows: %lld", n_traj_rows_for_pip));

    if (detector_filter > 0) {
        lat.DrawLatex(0.14, 0.78, Form("Detector filter: %d", detector_filter));
    }
    else {
        lat.DrawLatex(0.14, 0.78, "Detector filter: none");
    }

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}

void plot_RECTraj_edge_pip(const char* input_hipo,
                           const char* output_png = "",
                           int require_exactly_one_rec_pip = 1,
                           int detector_filter = -1,
                           int layer_filter = -1,
                           int nBins = 120,
                           double edgeMin = -100.0,
                           double edgeMax = 100.0,
                           bool logy = false)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECTraj_edge_pip_" + base_name_no_ext(input_hipo);

        if (detector_filter > 0) {
            outname += Form("_det%d", detector_filter);
        }

        if (layer_filter > 0) {
            outname += Form("_layer%d", layer_filter);
        }

        outname += ".png";
    }
    else {
        std::string s = outname;
        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recTraj(factory.getSchema("REC::Traj"));

    std::string title = "REC::Traj edge for #pi^{+}-associated rows";

    if (detector_filter > 0) {
        title += Form(", detector = %d", detector_filter);
    }

    if (layer_filter > 0) {
        title += Form(", layer = %d", layer_filter);
    }

    title += ";REC::Traj.edge (cm);Counts";

    TH1D* h_edge = new TH1D(
        "h_traj_edge",
        title.c_str(),
        nBins, edgeMin, edgeMax
    );

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_traj_for_pip = 0;

    long long n_traj_rows_total = 0;
    long long n_traj_rows_for_pip = 0;
    long long n_traj_rows_after_detector_filter = 0;
    long long n_traj_rows_after_layer_filter = 0;

    long long n_edge_filled = 0;
    long long n_edge_undefined = 0;
    long long n_edge_negative = 0;
    long long n_edge_positive = 0;
    long long n_edge_near_zero = 0;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recTraj);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_traj_for_this_pip = false;

        const int n_traj = recTraj.getRows();
        n_traj_rows_total += n_traj;

        for (int j = 0; j < n_traj; j++) {
            const int pidx = recTraj.getShort("pindex", j);

            // Critical association:
            // REC::Traj row belongs to pi+ iff REC::Traj.pindex == pion REC::Particle row.
            if (pidx != pip_pindex) continue;

            found_traj_for_this_pip = true;
            n_traj_rows_for_pip++;

            const int det = recTraj.getByte("detector", j);

            if (detector_filter > 0 && det != detector_filter) continue;
            n_traj_rows_after_detector_filter++;

            const int layer = recTraj.getByte("layer", j);

            if (layer_filter > 0 && layer != layer_filter) continue;
            n_traj_rows_after_layer_filter++;

            const double edge = recTraj.getFloat("edge", j);

            // Undefined boundary values are typically around -99 cm.
            if (edge <= -98.0) {
                n_edge_undefined++;
                continue;
            }

            h_edge->Fill(edge);
            n_edge_filled++;

            if (edge < 0.0) n_edge_negative++;
            if (edge > 0.0) n_edge_positive++;
            if (std::abs(edge) < 2.0) n_edge_near_zero++;
        }

        if (!found_traj_for_this_pip) {
            n_events_no_traj_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                         " << n_events << "\n"
              << "  Events with exactly one REC pi+:      " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:           " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:     " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                          " << n_events_used << "\n"
              << "  Total REC::Traj rows:                 " << n_traj_rows_total << "\n"
              << "  Pion-associated REC::Traj rows:       " << n_traj_rows_for_pip << "\n"
              << "  Rows after detector filter:           " << n_traj_rows_after_detector_filter << "\n"
              << "  Rows after layer filter:              " << n_traj_rows_after_layer_filter << "\n"
              << "  Filled edge entries:                  " << n_edge_filled << "\n"
              << "  Undefined edge entries skipped:       " << n_edge_undefined << "\n"
              << "  edge < 0 entries:                     " << n_edge_negative << "\n"
              << "  edge > 0 entries:                     " << n_edge_positive << "\n"
              << "  |edge| < 2 cm entries:                " << n_edge_near_zero << "\n"
              << "  Events with no pion REC::Traj row:    " << n_events_no_traj_for_pip << "\n";

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c_traj_edge", "REC::Traj edge for pion", 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    if (logy) {
        c->SetLogy();
    }

    h_edge->SetLineWidth(2);
    h_edge->SetMinimum(0.0);
    h_edge->Draw("hist");

    TLine* zero = new TLine(0.0, 0.0, 0.0, h_edge->GetMaximum() * 1.05);
    zero->SetLineStyle(2);
    zero->SetLineWidth(2);
    zero->Draw("same");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(0.16, 0.82, Form("File: %s", base_name_no_ext(input_hipo).c_str()));
    lat.DrawLatex(0.16, 0.77, Form("Events used: %d", n_events_used));
    lat.DrawLatex(0.16, 0.72, Form("#pi^{+}-associated REC::Traj rows: %lld", n_traj_rows_for_pip));

    if (detector_filter > 0) {
        lat.DrawLatex(0.16, 0.67, Form("Detector filter: %d", detector_filter));
    }
    else {
        lat.DrawLatex(0.16, 0.67, "Detector filter: none");
    }

    if (layer_filter > 0) {
        lat.DrawLatex(0.16, 0.62, Form("Layer filter: %d", layer_filter));
    }
    else {
        lat.DrawLatex(0.16, 0.62, "Layer filter: none");
    }

    lat.DrawLatex(0.55, 0.78, "edge > 0: inside contour");
    lat.DrawLatex(0.55, 0.73, "edge #approx 0: near detector edge");
    lat.DrawLatex(0.55, 0.68, "edge < 0: outside contour");

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}

/////////////////////////////// REC::Scintillator plots////////////////////////////////////////

void plot_RECScintillator_detector_pip(const char* input_hipo,
                                       const char* output_png = "",
                                       int require_exactly_one_rec_pip = 1,
                                       int nBins = 31,
                                       double xMin = -0.5,
                                       double xMax = 30.5)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECScintillator_detector_pip_" + base_name_no_ext(input_hipo) + ".png";
    }
    else {
        std::string s = outname;
        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recScint(factory.getSchema("REC::Scintillator"));

    TH1D* h_det = new TH1D(
        "h_scint_detector",
        "REC::Scintillator detector for #pi^{+}-associated rows;REC::Scintillator.detector;Counts",
        nBins, xMin, xMax
    );

    // Label integer detector-code bins
    for (int det = int(xMin + 0.5); det <= int(xMax - 0.5); det++) {
        if (det < 0) continue;
        int bin = h_det->GetXaxis()->FindBin(det);
        h_det->GetXaxis()->SetBinLabel(bin, Form("%d", det));
    }

    h_det->GetXaxis()->LabelsOption("h");
    h_det->GetXaxis()->SetLabelSize(0.030);

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_scint_for_pip = 0;

    long long n_scint_rows_total = 0;
    long long n_scint_rows_for_pip = 0;

    std::map<int, long long> detector_counts;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recScint);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_scint_for_this_pip = false;

        const int n_scint = recScint.getRows();
        n_scint_rows_total += n_scint;

        for (int j = 0; j < n_scint; j++) {
            const int pidx = recScint.getShort("pindex", j);

            // Critical association:
            // REC::Scintillator row belongs to pi+ iff pindex equals pion REC::Particle row.
            if (pidx != pip_pindex) continue;

            found_scint_for_this_pip = true;
            n_scint_rows_for_pip++;

            const int det = recScint.getByte("detector", j);

            h_det->Fill(det);
            detector_counts[det]++;
        }

        if (!found_scint_for_this_pip) {
            n_events_no_scint_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                              " << n_events << "\n"
              << "  Events with exactly one REC pi+:           " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:                " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:          " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                               " << n_events_used << "\n"
              << "  Total REC::Scintillator rows:              " << n_scint_rows_total << "\n"
              << "  Pion-associated REC::Scintillator rows:    " << n_scint_rows_for_pip << "\n"
              << "  Events with no pion REC::Scintillator row: " << n_events_no_scint_for_pip << "\n";

    std::cout << "\nREC::Scintillator.detector counts for pion-associated rows:\n";

    for (const auto& kv : detector_counts) {
        std::cout << "  detector = " << std::setw(4) << kv.first
                  << "  count = " << std::setw(8) << kv.second;

        if (kv.first == 3)  std::cout << "  CND";
        if (kv.first == 4)  std::cout << "  CTOF";
        if (kv.first == 12) std::cout << "  FTOF";

        std::cout << "\n";
    }

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c_scint_detector", "REC::Scintillator detector for pion", 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.14);

    h_det->SetLineWidth(2);
    h_det->SetMinimum(0.0);
    h_det->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(0.16, 0.82, Form("File: %s", base_name_no_ext(input_hipo).c_str()));
    lat.DrawLatex(0.16, 0.77, Form("Events used: %d", n_events_used));
    lat.DrawLatex(0.16, 0.72, Form("#pi^{+}-associated REC::Scintillator rows: %lld", n_scint_rows_for_pip));

    lat.DrawLatex(0.55, 0.78, "Useful detector IDs:");
    lat.DrawLatex(0.55, 0.73, "3 = CND");
    lat.DrawLatex(0.55, 0.68, "4 = CTOF");
    lat.DrawLatex(0.55, 0.63, "12 = FTOF");

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}

void plot_RECScintillator_layer_pip(const char* input_hipo,
                                    const char* output_png = "",
                                    int require_exactly_one_rec_pip = 1,
                                    int detector_filter = 12,
                                    int nBins = 10,
                                    double xMin = 0.5,
                                    double xMax = 10.5)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECScintillator_layer_pip_" + base_name_no_ext(input_hipo);

        if (detector_filter > 0) {
            outname += Form("_det%d", detector_filter);
        }

        outname += ".png";
    }
    else {
        std::string s = outname;
        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recScint(factory.getSchema("REC::Scintillator"));

    std::string title = "REC::Scintillator layer for #pi^{+}-associated rows";

    if (detector_filter > 0) {
        title += Form(", detector = %d", detector_filter);
    }

    title += ";REC::Scintillator.layer;Counts";

    TH1D* h_layer = new TH1D(
        "h_scint_layer",
        title.c_str(),
        nBins, xMin, xMax
    );

    for (int layer = int(xMin + 0.5); layer <= int(xMax - 0.5); layer++) {
        int bin = h_layer->GetXaxis()->FindBin(layer);
        h_layer->GetXaxis()->SetBinLabel(bin, Form("%d", layer));
    }

    h_layer->GetXaxis()->LabelsOption("h");
    h_layer->GetXaxis()->SetLabelSize(0.035);

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_scint_for_pip = 0;

    long long n_scint_rows_total = 0;
    long long n_scint_rows_for_pip = 0;
    long long n_scint_rows_after_detector_filter = 0;

    std::map<int, long long> layer_counts;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recScint);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_scint_for_this_pip = false;

        const int n_scint = recScint.getRows();
        n_scint_rows_total += n_scint;

        for (int j = 0; j < n_scint; j++) {
            const int pidx = recScint.getShort("pindex", j);

            // Critical association:
            // REC::Scintillator row belongs to pi+ iff pindex equals pion REC::Particle row.
            if (pidx != pip_pindex) continue;

            found_scint_for_this_pip = true;
            n_scint_rows_for_pip++;

            const int det = recScint.getByte("detector", j);

            if (detector_filter > 0 && det != detector_filter) continue;

            const int layer = recScint.getByte("layer", j);

            h_layer->Fill(layer);
            layer_counts[layer]++;

            n_scint_rows_after_detector_filter++;
        }

        if (!found_scint_for_this_pip) {
            n_events_no_scint_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                              " << n_events << "\n"
              << "  Events with exactly one REC pi+:           " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:                " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:          " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                               " << n_events_used << "\n"
              << "  Total REC::Scintillator rows:              " << n_scint_rows_total << "\n"
              << "  Pion-associated REC::Scintillator rows:    " << n_scint_rows_for_pip << "\n"
              << "  Rows after detector filter:                " << n_scint_rows_after_detector_filter << "\n"
              << "  Events with no pion REC::Scintillator row: " << n_events_no_scint_for_pip << "\n";

    std::cout << "\nREC::Scintillator.layer counts for pion-associated rows";
    if (detector_filter > 0) {
        std::cout << " with detector = " << detector_filter;
    }
    std::cout << ":\n";

    for (const auto& kv : layer_counts) {
        std::cout << "  layer = " << std::setw(3) << kv.first
                  << "  count = " << std::setw(8) << kv.second;

        if (detector_filter == 12) {
            if (kv.first == 1) std::cout << "  FTOF1A";
            if (kv.first == 2) std::cout << "  FTOF1B";
            if (kv.first == 3) std::cout << "  FTOF2";
        }

        std::cout << "\n";
    }

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c_scint_layer", "REC::Scintillator layer for pion", 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.14);

    h_layer->SetLineWidth(2);
    h_layer->SetMinimum(0.0);
    h_layer->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(0.16, 0.82, Form("File: %s", base_name_no_ext(input_hipo).c_str()));
    lat.DrawLatex(0.16, 0.77, Form("Events used: %d", n_events_used));
    lat.DrawLatex(0.16, 0.72, Form("#pi^{+}-associated REC::Scintillator rows: %lld", n_scint_rows_for_pip));

    if (detector_filter > 0) {
        lat.DrawLatex(0.16, 0.67, Form("Detector filter: %d", detector_filter));
    }
    else {
        lat.DrawLatex(0.16, 0.67, "Detector filter: none");
    }

    if (detector_filter == 12) {
        lat.DrawLatex(0.55, 0.78, "FTOF layers:");
        lat.DrawLatex(0.55, 0.73, "1 = FTOF1A");
        lat.DrawLatex(0.55, 0.68, "2 = FTOF1B");
        lat.DrawLatex(0.55, 0.63, "3 = FTOF2");
    }

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}

void plot_RECScintillator_component_pip(const char* input_hipo,
                                        const char* output_png = "",
                                        int require_exactly_one_rec_pip = 1,
                                        int detector_filter = 12,
                                        int layer_filter = -1,
                                        int nBins = 101,
                                        double xMin = -0.5,
                                        double xMax = 100.5)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECScintillator_component_pip_" + base_name_no_ext(input_hipo);

        if (detector_filter > 0) outname += Form("_det%d", detector_filter);
        if (layer_filter > 0) outname += Form("_layer%d", layer_filter);

        outname += ".png";
    }
    else {
        std::string s = outname;
        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recScint(factory.getSchema("REC::Scintillator"));

    std::string title = "REC::Scintillator component for #pi^{+}-associated rows";

    if (detector_filter > 0) title += Form(", detector = %d", detector_filter);
    if (layer_filter > 0)    title += Form(", layer = %d", layer_filter);

    title += ";REC::Scintillator.component;Counts";

    TH1D* h_comp = new TH1D(
        "h_scint_component",
        title.c_str(),
        nBins, xMin, xMax
    );

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_scint_for_pip = 0;

    long long n_scint_rows_total = 0;
    long long n_scint_rows_for_pip = 0;
    long long n_scint_rows_after_detector_filter = 0;
    long long n_scint_rows_after_layer_filter = 0;
    long long n_component_filled = 0;

    std::map<int, long long> component_counts;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recScint);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_scint_for_this_pip = false;

        const int n_scint = recScint.getRows();
        n_scint_rows_total += n_scint;

        for (int j = 0; j < n_scint; j++) {
            const int pidx = recScint.getShort("pindex", j);

            // Critical association:
            // REC::Scintillator row belongs to pi+ iff pindex equals pion REC::Particle row.
            if (pidx != pip_pindex) continue;

            found_scint_for_this_pip = true;
            n_scint_rows_for_pip++;

            const int det = recScint.getByte("detector", j);

            if (detector_filter > 0 && det != detector_filter) continue;
            n_scint_rows_after_detector_filter++;

            const int layer = recScint.getByte("layer", j);

            if (layer_filter > 0 && layer != layer_filter) continue;
            n_scint_rows_after_layer_filter++;

            const int component = recScint.getShort("component", j);

            h_comp->Fill(component);
            component_counts[component]++;
            n_component_filled++;
        }

        if (!found_scint_for_this_pip) {
            n_events_no_scint_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                              " << n_events << "\n"
              << "  Events with exactly one REC pi+:           " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:                " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:          " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                               " << n_events_used << "\n"
              << "  Total REC::Scintillator rows:              " << n_scint_rows_total << "\n"
              << "  Pion-associated REC::Scintillator rows:    " << n_scint_rows_for_pip << "\n"
              << "  Rows after detector filter:                " << n_scint_rows_after_detector_filter << "\n"
              << "  Rows after layer filter:                   " << n_scint_rows_after_layer_filter << "\n"
              << "  Filled component entries:                  " << n_component_filled << "\n"
              << "  Events with no pion REC::Scintillator row: " << n_events_no_scint_for_pip << "\n";

    std::cout << "\nREC::Scintillator.component counts for pion-associated rows";
    if (detector_filter > 0) std::cout << " with detector = " << detector_filter;
    if (layer_filter > 0)    std::cout << ", layer = " << layer_filter;
    std::cout << ":\n";

    for (const auto& kv : component_counts) {
        std::cout << "  component = " << std::setw(4) << kv.first
                  << "  count = " << std::setw(8) << kv.second
                  << "\n";
    }

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c_scint_component", "REC::Scintillator component for pion", 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    h_comp->SetLineWidth(2);
    h_comp->SetMinimum(0.0);
    h_comp->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(0.16, 0.82, Form("File: %s", base_name_no_ext(input_hipo).c_str()));
    lat.DrawLatex(0.16, 0.77, Form("Events used: %d", n_events_used));
    lat.DrawLatex(0.16, 0.72, Form("#pi^{+}-associated REC::Scintillator rows: %lld", n_scint_rows_for_pip));

    if (detector_filter > 0) {
        lat.DrawLatex(0.16, 0.67, Form("Detector filter: %d", detector_filter));
    }
    else {
        lat.DrawLatex(0.16, 0.67, "Detector filter: none");
    }

    if (layer_filter > 0) {
        lat.DrawLatex(0.16, 0.62, Form("Layer filter: %d", layer_filter));
    }
    else {
        lat.DrawLatex(0.16, 0.62, "Layer filter: none");
    }

    if (detector_filter == 12) {
        lat.DrawLatex(0.55, 0.78, "FTOF layers:");
        lat.DrawLatex(0.55, 0.73, "1 = FTOF1A");
        lat.DrawLatex(0.55, 0.68, "2 = FTOF1B");
        lat.DrawLatex(0.55, 0.63, "3 = FTOF2");
    }

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}

void plot_RECScintillator_energy_pip(const char* input_hipo,
                                     const char* output_png = "",
                                     int require_exactly_one_rec_pip = 1,
                                     int detector_filter = 12,
                                     int layer_filter = -1,
                                     int nBins = 100,
                                     double xMin = 0.0,
                                     double xMax = 20.0,
                                     bool logy = false)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECScintillator_energy_pip_" + base_name_no_ext(input_hipo);

        if (detector_filter > 0) outname += Form("_det%d", detector_filter);
        if (layer_filter > 0)    outname += Form("_layer%d", layer_filter);

        outname += ".png";
    }
    else {
        std::string s = outname;
        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recScint(factory.getSchema("REC::Scintillator"));

    std::string title = "REC::Scintillator energy for #pi^{+}-associated rows";

    if (detector_filter > 0) title += Form(", detector = %d", detector_filter);
    if (layer_filter > 0)    title += Form(", layer = %d", layer_filter);

    title += ";REC::Scintillator.energy;Counts";

    TH1D* h_energy = new TH1D(
        "h_scint_energy",
        title.c_str(),
        nBins, xMin, xMax
    );

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_scint_for_pip = 0;

    long long n_scint_rows_total = 0;
    long long n_scint_rows_for_pip = 0;
    long long n_scint_rows_after_detector_filter = 0;
    long long n_scint_rows_after_layer_filter = 0;
    long long n_energy_filled = 0;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recScint);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_scint_for_this_pip = false;

        const int n_scint = recScint.getRows();
        n_scint_rows_total += n_scint;

        for (int j = 0; j < n_scint; j++) {
            const int pidx = recScint.getShort("pindex", j);

            // REC::Scintillator row belongs to pi+ iff pindex equals pion REC::Particle row.
            if (pidx != pip_pindex) continue;

            found_scint_for_this_pip = true;
            n_scint_rows_for_pip++;

            const int det = recScint.getByte("detector", j);

            if (detector_filter > 0 && det != detector_filter) continue;
            n_scint_rows_after_detector_filter++;

            const int layer = recScint.getByte("layer", j);

            if (layer_filter > 0 && layer != layer_filter) continue;
            n_scint_rows_after_layer_filter++;

            const double energy = recScint.getFloat("energy", j);

            h_energy->Fill(energy);
            n_energy_filled++;
        }

        if (!found_scint_for_this_pip) {
            n_events_no_scint_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                              " << n_events << "\n"
              << "  Events with exactly one REC pi+:           " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:                " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:          " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                               " << n_events_used << "\n"
              << "  Total REC::Scintillator rows:              " << n_scint_rows_total << "\n"
              << "  Pion-associated REC::Scintillator rows:    " << n_scint_rows_for_pip << "\n"
              << "  Rows after detector filter:                " << n_scint_rows_after_detector_filter << "\n"
              << "  Rows after layer filter:                   " << n_scint_rows_after_layer_filter << "\n"
              << "  Filled energy entries:                     " << n_energy_filled << "\n"
              << "  Events with no pion REC::Scintillator row: " << n_events_no_scint_for_pip << "\n";

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c_scint_energy", "REC::Scintillator energy for pion", 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    if (logy) {
        c->SetLogy();
    }

    h_energy->SetLineWidth(2);
    h_energy->SetMinimum(0.0);
    h_energy->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(0.16, 0.82, Form("File: %s", base_name_no_ext(input_hipo).c_str()));
    lat.DrawLatex(0.16, 0.77, Form("Events used: %d", n_events_used));
    lat.DrawLatex(0.16, 0.72, Form("#pi^{+}-associated REC::Scintillator rows: %lld", n_scint_rows_for_pip));

    if (detector_filter > 0) lat.DrawLatex(0.16, 0.67, Form("Detector filter: %d", detector_filter));
    else                     lat.DrawLatex(0.16, 0.67, "Detector filter: none");

    if (layer_filter > 0) lat.DrawLatex(0.16, 0.62, Form("Layer filter: %d", layer_filter));
    else                  lat.DrawLatex(0.16, 0.62, "Layer filter: none");

    if (detector_filter == 12) {
        lat.DrawLatex(0.55, 0.78, "FTOF layers:");
        lat.DrawLatex(0.55, 0.73, "1 = FTOF1A");
        lat.DrawLatex(0.55, 0.68, "2 = FTOF1B");
        lat.DrawLatex(0.55, 0.63, "3 = FTOF2");
    }

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}


/////////////////////////////////////////////REC:: Calorimeter plots///////////////////////////////////////////

void plot_RECCalorimeter_detector_pip(const char* input_hipo,
                                      const char* output_png = "",
                                      int require_exactly_one_rec_pip = 1,
                                      int nBins = 31,
                                      double xMin = -0.5,
                                      double xMax = 30.5)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECCalorimeter_detector_pip_" + base_name_no_ext(input_hipo) + ".png";
    }
    else {
        std::string s = outname;
        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recCal(factory.getSchema("REC::Calorimeter"));

    TH1D* h_det = new TH1D(
        "h_cal_detector",
        "REC::Calorimeter detector for #pi^{+}-associated rows;REC::Calorimeter.detector;Counts",
        nBins, xMin, xMax
    );

    for (int det = int(xMin + 0.5); det <= int(xMax - 0.5); det++) {
        if (det < 0) continue;
        int bin = h_det->GetXaxis()->FindBin(det);
        h_det->GetXaxis()->SetBinLabel(bin, Form("%d", det));
    }

    h_det->GetXaxis()->LabelsOption("h");
    h_det->GetXaxis()->SetLabelSize(0.030);

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_cal_for_pip = 0;

    long long n_cal_rows_total = 0;
    long long n_cal_rows_for_pip = 0;

    std::map<int, long long> detector_counts;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recCal);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_cal_for_this_pip = false;

        const int n_cal = recCal.getRows();
        n_cal_rows_total += n_cal;

        for (int j = 0; j < n_cal; j++) {
            const int pidx = recCal.getShort("pindex", j);

            // Critical association:
            // REC::Calorimeter row belongs to pi+ iff pindex equals pion REC::Particle row.
            if (pidx != pip_pindex) continue;

            found_cal_for_this_pip = true;
            n_cal_rows_for_pip++;

            const int det = recCal.getByte("detector", j);

            h_det->Fill(det);
            detector_counts[det]++;
        }

        if (!found_cal_for_this_pip) {
            n_events_no_cal_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                              " << n_events << "\n"
              << "  Events with exactly one REC pi+:           " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:                " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:          " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                               " << n_events_used << "\n"
              << "  Total REC::Calorimeter rows:               " << n_cal_rows_total << "\n"
              << "  Pion-associated REC::Calorimeter rows:     " << n_cal_rows_for_pip << "\n"
              << "  Events with no pion REC::Calorimeter row:  " << n_events_no_cal_for_pip << "\n";

    std::cout << "\nREC::Calorimeter.detector counts for pion-associated rows:\n";

    for (const auto& kv : detector_counts) {
        std::cout << "  detector = " << std::setw(4) << kv.first
                  << "  count = " << std::setw(8) << kv.second;

        if (kv.first == 7) std::cout << "  ECAL";

        std::cout << "\n";
    }

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c_cal_detector", "REC::Calorimeter detector for pion", 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.14);

    h_det->SetLineWidth(2);
    h_det->SetMinimum(0.0);
    h_det->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(0.16, 0.82, Form("File: %s", base_name_no_ext(input_hipo).c_str()));
    lat.DrawLatex(0.16, 0.77, Form("Events used: %d", n_events_used));
    lat.DrawLatex(0.16, 0.72, Form("#pi^{+}-associated REC::Calorimeter rows: %lld", n_cal_rows_for_pip));

    lat.DrawLatex(0.55, 0.78, "Useful detector IDs:");
    lat.DrawLatex(0.55, 0.73, "7 = ECAL");

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}

void plot_RECCalorimeter_layer_pip(const char* input_hipo,
                                   const char* output_png = "",
                                   int require_exactly_one_rec_pip = 1,
                                   int detector_filter = 7,
                                   int nBins = 10,
                                   double xMin = 0.5,
                                   double xMax = 10.5)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECCalorimeter_layer_pip_" + base_name_no_ext(input_hipo);

        if (detector_filter > 0) outname += Form("_det%d", detector_filter);

        outname += ".png";
    }
    else {
        std::string s = outname;
        if (s.find("/") == std::string::npos) outname = "plots/" + s;
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recCal(factory.getSchema("REC::Calorimeter"));

    std::string title = "REC::Calorimeter layer for #pi^{+}-associated rows";

    if (detector_filter > 0) title += Form(", detector = %d", detector_filter);

    title += ";REC::Calorimeter.layer;Counts";

    TH1D* h_layer = new TH1D("h_cal_layer", title.c_str(), nBins, xMin, xMax);

    for (int layer = int(xMin + 0.5); layer <= int(xMax - 0.5); layer++) {
        int bin = h_layer->GetXaxis()->FindBin(layer);
        h_layer->GetXaxis()->SetBinLabel(bin, Form("%d", layer));
    }

    h_layer->GetXaxis()->LabelsOption("h");
    h_layer->GetXaxis()->SetLabelSize(0.035);

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_cal_for_pip = 0;

    long long n_cal_rows_total = 0;
    long long n_cal_rows_for_pip = 0;
    long long n_cal_rows_after_detector_filter = 0;

    std::map<int, long long> layer_counts;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recCal);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) n_events_with_one_pip++;

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_cal_for_this_pip = false;

        const int n_cal = recCal.getRows();
        n_cal_rows_total += n_cal;

        for (int j = 0; j < n_cal; j++) {
            const int pidx = recCal.getShort("pindex", j);

            // REC::Calorimeter row belongs to pi+ iff pindex equals pion REC::Particle row.
            if (pidx != pip_pindex) continue;

            found_cal_for_this_pip = true;
            n_cal_rows_for_pip++;

            const int det = recCal.getByte("detector", j);

            if (detector_filter > 0 && det != detector_filter) continue;

            const int layer = recCal.getByte("layer", j);

            h_layer->Fill(layer);
            layer_counts[layer]++;
            n_cal_rows_after_detector_filter++;
        }

        if (!found_cal_for_this_pip) n_events_no_cal_for_pip++;

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                              " << n_events << "\n"
              << "  Events with exactly one REC pi+:           " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:                " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:          " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                               " << n_events_used << "\n"
              << "  Total REC::Calorimeter rows:               " << n_cal_rows_total << "\n"
              << "  Pion-associated REC::Calorimeter rows:     " << n_cal_rows_for_pip << "\n"
              << "  Rows after detector filter:                " << n_cal_rows_after_detector_filter << "\n"
              << "  Events with no pion REC::Calorimeter row:  " << n_events_no_cal_for_pip << "\n";

    std::cout << "\nREC::Calorimeter.layer counts for pion-associated rows";
    if (detector_filter > 0) std::cout << " with detector = " << detector_filter;
    std::cout << ":\n";

    for (const auto& kv : layer_counts) {
        std::cout << "  layer = " << std::setw(4) << kv.first
                  << "  count = " << std::setw(8) << kv.second;

        if (detector_filter == 7) {
            if (kv.first == 1) std::cout << "  PCAL";
            if (kv.first == 4) std::cout << "  EC inner";
            if (kv.first == 7) std::cout << "  EC outer";
        }

        std::cout << "\n";
    }

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas("c_cal_layer", "REC::Calorimeter layer for pion", 900, 700);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.14);

    h_layer->SetLineWidth(2);
    h_layer->SetMinimum(0.0);
    h_layer->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(0.16, 0.82, Form("File: %s", base_name_no_ext(input_hipo).c_str()));
    lat.DrawLatex(0.16, 0.77, Form("Events used: %d", n_events_used));
    lat.DrawLatex(0.16, 0.72, Form("#pi^{+}-associated REC::Calorimeter rows: %lld", n_cal_rows_for_pip));

    if (detector_filter > 0) lat.DrawLatex(0.16, 0.67, Form("Detector filter: %d", detector_filter));
    else lat.DrawLatex(0.16, 0.67, "Detector filter: none");

    if (detector_filter == 7) {
        lat.DrawLatex(0.55, 0.78, "ECAL layers:");
        lat.DrawLatex(0.55, 0.73, "1 = PCAL");
        lat.DrawLatex(0.55, 0.68, "4 = EC inner");
        lat.DrawLatex(0.55, 0.63, "7 = EC outer");
    }

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}

void plot_deltaVz_pip(const char* input_hipo,
                      const char* output_png = "",
                      int require_exactly_one_rec_pip = 1,
                      int require_exactly_one_gen_pip = 1,
                      int nBins = 160,
                      double xMin = -20.0,
                      double xMax = 20.0,
                      bool logy = false)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/deltaVz_pip_" + base_name_no_ext(input_hipo) + ".png";
    }
    else {
        std::string s = outname;

        // If only a filename is supplied, save it inside plots/.
        // If a path is supplied, preserve that path.
        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank mcLund(factory.getSchema("MC::Lund"));

    TH1D* h_delta_vz = new TH1D(
        "h_delta_vz",
        "#Delta v_{z} for #pi^{+};"
        "#Delta v_{z} = v_{z,rec} - v_{z,gen} (cm);"
        "Counts",
        nBins, xMin, xMax
    );

    long long n_events = 0;
    long long n_events_used = 0;

    long long n_events_skipped_no_rec_pip = 0;
    long long n_events_skipped_multi_rec_pip = 0;

    long long n_events_skipped_no_gen_pip = 0;
    long long n_events_skipped_multi_gen_pip = 0;

    long long n_delta_vz_filled = 0;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(mcLund);

        n_events++;

        // ------------------------------------------------------------
        // Find reconstructed pi+
        // ------------------------------------------------------------
        int n_rec_pip = 0;
        int rec_pip_index = -1;

        for (int i = 0; i < recParticle.getRows(); i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                rec_pip_index = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_rec_pip++;
            continue;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_rec_pip++;
            continue;
        }

        // ------------------------------------------------------------
        // Find generated pi+
        // ------------------------------------------------------------
        int n_gen_pip = 0;
        int gen_pip_index = -1;

        for (int i = 0; i < mcLund.getRows(); i++) {
            const int pid = mcLund.getInt("pid", i);

            if (pid == 211) {
                n_gen_pip++;
                gen_pip_index = i;
            }
        }

        if (n_gen_pip == 0) {
            n_events_skipped_no_gen_pip++;
            continue;
        }

        if (require_exactly_one_gen_pip && n_gen_pip != 1) {
            n_events_skipped_multi_gen_pip++;
            continue;
        }

        // ------------------------------------------------------------
        // Calculate delta vz
        // ------------------------------------------------------------
        const double vz_rec = recParticle.getFloat("vz", rec_pip_index);
        const double vz_gen = mcLund.getFloat("vz", gen_pip_index);

        const double delta_vz = vz_rec - vz_gen;

        h_delta_vz->Fill(delta_vz);

        n_delta_vz_filled++;
        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                           " << n_events << "\n"
              << "  Events used:                           " << n_events_used << "\n"
              << "  Events skipped, no REC pi+:            " << n_events_skipped_no_rec_pip << "\n"
              << "  Events skipped, multiple REC pi+:      " << n_events_skipped_multi_rec_pip << "\n"
              << "  Events skipped, no MC::Lund pi+:       " << n_events_skipped_no_gen_pip << "\n"
              << "  Events skipped, multiple MC::Lund pi+: " << n_events_skipped_multi_gen_pip << "\n"
              << "  Filled delta-vz entries:               " << n_delta_vz_filled << "\n";

    std::cout << "\nDistribution summary:\n"
              << "  Mean delta vz:  " << h_delta_vz->GetMean() << " cm\n"
              << "  RMS delta vz:   " << h_delta_vz->GetRMS() << " cm\n";

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas(
        "c_delta_vz",
        "Delta vz for pion",
        900, 700
    );

    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    if (logy) {
        c->SetLogy();
    }

    h_delta_vz->SetLineWidth(2);
    h_delta_vz->SetMinimum(0.0);
    h_delta_vz->Draw("hist");

    TLine* zero = new TLine(
        0.0,
        0.0,
        0.0,
        h_delta_vz->GetMaximum() * 1.05
    );

    zero->SetLineStyle(2);
    zero->SetLineWidth(2);
    zero->Draw("same");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(
        0.16, 0.82,
        Form("File: %s", base_name_no_ext(input_hipo).c_str())
    );

    lat.DrawLatex(
        0.16, 0.77,
        Form("Events used: %lld", n_events_used)
    );

    lat.DrawLatex(
        0.16, 0.72,
        Form("Mean #Delta v_{z}: %.3f cm", h_delta_vz->GetMean())
    );

    lat.DrawLatex(
        0.16, 0.67,
        Form("RMS #Delta v_{z}: %.3f cm", h_delta_vz->GetRMS())
    );

    lat.DrawLatex(
        0.16, 0.62,
        "#Delta v_{z} = v_{z,rec}(#pi^{+}) - v_{z,gen}(#pi^{+})"
    );

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}

void plot_RECTrack_status_pip(const char* input_hipo,
                              const char* output_png = "",
                              int require_exactly_one_rec_pip = 1,
                              int detector_filter = -1,
                              int nBins = 21,
                              double xMin = -10.5,
                              double xMax = 10.5)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/RECTrack_status_pip_" + base_name_no_ext(input_hipo);

        if (detector_filter > 0) {
            outname += Form("_det%d", detector_filter);
        }

        outname += ".png";
    }
    else {
        std::string s = outname;

        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank recTrack(factory.getSchema("REC::Track"));

    std::string title = "REC::Track status for #pi^{+}-associated rows";

    if (detector_filter > 0) {
        title += Form(", detector = %d", detector_filter);
    }

    title += ";REC::Track.status;Counts";

    TH1D* h_status = new TH1D(
        "h_track_status",
        title.c_str(),
        nBins, xMin, xMax
    );

    for (int status = int(xMin + 0.5); status <= int(xMax - 0.5); status++) {
        const int bin = h_status->GetXaxis()->FindBin(status);
        h_status->GetXaxis()->SetBinLabel(bin, Form("%d", status));
    }

    h_status->GetXaxis()->LabelsOption("h");
    h_status->GetXaxis()->SetLabelSize(0.030);

    int n_events = 0;
    int n_events_with_one_pip = 0;
    int n_events_used = 0;
    int n_events_skipped_no_pip = 0;
    int n_events_skipped_multi_pip = 0;
    int n_events_no_track_for_pip = 0;

    long long n_track_rows_total = 0;
    long long n_track_rows_for_pip = 0;
    long long n_track_rows_after_detector_filter = 0;
    long long n_status_filled = 0;

    std::map<int, long long> status_counts;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(recTrack);

        n_events++;

        const int n_rec_particles = recParticle.getRows();

        int n_rec_pip = 0;
        int pip_pindex = -1;

        for (int i = 0; i < n_rec_particles; i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                pip_pindex = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_pip++;
            continue;
        }

        if (n_rec_pip == 1) {
            n_events_with_one_pip++;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_pip++;
            continue;
        }

        bool found_track_for_this_pip = false;

        const int n_tracks = recTrack.getRows();
        n_track_rows_total += n_tracks;

        for (int j = 0; j < n_tracks; j++) {
            const int pidx = recTrack.getShort("pindex", j);

            // REC::Track row belongs to the reconstructed pi+ iff
            // REC::Track.pindex equals the REC::Particle row index of the pi+.
            if (pidx != pip_pindex) continue;

            found_track_for_this_pip = true;
            n_track_rows_for_pip++;

            const int det = recTrack.getByte("detector", j);

            if (detector_filter > 0 && det != detector_filter) continue;
            n_track_rows_after_detector_filter++;

            const int status = recTrack.getShort("status", j);

            h_status->Fill(status);
            status_counts[status]++;
            n_status_filled++;
        }

        if (!found_track_for_this_pip) {
            n_events_no_track_for_pip++;
        }

        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                         " << n_events << "\n"
              << "  Events with exactly one REC pi+:     " << n_events_with_one_pip << "\n"
              << "  Events skipped, no REC pi+:          " << n_events_skipped_no_pip << "\n"
              << "  Events skipped, multiple REC pi+:    " << n_events_skipped_multi_pip << "\n"
              << "  Events used:                         " << n_events_used << "\n"
              << "  Total REC::Track rows:               " << n_track_rows_total << "\n"
              << "  Pion-associated REC::Track rows:     " << n_track_rows_for_pip << "\n"
              << "  Rows after detector filter:          " << n_track_rows_after_detector_filter << "\n"
              << "  Filled status entries:               " << n_status_filled << "\n"
              << "  Events with no pion REC::Track row:  " << n_events_no_track_for_pip << "\n";

    std::cout << "\nREC::Track.status counts for pion-associated rows";

    if (detector_filter > 0) {
        std::cout << " with detector = " << detector_filter;
    }

    std::cout << ":\n";

    for (const auto& kv : status_counts) {
        std::cout << "  status = " << std::setw(4) << kv.first
                  << "  count = " << std::setw(8) << kv.second
                  << "\n";
    }

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas(
        "c_track_status",
        "REC::Track status for pion",
        900, 700
    );

    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.14);

    h_status->SetLineWidth(2);
    h_status->SetMinimum(0.0);
    h_status->Draw("hist");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(
        0.16, 0.82,
        Form("File: %s", base_name_no_ext(input_hipo).c_str())
    );

    lat.DrawLatex(
        0.16, 0.77,
        Form("Events used: %d", n_events_used)
    );

    lat.DrawLatex(
        0.16, 0.72,
        Form("#pi^{+}-associated REC::Track rows: %lld", n_track_rows_for_pip)
    );

    if (detector_filter > 0) {
        lat.DrawLatex(
            0.16, 0.67,
            Form("Detector filter: %d", detector_filter)
        );
    }
    else {
        lat.DrawLatex(0.16, 0.67, "Detector filter: none");
    }

    lat.DrawLatex(0.58, 0.77, "REC::Track detector IDs:");
    lat.DrawLatex(0.58, 0.72, "5 = CVT / CD-like");
    lat.DrawLatex(0.58, 0.67, "6 = DC / FD-like");

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}

void plot_deltaTheta_pip(const char* input_hipo,
                         const char* output_png = "",
                         int require_exactly_one_rec_pip = 1,
                         int require_exactly_one_gen_pip = 1,
                         int nBins = 160,
                         double xMin = -10.0,
                         double xMax = 10.0,
                         bool logy = false)
{
    std::cout << "Input HIPO: " << input_hipo << "\n";

    gSystem->mkdir("plots", true);

    std::string outname = output_png;

    if (outname.empty()) {
        outname = "plots/deltaTheta_pip_" + base_name_no_ext(input_hipo) + ".png";
    }
    else {
        std::string s = outname;

        if (s.find("/") == std::string::npos) {
            outname = "plots/" + s;
        }
    }

    std::cout << "Output PNG: " << outname << "\n";

    hipo::reader reader;
    reader.open(input_hipo);

    hipo::dictionary factory;
    reader.readDictionary(factory);

    hipo::event event;

    hipo::bank recParticle(factory.getSchema("REC::Particle"));
    hipo::bank mcLund(factory.getSchema("MC::Lund"));

    TH1D* h_delta_theta = new TH1D(
        "h_delta_theta",
        "#Delta #theta for #pi^{+};"
        "#Delta #theta = #theta_{rec} - #theta_{gen} (deg);"
        "Counts",
        nBins, xMin, xMax
    );

    long long n_events = 0;
    long long n_events_used = 0;

    long long n_events_skipped_no_rec_pip = 0;
    long long n_events_skipped_multi_rec_pip = 0;

    long long n_events_skipped_no_gen_pip = 0;
    long long n_events_skipped_multi_gen_pip = 0;

    long long n_events_skipped_bad_rec_p = 0;
    long long n_events_skipped_bad_gen_p = 0;

    long long n_delta_theta_filled = 0;

    while (reader.next()) {
        reader.read(event);

        event.getStructure(recParticle);
        event.getStructure(mcLund);

        n_events++;

        // ------------------------------------------------------------
        // Find reconstructed pi+
        // ------------------------------------------------------------
        int n_rec_pip = 0;
        int rec_pip_index = -1;

        for (int i = 0; i < recParticle.getRows(); i++) {
            const int pid = recParticle.getInt("pid", i);

            if (pid == 211) {
                n_rec_pip++;
                rec_pip_index = i;
            }
        }

        if (n_rec_pip == 0) {
            n_events_skipped_no_rec_pip++;
            continue;
        }

        if (require_exactly_one_rec_pip && n_rec_pip != 1) {
            n_events_skipped_multi_rec_pip++;
            continue;
        }

        // ------------------------------------------------------------
        // Find generated pi+
        // ------------------------------------------------------------
        int n_gen_pip = 0;
        int gen_pip_index = -1;

        for (int i = 0; i < mcLund.getRows(); i++) {
            const int pid = mcLund.getInt("pid", i);

            if (pid == 211) {
                n_gen_pip++;
                gen_pip_index = i;
            }
        }

        if (n_gen_pip == 0) {
            n_events_skipped_no_gen_pip++;
            continue;
        }

        if (require_exactly_one_gen_pip && n_gen_pip != 1) {
            n_events_skipped_multi_gen_pip++;
            continue;
        }

        // ------------------------------------------------------------
        // Reconstructed theta
        // ------------------------------------------------------------
        const double px_rec = recParticle.getFloat("px", rec_pip_index);
        const double py_rec = recParticle.getFloat("py", rec_pip_index);
        const double pz_rec = recParticle.getFloat("pz", rec_pip_index);

        const double p_rec = std::sqrt(
            px_rec * px_rec +
            py_rec * py_rec +
            pz_rec * pz_rec
        );

        if (p_rec <= 0.0 || !std::isfinite(p_rec)) {
            n_events_skipped_bad_rec_p++;
            continue;
        }

        double cos_theta_rec = pz_rec / p_rec;

        if (cos_theta_rec > 1.0) cos_theta_rec = 1.0;
        if (cos_theta_rec < -1.0) cos_theta_rec = -1.0;

        const double theta_rec =
            std::acos(cos_theta_rec) * 180.0 / M_PI;

        // ------------------------------------------------------------
        // Generated theta
        // ------------------------------------------------------------
        const double px_gen = mcLund.getFloat("px", gen_pip_index);
        const double py_gen = mcLund.getFloat("py", gen_pip_index);
        const double pz_gen = mcLund.getFloat("pz", gen_pip_index);

        const double p_gen = std::sqrt(
            px_gen * px_gen +
            py_gen * py_gen +
            pz_gen * pz_gen
        );

        if (p_gen <= 0.0 || !std::isfinite(p_gen)) {
            n_events_skipped_bad_gen_p++;
            continue;
        }

        double cos_theta_gen = pz_gen / p_gen;

        if (cos_theta_gen > 1.0) cos_theta_gen = 1.0;
        if (cos_theta_gen < -1.0) cos_theta_gen = -1.0;

        const double theta_gen =
            std::acos(cos_theta_gen) * 180.0 / M_PI;

        // ------------------------------------------------------------
        // Delta theta
        // ------------------------------------------------------------
        const double delta_theta = theta_rec - theta_gen;

        h_delta_theta->Fill(delta_theta);

        n_delta_theta_filled++;
        n_events_used++;
    }

    std::cout << "\nCut flow / sanity check:\n"
              << "  Events read:                           " << n_events << "\n"
              << "  Events used:                           " << n_events_used << "\n"
              << "  Events skipped, no REC pi+:            " << n_events_skipped_no_rec_pip << "\n"
              << "  Events skipped, multiple REC pi+:      " << n_events_skipped_multi_rec_pip << "\n"
              << "  Events skipped, no MC::Lund pi+:       " << n_events_skipped_no_gen_pip << "\n"
              << "  Events skipped, multiple MC::Lund pi+: " << n_events_skipped_multi_gen_pip << "\n"
              << "  Events skipped, bad REC momentum:      " << n_events_skipped_bad_rec_p << "\n"
              << "  Events skipped, bad GEN momentum:      " << n_events_skipped_bad_gen_p << "\n"
              << "  Filled delta-theta entries:            " << n_delta_theta_filled << "\n";

    std::cout << "\nDistribution summary:\n"
              << "  Mean delta theta:  " << h_delta_theta->GetMean() << " deg\n"
              << "  RMS delta theta:   " << h_delta_theta->GetRMS() << " deg\n";

    gStyle->SetOptStat(1110);

    TCanvas* c = new TCanvas(
        "c_delta_theta",
        "Delta theta for pion",
        900, 700
    );

    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    if (logy) {
        c->SetLogy();
    }

    h_delta_theta->SetLineWidth(2);
    h_delta_theta->SetMinimum(0.0);
    h_delta_theta->Draw("hist");

    TLine* zero = new TLine(
        0.0,
        0.0,
        0.0,
        h_delta_theta->GetMaximum() * 1.05
    );

    zero->SetLineStyle(2);
    zero->SetLineWidth(2);
    zero->Draw("same");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.030);

    lat.DrawLatex(
        0.16, 0.82,
        Form("File: %s", base_name_no_ext(input_hipo).c_str())
    );

    lat.DrawLatex(
        0.16, 0.77,
        Form("Events used: %lld", n_events_used)
    );

    lat.DrawLatex(
        0.16, 0.72,
        Form("Mean #Delta #theta: %.4f deg", h_delta_theta->GetMean())
    );

    lat.DrawLatex(
        0.16, 0.67,
        Form("RMS #Delta #theta: %.4f deg", h_delta_theta->GetRMS())
    );

    lat.DrawLatex(
        0.16, 0.62,
        "#Delta #theta = #theta_{rec}(#pi^{+}) - #theta_{gen}(#pi^{+})"
    );

    c->SaveAs(outname.c_str());

    std::cout << "\nSaved: " << outname << "\n";
}