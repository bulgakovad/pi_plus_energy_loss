// clas12root -q -b hipo2root.c  --in=<filelist.dat> [--out=output.root]
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <TSystem.h>

#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TBenchmark.h>
#include "clas12reader.h"

using namespace clas12;

struct Args {
  TString inList;
  TString outRoot;   // optional
};

static Args parse_args() {
  Args a;
  for (Int_t i = 1; i < gApplication->Argc(); ++i) {
    TString opt = gApplication->Argv(i);
    if (opt.BeginsWith("--in=")) {
      a.inList = opt(5, opt.Length() - 5);     // everything after "--in="
    } else if (opt.BeginsWith("--out=")) {
      a.outRoot = opt(6, opt.Length() - 6);    // everything after "--out="
    } else if (opt.EndsWith(".dat") || opt.EndsWith(".txt")) {
      // allow bare list file as a convenience
      a.inList = opt;
    }
  }
  if (a.inList.IsNull()) {
    std::cerr << "Usage: clas12root -q -b hipo2root.c --in=<filelist.dat> [--out=output.root]\n";
    gSystem->Exit(1);
  }
  if (a.outRoot.IsNull()) {
    // default: <listBasename>_no_edge.root in CWD
    TString base = gSystem->BaseName(a.inList);
    base.ReplaceAll(".dat", "");
    base.ReplaceAll(".txt", "");
    a.outRoot = Form("../../data/%s.root", base.Data());
  }
  return a;
}

void ProcessHipo(const Args& args);

void hipo2root() {
  auto args = parse_args();
  ProcessHipo(args);
}

void ProcessHipo(const Args& args) {
  // --- Output file / tree
  TFile outFile(args.outRoot, "RECREATE");
  if (outFile.IsZombie()) {
    std::cerr << "ERROR: cannot create output ROOT file: " << args.outRoot << "\n";
    gSystem->Exit(2);
  }
  TTree out_tree("out_tree", "out_tree");

  // --- Branch variables
  float px_piplus_gen, py_piplus_gen, pz_piplus_gen, p_piplus_gen;
  float px_piplus_rec, py_piplus_rec, pz_piplus_rec, p_piplus_rec;
  float vx_piplus, vy_piplus, vz_piplus;
  int   pid_piplus, status_piplus, sector_piplus;

  float px_electron_gen, py_electron_gen, pz_electron_gen, p_electron_gen;
  float px_electron_rec, py_electron_rec, pz_electron_rec, p_electron_rec;
  int   pid_electron, status_electron;

  float edge1_electron, edge2_electron, edge3_electron;
  float edge1_piplus,   edge2_piplus,   edge3_piplus;

  float x1_piplus,   y1_piplus,   z1_piplus;
  float x1_electron, y1_electron, z1_electron;

  // --- Branches
  out_tree.Branch("px_piplus_gen", &px_piplus_gen);
  out_tree.Branch("py_piplus_gen", &py_piplus_gen);
  out_tree.Branch("pz_piplus_gen", &pz_piplus_gen);
  out_tree.Branch("px_piplus_rec", &px_piplus_rec);
  out_tree.Branch("py_piplus_rec", &py_piplus_rec);
  out_tree.Branch("pz_piplus_rec", &pz_piplus_rec);
  out_tree.Branch("p_piplus_gen", &p_piplus_gen);
  out_tree.Branch("p_piplus_rec", &p_piplus_rec);
  out_tree.Branch("vx_piplus", &vx_piplus);
  out_tree.Branch("vy_piplus", &vy_piplus);
  out_tree.Branch("vz_piplus", &vz_piplus);
  out_tree.Branch("pid_piplus", &pid_piplus);
  out_tree.Branch("status_piplus", &status_piplus);
  out_tree.Branch("sector_piplus", &sector_piplus);

  out_tree.Branch("px_electron_gen", &px_electron_gen);
  out_tree.Branch("py_electron_gen", &py_electron_gen);
  out_tree.Branch("pz_electron_gen", &pz_electron_gen);
  out_tree.Branch("p_electron_gen", &p_electron_gen);
  out_tree.Branch("px_electron_rec", &px_electron_rec);
  out_tree.Branch("py_electron_rec", &py_electron_rec);
  out_tree.Branch("pz_electron_rec", &pz_electron_rec);
  out_tree.Branch("p_electron_rec", &p_electron_rec);
  out_tree.Branch("pid_electron", &pid_electron);
  out_tree.Branch("status_electron", &status_electron);

  out_tree.Branch("edge1_electron", &edge1_electron);
  out_tree.Branch("edge2_electron", &edge2_electron);
  out_tree.Branch("edge3_electron", &edge3_electron);
  out_tree.Branch("edge1_piplus", &edge1_piplus);
  out_tree.Branch("edge2_piplus", &edge2_piplus);
  out_tree.Branch("edge3_piplus", &edge3_piplus);

  out_tree.Branch("x1_piplus", &x1_piplus);
  out_tree.Branch("y1_piplus", &y1_piplus);
  out_tree.Branch("z1_piplus", &z1_piplus);

  out_tree.Branch("x1_electron", &x1_electron);
  out_tree.Branch("y1_electron", &y1_electron);
  out_tree.Branch("z1_electron", &z1_electron);

  // --- Read file list
  std::ifstream flist(args.inList.Data());
  if (!flist.is_open()) {
    std::cerr << "ERROR: cannot open list file: " << args.inList << "\n";
    gSystem->Exit(3);
  }
  std::vector<std::string> data;
  for (std::string s; std::getline(flist, s);) if (!s.empty()) data.push_back(s);

  gBenchmark->Start("timer");
  std::cout << "Reading HIPO…\n";

  Long64_t total_events = 0;
  Long64_t total_events_check = 0;
  Long64_t events_kept = 0;

  for (const auto& filePath : data) {
    hipo::reader reader;
    reader.open(filePath.c_str());
    if (!reader.is_open()) {
      std::cerr << "WARNING: cannot open " << filePath << " (skipping)\n";
      continue;
    }
    hipo::dictionary dict; reader.readDictionary(dict);
    if (!dict.hasSchema("REC::Particle") || !dict.hasSchema("MC::Particle")) {
      std::cerr << "WARNING: required banks missing in " << filePath << " (skipping)\n";
      continue;
    }

    hipo::event event;
    hipo::bank  REC_particle(dict.getSchema("REC::Particle"));
    hipo::bank  MC_particle (dict.getSchema("MC::Particle"));
    hipo::bank  REC_track   (dict.getSchema("REC::Track"));
    hipo::bank  REC_traj    (dict.getSchema("REC::Traj"));

    Long64_t n_this_file = 0;

    while (reader.next()) {
      reader.read(event);               // <-- read every event (don’t skip the first)
      ++total_events;
      ++n_this_file;

      event.getStructure(REC_particle);
      event.getStructure(MC_particle);
      event.getStructure(REC_track);
      event.getStructure(REC_traj);

      const int Nrec = REC_particle.getRows();
      const int Nmc  = MC_particle.getRows();
      if (Nrec <= 0 || Nmc <= 0) continue;

      // --- pick one pi+/electron in REC
      int idx_piplus_rec = -1, idx_e_rec = -1;
      for (int i = 0; i < Nrec; ++i) {
        const int pid = REC_particle.getInt("pid", i);
        if (pid == 211 && idx_piplus_rec < 0) idx_piplus_rec = i;
        else if (pid == 11 && idx_e_rec < 0) idx_e_rec = i;
        if (idx_piplus_rec >= 0 && idx_e_rec >= 0) break;
      }
      if (idx_piplus_rec < 0 || idx_e_rec < 0) continue;

      // --- pick one pi+/electron in MC
      int idx_piplus_mc = -1, idx_e_mc = -1;
      for (int i = 0; i < Nmc; ++i) {
        const int pid = MC_particle.getInt("pid", i);
        if (pid == 211 && idx_piplus_mc < 0) idx_piplus_mc = i;
        else if (pid == 11 && idx_e_mc < 0) idx_e_mc = i;
        if (idx_piplus_mc >= 0 && idx_e_mc >= 0) break;
      }
      if (idx_piplus_mc < 0 || idx_e_mc < 0) continue;

      // --- reset all per-event sentinels (prevents carry-over)
      edge1_electron = edge2_electron = edge3_electron = -1.f;
      edge1_piplus   = edge2_piplus   = edge3_piplus   = -1.f;
      x1_piplus = y1_piplus = z1_piplus = -1000.f;
      x1_electron = y1_electron = z1_electron = -1000.f;
      sector_piplus = -1;

      // --- fill pi+
      px_piplus_gen = MC_particle.getFloat("px", idx_piplus_mc);
      py_piplus_gen = MC_particle.getFloat("py", idx_piplus_mc);
      pz_piplus_gen = MC_particle.getFloat("pz", idx_piplus_mc);
      p_piplus_gen = std::sqrt(px_piplus_gen*px_piplus_gen + py_piplus_gen*py_piplus_gen + pz_piplus_gen*pz_piplus_gen);

      px_piplus_rec = REC_particle.getFloat("px", idx_piplus_rec);
      py_piplus_rec = REC_particle.getFloat("py", idx_piplus_rec);
      pz_piplus_rec = REC_particle.getFloat("pz", idx_piplus_rec);
      p_piplus_rec = std::sqrt(px_piplus_rec*px_piplus_rec + py_piplus_rec*py_piplus_rec + pz_piplus_rec*pz_piplus_rec);

      vx_piplus = REC_particle.getFloat("vx", idx_piplus_rec);
      vy_piplus = REC_particle.getFloat("vy", idx_piplus_rec);
      vz_piplus = REC_particle.getFloat("vz", idx_piplus_rec);

      pid_piplus    = REC_particle.getInt("pid",    idx_piplus_rec);
      status_piplus = REC_particle.getInt("status", idx_piplus_rec);

      // sector from REC::Track
      for (int i = 0, Nt = REC_track.getRows(); i < Nt; ++i) {
        if (REC_track.getInt("pindex", i) == idx_piplus_rec) {
          sector_piplus = REC_track.getInt("sector", i);
          break;
        }
      }

      // --- fill electron
      px_electron_gen = MC_particle.getFloat("px", idx_e_mc);
      py_electron_gen = MC_particle.getFloat("py", idx_e_mc);
      pz_electron_gen = MC_particle.getFloat("pz", idx_e_mc);
      p_electron_gen = std::sqrt(px_electron_gen*px_electron_gen + py_electron_gen*py_electron_gen + pz_electron_gen*pz_electron_gen);

      px_electron_rec = REC_particle.getFloat("px", idx_e_rec);
      py_electron_rec = REC_particle.getFloat("py", idx_e_rec);
      pz_electron_rec = REC_particle.getFloat("pz", idx_e_rec);
      p_electron_rec = std::sqrt(px_electron_rec*px_electron_rec + py_electron_rec*py_electron_rec + pz_electron_rec*pz_electron_rec);

      pid_electron    = REC_particle.getInt("pid",    idx_e_rec);
      status_electron = REC_particle.getInt("status", idx_e_rec);

      // --- DC traj (detector 6), grab layer 6/18/36 and DC1 xyz
      for (int i = 0, Ntj = REC_traj.getRows(); i < Ntj; ++i) {
        if (REC_traj.getInt("detector", i) != 6) continue;
        const int pidx  = REC_traj.getInt("pindex", i);
        const int layer = REC_traj.getInt("layer",  i);
        const float edge = REC_traj.getFloat("edge", i);

        if (pidx == idx_e_rec) {
          if (layer == 6)  { edge1_electron = edge; x1_electron = REC_traj.getFloat("x", i); y1_electron = REC_traj.getFloat("y", i); z1_electron = REC_traj.getFloat("z", i); }
          if (layer == 18) edge2_electron = edge;
          if (layer == 36) edge3_electron = edge;
        } else if (pidx == idx_piplus_rec) {
          if (layer == 6)  { edge1_piplus   = edge; x1_piplus   = REC_traj.getFloat("x", i); y1_piplus   = REC_traj.getFloat("y", i); z1_piplus   = REC_traj.getFloat("z", i); }
          if (layer == 18) edge2_piplus   = edge;
          if (layer == 36) edge3_piplus   = edge;
        }
      }

      out_tree.Fill();
      ++events_kept;
    } // while events

    std::cout << filePath << " : " << n_this_file << " events\n";
    total_events_check += n_this_file;
  } // for files

  outFile.Write();
  outFile.Close();

  std::cout << "Wrote data into: " << args.outRoot << "\n";
  std::cout << "Total events (streamed)    : " << total_events << "\n";
  std::cout << "Total events (per-file sum): " << total_events_check << "\n";
  std::cout << "Events kept (rec+mc e&pi+) : " << events_kept << "\n";
  gBenchmark->Show("timer");
}
