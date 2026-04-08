// clas12root -q -b h2r_fixed.c --in=<filelist.dat> [--out=output.root]
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <initializer_list>
#include <cctype>
#include <sys/types.h>   // CHANGED: for fork/wait
#include <sys/wait.h>    // CHANGED: for fork/wait
#include <unistd.h>      // CHANGED: for fork/getpid/_exit
#include <TSystem.h>

#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TBenchmark.h>
#include <TFileMerger.h> // CHANGED: merge temp ROOT files
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
      a.inList = opt(5, opt.Length() - 5);
    } else if (opt.BeginsWith("--out=")) {
      a.outRoot = opt(6, opt.Length() - 6);
    } else if (opt.EndsWith(".dat") || opt.EndsWith(".txt")) {
      a.inList = opt;
    }
  }

  if (a.inList.IsNull()) {
    std::cerr << "Usage: clas12root -q -b hipo2root.c --in=<filelist.dat> [--out=output.root]\n";
    gSystem->Exit(1);
  }

  if (a.outRoot.IsNull()) {
    TString base = gSystem->BaseName(a.inList);
    base.ReplaceAll(".dat", "");
    base.ReplaceAll(".txt", "");
    a.outRoot = Form("../../data/%s.root", base.Data());
  }

  return a;
}

// ========================= CHANGED: helpers =========================
static void trim_in_place(std::string& s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
    [](unsigned char ch) { return !std::isspace(ch); }));
  s.erase(std::find_if(s.rbegin(), s.rend(),
    [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}

static bool has_entry(const hipo::schema& sch, const char* name) {
  return sch.getEntryOrder(name) >= 0;
}

static bool schema_has_all_entries(const hipo::schema& sch,
                                   std::initializer_list<const char*> names) {
  for (const char* name : names) {
    if (!has_entry(sch, name)) return false;
  }
  return true;
}

static bool read_info_file(const std::string& infoPath,
                           Long64_t& n_this_file,
                           Long64_t& events_kept) {
  std::ifstream fin(infoPath);
  if (!fin.is_open()) return false;
  fin >> n_this_file >> events_kept;
  return static_cast<bool>(fin);
}
// ======================= END CHANGED: helpers =======================

void ProcessHipo(const Args& args);

// ========================= CHANGED: process ONE file in child =========================
static bool ProcessOneHipoFile(const std::string& filePath,
                               const std::string& outPath,
                               const std::string& infoPath) {
  // --- Output file / tree for this one HIPO file
  TFile outFile(outPath.c_str(), "RECREATE");
  if (outFile.IsZombie()) {
    std::cerr << "WARNING: cannot create temp ROOT file for " << filePath << "\n";
    return false;
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

  // --- Open HIPO
  hipo::reader reader;
  reader.open(filePath.c_str());
  if (!reader.is_open()) {
    std::cerr << "WARNING: cannot open " << filePath << " (skipping)\n";
    return false;
  }

  hipo::dictionary dict;
  reader.readDictionary(dict);

  // --- Guard required banks
  const bool has_REC_Particle = dict.hasSchema("REC::Particle");
  const bool has_MC_Particle  = dict.hasSchema("MC::Particle");
  const bool has_REC_Track    = dict.hasSchema("REC::Track");
  const bool has_REC_Traj     = dict.hasSchema("REC::Traj");

  if (!has_REC_Particle || !has_MC_Particle || !has_REC_Track || !has_REC_Traj) {
    std::cerr << "WARNING: required banks missing in " << filePath << " (skipping)\n";
    std::cerr << "         REC::Particle=" << has_REC_Particle
              << "  MC::Particle="        << has_MC_Particle
              << "  REC::Track="          << has_REC_Track
              << "  REC::Traj="           << has_REC_Traj << "\n";
    return false;
  }

  const auto& sch_REC_particle = dict.getSchema("REC::Particle");
  const auto& sch_MC_particle  = dict.getSchema("MC::Particle");
  const auto& sch_REC_track    = dict.getSchema("REC::Track");
  const auto& sch_REC_traj     = dict.getSchema("REC::Traj");

  // --- Guard required fields
  const bool rec_particle_ok = schema_has_all_entries(
    sch_REC_particle, {"pid","px","py","pz","vx","vy","vz","status"}
  );
  const bool mc_particle_ok = schema_has_all_entries(
    sch_MC_particle, {"pid","px","py","pz"}
  );
  const bool rec_track_ok = schema_has_all_entries(
    sch_REC_track, {"pindex","sector"}
  );
  const bool rec_traj_ok = schema_has_all_entries(
    sch_REC_traj, {"detector","pindex","layer","edge","x","y","z"}
  );

  if (!rec_particle_ok || !mc_particle_ok || !rec_track_ok || !rec_traj_ok) {
    std::cerr << "WARNING: required bank fields missing in " << filePath << " (skipping)\n";
    std::cerr << "         REC::Particle fields ok = " << rec_particle_ok << "\n";
    std::cerr << "         MC::Particle  fields ok = " << mc_particle_ok  << "\n";
    std::cerr << "         REC::Track    fields ok = " << rec_track_ok    << "\n";
    std::cerr << "         REC::Traj     fields ok = " << rec_traj_ok     << "\n";
    return false;
  }

  hipo::event event;
  hipo::bank  REC_particle(sch_REC_particle);
  hipo::bank  MC_particle (sch_MC_particle);
  hipo::bank  REC_track   (sch_REC_track);
  hipo::bank  REC_traj    (sch_REC_traj);

  Long64_t n_this_file = 0;
  Long64_t events_kept = 0;

  // --- SAME event logic as before
  while (reader.next()) {
    reader.read(event);   // if a bad event aborts, only child dies
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

    // --- reset all per-event sentinels
    edge1_electron = edge2_electron = edge3_electron = -1.f;
    edge1_piplus   = edge2_piplus   = edge3_piplus   = -1.f;
    x1_piplus = y1_piplus = z1_piplus = -1000.f;
    x1_electron = y1_electron = z1_electron = -1000.f;
    sector_piplus = -1;

    // --- fill pi+
    px_piplus_gen = MC_particle.getFloat("px", idx_piplus_mc);
    py_piplus_gen = MC_particle.getFloat("py", idx_piplus_mc);
    pz_piplus_gen = MC_particle.getFloat("pz", idx_piplus_mc);
    p_piplus_gen = std::sqrt(px_piplus_gen*px_piplus_gen +
                             py_piplus_gen*py_piplus_gen +
                             pz_piplus_gen*pz_piplus_gen);

    px_piplus_rec = REC_particle.getFloat("px", idx_piplus_rec);
    py_piplus_rec = REC_particle.getFloat("py", idx_piplus_rec);
    pz_piplus_rec = REC_particle.getFloat("pz", idx_piplus_rec);
    p_piplus_rec = std::sqrt(px_piplus_rec*px_piplus_rec +
                             py_piplus_rec*py_piplus_rec +
                             pz_piplus_rec*pz_piplus_rec);

    vx_piplus = REC_particle.getFloat("vx", idx_piplus_rec);
    vy_piplus = REC_particle.getFloat("vy", idx_piplus_rec);
    vz_piplus = REC_particle.getFloat("vz", idx_piplus_rec);

    pid_piplus    = REC_particle.getInt("pid",    idx_piplus_rec);
    status_piplus = REC_particle.getInt("status", idx_piplus_rec);

    // --- sector from REC::Track
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
    p_electron_gen = std::sqrt(px_electron_gen*px_electron_gen +
                               py_electron_gen*py_electron_gen +
                               pz_electron_gen*pz_electron_gen);

    px_electron_rec = REC_particle.getFloat("px", idx_e_rec);
    py_electron_rec = REC_particle.getFloat("py", idx_e_rec);
    pz_electron_rec = REC_particle.getFloat("pz", idx_e_rec);
    p_electron_rec = std::sqrt(px_electron_rec*px_electron_rec +
                               py_electron_rec*py_electron_rec +
                               pz_electron_rec*pz_electron_rec);

    pid_electron    = REC_particle.getInt("pid",    idx_e_rec);
    status_electron = REC_particle.getInt("status", idx_e_rec);

    // --- DC traj
    for (int i = 0, Ntj = REC_traj.getRows(); i < Ntj; ++i) {
      if (REC_traj.getInt("detector", i) != 6) continue;
      const int pidx   = REC_traj.getInt("pindex", i);
      const int layer  = REC_traj.getInt("layer",  i);
      const float edge = REC_traj.getFloat("edge", i);

      if (pidx == idx_e_rec) {
        if (layer == 6)  {
          edge1_electron = edge;
          x1_electron = REC_traj.getFloat("x", i);
          y1_electron = REC_traj.getFloat("y", i);
          z1_electron = REC_traj.getFloat("z", i);
        }
        if (layer == 18) edge2_electron = edge;
        if (layer == 36) edge3_electron = edge;
      } else if (pidx == idx_piplus_rec) {
        if (layer == 6)  {
          edge1_piplus = edge;
          x1_piplus = REC_traj.getFloat("x", i);
          y1_piplus = REC_traj.getFloat("y", i);
          z1_piplus = REC_traj.getFloat("z", i);
        }
        if (layer == 18) edge2_piplus = edge;
        if (layer == 36) edge3_piplus = edge;
      }
    }

    out_tree.Fill();
    ++events_kept;
  }

  outFile.Write();
  outFile.Close();

  std::ofstream info(infoPath.c_str());
  if (info.is_open()) {
    info << n_this_file << " " << events_kept << "\n";
  }

  std::cout << filePath << " : " << n_this_file << " events\n";
  return true;
}
// ======================= END CHANGED: process ONE file in child =======================

void h2r_fixed() {
  auto args = parse_args();
  ProcessHipo(args);
}

void ProcessHipo(const Args& args) {
  // --- Read file list
  std::ifstream flist(args.inList.Data());
  if (!flist.is_open()) {
    std::cerr << "ERROR: cannot open list file: " << args.inList << "\n";
    gSystem->Exit(3);
  }

  std::vector<std::string> data;
  for (std::string s; std::getline(flist, s); ) {
    trim_in_place(s);
    if (!s.empty()) data.push_back(s);
  }

  gBenchmark->Start("timer");
  std::cout << "Reading HIPO...\n";

  // CHANGED: parent now supervises child processes and skips bad files
  TString skippedLog = args.outRoot;
  if (skippedLog.EndsWith(".root")) skippedLog.ReplaceAll(".root", "_skipped.txt");
  else skippedLog += "_skipped.txt";

  std::ofstream skipped(skippedLog.Data());
  std::vector<std::string> goodTmpRoots;

  Long64_t total_events_check = 0;
  Long64_t events_kept_total  = 0;
  int files_ok = 0;
  int files_skipped = 0;

  TString tempDir = gSystem->TempDirectory();
  const int parent_pid = static_cast<int>(getpid());

  for (size_t i = 0; i < data.size(); ++i) {
    const std::string& filePath = data[i];
    std::cout << "Opening: " << filePath << std::endl;

    TString tmpRoot = Form("%s/hipo2root_tmp_%d_%06zu.root",
                           tempDir.Data(), parent_pid, i);
    TString tmpInfo = Form("%s/hipo2root_tmp_%d_%06zu.info",
                           tempDir.Data(), parent_pid, i);

    pid_t pid = fork();

    if (pid < 0) {
      std::cerr << "WARNING: fork() failed for " << filePath << " (skipping)\n";
      if (skipped.is_open()) skipped << filePath << "    fork_failed\n";
      ++files_skipped;
      gSystem->Unlink(tmpRoot.Data());
      gSystem->Unlink(tmpInfo.Data());
      continue;
    }

    if (pid == 0) {
      bool ok = ProcessOneHipoFile(filePath, tmpRoot.Data(), tmpInfo.Data());
      _exit(ok ? 0 : 100); // child exits cleanly; parent continues
    }

    int status = 0;
    waitpid(pid, &status, 0);

    if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
      Long64_t n_this_file = 0;
      Long64_t events_kept = 0;
      if (read_info_file(tmpInfo.Data(), n_this_file, events_kept)) {
        total_events_check += n_this_file;
        events_kept_total  += events_kept;
      }
      goodTmpRoots.push_back(tmpRoot.Data());
      ++files_ok;
    } else {
      ++files_skipped;
      if (skipped.is_open()) {
        if (WIFSIGNALED(status)) {
          skipped << filePath << "    signal=" << WTERMSIG(status) << "\n";
        } else if (WIFEXITED(status)) {
          skipped << filePath << "    exit_code=" << WEXITSTATUS(status) << "\n";
        } else {
          skipped << filePath << "    unknown_failure\n";
        }
      }
      std::cerr << "WARNING: skipping bad file " << filePath << "\n";
      gSystem->Unlink(tmpRoot.Data());
    }

    gSystem->Unlink(tmpInfo.Data());
  }

  // --- Merge successful temp ROOT files
  if (goodTmpRoots.empty()) {
    std::cerr << "ERROR: no good files were processed.\n";
    TFile outFile(args.outRoot, "RECREATE");
    outFile.Close();
    std::cout << "Created empty output: " << args.outRoot << "\n";
    std::cout << "Skipped-file log: " << skippedLog << "\n";
    gBenchmark->Show("timer");
    return;
  }

  TFileMerger merger;
  merger.OutputFile(args.outRoot.Data(), "RECREATE");

  for (const auto& tmp : goodTmpRoots) {
    merger.AddFile(tmp.c_str());
  }

  Bool_t merge_ok = merger.Merge();
  if (!merge_ok) {
    std::cerr << "ERROR: merge failed for temp ROOT files.\n";
    gSystem->Exit(4);
  }

  for (const auto& tmp : goodTmpRoots) {
    gSystem->Unlink(tmp.c_str());
  }

  std::cout << "Wrote data into: " << args.outRoot << "\n";
  std::cout << "Total events (per-file sum): " << total_events_check << "\n";
  std::cout << "Events kept (rec+mc e&pi+) : " << events_kept_total << "\n";
  std::cout << "Files OK                   : " << files_ok << "\n";
  std::cout << "Files skipped              : " << files_skipped << "\n";
  std::cout << "Skipped-file log           : " << skippedLog << "\n";
  gBenchmark->Show("timer");
}