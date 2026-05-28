// clas12root -q -b h2r_multypion.c --in=<filelist.dat> [--out=output.root]
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
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <TSystem.h>

#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TBenchmark.h>
#include <TFileMerger.h>
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
    std::cerr << "Usage: clas12root -q -b h2r_multypion.c --in=<filelist.dat> [--out=output.root]\n";
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

// ========================= helpers =========================
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
                           Long64_t& events_kept,
                           Long64_t& events_with_zero_rec_piplus,
                           Long64_t& events_with_one_rec_piplus,
                           Long64_t& events_with_gt1_rec_piplus,
                           Long64_t& total_rec_piplus_seen) {
  std::ifstream fin(infoPath);
  if (!fin.is_open()) return false;
  fin >> n_this_file
      >> events_kept
      >> events_with_zero_rec_piplus
      >> events_with_one_rec_piplus
      >> events_with_gt1_rec_piplus
      >> total_rec_piplus_seen;
  return static_cast<bool>(fin);
}
// ======================= END helpers =======================

void ProcessHipo(const Args& args);

// ========================= process ONE file in child =========================
static bool ProcessOneHipoFile(const std::string& filePath,
                               const std::string& outPath,
                               const std::string& infoPath) {
  TFile outFile(outPath.c_str(), "RECREATE");
  if (outFile.IsZombie()) {
    std::cerr << "WARNING: cannot create temp ROOT file for " << filePath << "\n";
    return false;
  }
  TTree out_tree("out_tree", "out_tree");

  // --- Legacy scalar branches (first REC pi+ only, for backward compatibility)
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

  // --- Multipion info per event
  int n_rec_piplus;

  std::vector<float> px_piplus_rec_all;
  std::vector<float> py_piplus_rec_all;
  std::vector<float> pz_piplus_rec_all;
  std::vector<float> p_piplus_rec_all;
  std::vector<float> vx_piplus_all;
  std::vector<float> vy_piplus_all;
  std::vector<float> vz_piplus_all;
  std::vector<int>   pid_piplus_all;
  std::vector<int>   status_piplus_all;
  std::vector<int>   sector_piplus_all;

  std::vector<float> edge1_piplus_all;
  std::vector<float> edge2_piplus_all;
  std::vector<float> edge3_piplus_all;

  std::vector<float> x1_piplus_all;
  std::vector<float> y1_piplus_all;
  std::vector<float> z1_piplus_all;

  // --- Legacy branches
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

  // --- New vector branches
  out_tree.Branch("n_rec_piplus", &n_rec_piplus);
  out_tree.Branch("px_piplus_rec_all", &px_piplus_rec_all);
  out_tree.Branch("py_piplus_rec_all", &py_piplus_rec_all);
  out_tree.Branch("pz_piplus_rec_all", &pz_piplus_rec_all);
  out_tree.Branch("p_piplus_rec_all",  &p_piplus_rec_all);
  out_tree.Branch("vx_piplus_all", &vx_piplus_all);
  out_tree.Branch("vy_piplus_all", &vy_piplus_all);
  out_tree.Branch("vz_piplus_all", &vz_piplus_all);
  out_tree.Branch("pid_piplus_all", &pid_piplus_all);
  out_tree.Branch("status_piplus_all", &status_piplus_all);
  out_tree.Branch("sector_piplus_all", &sector_piplus_all);
  out_tree.Branch("edge1_piplus_all", &edge1_piplus_all);
  out_tree.Branch("edge2_piplus_all", &edge2_piplus_all);
  out_tree.Branch("edge3_piplus_all", &edge3_piplus_all);
  out_tree.Branch("x1_piplus_all", &x1_piplus_all);
  out_tree.Branch("y1_piplus_all", &y1_piplus_all);
  out_tree.Branch("z1_piplus_all", &z1_piplus_all);

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
    std::cerr << "  REC::Particle=" << has_REC_Particle
              << "  MC::Particle="  << has_MC_Particle
              << "  REC::Track="    << has_REC_Track
              << "  REC::Traj="     << has_REC_Traj << "\n";
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

  // --- Multiplicity counters
  Long64_t events_with_zero_rec_piplus = 0;
  Long64_t events_with_one_rec_piplus  = 0;
  Long64_t events_with_gt1_rec_piplus  = 0;
  Long64_t total_rec_piplus_seen       = 0;

  while (reader.next()) {
    reader.read(event);
    ++n_this_file;

    event.getStructure(REC_particle);
    event.getStructure(MC_particle);
    event.getStructure(REC_track);
    event.getStructure(REC_traj);

    const int Nrec = REC_particle.getRows();
    const int Nmc  = MC_particle.getRows();

    // ========================= single source of truth for REC multiplicities =========================
    std::vector<int> rec_piplus_indices;
    rec_piplus_indices.reserve(std::max(0, Nrec));

    int idx_e_rec = -1;
    int n_trigger_electrons = 0;
    int n_rec_electrons_total = 0;

    for (int i = 0; i < Nrec; ++i) {
      const int pid    = REC_particle.getInt("pid", i);
      const int status = REC_particle.getInt("status", i);

      if (pid == 211) {
        rec_piplus_indices.push_back(i);
      }

      if (pid == 11) {
      ++n_rec_electrons_total;

      if (status < 0) {
        ++n_trigger_electrons;
        if (idx_e_rec < 0) idx_e_rec = i;
      }
}




    }

    n_rec_piplus = static_cast<int>(rec_piplus_indices.size());

    if (n_rec_piplus == 0) ++events_with_zero_rec_piplus;
    else if (n_rec_piplus == 1) ++events_with_one_rec_piplus;
    else ++events_with_gt1_rec_piplus;

    total_rec_piplus_seen += n_rec_piplus;
    // ======================= end single source of truth for REC multiplicities =======================

    // only after that do physics-selection cuts
    if (Nrec <= 0 || Nmc <= 0) continue;
    if (rec_piplus_indices.empty()) continue;
    // exactly one REC electron total, and it must be the trigger electron
    if (n_rec_electrons_total != 1) continue;
    if (n_trigger_electrons != 1 || idx_e_rec < 0) continue;

    // --- pick one pi+/electron in MC
    int idx_piplus_mc = -1, idx_e_mc = -1;
    for (int i = 0; i < Nmc; ++i) {
      const int pid = MC_particle.getInt("pid", i);
      if (pid == 211 && idx_piplus_mc < 0) idx_piplus_mc = i;
      else if (pid == 11 && idx_e_mc < 0) idx_e_mc = i;
      if (idx_piplus_mc >= 0 && idx_e_mc >= 0) break;
    }
    if (idx_piplus_mc < 0 || idx_e_mc < 0) continue;

    // --- reset legacy sentinels
    edge1_electron = edge2_electron = edge3_electron = -1.f;
    edge1_piplus   = edge2_piplus   = edge3_piplus   = -1.f;
    x1_piplus = y1_piplus = z1_piplus = -1000.f;
    x1_electron = y1_electron = z1_electron = -1000.f;
    sector_piplus = -1;

    // --- clear vector branches each event
    px_piplus_rec_all.clear();
    py_piplus_rec_all.clear();
    pz_piplus_rec_all.clear();
    p_piplus_rec_all.clear();
    vx_piplus_all.clear();
    vy_piplus_all.clear();
    vz_piplus_all.clear();
    pid_piplus_all.clear();
    status_piplus_all.clear();
    sector_piplus_all.clear();
    edge1_piplus_all.clear();
    edge2_piplus_all.clear();
    edge3_piplus_all.clear();
    x1_piplus_all.clear();
    y1_piplus_all.clear();
    z1_piplus_all.clear();

    // --- fill generated pi+
    px_piplus_gen = MC_particle.getFloat("px", idx_piplus_mc);
    py_piplus_gen = MC_particle.getFloat("py", idx_piplus_mc);
    pz_piplus_gen = MC_particle.getFloat("pz", idx_piplus_mc);
    p_piplus_gen = std::sqrt(px_piplus_gen*px_piplus_gen +
                             py_piplus_gen*py_piplus_gen +
                             pz_piplus_gen*pz_piplus_gen);

    // --- loop over ALL REC pi+
    for (size_t ip = 0; ip < rec_piplus_indices.size(); ++ip) {
      const int idx_piplus_rec = rec_piplus_indices[ip];

      const float px = REC_particle.getFloat("px", idx_piplus_rec);
      const float py = REC_particle.getFloat("py", idx_piplus_rec);
      const float pz = REC_particle.getFloat("pz", idx_piplus_rec);
      const float p  = std::sqrt(px*px + py*py + pz*pz);

      const float vx = REC_particle.getFloat("vx", idx_piplus_rec);
      const float vy = REC_particle.getFloat("vy", idx_piplus_rec);
      const float vz = REC_particle.getFloat("vz", idx_piplus_rec);

      const int pid    = REC_particle.getInt("pid",    idx_piplus_rec);
      const int status = REC_particle.getInt("status", idx_piplus_rec);

      int   sector = -1;
      float e1 = -1.f, e2 = -1.f, e3 = -1.f;
      float x1 = -1000.f, y1 = -1000.f, z1 = -1000.f;

      for (int i = 0, Nt = REC_track.getRows(); i < Nt; ++i) {
        if (REC_track.getInt("pindex", i) == idx_piplus_rec) {
          sector = REC_track.getInt("sector", i);
          break;
        }
      }

      for (int i = 0, Ntj = REC_traj.getRows(); i < Ntj; ++i) {
        if (REC_traj.getInt("detector", i) != 6) continue;
        if (REC_traj.getInt("pindex", i) != idx_piplus_rec) continue;

        const int layer = REC_traj.getInt("layer", i);
        const float edge = REC_traj.getFloat("edge", i);

        if (layer == 6) {
          e1 = edge;
          x1 = REC_traj.getFloat("x", i);
          y1 = REC_traj.getFloat("y", i);
          z1 = REC_traj.getFloat("z", i);
        }
        if (layer == 18) e2 = edge;
        if (layer == 36) e3 = edge;
      }

      px_piplus_rec_all.push_back(px);
      py_piplus_rec_all.push_back(py);
      pz_piplus_rec_all.push_back(pz);
      p_piplus_rec_all.push_back(p);
      vx_piplus_all.push_back(vx);
      vy_piplus_all.push_back(vy);
      vz_piplus_all.push_back(vz);
      pid_piplus_all.push_back(pid);
      status_piplus_all.push_back(status);
      sector_piplus_all.push_back(sector);

      edge1_piplus_all.push_back(e1);
      edge2_piplus_all.push_back(e2);
      edge3_piplus_all.push_back(e3);

      x1_piplus_all.push_back(x1);
      y1_piplus_all.push_back(y1);
      z1_piplus_all.push_back(z1);

      // keep old scalar branches equal to the FIRST reconstructed pi+
      if (ip == 0) {
        px_piplus_rec = px;
        py_piplus_rec = py;
        pz_piplus_rec = pz;
        p_piplus_rec  = p;
        vx_piplus = vx;
        vy_piplus = vy;
        vz_piplus = vz;
        pid_piplus = pid;
        status_piplus = status;
        sector_piplus = sector;
        edge1_piplus = e1;
        edge2_piplus = e2;
        edge3_piplus = e3;
        x1_piplus = x1;
        y1_piplus = y1;
        z1_piplus = z1;
      }
    }

    // --- fill electron (unique trigger electron: pid==11 and status<0)
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

    pid_electron    = REC_particle.getInt("pid", idx_e_rec);
    status_electron = REC_particle.getInt("status", idx_e_rec);

    // --- electron DC traj
    for (int i = 0, Ntj = REC_traj.getRows(); i < Ntj; ++i) {
      if (REC_traj.getInt("detector", i) != 6) continue;
      if (REC_traj.getInt("pindex", i) != idx_e_rec) continue;

      const int layer  = REC_traj.getInt("layer", i);
      const float edge = REC_traj.getFloat("edge", i);

      if (layer == 6) {
        edge1_electron = edge;
        x1_electron = REC_traj.getFloat("x", i);
        y1_electron = REC_traj.getFloat("y", i);
        z1_electron = REC_traj.getFloat("z", i);
      }
      if (layer == 18) edge2_electron = edge;
      if (layer == 36) edge3_electron = edge;
    }

    out_tree.Fill();
    ++events_kept;
  }

  outFile.Write();
  outFile.Close();

  std::ofstream info(infoPath.c_str());
  if (info.is_open()) {
    info << n_this_file << " "
         << events_kept << " "
         << events_with_zero_rec_piplus << " "
         << events_with_one_rec_piplus << " "
         << events_with_gt1_rec_piplus << " "
         << total_rec_piplus_seen << "\n";
  }

  std::cout << filePath << " : " << n_this_file << " events\n";
  return true;
}
// ======================= END process ONE file in child =======================

void h2r_multypion() {
  auto args = parse_args();
  ProcessHipo(args);
}

void ProcessHipo(const Args& args) {
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

  TString skippedLog = args.outRoot;
  if (skippedLog.EndsWith(".root")) skippedLog.ReplaceAll(".root", "_skipped.txt");
  else skippedLog += "_skipped.txt";

  std::ofstream skipped(skippedLog.Data());
  std::vector<std::string> goodTmpRoots;

  Long64_t total_events_check = 0;
  Long64_t events_kept_total  = 0;

  Long64_t total_events_with_zero_rec_piplus = 0;
  Long64_t total_events_with_one_rec_piplus  = 0;
  Long64_t total_events_with_gt1_rec_piplus  = 0;
  Long64_t total_rec_piplus_seen             = 0;

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
      _exit(ok ? 0 : 100);
    }

    int status = 0;
    waitpid(pid, &status, 0);

    if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
      Long64_t n_this_file = 0;
      Long64_t events_kept = 0;
      Long64_t events_with_zero_rec_piplus = 0;
      Long64_t events_with_one_rec_piplus  = 0;
      Long64_t events_with_gt1_rec_piplus  = 0;
      Long64_t rec_piplus_seen_this_file   = 0;

      if (read_info_file(tmpInfo.Data(),
                         n_this_file,
                         events_kept,
                         events_with_zero_rec_piplus,
                         events_with_one_rec_piplus,
                         events_with_gt1_rec_piplus,
                         rec_piplus_seen_this_file)) {
        total_events_check += n_this_file;
        events_kept_total  += events_kept;
        total_events_with_zero_rec_piplus += events_with_zero_rec_piplus;
        total_events_with_one_rec_piplus  += events_with_one_rec_piplus;
        total_events_with_gt1_rec_piplus  += events_with_gt1_rec_piplus;
        total_rec_piplus_seen             += rec_piplus_seen_this_file;
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
  std::cout << "Events kept (exactly 1 REC electron total, trigger e + >=1 rec pi+ + mc e/pi+): " << events_kept_total << "\n";
  std::cout << "Files OK                   : " << files_ok << "\n";
  std::cout << "Files skipped              : " << files_skipped << "\n";
  std::cout << "Skipped-file log           : " << skippedLog << "\n";
  std::cout << "REC pi+ multiplicity summary over decoded events:\n";
  std::cout << "  Events with exactly 1 REC pi+  : " << total_events_with_one_rec_piplus << "\n";
  std::cout << "  Events with >1 REC pi+         : " << total_events_with_gt1_rec_piplus << "\n";
  std::cout << "  Events with 0 REC pi+          : " << total_events_with_zero_rec_piplus << "\n";
  std::cout << "  Total REC pi+ seen             : " << total_rec_piplus_seen << "\n";
  gBenchmark->Show("timer");
}