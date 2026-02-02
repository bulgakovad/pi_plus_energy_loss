#!/usr/bin/env python3
import os
import numpy as np
import multiprocessing as mp
import ROOT

# --- constants (GeV, radians) ---
M_E  = 0.000511
M_PI = 0.13957039
DEG2RAD = np.pi / 180.0

def _seed_this_process():
    """
    Unique RNG seeds per worker so parallel jobs don't clone events.
    Uses PID + time so seeds differ even if processes respawn quickly.
    """
    pid = os.getpid()
    # ROOT seed expects an integer; mix pid with nanoseconds
    seed = (pid ^ (time_ns() & 0xFFFFFFFF)) & 0x7FFFFFFF
    ROOT.gRandom.SetSeed(int(seed))
    np.random.seed(int(seed) & 0xFFFFFFFF)

def time_ns():
    # lightweight local helper (avoids importing time globally if you prefer)
    import time
    return time.time_ns()

def generate_lund_file(file_index, output_folder, num_events_per_file, Ebeam):
    path = os.path.join(output_folder, f"piplus_electron_lund_{file_index}.lund")
    lines = []

    for _ in range(num_events_per_file):
        # ---- header: 10 fields ----
        # Npart  A  Z  tgtPol  spinZ  beamPID  Ebeam  struckNucleon  processID  weight
        Npart = 2
        A, Z = 1, 1
        tgtPol = 0.0
        spinZ  = 0.0
        beamPID = 11
        struck = 2212
        processID = 0
        weight = ROOT.gRandom.Uniform(0.2, 2.0)

        lines.append(f"{Npart} {A} {Z} {tgtPol} {spinZ} {beamPID} {Ebeam:.6f} {struck} {processID} {weight}\n")

        # ---- trigger electron (random toy kinematics) ----
        # keep it strictly below Ebeam to avoid edge cases
        p_e  = ROOT.gRandom.Uniform(1.5, max(1.501, Ebeam - 1e-3))
        th_e = ROOT.gRandom.Uniform(5.0, 25.0) * DEG2RAD
        ph_e = ROOT.gRandom.Uniform(-180.0, 180.0) * DEG2RAD
        px_e = p_e * np.sin(th_e) * np.cos(ph_e)
        py_e = p_e * np.sin(th_e) * np.sin(ph_e)
        pz_e = p_e * np.cos(th_e)
        E_e  = np.sqrt(px_e*px_e + py_e*py_e + pz_e*pz_e + M_E*M_E)

        # shared vertex (cm) — vx=vy=0, vz smeared
        vz = np.random.normal(-3.0, 2.5)

        # index  lifetime(ns)  type  pid   parent firstDau  px        py        pz        E         m        vx       vy       vz
        lines.append(
            f"1 -1 1 11   0      0        "
            f"{px_e:.6f} {py_e:.6f} {pz_e:.6f} {E_e:.6f} {M_E:.6f} "
            f"0.000000 0.000000 {vz:.6f}\n"
        )

        # ---- pi+ (uncorrelated toy hadron) ----
        p_pi  = ROOT.gRandom.Uniform(0.1, max(0.101, Ebeam - 1e-3))
        th_pi = ROOT.gRandom.Uniform(3.0, 45.0) * DEG2RAD
        ph_pi = ROOT.gRandom.Uniform(-180.0, 180.0) * DEG2RAD
        px_pi = p_pi * np.sin(th_pi) * np.cos(ph_pi)
        py_pi = p_pi * np.sin(th_pi) * np.sin(ph_pi)
        pz_pi = p_pi * np.cos(th_pi)
        E_pi  = np.sqrt(px_pi*px_pi + py_pi*py_pi + pz_pi*pz_pi + M_PI*M_PI)

        lines.append(
            f"2 -1 1 211  0      0        "
            f"{px_pi:.6f} {py_pi:.6f} {pz_pi:.6f} {E_pi:.6f} {M_PI:.6f} "
            f"0.000000 0.000000 {vz:.6f}\n"
        )

    with open(path, "w") as f:
        f.writelines(lines)

def generate_lund_files(output_folder, num_files, num_events_per_file, Ebeam):
    os.makedirs(output_folder, exist_ok=True)

    tasks = [(i, output_folder, num_events_per_file, Ebeam) for i in range(1, num_files + 1)]

    # Seed once per worker (IMPORTANT): avoids identical sequences for multiple files in the same worker
    with mp.Pool(processes=mp.cpu_count(), initializer=_seed_this_process) as pool:
        pool.starmap(generate_lund_file, tasks)

    print(f"Generated {num_files} LUND files with {num_events_per_file} events each at {output_folder}")

if __name__ == "__main__":
    outdir = "/lustre24/expphy/volatile/clas12/bulgakov/lund_files/piplus_electron_lund_multithread_2"

    # RGK common beam energies: 6.535 or 7.546 (pick what you need)
    Ebeam = 6.535

    num_files = 1000
    num_events_per_file = 5000
    generate_lund_files(outdir, num_files, num_events_per_file, Ebeam)
