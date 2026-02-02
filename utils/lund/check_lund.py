#!/usr/bin/env python3
import ROOT
import os
import sys
import glob
import math

# ----------------- Constants -----------------
TARGET_MASS = 0.9382720813  # GeV (proton target mass for x_Bj)
DEFAULT_BEAM_ENERGY = 6.535  # GeV (fallback if header doesn't provide Ebeam)

# ----------------- Globals -------------------
event_count = 0

# ----------------- Histograms ----------------
h2_q2_xbj = ROOT.TH2D(
    "h2_q2_xbj",
    "Q^{2} vs x_{Bj};x_{Bj};Q^{2} [GeV^{2}]",
    100, 0.0, 3.0, 100, 0.0, 12.0
)

h1_q2 = ROOT.TH1D(
    "h1_q2",
    "Q^{2} Distribution;Q^{2} [GeV^{2}];Counts",
    100, 0.0, 15.0
)

h1_xbj = ROOT.TH1D(
    "h1_xbj",
    "x_{Bj} Distribution;x_{Bj};Counts",
    100, 0.0, 3.0
)

# Q2 vs pi+ momentum magnitude
h2_q2_pip = ROOT.TH2D(
    "h2_q2_pip",
    "Q^{2} vs #pi^{+} Momentum Magnitude;|p_{#pi^{+}}| [GeV];Q^{2} [GeV^{2}]",
    100, 0.0, 10.0, 100, 0.0, 12.0
)

# --- Electron histograms ---
h2_e_theta_p = ROOT.TH2D(
    "h2_e_theta_p",
    "Electron #theta vs |p|;|p_{e}| [GeV];#theta_{e} [deg]",
    100, 0, 10, 100, 0, 180
)
h2_e_phi_p = ROOT.TH2D(
    "h2_e_phi_p",
    "Electron #phi vs |p|;|p_{e}| [GeV];#phi_{e} [deg]",
    100, 0, 10, 100, -180, 180
)
h2_e_theta_phi = ROOT.TH2D(
    "h2_e_theta_phi",
    "Electron #theta vs #phi;#phi_{e} [deg];#theta_{e} [deg]",
    100, -180, 180, 100, 0, 180
)

# --- Pi+ histograms ---
h2_pip_theta_p = ROOT.TH2D(
    "h2_pip_theta_p",
    "#pi^{+} #theta vs |p|;|p_{#pi^{+}}| [GeV];#theta_{#pi^{+}} [deg]",
    100, 0, 10, 100, 0, 180
)
h2_pip_phi_p = ROOT.TH2D(
    "h2_pip_phi_p",
    "#pi^{+} #phi vs |p|;|p_{#pi^{+}}| [GeV];#phi_{#pi^{+}} [deg]",
    100, 0, 10, 100, -180, 180
)
h2_pip_theta_phi = ROOT.TH2D(
    "h2_pip_theta_phi",
    "#pi^{+} #theta vs #phi;#phi_{#pi^{+}} [deg];#theta_{#pi^{+}} [deg]",
    100, -180, 180, 100, 0, 180
)


def _parse_int(token, default=None):
    try:
        return int(token)
    except Exception:
        return default


def _parse_float(token, default=None):
    try:
        return float(token)
    except Exception:
        return default


def process_file(filepath):
    """
    Read a single LUND file and fill histograms.
    Assumes event header line followed by exactly N particle lines.
    Blank lines are ignored.
    """
    global event_count

    with open(filepath, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 0
    n_lines = len(lines)
    while i < n_lines:
        # --- Header ---
        header_tokens = lines[i].split()
        n_particles = _parse_int(header_tokens[0], default=None)
        if n_particles is None:
            i += 1
            continue

        # Ebeam is token[6] in your generator header:
        # Npart A Z tgtPol spinZ beamPID Ebeam struckNucleon processID weight
        ebeam = DEFAULT_BEAM_ENERGY
        if len(header_tokens) >= 7:
            ebeam = _parse_float(header_tokens[6], default=DEFAULT_BEAM_ENERGY)

        beam_vec = ROOT.TLorentzVector(0.0, 0.0, ebeam, ebeam)

        # Sanity: ensure we have enough lines left for the event
        if i + 1 + n_particles > n_lines:
            break

        # --- Particle block ---
        event_lines = lines[i + 1: i + 1 + n_particles]
        i += 1 + n_particles

        scattered_electron = None
        pip_momentum_mag = None

        for line in event_lines:
            parts = line.split()
            if len(parts) < 14:
                continue

            # LUND particle fields (indices):
            # 0:index 1:lifetime 2:type 3:pid 4:parent 5:firstDau 6:px 7:py 8:pz 9:E 10:m 11:vx 12:vy 13:vz
            pid = _parse_int(parts[3], default=0)

            px = _parse_float(parts[6], default=0.0)
            py = _parse_float(parts[7], default=0.0)
            pz = _parse_float(parts[8], default=0.0)
            E  = _parse_float(parts[9], default=0.0)

            p_mag = math.sqrt(px*px + py*py + pz*pz)
            if p_mag > 0.0:
                cos_th = max(-1.0, min(1.0, pz / p_mag))
                theta = math.degrees(math.acos(cos_th))
            else:
                theta = 0.0
            phi = math.degrees(math.atan2(py, px))

            if pid == 11:
                scattered_electron = ROOT.TLorentzVector(px, py, pz, E)
                h2_e_theta_p.Fill(p_mag, theta)
                h2_e_phi_p.Fill(p_mag, phi)
                h2_e_theta_phi.Fill(phi, theta)

            elif pid == 211:
                pip_momentum_mag = p_mag
                h2_pip_theta_p.Fill(p_mag, theta)
                h2_pip_phi_p.Fill(p_mag, phi)
                h2_pip_theta_phi.Fill(phi, theta)

        # --- Kinematics from scattered electron ---
        if scattered_electron:
            q = beam_vec - scattered_electron
            Q2 = -q.Mag2()
            nu = q.E()

            if Q2 > 0.0 and nu > 0.0:
                xbj = Q2 / (2.0 * TARGET_MASS * nu)
                if xbj > 0.0:
                    h2_q2_xbj.Fill(xbj, Q2)
                    h1_q2.Fill(Q2)
                    h1_xbj.Fill(xbj)
                    if pip_momentum_mag is not None:
                        h2_q2_pip.Fill(pip_momentum_mag, Q2)
                    event_count += 1


def draw_with_text(hist, outname, text_x=0.2, text_y=0.85):
    """
    Draws a histogram and stamps N = total processed events.
    Uses COLZ for TH2; default draw for TH1.
    """
    outfolder = "output_plots"
    os.makedirs(outfolder, exist_ok=True)
    outname = os.path.join(outfolder, outname)
    
    
    c = ROOT.TCanvas("c_" + hist.GetName(), "", 900, 700)
    draw_opt = "COLZ" if hist.InheritsFrom("TH2") else ""
    hist.Draw(draw_opt)

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.04)
    lat.DrawLatex(text_x, text_y, f"N = {event_count}")

    c.SaveAs(outname)
    c.Close()


def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python3 check_lund.py /path/to/file.lund [prefix]")
        print("  python3 check_lund.py '/path/to/dir/*.lund' [prefix]")
        print("  python3 check_lund.py /path/to/dir [prefix]")
        sys.exit(1)

    pattern = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) >= 3 else "piplus"

    # Build a list of files: support directory, glob, or single file
    paths = []
    if os.path.isdir(pattern):
        paths = sorted(glob.glob(os.path.join(pattern, "*.lund")))
    else:
        matched = glob.glob(pattern)
        if matched:
            paths = sorted(p for p in matched if os.path.isfile(p))
        elif os.path.isfile(pattern):
            paths = [pattern]

    if not paths:
        print(f"No LUND files matched: {pattern}")
        sys.exit(1)

    print(f"Processing {len(paths)} file(s)")
    for p in paths:
        print(f"  -> {p}")
        process_file(p)

    # Draw/save plots once after all files are processed
    draw_with_text(h2_q2_xbj, f"{prefix}_Q2_vs_xbj.pdf")
    draw_with_text(h1_q2,    f"{prefix}_Q2_hist.pdf")
    draw_with_text(h1_xbj,   f"{prefix}_xbj_hist.pdf")
    draw_with_text(h2_q2_pip, f"{prefix}_Q2_vs_piplus_pmag.pdf")

    # Electron plots
    draw_with_text(h2_e_theta_p,   f"{prefix}_electron_theta_vs_p.pdf")
    draw_with_text(h2_e_phi_p,     f"{prefix}_electron_phi_vs_p.pdf")
    draw_with_text(h2_e_theta_phi, f"{prefix}_electron_theta_vs_phi.pdf")

    # Pi+ plots
    draw_with_text(h2_pip_theta_p,   f"{prefix}_piplus_theta_vs_p.pdf")
    draw_with_text(h2_pip_phi_p,     f"{prefix}_piplus_phi_vs_p.pdf")
    draw_with_text(h2_pip_theta_phi, f"{prefix}_piplus_theta_vs_phi.pdf")

    print("Done. Wrote:")
    print(f"{prefix}_Q2_vs_xbj.pdf")
    print(f"{prefix}_Q2_hist.pdf")
    print(f"{prefix}_xbj_hist.pdf")
    print(f"{prefix}_Q2_vs_piplus_pmag.pdf")
    print(f"{prefix}_electron_theta_vs_p.pdf")
    print(f"{prefix}_electron_phi_vs_p.pdf")
    print(f"{prefix}_electron_theta_vs_phi.pdf")
    print(f"{prefix}_piplus_theta_vs_p.pdf")
    print(f"{prefix}_piplus_phi_vs_p.pdf")
    print(f"{prefix}_piplus_theta_vs_phi.pdf")


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    main()
