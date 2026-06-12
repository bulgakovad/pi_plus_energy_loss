#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt


def momentum(px, py, pz):
    return math.sqrt(px * px + py * py + pz * pz)


def theta_deg(px, py, pz):
    p = momentum(px, py, pz)
    if p == 0:
        return float("nan")

    cos_theta = pz / p
    cos_theta = max(-1.0, min(1.0, cos_theta))

    return math.degrees(math.acos(cos_theta))


def read_lund_kinematics(lund_file):
    """
    Reads pi+ and electron kinematics from a LUND file.

    Assumed particle-line format:

        index charge type pid parent daughter px py pz E mass vx vy vz

    Therefore:
        pid = column 4 -> cols[3]
        px  = column 7 -> cols[6]
        py  = column 8 -> cols[7]
        pz  = column 9 -> cols[8]
    """

    data = {
        "piplus": {
            "p": [],
            "theta": [],
        },
        "electron": {
            "p": [],
            "theta": [],
        },
    }

    with open(lund_file, "r") as f:
        while True:
            header = f.readline()

            if not header:
                break

            if not header.strip():
                continue

            header_cols = header.split()

            try:
                n_particles = int(header_cols[0])
            except ValueError:
                raise ValueError(f"Bad LUND header line:\n{header}")

            for _ in range(n_particles):
                line = f.readline()

                if not line:
                    raise ValueError("Unexpected end of file while reading particles.")

                cols = line.split()

                if len(cols) < 9:
                    continue

                try:
                    pid = int(cols[3])
                    px = float(cols[6])
                    py = float(cols[7])
                    pz = float(cols[8])
                except ValueError:
                    continue

                p = momentum(px, py, pz)
                theta = theta_deg(px, py, pz)

                if pid == 211:
                    data["piplus"]["p"].append(p)
                    data["piplus"]["theta"].append(theta)

                elif pid == 11:
                    data["electron"]["p"].append(p)
                    data["electron"]["theta"].append(theta)

    return data


def print_summary(label, p_values, theta_values):
    print(f"{label}:")
    print(f"  entries: {len(p_values)}")

    if len(p_values) == 0:
        print("  no entries found")
        print()
        return

    print(f"  p range:     {min(p_values):.6f} to {max(p_values):.6f} GeV")
    print(f"  theta range: {min(theta_values):.6f} to {max(theta_values):.6f} deg")
    print()


def add_entry_text(ax, n_entries):
    ax.text(
        0.95,
        0.95,
        f"N = {n_entries}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Plot pi+ and electron momentum/theta distributions from a LUND file."
    )

    parser.add_argument(
        "lund_file",
        type=Path,
        help="Input LUND file",
    )

    parser.add_argument(
        "--out",
        type=Path,
        default=Path("lund_piplus_electron_sanity_check.png"),
        help="Output plot filename",
    )

    parser.add_argument(
        "--p-min",
        type=float,
        default=0.8,
        help="Lower edge of expected pi+ momentum bin",
    )

    parser.add_argument(
        "--p-max",
        type=float,
        default=1.4,
        help="Upper edge of expected pi+ momentum bin",
    )

    parser.add_argument(
        "--theta-min",
        type=float,
        default=30.0,
        help="Lower edge of expected pi+ theta bin in degrees",
    )

    parser.add_argument(
        "--theta-max",
        type=float,
        default=47.0,
        help="Upper edge of expected pi+ theta bin in degrees",
    )

    parser.add_argument(
        "--bins",
        type=int,
        default=50,
        help="Number of histogram bins",
    )

    args = parser.parse_args()

    if not args.lund_file.exists():
        raise FileNotFoundError(f"Input file does not exist: {args.lund_file}")

    data = read_lund_kinematics(args.lund_file)

    piplus_p = data["piplus"]["p"]
    piplus_theta = data["piplus"]["theta"]

    electron_p = data["electron"]["p"]
    electron_theta = data["electron"]["theta"]

    print(f"Read file: {args.lund_file}")
    print()

    print_summary("pi+", piplus_p, piplus_theta)
    print_summary("electron", electron_p, electron_theta)

    if len(piplus_p) == 0:
        print("WARNING: no pi+ particles with pid=211 found.")

    if len(electron_p) == 0:
        print("WARNING: no electrons with pid=11 found.")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # Top left: pi+ momentum
    axes[0, 0].hist(piplus_p, bins=args.bins)
    axes[0, 0].axvline(args.p_min, linestyle="--", linewidth=2)
    axes[0, 0].axvline(args.p_max, linestyle="--", linewidth=2)
    axes[0, 0].set_xlabel(r"$p_{\pi^+}$ [GeV]")
    axes[0, 0].set_ylabel("Counts")
    axes[0, 0].set_title(r"$\pi^+$ momentum distribution")
    axes[0, 0].grid(alpha=0.3)
    add_entry_text(axes[0, 0], len(piplus_p))

    # Top right: electron momentum
    axes[0, 1].hist(electron_p, bins=args.bins)
    axes[0, 1].set_xlabel(r"$p_{e^-}$ [GeV]")
    axes[0, 1].set_ylabel("Counts")
    axes[0, 1].set_title(r"$e^-$ momentum distribution")
    axes[0, 1].grid(alpha=0.3)
    add_entry_text(axes[0, 1], len(electron_p))

    # Bottom left: pi+ theta
    axes[1, 0].hist(piplus_theta, bins=args.bins)
    axes[1, 0].axvline(args.theta_min, linestyle="--", linewidth=2)
    axes[1, 0].axvline(args.theta_max, linestyle="--", linewidth=2)
    axes[1, 0].set_xlabel(r"$\theta_{\pi^+}$ [deg]")
    axes[1, 0].set_ylabel("Counts")
    axes[1, 0].set_title(r"$\pi^+$ polar angle distribution")
    axes[1, 0].grid(alpha=0.3)
    add_entry_text(axes[1, 0], len(piplus_theta))

    # Bottom right: electron theta
    axes[1, 1].hist(electron_theta, bins=args.bins)
    axes[1, 1].set_xlabel(r"$\theta_{e^-}$ [deg]")
    axes[1, 1].set_ylabel("Counts")
    axes[1, 1].set_title(r"$e^-$ polar angle distribution")
    axes[1, 1].grid(alpha=0.3)
    add_entry_text(axes[1, 1], len(electron_theta))

    fig.suptitle(f"LUND sanity check: {args.lund_file.name}")
    fig.tight_layout()

    fig.savefig(args.out, dpi=200)
    print(f"Saved plot to: {args.out.resolve()}")

    plt.show()


if __name__ == "__main__":
    main()