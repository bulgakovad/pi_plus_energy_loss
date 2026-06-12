#!/usr/bin/env python3
"""
Select ALL LUND events where the generated pi+ satisfies:

    p     in [1.0, 1.2] GeV
    theta in [38, 39] degrees

The script scans all .lund files in a given input folder and writes all selected
events into a new LUND file in the folder where the script is run.

Usage:
    python select_all_lund_piplus_bin.py /path/to/lund_folder

Optional:
    python select_all_lund_piplus_bin.py /path/to/lund_folder --out all_selected.lund
"""

import argparse
import math
from pathlib import Path


def momentum(px, py, pz):
    return math.sqrt(px * px + py * py + pz * pz)


def theta_deg(px, py, pz):
    p = momentum(px, py, pz)
    if p == 0:
        return float("nan")

    cos_theta = pz / p
    cos_theta = max(-1.0, min(1.0, cos_theta))

    return math.degrees(math.acos(cos_theta))


def parse_lund_events(filename):
    """
    Yield one full LUND event as a list of lines.

    Assumes each event starts with a header line whose first column is the
    number of particles in the event.
    """

    with open(filename, "r") as f:
        while True:
            header = f.readline()

            if not header:
                break

            if not header.strip():
                continue

            cols = header.split()

            try:
                n_particles = int(cols[0])
            except ValueError:
                raise ValueError(f"Bad LUND header in {filename}:\n{header}")

            event_lines = [header]

            for _ in range(n_particles):
                line = f.readline()

                if not line:
                    raise ValueError(
                        f"Unexpected end of file in {filename} after header:\n{header}"
                    )

                event_lines.append(line)

            yield event_lines


def event_has_good_piplus(
    event_lines,
    p_min=1.0,
    p_max=1.2,
    theta_min=38.0,
    theta_max=39.0,
):
    """
    Return True if event contains pi+ with p/theta inside the requested bin.

    Particle-line format assumed:

        index charge type pid parent daughter px py pz E mass vx vy vz

    Therefore:
        pid = cols[3]
        px  = cols[6]
        py  = cols[7]
        pz  = cols[8]
    """

    for line in event_lines[1:]:
        cols = line.split()

        if len(cols) < 9:
            continue

        try:
            pid = int(cols[3])
        except ValueError:
            continue

        if pid != 211:
            continue

        px = float(cols[6])
        py = float(cols[7])
        pz = float(cols[8])

        p = momentum(px, py, pz)
        theta = theta_deg(px, py, pz)

        if p_min <= p <= p_max and theta_min <= theta <= theta_max:
            return True

    return False


def main():
    parser = argparse.ArgumentParser(
        description="Select all LUND events with generated pi+ in a given p/theta bin."
    )

    parser.add_argument(
        "input_folder",
        type=Path,
        help="Folder containing .lund files",
    )

    parser.add_argument(
        "--out",
        type=Path,
        default=Path("all_piplus_p1p0_1p2_theta38_39.lund"),
        help="Output LUND file name.",
    )

    parser.add_argument("--p-min", type=float, default=1.0)
    parser.add_argument("--p-max", type=float, default=1.2)
    parser.add_argument("--theta-min", type=float, default=38.0)
    parser.add_argument("--theta-max", type=float, default=39.0)
    parser.add_argument("--N_evnts",type=int,default=-1,help="Maximum number of selected events to write. Use -1 to select all events.",)

    args = parser.parse_args()

    input_folder = args.input_folder

    if not input_folder.exists():
        raise FileNotFoundError(f"Input folder does not exist: {input_folder}")

    if not input_folder.is_dir():
        raise NotADirectoryError(f"Input path is not a folder: {input_folder}")

    lund_files = sorted(input_folder.glob("*.lund"))

    if not lund_files:
        raise FileNotFoundError(f"No .lund files found in {input_folder}")


    n_seen = 0
    n_selected = 0
    reached_limit = False

    with open(args.out, "w") as fout:
        for lund_file in lund_files:
            print(f"Reading {lund_file}")

            file_seen = 0
            file_selected = 0

            for event in parse_lund_events(lund_file):
                n_seen += 1
                file_seen += 1

                if event_has_good_piplus(
                    event,
                    p_min=args.p_min,
                    p_max=args.p_max,
                    theta_min=args.theta_min,
                    theta_max=args.theta_max,
                ):
                    fout.writelines(event)
                    n_selected += 1
                    file_selected += 1

                    # Stop after writing the requested number of selected events.
                    # If args.N_evnts <= 0, keep all matching events.
                    if args.N_evnts > 0 and n_selected >= args.N_evnts:
                        reached_limit = True
                        break

            print(f"  events checked:  {file_seen}")
            print(f"  events selected: {file_selected}")

            if reached_limit:
                break

    print()
    print("Done.")
    print(f"Total events checked:    {n_seen}")
    print(f"Total events selected:   {n_selected}")
    print(f"Output file:             {args.out.resolve()}")

    if n_selected == 0:
        print()
        print("Warning: no matching events found. Either the bin is empty, or the column parsing is wrong.")


if __name__ == "__main__":
    main()