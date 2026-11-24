"""Reconstruct ascending differential phase for the two T1 patch pairs."""
import csv
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def connected_dates(edges):
    neighbors = {}
    for first, last in edges:
        neighbors.setdefault(first, set()).add(last)
        neighbors.setdefault(last, set()).add(first)
    groups = []
    while neighbors:
        pending = [next(iter(neighbors))]
        group = set()
        while pending:
            day = pending.pop()
            if day not in group:
                group.add(day)
                pending.extend(neighbors.pop(day) - group)
        groups.append(sorted(group))
    return min(groups, key=lambda g: (-len(g), g[0])) if groups else []


def invert(edges, threshold):
    days = sorted({day for pair in edges for day in pair})
    passed, failed = set(), set()
    for a, b, c in combinations(days, 3):
        triangle = [(a, b), (b, c), (a, c)]
        if all(pair in edges for pair in triangle):
            residual = edges[a, b] + edges[b, c] - edges[a, c]
            bucket = passed if abs(residual) <= threshold else failed
            bucket.update(triangle)
    reliable = passed - failed
    dates = connected_dates(reliable)
    if not dates:
        return []
    index = {day: i for i, day in enumerate(dates)}
    pairs = sorted((a, b) for a, b in reliable if a in index and b in index)
    design = np.zeros((len(pairs), len(dates)))
    for row, (a, b) in enumerate(pairs):
        design[row, index[a]] = -1
        design[row, index[b]] = 1
    values = [edges[pair] for pair in pairs]
    solution = np.linalg.lstsq(design[:, 1:], values, rcond=None)[0]
    return zip(dates, np.r_[0.0, solution])


def main():
    output = ROOT / "results"
    output.mkdir(exist_ok=True)
    for period in ("pilot", "full"):
        with (ROOT / "marvdasht_revision" / f"{period}_phase_samples.csv").open() as src:
            samples = list(csv.DictReader(src))
        with (output / f"{period}_closure_filtered_phase.csv").open("w", newline="") as dest:
            writer = csv.writer(dest)
            writer.writerow(["site", "closure_threshold_rad", "date", "phase_rad"])
            for site in ("T1_SE6", "T1_SE3"):
                local, control = site + "_local", site + "_control"
                edges = {tuple(row["pair"].split("_")): float(row[local]) - float(row[control])
                         for row in samples if row[local] and row[control]}
                for threshold in (0.3, 1.0, np.pi):
                    for day, phase in invert(edges, threshold):
                        writer.writerow([site, threshold, day, phase])


if __name__ == "__main__":
    main()
