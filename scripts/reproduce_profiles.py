"""Rebuild Fault B profiles and test sensitivity to spatial covariance."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
import rasterio
from scipy.linalg import cholesky, solve_triangular
from scipy.ndimage import map_coordinates
from shapely.geometry import LineString

DATA = Path(__file__).resolve().parents[1]
OUT = DATA / "results"
project = Transformer.from_crs(4326, 32639, always_xy=True)
unproject = Transformer.from_crs(32639, 4326, always_xy=True)


def read_raster(name):
    with rasterio.open(DATA / "rasters" / name) as src:
        return src.read(1, masked=True).filled(np.nan), src.transform, tuple(src.bounds)


corrected, transform, bounds = read_raster("vel_los_corrected_mm_yr.tif")
uncorrected, raw_transform, _ = read_raster("vel_los_uncorrected_mm_yr.tif")
assert transform == raw_transform
faults = json.loads((DATA / "faults/FaultA_FaultB.geojson").read_text())
trace = np.array(next(f["geometry"]["coordinates"] for f in faults["features"]
                      if f["properties"]["name"] == "Fault B"))
line = LineString(np.column_stack(project.transform(*trace.T)))


def sample(raster, points):
    lon, lat = unproject.transform(points[..., 0], points[..., 1])
    col, row = ~transform * (lon, lat)
    return map_coordinates(raster, [row - .5, col - .5], order=1,
                           mode="constant", cval=np.nan, prefilter=False)


def average(samples, mask):
    count = mask.sum(axis=0)
    return np.divide(np.where(mask, samples, 0).sum(axis=0), count,
                     out=np.full(count.shape, np.nan), where=count > 0)


def profile(center, normal, x, offsets):
    tangent = np.array([-normal[1], normal[0]])
    points = center + 1000 * (x[None, :, None] * normal
                              + np.asarray(offsets)[:, None, None] * tangent)
    final = sample(corrected, points)
    raw = sample(uncorrected, points)
    mask = np.isfinite(final) & np.isfinite(raw)
    return average(final, mask), average(raw, mask), mask.sum(axis=0)


def fit_models(x, y, correlation_km):
    # Same covariance and observations for all candidates in each comparison.
    distance = np.abs(x[:, None] - x)
    cov = np.eye(len(x)) if not correlation_km else (
        .9 * np.exp(-distance / correlation_km) + .1 * np.eye(len(x)))
    whiten = lambda values: solve_triangular(cholesky(cov, lower=True), values, lower=True)
    yw = whiten(y)
    baseline = np.array([np.ones_like(x), x / 25]).T
    edge = (x > 0).astype(float)
    smooth = [np.arctan(x / d) * 2 / np.pi for d in (0.5, 1, 2, 4, 8, 16, 32)]
    # Grid search prevents single-start failures. Grid bounds are diagnostic choices.
    centers = np.arange(x.min() + 2, x.max() - 1, 1)
    widths = np.geomspace(1, 30, 21)
    bowl = np.exp(-.5 * ((x[:, None] - np.repeat(centers, len(widths))) /
                         np.tile(widths, len(centers))) ** 2)
    kernels = np.column_stack([baseline, edge, *smooth, bowl])
    kw = whiten(kernels)
    scores = {}
    fitted = {}

    def evaluate(name, candidates, k, bowl_model=False):
        best = None
        for cols in candidates:
            design = kw[:, cols]
            p = np.linalg.lstsq(design, yw, rcond=None)[0]
            if bowl_model and p[-1] >= 0:
                continue
            rss = np.sum((yw - design @ p) ** 2)
            score = len(x) * np.log(rss / len(x)) + 2 * k
            if best is None or score < best[0]:
                best = (score, cols, p)
        if best is not None:
            scores[name] = best[0]
            fitted[name] = kernels[:, best[1]] @ best[2]

    evaluate("trend", [[0, 1]], 3)
    evaluate("step", [[0, 1, 2]], 4)
    evaluate("smooth", ([0, 1, j] for j in range(3, 10)), 5)
    evaluate("bowl", ([0, 1, j] for j in range(10, kernels.shape[1])), 6, True)
    evaluate("joint", ([0, 1, 2, j] for j in range(10, kernels.shape[1])), 7, True)
    return scores, fitted


def main():
    OUT.mkdir(exist_ok=True)
    x = np.arange(-50., 51.)
    summaries, rows, scores = [], [], []
    fig, axes = plt.subplots(4, 1, figsize=(8, 8), sharex=True)
    map_fig, map_ax = plt.subplots(figsize=(8, 5), layout="constrained")
    west, south, east, north = bounds
    im = map_ax.imshow(corrected, extent=(west, east, south, north),
                       vmin=-75, vmax=10, cmap="RdBu", interpolation="nearest")
    map_ax.plot(*trace.T, "k-", label="Mapped lineament B")
    for i, ax in enumerate(axes, 1):
        source = np.genfromtxt(DATA / f"transects/FaultB_T{i}_swath12km.csv",
                               delimiter=",", names=True)
        xy = np.column_stack(project.transform(source["lon"], source["lat"]))
        located = np.isfinite(xy).all(axis=1)
        intercept, direction = np.linalg.lstsq(
            np.column_stack([np.ones(located.sum()), source["distance_km"][located]]),
            xy[located], rcond=None)[0]
        normal = direction / np.linalg.norm(direction)
        final, raw, count = profile(intercept, normal, x, [-6, -3, 0, 3, 6])
        center_lon, center_lat = unproject.transform(*intercept)
        ll = unproject.transform(*(intercept + 1000 * x[:, None] * normal).T)
        map_ax.plot(*ll, color="#762087", lw=.8)
        map_ax.text(center_lon, center_lat, f"T{i}", color="#762087")
        rmse = np.sqrt(np.nanmean((source["los_vel_mm_yr"] - final[25:76]) ** 2))
        row = dict(transect=f"T{i}", rmse=float(rmse))
        for label, values in [("corrected", final), ("uncorrected", raw)]:
            good = count >= 3
            j = np.argmin(np.where(good, values, np.inf))
            row[label] = [float(x[j]), float(values[j]), int(count[j])]
            for half in (25, 50):
                use = good & (abs(x) <= half)
                if use.sum() < 30 or min((x[use] < 0).sum(), (x[use] > 0).sum()) < 10:
                    continue
                for length in (0, 1, 3, 6):
                    s, curves = fit_models(x[use], values[use], length)
                    scores.append(dict(transect=f"T{i}", field=label, half_km=half,
                                       correlation_km=length, n=int(use.sum()), **s))
                    if label == "corrected" and half == 50 and length == 3:
                        for name, color in [("step", "#b66f23"), ("smooth", "#518341"), ("bowl", "#007e89"),
                                            ("joint", "#8548a4")]:
                            if name in curves:
                                curve = np.full(x.shape, np.nan)
                                curve[use] = curves[name]
                                ax.plot(x, curve, color=color, lw=1, label=name)
        ax.plot(x, np.where(count >= 3, final, np.nan), "k", lw=1.4, label="Corrected")
        ax.plot(x, np.where(count >= 3, raw, np.nan), color=".60", lw=.8, label="Uncorrected")
        ax.scatter(x[count < 3], final[count < 3], s=7, color=".7", label="1–2 lines")
        for v in (-25, 0, 25):
            ax.axvline(v, color=".45", lw=.6, ls=":" if v else "--")
        ax.set_ylabel(f"T{i}\nmm/yr")
        summaries.append(row)
        rows.extend(dict(transect=f"T{i}", km=a, corrected=b, uncorrected=c, valid_of_5=d)
                    for a, b, c, d in zip(x, final, raw, count))
    axes[0].legend(loc="lower left", fontsize=7, ncol=3)
    axes[-1].set_xlabel("Signed distance from profile origin (km)")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_profiles.png", dpi=260)
    plt.close(fig)
    map_ax.set(xlim=(52.2, 53.85), ylim=(29.0, 30.4), xlabel="Longitude (°E)",
               ylabel="Latitude (°N)", aspect=1 / np.cos(np.deg2rad(29.7)))
    map_ax.legend(loc="upper right")
    map_fig.colorbar(im, ax=map_ax, label="Corrected LOS velocity (mm/yr)")
    map_fig.savefig(OUT / "Figure_map.png", dpi=260)
    plt.close(map_fig)
    region, panels = plt.subplots(2, 1, figsize=(7, 9), layout="constrained")
    all_lines = [np.array(f["geometry"]["coordinates"]) for f in faults["features"]]
    for ax, raster, label in zip(panels, [uncorrected, corrected], ["Uncorrected", "Corrected"]):
        im = ax.imshow(raster, extent=(west, east, south, north), vmin=-65, vmax=8,
                       cmap="RdBu", interpolation="nearest")
        for coords, style in zip(all_lines, ["-", "--"]):
            ax.plot(*coords.T, color="k", ls=style, lw=.8)
        for i in range(1, 5):
            a = np.genfromtxt(DATA / f"transects/FaultA_T{i}_swath12km.csv",
                              delimiter=",", names=True)
            ax.plot(a["lon"], a["lat"], color="#552064", lw=.8)
            j = int(np.argmin(abs(a["distance_km"])))
            ax.annotate(f"A{i}", (a["lon"][j], a["lat"][j]), xytext=(2, 2),
                        textcoords="offset points", fontsize=8)
        ax.set(xlim=(west, east), ylim=(south, north), title=label, ylabel="Latitude (°N)",
               aspect=1 / np.cos(np.deg2rad(30.4)))
    panels[-1].set_xlabel("Longitude (°E)")
    region.colorbar(im, ax=panels, shrink=.55, label="LOS velocity (mm/yr)")
    region.savefig(OUT / "Figure_regional.png", dpi=260)
    plt.close(region)

    # Map-scale association diagnostic using local normals, with explicit gaps.
    along_rows = []
    for d in np.arange(2000, line.length - 2000, 2000):
        center = np.array(line.interpolate(d).coords[0])
        tangent = np.array(line.interpolate(d + 1000).coords[0]) - line.interpolate(d - 1000).coords[0]
        normal = np.array([-tangent[1], tangent[0]]) / np.linalg.norm(tangent)
        xx = np.r_[np.arange(-5., 0.), np.arange(1., 6.)]
        final, raw, count = profile(center, normal, xx, [-3, -1.5, 0, 1.5, 3])
        valid = count >= 3
        a, b = valid & (xx < 0), valid & (xx > 0)
        if min(a.sum(), b.sum()) >= 3:
            lon, lat = unproject.transform(*center)
            along_rows.append(dict(chainage_km=d / 1000, lon=lon, lat=lat,
                                   corrected_contrast=np.mean(final[b]) - np.mean(final[a]),
                                   uncorrected_contrast=np.mean(raw[b]) - np.mean(raw[a])))
    for name, data in [("profiles", rows), ("model_sensitivity", scores),
                       ("along_trace", along_rows)]:
        with (OUT / f"{name}.csv").open("w", newline="") as dest:
            writer = csv.DictWriter(dest, sorted(set().union(*(r.keys() for r in data))))
            writer.writeheader()
            writer.writerows(data)
    fig, ax = plt.subplots(figsize=(8, 3.1), layout="constrained")
    chain = [r["chainage_km"] for r in along_rows]
    ax.axhline(0, lw=.7, color=".5")
    for group, color, marker in [("uncorrected", ".5", "o"), ("corrected", "#126783", ".")]:
        ax.plot(chain, [r[f"{group}_contrast"] for r in along_rows], marker=marker,
                ls="none", color=color, fillstyle="none" if marker == "o" else "full",
                label=group.capitalize())
    ax.set(xlabel="Distance along mapped lineament from NW endpoint (km)",
           ylabel="NE minus SW LOS\nvelocity (mm/yr)")
    ax.legend(frameon=False)
    fig.savefig(OUT / "Figure_contrast.png", dpi=260)
    plt.close(fig)
    (OUT / "profile_summary.json").write_text(json.dumps(summaries, indent=2))
    time_file = next((DATA / "timeseries").glob("anomaly_center*.csv"))
    data = np.genfromtxt(time_file, delimiter=",", names=True)
    dates = np.array([f"{int(v):08d}" for v in data["date"]])
    dates = np.array([f"{d[:4]}-{d[4:6]}-{d[6:]}" for d in dates], dtype="datetime64[D]")
    y = data["cum_los_mm_corrected"]
    valid = np.isfinite(y)
    dates, y = dates[valid], y[valid]
    t = (dates - dates[0]).astype(float) / 365.2425
    design = np.column_stack([np.ones_like(t), t, np.sin(2*np.pi*t), np.cos(2*np.pi*t)])
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    resid = y - design @ beta
    bread = np.linalg.inv(design.T @ design)
    score = resid[:, None] * design
    fit_summary = dict(n=len(t), slope=float(beta[1]), annual_amplitude=float(np.hypot(*beta[2:])))
    for days in (90, 180, 365):
        taper = np.maximum(1 - abs(t[:, None] - t) * 365.2425 / days, 0)
        cov = bread @ score.T @ taper @ score @ bread
        test = beta[2:] @ np.linalg.solve(cov[2:, 2:], beta[2:])
        fit_summary[str(days)] = dict(slope_se=float(np.sqrt(cov[1, 1])),
                                      annual_wald_p=float(np.exp(-test/2)))
    (OUT / "temporal_sensitivity.json").write_text(json.dumps(fit_summary, indent=2))
    fig, ax = plt.subplots(figsize=(8, 3.1), layout="constrained")
    ax.scatter(dates, y, s=7, color=".55", label="Corrected LOS")
    ax.plot(dates, design @ beta, color="#126783", lw=1, label="Trend + annual")
    ax.set(xlabel="Acquisition date", ylabel="Relative LOS (mm)")
    ax.legend(frameon=False)
    fig.savefig(OUT / "Figure_timeseries.png", dpi=260)
    plt.close(fig)


if __name__ == "__main__":
    main()
