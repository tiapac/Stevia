#!/usr/bin/env python3
"""Plot all stellar tracks listed in stellar_evolution_output.txt."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DATA_COLUMNS = (
    "time_myr",
    "mass_ms",
    "teff_k",
    "lum_lsun",
    "mdot_msun_per_yr",
    "vterm_kms",
    "gamma_ed",
)


def read_track(track_path: Path) -> dict[str, np.ndarray]:
    """Parse the simple text format emitted by main.f90."""

    lines = [line.strip() for line in track_path.read_text().splitlines() if line.strip()]
    data_rows: list[list[float]] = []
    mass_ini = lifetime = np.nan
    nsteps = None

    for line in lines:
        if line.lstrip().startswith("#"):
            continue
        values = [float(tok) for tok in line.split()]
        if np.isnan(mass_ini):
            mass_ini = values[0]
            lifetime = values[1] if len(values) > 1 else np.nan
            continue
        if nsteps is None:
            nsteps = int(values[0])
            continue
        data_rows.append(values)
    
    data = np.asarray(data_rows, dtype=float)
    if data.shape[1] != len(DATA_COLUMNS):
        raise ValueError(
            f"Expected {len(DATA_COLUMNS)} columns in {track_path}, found {data.shape[1]}"
        )
    
    return {
        "mass_ini": mass_ini,
        "lifetime": lifetime,
        "nsteps": nsteps if nsteps is not None else data.shape[0],
        "columns": DATA_COLUMNS,
        "data": data,
        "path": track_path,
    }


def parse_mass_list(summary_file: Path) -> list[float]:
    first_line = summary_file.read_text().splitlines()[0]
    tokens = first_line.replace(",", " ").split()
    masses = []
    for token in tokens:
        try:
            masses.append(float(token))
        except ValueError:
            continue
    return masses

def plot_all_tracks(tracks: list[dict[str, np.ndarray]], summary_file: Path) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
    axes = axes.flatten()
    metrics = [
        (1, "Current Mass [M$_{\\odot}$]"),
        (2, "Effective Temperature [K]"),
        (3, "Luminosity [L$_{\\odot}$]"),
        (4, "Mass-Loss [M$_{\\odot}$/yr]"),
        (5, "Terminal Velocity [km/s]"),
        (6, "Gamma Eddington"),
    ]

    cmap = plt.get_cmap("Blues")
    colors = cmap(np.linspace(0.3, 0.9, len(tracks)))

    for track, color in zip(tracks, colors):
        time = track["data"][:, 0]
        lifetime = track["lifetime"] if not np.isnan(track["lifetime"]) else time[-1]
        t_norm = time / lifetime if lifetime else time
        label = f"{track['mass_ini']:.1f} M$_{{\\odot}}$"
        for ax, (idx, ylabel) in zip(axes, metrics):
            ax.plot(t_norm, track["data"][:, idx], label=label, color=color)
            ax.set_ylabel(ylabel)
            if idx in [3, 4, ]:
                ax.set_yscale("log")

    for ax in axes[3:]:
        ax.set_xlabel("Normalized Time (t / t_SN)")
    axes[0].legend(loc="best", ncol=2)
    
    fig.suptitle("Stellar Evolution Tracks (normalized time)")
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    out_path = summary_file.parent / f"{summary_file.stem}_normalized.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def main() -> None:
    summary = Path("stellar_evolution_output/stellar_evolution_output.txt")
    masses = parse_mass_list(summary)


    tracks = []
    for mass in masses:
        # summary_file.parent
        track_path = summary.parent / f"stellar_evolution_mass_{mass:0.2f}.txt"        
        tracks.append(read_track(track_path))


    out_fig = plot_all_tracks(tracks, summary)
    print(f"Saved {out_fig}")


if __name__ == "__main__":
    main()
