# `vaspout_h5.py` — README

**Purpose:** `vaspout_h5.py` is a lightweight Python tool for **post-processing VASP `vaspout.h5` files**. It extracts band structure, density of states, and magnetic-moment information, and produces publication-ready plots (PNG) and optional gnuplot-ready data (`.txt`). It can be used interactively from Python or from the command line, and there is a helper `mass_process.sh` for parallel batch runs.

---

## Table of contents

* [Quick facts](#quick-facts)
* [Requirements](#requirements)
* [Files in this package](#files-in-this-package)
* [Command-line usage (examples)](#command-line-usage-examples)
* [Python API (examples)](#python-api-examples)
* [Detailed options for `plot_BS` / `BandStructurePlot`](#detailed-options-for-plot_bs--bandstructureplot)
* [`description` / `kpaths` / `bands` syntax (full spec)](#description--kpaths--bands-syntax-full-spec)

  * [Projection `description` string grammar](#projection-description-string-grammar)
  * [K-paths (`--kpaths`) syntax](#k-paths---kpaths-syntax)
  * [Band selection (`--bands`) syntax](#band-selection---bands-syntax)
  * [Valid orbital specifiers](#valid-orbital-specifiers)
* [Outputs and filenames](#outputs-and-filenames)
* [Batch processing with `mass_process.sh`](#batch-processing-with-mass_processsh)

  * [Usage and behavior](#usage-and-behavior)
  * [Notes & examples](#notes--examples)
  * [Optional: process only directories containing `vaspout.h5` (one-liner)](#optional-process-only-directories-containing-vaspouth5-one-liner)
* [MAGMOM extraction](#magmom-extraction)
* [Troubleshooting & tips](#troubleshooting--tips)

---

## Quick facts

* **Input:** `vaspout.h5` (VASP HDF5 output)
* **Default plot output:** **PNG**
* **Gnuplot data:** optional `.txt` (when `--gnuplot` is used)
* **Interfaces:** CLI + Python API
* **Intended use:** research post-processing (publication figures)

---

## Requirements

* Python 3.8+
* Python packages:

  * `h5py`
  * `numpy`
  * `matplotlib`
* Local modules (must be available on `PYTHONPATH` or in the same directory):

  * `Plotting` (must export `BandStructurePlot`, `DensityOfStatesPlot`)
  * `InformationCollector`

---

## Files in this package

* `vaspout_h5.py` — main module and CLI.
* `mass_process.sh` — batch-run helper (included below).
* `Plotting.py`, `InformationCollector.py` — required modules expected by `vaspout_h5.py` (not included here).

---

## Command-line usage (examples)

Basic pattern:

```bash
python vaspout_h5.py <calc_type> "<description>" [--kpaths "<kpaths>"] [--bands "<bands>"] [other options]
```

Examples:

```bash
# Band structure (projected)
python vaspout_h5.py band "1,2 Mo dxy" --kpaths "1 2 3 -4" --bands "5-15" --E0 fermi

# DOS (total)
python vaspout_h5.py dos "total" --E0 0 --gnuplot
```

Helpful CLI options:

* `calc_type` — `band` or `dos` (required).
* `description` — projection string (see grammar below).
* `--kpaths` — ordering of k-path segments (BS only).
* `--bands` — select specific bands/ranges.
* `--E0` — energy reference (float) or `"fermi"`.
* `--gnuplot` — export plot data to `.txt`.
* `--bandnums` — display band indices on BS plots.
* `--color`, `--mult`, `--alpha` — visual controls.
* `--folder` — override output folder (default: folder containing `vaspout.h5`).

---

## Python API (examples)

```python
from vaspout_h5 import vaspout_h5, plot_BS_DOS

# Load (defaults to 'vaspout.h5' if omitted)
BS = vaspout_h5("path/to/BS/vaspout.h5")
DOS = vaspout_h5("path/to/DOS/vaspout.h5")

# Extract MAGMOM to POSCAR
DOS.extract_magmom()

# Plot band structure (projected)
BS.plot_BS(description="1,2 Mo dxy", kpaths="1 2 3 -4", bands="5-15", E0="fermi", save=True)

# Plot DOS (total)
DOS.plot_DOS(description="total", E0=0, save=True, gnuplot=True)

# Combined BS + DOS figure (utility function)
plot_BS_DOS(BS, DOS, description="1,2 Mo dxy")
```

---

## Detailed options for `plot_BS` / `BandStructurePlot`

The `plot_BS` wrapper calls `BandStructurePlot(ic, **kwargs)`. Below is a complete, function-level description of the keyword arguments supported by `plot_BS` and how they affect plotting. (This mirrors the implementation in `BandStructurePlot`.)

### Signature (relevant kwargs)

```py
BandStructurePlot(ic: InformationCollector,
                  min_diff=None,
                  ax=None,
                  E0=None,
                  gnuplot=False,
                  folder='',
                  save=False,
                  bandnums=False,
                  color=None,
                  mult=1,
                  alpha=0.1,
                  linestyle='-',
                  **kwargs)
```

### Parameter details

* `ic` (InformationCollector)
  The populated information container created by `vaspout_h5.read_BS()` (contains Emat, Pmat, K-points, Efermi, Evalmax, Egap, ions, num_ions, etc.).

* `min_diff` (None | float) — *default: None*
  Threshold used when grouping consecutive bands for labeling (`bandnums=True`). If `None`, it defaults to `(Emax - Emin) / 100`, where `Emax`/`Emin` are the max/min values of the band energies being plotted. Bands closer than `min_diff` are considered part of the same group for number-range annotation.

* `ax` (matplotlib.axes or None) — *default: None*
  Axis to draw onto. If `None`, uses `plt.gca()`.

* `E0` (float | 'fermi' | None) — *default: None*
  Energy reference subtracted from the band energies before plotting:

  * `None`: subtract **Evalmax** (valence-band maximum) — so 0 is VB maximum.
  * `'fermi'`: subtract **ic.Efer** (Fermi level).
  * numeric: subtract that numeric value.
    This is the same logic used by `vaspout_h5.plot_BS` before calling the plotting routine.

* `gnuplot` (bool) — *default: False*
  If `True`, the function calls `BandToGnuplot(...)` to write out gnuplot-friendly data files. The gnuplot output routine receives (Pmat, Emat, Kpts, Dvec, projection, bands, description, kpaths, folder). Files are written to `folder` (or `ic.folder` if `folder==''`).

* `folder` (str) — *default: ''*
  Directory to which outputs (gnuplot data and saved images) are written. If `''`, the folder from `ic.folder` (set from the `vaspout.h5` path) is used.

* `save` (bool) — *default: False*
  If `True`, the function saves the figure immediately with a base filename derived from the `description`: `folder + description.replace(' ', '_')`. (Note: Matplotlib determines final extension/format; defaults are PNG/Matplotlib config.)

* `bandnums` (bool) — *default: False*
  If `True`, the routine annotates band numbers (or ranges) at the right edge of the plot. Grouping is controlled by `min_diff`. The algorithm compares final k-point energies (`Emat[-1, :]`) to detect separations between bands.

* `color` (int | float | str | None) — *default: None*
  Color specifier passed to Matplotlib plotting calls.

* `mult` (float) — *default: 1*
  Multiplier applied to marker sizes for projected points (dot area is `abs(Pmat) * 100 * mult`).

* `alpha` (float) — *default: 0.1*
  Alpha (opacity) for the scatter markers showing projections. Range `[0,1]`.

* `linestyle` (str or dash tuple) — *default: '-'*
  Line style used when plotting "clean" band plots (see `description` semantics below). Accepts any Matplotlib linestyle specifier (e.g., `'-'`, `'--'`, `(0, (1,1))`, etc).

### Plotting behavior & notes

* **Energy shifting** — energies (`Emat`) are shifted by subtracting:

  * `Evalmax` if `E0 is None`
  * `ic.Efer` if `E0 == 'fermi'`
  * the numeric `E0` otherwise

* **Projection handling** — `ProjectionPainter(ions, num_ions, description)` is used to build the projection specification. The code attempts to index `Emat[..., projection['direction']]` to choose the directional component; if that fails, it falls back to `Emat[...,0]`.

* **Band selection** — before plotting, energies and projectors are sliced using the band-selection produced by `BandWriter(ic.bands)` (this is applied as `Emat = Emat[:, bands]` and `Pmat = Pmat[:, bands]`).

* **'clean' description special-case**
  If the `description` string contains `'clean'`:

  * The routine plots plain band lines using `linestyle` and `color`.
    If the `description` does not contain `'clean'`:
  * The routine plots the band lines in black (`':'`) and overlays a scatter of projected weights:

    ```py
    ax.scatter(..., s = abs(Pmat)*100*mult, alpha=alpha, label=description, c=color)
    ```

  This lets you choose between a plain band plot (`"clean"`) or a projected dot plot.

* **K-line vertical ticks & energy references**

  * Vertical lines (`vlines`) are drawn at k-path tick positions (`Ktik`) using a dotted line.
  * A horizontal line is drawn at energy 0. If `E0 is None` (i.e., energies were shifted by Evalmax), an extra horizontal line at the calculated `Egap` is drawn and labeled.

* **Gnuplot export**
  If `gnuplot=True`, `BandToGnuplot(...)` is invoked with the processed matrices and projection—this writes plot data files to `folder`.

* **Saving**
  If `save=True`, the saved filename base is `folder + description.replace(' ', '_')`. If you need a specific filename or format (PDF, SVG), call `plot_BS`/`BandStructurePlot` from Python and use `plt.savefig("myfile.pdf")` explicitly.

* **Axis and limits**
  The function sets x-ticks from k-path ticks and labels, x-limits across the k-path, and y-limits to the min/max of the plotted energies.

### Error / edge cases

* If the projection indexing tries to access a non-existent direction, an `IndexError` is caught and the code falls back to the first component (`Emat[...,0]`).
* `mass_process.sh` and other wrappers may call `plot_BS` with CLI-parsed args; ensure `description` and `bands` are passed in the expected formats so band slicing and projections behave consistently.

---


## `description` / `kpaths` / `bands` syntax (full spec)

### Projection `description` string grammar

The `description` string tells the plotting layer which atoms and orbitals to sum (project) when producing projected plots.

**Examples:**

* `total` — total DOS
* `1,2 Mo x p dxy` — atoms 1 and 2 of species `Mo`, x-direction projection, orbitals `p` and `dxy`
* `3 O pz` — atom 3 of species `O`, orbital `pz`

**Components (order matters):**

1. **Atom indices (optional)** — comma-separated, 1-based indices. Example: `1,2,5`. If omitted, projection applies to all atoms of the specified species.
2. **Species (recommended)** — element symbol matching `ion_types` in `vaspout.h5` (e.g., `Mo`, `O`, `Te`).
   **Important:** summation across different species in a single description is **not supported** (e.g., `"Mo Te"` is invalid).
3. **Direction (optional)** — `x`, `y`, `z` for vector-like projectors.
4. **Orbitals (one or more)** — e.g., `s`, `p`, `px`, `dxy`, etc.

If the `description` is empty (`''`), plotting functions use internal defaults.

---

### K-paths (`--kpaths`) syntax

For line-mode band-structure calculations where k-paths are split into numbered segments:

* Provide a space-separated list of integers, where each integer references a path segment.
* A **negative** integer reverses the direction of that segment.
* Example: `"1 2 3 -4"` → traverse segments 1 → 2 → 3 → segment 4 reversed.

If omitted, the k-path order from `vaspout.h5` is used.

---

### Band selection (`--bands`) syntax

* Single band: `"5"`
* Range: `"3-10"` (inclusive)
* List: `"3,7,12"`
* Combination: `"3-6,9,12-14"`

Indices are 1-based.

---

### Valid orbital specifiers

Common (conventional) labels:

* `s`
* `p`, `px`, `py`, `pz`  (or `p` to sum all p orbitals)
* `d`, `dxy`, `dxz`, `dyz`, `dx2-y2` (or `dx2_y2` depending on plotting module), `dz2`

> If your `Plotting` implementation expects slightly different names (e.g., `dx2y2` vs `dx2-y2`), use the naming scheme your plotting module supports.

---

## Outputs and filenames

* **PNG image(s)** — default output (e.g., `band_plot` and `dos_plot`; Matplotlib default format is PNG).
* **Gnuplot data** — optional `.txt` generated if `--gnuplot` is used.
* **MAGMOM** — written by `extract_magmom()` as `MAGMOM` formatted for `INCAR`.
* Use `--folder` to override output directory (default is the folder containing `vaspout.h5`).

> For exact control of formats and filenames, call plotting functions from Python and use `plt.savefig("name.ext")`.

---

## Batch processing with `mass_process.sh`

`mass_process.sh` runs `vaspout_h5.py` in multiple directories in parallel and throttles the number of concurrent jobs.

### Usage

```
./mass_process.sh max_jobs parent_folder calc_folder [args...]
```

* `max_jobs` — maximum number of concurrent jobs (integer).
* `parent_folder` — directory under which the script searches.
* `calc_folder` — directory **name** to match (e.g., `calc`).
* `[args...]` — arguments passed to `vaspout_h5.py` (e.g., `band "1,2 Mo dxy" --E0 fermi`).

### Full `mass_process.sh` (verbatim)

```bash
#!/bin/bash
# Usage: ./mass_process.sh max_jobs parent_folder calc_folder [args...]
# Runs ./vaspout_h5.py inside every matching subdir, keeping at most max_jobs in parallel.

MAX_JOBS=$1
parent_folder=$2
calc_folder=$3
shift 3
script_args=("$@")

if [[ -z "$MAX_JOBS" || -z "$parent_folder" || -z "$calc_folder" ]]; then
    echo "Usage: $0 max_jobs parent_folder calc_folder [args...]"
    exit 1
fi

# Fixed script to execute inside each subdir
START_DIR="$(pwd)/${0%/*}"
SCRIPT="$START_DIR/vaspout_h5.py"

# Function to throttle jobs to MAX_JOBS
function throttle() {
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
        sleep 1
    done
}

# Find all subdirectories named $calc_folder and run in parallel
find "$parent_folder" -type d -name "$calc_folder" | while read -r dir; do
    throttle
    (
        echo ">>> Entering $dir"
        cd "$dir" || exit 1
        "python" "$SCRIPT" "${script_args[@]}"
    ) &
done

wait
echo "All jobs finished."
```

Source: `mass_process.sh`.&#x20;

### Important behavior details & recommendations

* The script uses `find "$parent_folder" -type d -name "$calc_folder"`, so it **matches directory names**, not directories that contain `vaspout.h5`. If a matched directory lacks `vaspout.h5`, `vaspout_h5.py` will likely error — the wrapper does not check for the file.
* `mass_process.sh` expects `vaspout_h5.py` to be located next to `mass_process.sh` (it computes `START_DIR` from the script location).
* The wrapper uses `python` explicitly; change it if you need another interpreter or virtualenv.
* To process only directories containing `vaspout.h5`, replace the `find` line with a command that finds `vaspout.h5` and returns its parent directories (example below).

---

### Optional: find only dirs that contain `vaspout.h5`

If you prefer to process only directories that actually contain `vaspout.h5`, use:

```bash
find /path/to/parent -type f -name vaspout.h5 -printf '%h\n' | sort -u | while read -r dir; do
  # run per-dir actions
done
```

---

## MAGMOM extraction

`vaspout_h5.py` includes `extract_magmom()`:

* Reads spin moments from `/intermediate/ion_dynamics/magnetism/spin_moments/values` in `vaspout.h5`.
* Handles ISPIN=1, ISPIN=2, and non-collinear outputs.
* Writes a `MAGMOM` file formatted for `INCAR`, arranged with aligned columns and `\` continuation where appropriate.

---

## Troubleshooting & tips

* **Missing modules (`Plotting`, `InformationCollector`)**: ensure these files are in the repo or on `PYTHONPATH`.
* **vaspout.h5 not found**: run `vaspout_h5.py` in the directory containing `vaspout.h5` or pass the correct path when instantiating `vaspout_h5(...)`.
* **Parallel jobs hang**: set a lower `max_jobs`.
* **Interpreter mismatch**: edit `mass_process.sh` (it uses `python`).
* **Different output format**: call plotting functions from Python and use `plt.savefig("name.pdf")`.

---
