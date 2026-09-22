# baikal_dist_calibr

Project for calibrating the light yield of the Baikal-GVD optical modules (OMs) as a function of distance to the muon track, based on MC data (atmospheric muons, 2020).

Overall pipeline: **list of MC ROOT files** (`configs/*`) → **`loadData.C`** turns raw BARS events into compact ROOT trees (`output/data/*.root`) → those trees are then read either by ROOT macros (`OMdistProfile.C`, `mcMuon.C`) to build "number of fired OMs vs distance" histograms/graphs and estimate ph.e., or by the Python notebook (`notebooks/data_research.ipynb`) for exploratory analysis via `uproot`/`pandas`.

## Setup and environment

- **`setup.sh`** — configures the ROOT/BARS environment (gcc-8.3.0 compiler, `thisroot.sh`, `BARS setenv`). Run it once at the start of the session: `. setup.sh` (via `source`, not as an executable, so the environment variables land in the current shell).
- **`.rootrc`** — ROOT config that adds the current directory and `helpers/`/`macros` to `MacroPath`, so macros in `src/` can include `../helpers/help_functions.C` via a relative path.

## File-list configs (`configs/`)

Plain text files listing paths to the source MC ROOT files (one path per line), passed to the macros as the `_filelist` argument:

- `onefile` — a single file, for quickly debugging a macro.
- `smalltest` — a small subset of one cluster's files, for fast runs.
- `cluster1` — all files of cluster 1.
- `twoclusters` — files from two clusters.
- `allfiles` — the full set of available MC files.

Pick a config depending on whether you need a quick test or full statistics.

## `src/` scripts (ROOT macros)

All are run from the project root inside ROOT after `setup.sh`, as `.X src/<name>.C("configs/<config>")`.

### `loadData.C`
Main preprocessing step: reads raw BARS MC events (`BEvent`, `BMCEvent`, `BRecoMuon`, `BGeomTel`) and writes three compact trees to `output/data/`:
- `mc_reco_event.root` (`eRecoTree`) — one entry per event: true (MC) and reconstructed track parameters (angles, energy, number of hits/strings, DEDX, Z-extent, etc.).
- `mc_reco_pulse.root` (`pRecoTree`) — one entry per pulse (fired OM): charge (`pulseLY`), time, distance from the track to the OM (`RecoDistToPoint_BH`), etc.
- `mc_gen_trk.root` (`tGenTree`) — one entry per MC muon (true angle/energy/delay).

Run:
```
.X src/loadData.C("configs/smalltest")
```
Arguments (see the function signature): `_saveEvent/_savePulse/_saveGentrk` — flags for which trees to write, `_minTD/_maxTD` — time window, `_minLY` — charge threshold (currently commented out in the code, i.e. not applied — filtering is deferred to the analysis stage in Python/downstream macros).

The output of this macro is the input for the `data_research.ipynb` notebook.

### `OMdistProfile.C`
Builds a "fraction of OMs with no signal vs distance from track" profile and derives a proxy photoelectron-count estimate from it (ph.e. ~ -ln(1 - signal_fraction), the zero-Poisson method) — overall, split by event groups (by default, source zenith angle), and split by different `pulseLY` cut thresholds. Also applies basic reconstruction-quality cuts (number of hits, strings, covariance-matrix status, angle, Z-extent) and a window on the expected-vs-reconstructed signal arrival time difference.

Run:
```
.X src/OMdistProfile.C("configs/cluster1")
```
Parameters: `_rmin/_rmax/_nsteps` — range and number of bins in distance (m), `_doGroupPlots/_doCutPlots` — toggle grouping by an event variable and by `pulseLY` thresholds, respectively. The grouping variable and its bin edges (zenith angle by default) can be changed by uncommenting one of the ready-made examples in the macro body (number of muons, number of strings, reconstructed theta).

Output: PDF plots in `output/figures/` (`totalnom_vs_dist.pdf`, `signom_vs_dist.pdf`, `phe_estimation*.pdf`) and all histograms/graphs in `output/data/OMdistProfileFile.root`.

### `mcMuon.C`
A more extensive exploratory macro over the same raw MC (without going through `loadData.C` first): builds histograms of the number of muons, true/reconstructed time, polar angle, pulse charge and time, plus 2D "charge vs distance to OM" histograms under several successive cuts (on number of muons, angle, time difference, charge) to visually assess the effect of each cut.

Run:
```
.X src/mcMuon.C("configs/onefile")
```
Main parameters — histogram ranges/bin counts and cut thresholds (`_minNM/_maxNM`, `_minPA/_maxPA`, `_minTD/_maxTD`, `_minTDEvsR/_maxTDEvsR`, `_minLY`, `_rmin/_rmax/_nsteps`, etc.), all with default values.

Output: PDF plots in `output/figures/` and histograms in `output/data/mcMuonFile.root`.

## Helper code (`helpers/`)

- **`help_functions.C`** — shared C++ helpers for the `src/` macros: `extractClusterId()` extracts the cluster number from a file name, `MakeLogGraph()` builds a `TGraphErrors` of `-ln(bin content)` with propagated errors (used for the ph.e. proxy estimate in `OMdistProfile.C`).
- **`analysis_help.py`** — Python helper for the notebook: `LY_vs_dist_ploter(df)` plots a 2D histogram of `pulseLY` vs `RecoDistToPoint_BH` (log color scale), annotated with the `theta` range and event count.

## Notebook `notebooks/data_research.ipynb`

Exploratory analysis of the data produced by `loadData.C`, using `uproot`/`pandas`/`matplotlib`/`seaborn`. Cell order:
1. Load `mc_reco_event.root` and `mc_reco_pulse.root` into `DataFrame`s via `uproot`.
2. `describe()` on both tables — a quick statistics overview.
3. `pairplot` and correlation matrix over key event variables (`nTrueMuons`, `RecoNHits`, `RecoNStrings`, `RecoDEDX`, `RecoTheta`, `RecoPhi`, `RecoZDist`).
4. `merge` the pulse table with the event table on `eventId` — producing `df_merged`, where every pulse carries its event's parameters.
5. Histograms of `RecoTheta`, `RecoDEDX`, `RecoZDist` distributions.
6. The final cell calls `analysis_help.LY_vs_dist_ploter(df_merged)` — the resulting charge-vs-distance plot.

Before running, make sure `output/data/` already contains the files generated by `loadData.C` (the notebook reads them via the relative path `../output/data/...`, so it must be run from `notebooks/`). The notebook itself appends `../helpers` to `sys.path` to import `analysis_help`.
