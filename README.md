# BGOTSCAN

Analysis of the BGO temperature scan of 2026-10-02: a BGO crystal read out by
64 SiPMs, held in a cryostat at temperatures between 128 K and 197 K, with the
SiPM bias voltage changed with temperature (and scanned at 185 K and 196 K).

## Data

The data are not in git. `2026_10_02_bgo_cryo_tempscan` is a symlink to the
local copy (`~/Data/2026_10_02_bgo_cryo_tempscan`), fetched from
`hyperion.sw.ehu.es:/data/lnavarro/cryostat_data/2026_10_02_bgo_cryo_tempscan`.
Inside:

- `data_cryo/*.h5`: one HDF5 file per run, written by WebDump (CAEN DT2745,
  Scope firmware). The name gives the crystal temperature (T3) in K, the SiPM
  bias (`39d5V` = 39.5 V) and the trigger threshold in ADC counts, e.g.
  `185K_39d5V_th19.h5`. Each file holds 5000 events of 64 channels x 3000
  samples at 125 MHz (8 ns), with 200 pre-trigger samples; about 1 GB each,
  LZF-compressed.
- `temperatures.txt`: the three thermometer readings per run (T1 on the
  aluminium case, T2 under it, T3 on the back face of the crystal).
- `README.txt`, `files.tex`, `analysis.ipynb`: Laura's notes and the
  first-look notebook (a copy is in `nb/`).

The file layout (datasets and attributes) is documented in
`docs/jscan_doc.tex`.

## Layout

- `src/jscan/`: Julia project `JScan`, the processing side (reads the runs,
  will hold the per-event passes). Tests in `src/jscan/test/`.
- `src/py/`: Python scripts for plots and quick looks. `wfview.py` is the
  waveform viewer; it reads the raw HDF5 files with h5py. Figures it saves go
  to `src/py/plots/` (not in git).
- `src/h5lzf/`: the LZF compression filter of h5py as an HDF5 plugin, needed
  for Julia to read the LZF-compressed files (h5py has it built in). See
  below.
- `docs/`: LaTeX documentation, one file per module or script
  (`jscan_doc.tex`, `wfview_doc.tex`).
- `nb/`: notebooks received from the authors, kept as received.

## Setup

Julia 1.12 and Python 3 with `h5py`, `numpy` and `matplotlib`.

1. Instantiate the Julia project:

       julia --project=src/jscan -e 'import Pkg; Pkg.instantiate()'

2. Build the LZF plugin (macOS; compiles against the HDF5 library that the
   Julia project uses, so rebuild it if that library is updated):

       sh src/h5lzf/build.sh

   This produces `src/h5lzf/libh5lzf.dylib` (not in git). `JScan.register_lzf()`
   tells HDF5 where it is; nothing else is needed.

## Tests

    julia --project=src/jscan -e 'import Pkg; Pkg.test()'

The tests use a small LZF-compressed fixture in `src/jscan/test/data/` with
the layout of a real run (made by `src/jscan/test/make_fixture.py`). One test
reads a real run and is skipped if the data symlink is absent.

## Viewer

    python3 src/py/wfview.py 128K_37V_th19 --events 0-3 --smooth 64
    python3 src/py/wfview.py 128K_37V_th19 --events 7 --channels
    python3 src/py/wfview.py 187K_38d5V_th19 --average --save 187K

See `docs/wfview_doc.tex` for the arguments and what each figure shows.
