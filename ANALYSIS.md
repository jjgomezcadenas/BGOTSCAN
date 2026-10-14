# Analysis log

Findings and decisions, dated, with the script or document each comes from.
The purpose of the analysis is to measure the relative light yield and the
energy resolution of the BGO crystal as a function of temperature and SiPM
overvoltage, from the 511 keV photopeak of the 22Na source.

## 2026-10-14

### Setup and data

- A 22Na source at a fixed position shoots 511 keV gammas into the crystal
  (and 1275 keV at lower rate). The photopeak of the charge spectrum of the
  selected events gives the light yield per run, in ADC·ns per keV, and its
  width the energy resolution. Only one of the two back-to-back 511 keV
  gammas can reach the crystal: no coincidence to exploit.
- 12 runs of 5000 events, 64 channels × 3000 samples at 125 MHz (24 µs,
  1.6 µs pre-trigger), self-trigger majority 2 at 19 to 20 ADC counts
  (`docs/jscan_doc.tex` for the file layout).
- Five runs at one nominal overvoltage: 37, 38, 39, 40, 41 V at T3 = 128.6,
  153.1, 164.8, 184.8, 195.8 K, a line of about 60 mV/K, the datasheet
  temperature coefficient of a Hamamatsu MPPC breakdown voltage. Two bias
  scans step down from the nominal point: 40, 39.5, 39, 38.5 V at 185 K and
  41, 40.5, 40, 39.5, 39 V at 196 K. The nominal overvoltage itself is not
  known (asked to Laura).
- Within each bias scan the crystal thermometer T3 rises between the first
  and the last run (184.8 to 187.4 K; 195.8 to 197.3 K) while the case
  thermometer T2 is steady or relaxing: the crystal lags the case. Use the
  measured T3 of each run, and expect a percent-level temperature effect
  across a bias scan on the light yield. (`runinfo.py`, `temperatures.txt`.)

### SiPMs, overvoltage, map and source (datasheet in `extDocs/`, Laura's answers of 2026-10-14)

- Hamamatsu S14160-6050HS, 64 of them at 6 mm pitch tiling the 48 × 48 mm²
  face of a 20 mm thick BGO crystal, Teflon wrapped. Datasheet at 25 °C:
  breakdown 38 V, recommended overvoltage 2.7 V, coefficient 34 mV/K, gain
  2.5e6 at 2.7 V and linear in overvoltage (about 0.9e6/V), PDE 50% at 2.7 V
  and saturated above 3 V, crosstalk 7%.
- Each SiPM goes through a shaper into the 50 Ω input of the digitizer
  (details to come). This is why the single-photoelectron pulse is 24 ns
  wide at half maximum (`template.py`), and the pulse height in counts
  depends on the shaper gain.
- The bias was set from the datasheet values and best guess, never from a
  measured breakdown. Extrapolating the 34 mV/K coefficient from 298 K to
  the crystal temperature T3, the overvoltage is 4.8, 4.9, 5.5, 5.8, 6.5 V
  along the scan (37 to 41 V) and 4.3 and 4.4 V at the lowest bias points
  (38.5 V at 187 K, 39 V at 197 K): no run is near breakdown, and the
  "nominal" overvoltage rises by 35% along the scan, which the light yield
  must be corrected for (gain, linear in overvoltage; PDE saturated). The
  SiPMs were probably warmer than the crystal, between T3 and T2 (5 to
  9 K above), which lowers these overvoltages by 0.2 to 0.3 V. The
  correction is measured from the bias scans: charge versus bias at 185
  and 196 K, straight lines in overvoltage whose zero is the breakdown
  voltage at that temperature.
- Temperature control: thermometers only measure; heating resistors set
  the temperature; runs were taken in a hurry, not fully stabilised. Use
  the measured T3 per run, with a few kelvin of uncertainty on the SiPM
  temperature.
- Channel map (`src/channel_map.csv`): row-major, x from -21.8 to +21.8 mm
  along the column, y from +21.8 mm (row 0) to -21.8 mm (row 7), 6.2 mm
  pitch inside each 4 × 4 matrix and 6.4 mm across the centre. Swapped
  cables are possible; the co-occurrence of channels in clusters over a
  run will test the map.
- Source: 22Na, 3 cm from the face, centred; the whole face is
  illuminated, so clusters should appear everywhere. Channels 48 to 55
  (y = -15.6 mm) and 11, 12 (y = +15.6 mm) light up together in the first
  events looked at, 31 mm apart: a map or trigger question to settle with
  the occupancy map.

### Acquisition and trigger (`src/py/runinfo.py`, `docs/runinfo_doc.tex`)

- Every run records its 5000 events in 40.7 s, 123 Hz, in bursts of 3 events
  1.3 ms apart separated by 22.8 ms gaps: the readout is saturated at about
  47 MB/s. The recorded rate is the same in all runs and carries no physics.
- The board's 24-bit trigger counter counts every trigger, recorded or not.
  It gives the true trigger rate: 1.9 kHz at 128 K and 37 V, 20 kHz at
  153 K, 66 kHz at 165 K, 132 kHz at 185 K and 40 V, 225 kHz at 196 K and
  41 V. At fixed bias (39 V) the rate falls with temperature (66, 43,
  2.7 kHz at 165, 187, 197 K), so the rate follows the overvoltage, not the
  dark-count temperature dependence: dark pulses of one or two
  photoelectrons crossing the threshold, squared by the majority-2 logic.
- The recorded events are a random sample of the triggers, 0.05% to 7%.
  The number of real events in a run is therefore of order
  5000 × (real rate / trigger rate), smallest at the warm end at nominal
  bias.

### Events (`src/py/selection.py`, `docs/selection_doc.tex`)

- Per channel the noise is white, 4 ADC counts rms per sample; single
  photoelectron pulses of about 15 counts and 100 to 200 ns are visible in
  the channels. There is no SPE run; these pulses are the handle for an SPE
  calibration and for the gain versus bias.
- The light of an event is concentrated in 5 to 12 channels, the same ones in
  most events (48 to 55, 11 and 12 at 197 K; 2 to 12 at 128 K): the fixed
  source position. The multiplicity and the channel map carry position
  information from the start.
- The common mode across channels is small, 3 to 4% of the per-channel
  variance, below one count rms. The per-sample median across channels
  leaks signal during a pulse (one or two counts, times 64 in the sum: 60 to
  80% of the charge lost). The estimate is therefore interpolated across the
  pulse window: median outside it, a straight line between 32-sample means
  at its ends inside it. This keeps the pulse and still removes the slow
  wander, which brings the control-window rms of the sum from 35 to 23
  counts.
- Baseline of the sum from the last 250 samples of the record, since the
  pre-trigger window often contains the start of the pulse (the trigger
  fires late on these slow, small pulses).
- Pulse decay about 1 µs at 197 K and several µs at 128 K (a 3 µs plateau
  seen in one event). The charge window, 1 µs before to 4 µs after the
  smoothed maximum, has to be set per temperature from the average pulse.
- Classes: REAL (charge above 5 noise widths and a rising edge), TAIL (charge
  above the cut, no rising edge: the pulse started before the record, only
  its tail was recorded; small pulses are affected more, which biases the
  spectrum), EMPTY, plus a PILE-UP flag when a second maximum outside the
  charge window is separated from the first by a valley of at least 5 noise
  widths. First 40 events of 197K_39V_th20: 40 REAL, 5 of them with pile-up.
  First 20 of 128K_37V_th19: 15 REAL, 1 with pile-up, 2 TAIL, 2 EMPTY.
- Zero suppression (`selection.py` steps 11 and 12): a channel is kept if
  its charge is above k times the per-channel charge noise (k = 3
  default, changeable at the prompt). In event 1 of 197K_39V_th20 nine
  channels above 8 sigma carry 109% of the 64-channel charge; the 30
  channels with negative charge (per-channel drift over the window) are
  worth -30%. Q_all (64 channels), Q_0 (positive channels only, a biased
  worst case) and Q_kept (above k) will all be computed and compared on
  the photopeak width. The kept sum has half the noise of the full sum.
  The full sum also shows an excursion of ±200 counts just before the
  trigger, absent from the kept channels: pickup at the trigger time in
  the channels without light, inside the charge window; to be looked at.
- Not done on purpose: no fixed veto around the main pulse for the pile-up
  search (to be judged from the time-difference distribution of the pass);
  Q/Q_ctrl is not printed (a single noise sample, meaningless on one
  event).

### Single photoelectron (`src/py/spe.py`, `src/py/template.py`, their docs)

- There is no SPE run; dark pulses were searched in the last 10 µs of the
  records. With 4 counts of noise per 8 ns sample and a pulse of about 15
  counts over three samples, a single photoelectron is a 4-sigma object,
  and a matched filter gains at most a factor 1.5.
- At 39 V and a 12-count threshold the candidates are noise maxima: their
  mean shape is a symmetric spike with the autocorrelation of the noise,
  and their number (0.65 per channel per 10 µs) matches 3-sigma noise,
  while the datasheet dark current scaled to 197 K predicts 0.06 real
  pulses in the same window. At 20 counts only 101 candidates in 300
  events.
- At 41 V and 20 counts, 1566 candidates in 300 events with a real pulse
  shape (rise over two samples, 24 ns at half maximum, tail to 50 ns),
  but sitting on an elevated, fluctuating level: at 6.5 V overvoltage
  the channels are busy with dark pulses, crosstalk and afterpulses.
- The refinement step of template.py (re-alignment on the filtered
  maximum) broadens the template by the alignment jitter; the first
  template, aligned on the raw maximum, is sharpened by the selection
  bias; the truth is in between. To be fixed or dropped.
- Conclusion for now: the gain per run comes from the bias scans, not
  from SPE spectra; the SPE needs a dark run with long records at high
  bias, or a lower-noise chain. The pulse shape at 41 V is kept as the
  template for a matched filter if one is needed.

### Consequences for the light yield

- Light yield versus temperature needs a gain and PDE correction, because
  the overvoltage is only nominally constant along the scan. The bias scans
  give the response versus overvoltage at two temperatures; the SPE pulses
  give the gain directly.
- The two lowest-bias runs (38.5 V at 187 K, 39 V at 197 K) have the
  cleanest samples (4 to 7% of the triggers are substantial) but the
  smallest signals; they are the place to develop the selection, not to
  measure yield.
- Decay time versus temperature needs only the pulse shape and is safe as
  long as the pulse is visible.

### Output of the pass (decided 2026-10-14)

One HDF5 file per run, written by Julia, read by Python: `events` (event
index, trigger id, timestamp, t_max, charge window bounds, baseline,
per-channel charge noise), `sum` (the corrected summed waveform per event,
3000 samples, float16: 30 MB per run), `qc` and `tc` (charge and time of
the maximum of every one of the 64 channels in the window, float32), and
metadata (run conditions, trigger rate, selection parameters, version).
Everything else (Q, Q_ctrl, rising edge, verdict, pile-up, Q_all, Q_0,
Q_kept, multiplicity, clusters, barycentre, average pulses, decay fits) is
computed offline in Python from these. Per-channel waveforms are not
written. The per-event decay fit is a second pass on `sum`, once the model
is chosen on the average pulse. Clustering: seed above 5 sigma, neighbours
on the map above 3 sigma, several clusters allowed (Compton).

### Next

- The first JScan pass over whole runs, with the definitions of
  `docs/selection_doc.tex`: per event the baseline, control rms, maximum and
  its time, rising edge, charge, multiplicity, verdict, written to a small
  table per run.
- Python plots of the distributions: the fraction of each class per run,
  the charge spectrum of the REAL events and its 511 keV photopeak, light
  yield and energy resolution versus T3 and bias, the average pulse and
  decay time per run, the time difference of the pile-up maxima.

### Open questions for Laura

1. Readout chain details: the shaper (type, gain, shaping time) and the
   50 Ω input; it sets what a single photoelectron looks like.
2. Trigger: the coincidence window of the majority-2 logic.
3. Later, when the analysis says what is needed: next data (a longer
   pre-trigger and record at the cold end, a dark run with long records at
   high bias for the SPE, stabilised temperatures).

Answered 2026-10-14: source (22Na, 3 cm, centred, Teflon wrapping); SiPM
model and bias rule (datasheet and best guess, no measured breakdown; SiPMs
probably warmer than the crystal); channel map (as in src/channel_map.csv,
swapped cables possible); temperature control (resistors, thermometers
only measure, not fully stabilised); trigger logic (majority 2); readout
saturation was not known.
