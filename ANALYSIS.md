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
- Not done on purpose: no fixed veto around the main pulse for the pile-up
  search (to be judged from the time-difference distribution of the pass);
  Q/Q_ctrl is not printed (a single noise sample, meaningless on one
  event).

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

1. SiPM model, breakdown voltage at room temperature, and the rule used to
   set the bias at each temperature (nominal overvoltage).
2. Channel map: which channel number sits where on the readout face.
3. Trigger: the coincidence window of the majority-2 logic, and whether the
   board trigger rate was seen during the runs.
4. Readout: was the saturation at 123 Hz known; can the throughput be raised
   or the trigger made more selective (threshold, majority, hardware sum).
5. Temperature control: which thermometer the controller regulates on, and
   the waiting time after a setpoint before the runs.
6. Next data: longer pre-trigger and record at LN2 (pulse longer than the
   4 µs window at 128 K; trigger later than 1.6 µs after the pulse start).
7. A dark run without the source at one or two temperatures, for the
   dark-count rate per channel and a clean SPE amplitude.

Answered: the source is 22Na at a fixed position (2026-10-14).
