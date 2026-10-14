# Analysis log

Findings and decisions, dated, with the script or document each comes from.
The purpose of the analysis is to measure the relative light yield and the
energy resolution of the BGO crystal as a function of temperature and SiPM
overvoltage, from the 511 keV photopeak of the 22Na source.

## 2026-10-15

### State of the energy resolution question (handoff)

- The 511 keV photopeak of 195K_41V_th20 is 45% FWHM in every charge
  estimator (sum, all channels, positive channels, zero-suppressed,
  equalized or not), with the photopeak fit started from the measured
  width (a fixed 8% start had produced spurious narrow fits, including a
  "27%" that was quoted and was wrong). Expected: 10 to 15% (840
  photoelectrons give 9% by themselves). So something is wrong, not the
  crystal.
- Per SiPM (events whose largest charge is in that SiPM): the photopeak
  position varies only +-10% across the face, but within one SiPM the
  peak is 38% wide; and the charge rises by 50% from events with light
  in <= 6 SiPMs to events with light in > 10 (multiplicity). An (x, y)
  map would gain little. Far from the dead block: 31%, next to it 43%.
- Noise: each biased channel carries a low-frequency noise (1/f-like
  below 0.5 MHz, partly coherent among biased channels, absent in the 12
  dead channels: figure src/py/plots/cs195w_wander.png). The per-channel
  charge noise over the 5 us window with the end-of-record baseline is
  2400 ADC ns (white: 700); the sum's is 66 000 (12% of the peak). Not
  the mains: folded at 50 Hz the baseline shows no harmonic (mains.py).
- The decisive clue: the scatter of the end-of-record baseline between
  events falls from 47 counts at 41 V (board trigger rate 225 kHz) to
  26, 23, 14 and 5.6 counts at 40.5, 40, 39.5 and 39 V (148, 97, 29,
  2.7 kHz), a factor 8 while the gain changes by 1.5: it follows the
  square root of the trigger rate. Hypothesis: every firing of the
  trigger logic, recorded or not (about 5 per 24 us record at 41 V),
  disturbs the biased channels, as the +-200-count excursion at t = 0
  shows for the firing that opens the record. One event looked at (step
  2 of selection.py, event 1 of 41 V) does not show a clear train of
  such excursions, so the mechanism is not confirmed; the scaling is.
- NEXT, decisive and cheap: pass of 197K_39V_th20 with the current
  settings (dead-channel common mode; its present pass file is from the
  interp days) and offline.py on it. If the photopeak is near 15% at
  39 V, the noise at high trigger rate is the culprit and the cure is a
  quieter trigger (hardware); if it stays at 45%, the width is elsewhere
  (then: look by eye at a narrow-spot and a wide-spot event, channel by
  channel, to see whether the 50% multiplicity effect is in the crystal
  or in our per-channel charge).
- Tried and discarded today (kept out of the code or as non-default
  options): "unlit channels" common mode (leaks a third of the charge),
  local baselines of the sum (post-pulse, line; the pre-trigger samples
  hold the trigger excursion, the record start holds the pulse rise of
  late triggers), per-channel line baseline in the pass (option
  --chbase line, worse), Q_kept at 3 sigma (loses variable light, 60%).
- Working rule adopted after today: one hypothesis, one test, the
  user's OK, one result; no chained trials.

## 2026-10-14

### Setup and data

- A 22Na source at a fixed position shoots 511 keV gammas into the crystal
  (and 1275 keV at lower rate). The photopeak of the charge spectrum of the
  selected events gives the light yield per run, in ADC·ns per keV, and its
  width the energy resolution. Only one of the two back-to-back 511 keV
  gammas can reach the crystal: no coincidence to exploit.
- 12 runs of 5000 events, 64 channels × 3000 samples at 125 MHz (24 µs,
  1.6 µs pre-trigger), self-trigger majority 2 at 19 to 20 ADC counts
  (`docs/manual/jscan_doc.tex` for the file layout).
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

### Thermal layout and the SiPM temperature (`docs/figs/thermometers.jpg`)

- The crystal hangs in an aluminium frame heated top (heater 1) and bottom
  (heater 2); T1 is on the top plate, T2 at the bottom of the frame, T3 at
  the centre of the back face between the Teflon-wrapped crystal and the
  SiPM array; the SiPM board is mounted on the frame. The crystal follows
  the frame with a lag (Teflon); the SiPMs sit between T3 and the frame
  temperature.
- The frame had a vertical gradient at the warm points: T1 = 181 K against
  T2 = 190 K in the 185 K group, 184 against 202 to 205 K in the 196 K
  group. If the array shares part of it, the bottom rows are warmer than
  the top rows: 20 K is 0.7 V of overvoltage, 10 to 15% of gain across
  the array in the warm runs, none at 128 K where the three thermometers
  agree within 3 K. Testable: the per-channel response (mean charge of a
  channel in the photopeak events that light it) against the row of the
  map, at 196 K and at 128 K. If present, equalization must be per run
  from the per-channel response; one more reason to keep the 64
  per-channel charges in the pass output.
- Operative SiPM temperature: T_SiPM = T3 + f (T2 - T3), with f between 0
  and 1 to be fixed from the data (the breakdown voltage from the bias
  scans at 185 and 196 K against the datasheet coefficient). Expected
  overvoltage per run, with Vbd(T) = 38 V + 34 mV/K (T - 298.15 K), for
  f = 0, 0.5 and 1:

| Run | Bias (V) | T3 (K) | (T2+T3)/2 (K) | T2 (K) | OV at T3 | OV at mid | OV at T2 |
|---|---|---|---|---|---|---|---|
| 128K_37V_th19 | 37.0 | 128.6 | 130.2 | 131.8 | 4.77 | 4.71 | 4.66 |
| 153K_38V_th19 | 38.0 | 153.1 | 155.3 | 157.5 | 4.93 | 4.86 | 4.78 |
| 164K_39V_th19 | 39.0 | 164.8 | 168.8 | 172.9 | 5.53 | 5.40 | 5.26 |
| 184K_40V_th19 | 40.0 | 184.8 | 187.5 | 190.1 | 5.85 | 5.76 | 5.67 |
| 185K_39d5V_th19 | 39.5 | 185.8 | 188.0 | 190.1 | 5.32 | 5.25 | 5.17 |
| 186K_39V_th19 | 39.0 | 186.6 | 188.4 | 190.2 | 4.79 | 4.73 | 4.67 |
| 187K_38d5V_th19 | 38.5 | 187.4 | 188.8 | 190.2 | 4.26 | 4.22 | 4.17 |
| 195K_41V_th20 | 41.0 | 195.8 | 200.4 | 205.0 | 6.48 | 6.32 | 6.17 |
| 196K_40d5V_th20 | 40.5 | 196.2 | 200.1 | 203.9 | 5.96 | 5.84 | 5.71 |
| 196K_40V_th20 | 40.0 | 196.6 | 199.8 | 203.1 | 5.45 | 5.34 | 5.23 |
| 196K_39d5V_th20 | 39.5 | 196.9 | 199.7 | 202.5 | 4.94 | 4.85 | 4.75 |
| 197K_39V_th20 | 39.0 | 197.3 | 199.6 | 201.9 | 4.43 | 4.35 | 4.27 |

  The spread between f = 0 and f = 1 is 0.1 V at 128 K and 0.2 to 0.3 V
  at the warm end, against a 35% rise of the overvoltage along the scan:
  the assumption on f is a second-order effect, the scan-wide rise is not.

### Acquisition and trigger (`src/py/runinfo.py`, `docs/manual/runinfo_doc.tex`)

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

### Events (`src/py/selection.py`, `docs/manual/selection_doc.tex`)

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
- Superseded below: with the full statistics of the pass the SPE peak is
  measurable per channel at 41 V (and only there); the pulse shape at
  41 V is the template of the matched filter.

### The 64 channels (`src/py/chanstat.py`, `docs/manual/chanstat_doc.tex`)

- Twelve channels are dead: 19 to 23, 25, 27, 33, 35, 41, 46, 47. Normal
  noise, no dark pulses, never lit in 500 events, charge noise of pure
  white noise. A block on the map covering 19% of the face; light falling
  on it is lost. To be masked in the pass and asked to Laura (known? SiPM,
  bias or cable?).
- Four channels have 30 to 50% more gain than the rest: 11, 12, 48, 49
  (five times the dark-pulse candidates above 20 counts, lit three times
  more often, 30 to 50% more charge when lit). The recurring cluster at
  these channels is a gain pattern, not a position. Equalization is
  mandatory and of this size.
- Channel 44 has a slow baseline wander (charge noise three times the
  others). Noise rms is uniform otherwise (3.4 to 4.3 counts).
- The lit fraction of the normal channels is flat across the face (7 to
  20%): the source illuminates the whole crystal; the map shows nothing
  wrong at this level.
- SPE template (`template.py`): the first template from 195K_41V_th20 at a
  20-count threshold is stable (identical from 300 to 2000 events), 24 ns
  at half maximum, matched-filter gain 1.5; written to
  `src/spe_template_195K_41V_th20.csv` (one file per run). The refined template of step 5 comes out with a
  longer tail in every variant (the filtered selection favours tailed
  pulses) and is kept as a diagnostic only.

### SPE with full statistics from the pass (`src/jscan`, `src/py/spefit.py`)

- The JScan pass reads a run in 27 s and writes 27 MB: per event the
  corrected summed waveform (Int16, scale 4), the charge and time of
  every channel in the window, and the SPE candidates of every channel
  (matched filter, 4.5 sigma). On 195K_41V_th20: 28 664 candidates, 0.09
  per channel and event; the Python prototypes are reproduced on a
  three-event fixture (225 of 225 comparisons).
- 195K_41V_th20 (overvoltage 6.5 V): a clear single-photoelectron peak
  in nearly every live channel, median 20.7 counts, cut on its low side
  by the 17-count threshold. 45 of 52 live channels fitted; normal
  channels within ±9% of one another, errors 1 to 4 counts. Ringers:
  11 at 28.2, 12 at 30.5, 49 at 31.4 counts (g = 1.36, 1.48, 1.52); 48
  bimodal (peaks near 22 and 33), no fit. Dead channels: 0.3 to 1.8 kHz
  of candidates, all at the threshold (the noise floor, a few percent of
  the live channels' 8 to 15 kHz).
- The 196 K bias scan (41, 40.5, 40, 39.5, 39 V, all with the 41 V
  template): the normal channels cannot be followed below 41 V, their
  true peak slides from 20.7 to about 14 counts, below the threshold, and
  the fit is pinned at 20. Channel 49, whose peak stays above the cut,
  gives 31.4, 29.2, 26.7 counts at 41, 40.5, 40 V: 4.7 counts per volt,
  zero at 34.3 V. The datasheet line gives 34.5 V at T3 and 34.9 V at T2:
  the breakdown at 196 K is where the datasheet puts it within half a
  volt; f (SiPM temperature between T3 and T2) is not resolvable this
  way.
- The ringer factor is constant with bias: channel 49 over the normal
  amplitude scaled by overvoltage gives 1.52, 1.53, 1.53, 1.48 at 41 to
  39.5 V. A lower breakdown would make it grow as the bias drops; a
  constant factor is a larger gain after the SiPM (the shaper channel).
  Hence one equalization table, from the 41 V run, applies to every run
  and temperature; channel 48 from its light response.
- Dark rate per live channel at 196 K: 12 kHz at 41 V, 1 kHz at 39 V
  (threshold effect and dark-count dependence together).

### Common mode, baseline and the noise of the sum (2026-10-14, evening)

- The interpolated common mode (median across channels outside the
  charge window, a line inside) distorted the pulse: the median outside
  the window contains the pulse tail, the tail was partly removed beyond
  the window (a kink in the average pulse at +4 us) and the line pulled
  the sum down along the window. It removed 20% of the charge. Replaced
  by the mean of the twelve dead channels, smoothed over 1 us, as the
  default of selection.py and of the pass: the average pulse is then a
  clean single exponential, 2.74 +- 0.01 us at 196 K, and the photopeak
  rose by 25%.
- The dead channels remove only what they share with the live ones. The
  sum of the live channels wanders by about 13 counts rms over 5 us,
  which the dead channels do not see (a wander that comes with the bias);
  the live channels without light cannot serve as monitor either (almost
  every channel carries some light in a 511 keV event). A local baseline
  for the sum does not help, the wander is on the scale of the window.
  Consequence: the charge of the sum has a noise of 66 000 ADC ns, 12% of
  the 511 keV peak, and is not the estimator for the resolution; the
  per-channel estimators that sum only the channels with light carry a
  tenth of that noise. The 23 counts rms that the median subtraction
  showed was a statistical artefact (a noisy estimate correlated with the
  sum), not a common mode removed.
- First resolution numbers (195K_41V_th20, REAL without pile-up): the
  511 keV photopeak in every estimator, the 1275 keV line at 2.5 times it;
  FWHM 41 to 53% with the dead-channel pass, 27 to 46% with the interp
  pass (whose line also removed the wander per channel). Far wider than
  the crystal should give: the dead block (19% of the face) and the
  position dependence are the suspects; the clustering and the position
  selection are next.
- Pile-up: a valley criterion alone flags the rise and the tail of the
  main pulse just outside the window; a veto of 2 us beyond the window
  removes most of it. Genuine random pile-up: about 9% of the events in
  24 us, a real event rate near 4 kHz.

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
  `docs/manual/selection_doc.tex`: per event the baseline, control rms, maximum and
  its time, rising edge, charge, multiplicity, verdict, written to a small
  table per run.
- Python plots of the distributions: the fraction of each class per run,
  the charge spectrum of the REAL events and its 511 keV photopeak, light
  yield and energy resolution versus T3 and bias, the average pulse and
  decay time per run, the time difference of the pile-up maxima.

### Open questions for Laura

1. The twelve dead channels (19 to 23, 25, 27, 33, 35, 41, 46, 47): known?
   SiPMs, bias or cables?
2. Readout chain details: the shaper (type, gain, shaping time) and the
   50 Ω input; it sets what a single photoelectron looks like.
3. Trigger: the coincidence window of the majority-2 logic.
4. Later, when the analysis says what is needed: next data (a longer
   pre-trigger and record at the cold end, a dark run with long records at
   high bias for the SPE, stabilised temperatures).

Answered 2026-10-14: source (22Na, 3 cm, centred, Teflon wrapping); SiPM
model and bias rule (datasheet and best guess, no measured breakdown; SiPMs
probably warmer than the crystal); channel map (as in src/channel_map.csv,
swapped cables possible); temperature control (resistors, thermometers
only measure, not fully stabilised); trigger logic (majority 2); readout
saturation was not known.
