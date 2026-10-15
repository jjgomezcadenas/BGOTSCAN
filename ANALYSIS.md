# Analysis log

Findings and decisions, dated, with the script or document each comes from.
The purpose of the analysis is to measure the relative light yield and the
energy resolution of the BGO crystal as a function of temperature and SiPM
overvoltage, from the 511 keV photopeak of the 22Na source.

## 2026-10-15

### The seed factor refined, the unmasked sum, and what the control charges measure (`src/py/offline.py`, `src/py/q0corr.py`, `docs/q0_correction.tex`)

- The seed-channel correction written out, step by step, in
  docs/q0_correction.tex (figures docs/figs/q0corr_195K_*.png by
  q0corr.py): the factor is the photopeak position of the events whose
  seed is channel k, relative to the median over channels, and the
  charge of the event is divided by the factor of its seed. On
  195K_41V_th20 the positions spread by 0.109 rms (0.81 to 1.26) before
  and 0.035 after, and the line goes from 38.2% to 20.7%. The 1275 keV
  line, which the factor never used, goes from a smeared bump to a
  resolved peak: the factor is a collection factor.
- Two more estimators (commit 608b972): Q_0_posf, the same correction
  with the factor from a Gaussian fit of each seed-channel population
  (100+ events) or the parabola-interpolated mode (40 to 99), and
  Q_all_posf, that factor measured on and applied to the unmasked sum
  of all live channels. The fitted factor changes nothing: 20.5%
  against 20.7% at 41 V, 19.6% against 19.5% at 40.5 V. The "8% FWHM
  granularity term" of the first version of the note was wrong (the
  0.035 after correction was the granularity of the re-measurement);
  struck from the note. The unmasked sum is 1.5 points wider (22.2 and
  21.6%), linear in the charges, without the pedestal of the mask (the
  masked peak is 2.8% higher), and its noise term from the control sums
  is exact.
- The noise term from the control sums over all REAL events is wrong:
  at 40.5 V it gives 22.1% for Q_all_posf, a line that is 21.6% wide;
  a linear estimator cannot be narrower than its noise. The fits are
  not at fault (checked by eye: the Gaussian follows the whole peak).
  Two tests added to offline.py (columns nz.EMPTY / after.E, figure
  ctrlvsq):
  - The control sum of Q_all_eq against the charge of the event, in
    twelve bins: the mean is flat (slope +0.20% at 41 V, -0.04% at
    40.5 V): the pulse does not feed the tail; the afterpulsing
    hypothesis is not supported, whatever feeds the tail is below 0.2%
    of the charge. The rms is flat at 40 000 to 43 000 ADC ns for all
    events below 0.7e6 and grows above, 55 000 to 115 000 across the
    1275 keV region, in both runs; unexplained. The all-REAL term
    (46 000) was inflated by the high-charge events. For the 511 keV
    line the term is the one at the photopeak: 42 000 +- 900 at 41 V,
    41 000 +- 900 at 40.5 V.
  - The EMPTY events as a quiet noise sample: useless. Those with a
    valid control window (154 and 175) give an rms of 550 000 and a
    mean of -79 000 +- 44 000: records with large structures the pulse
    search did not call a pulse. The class EMPTY is not "quiet record".
- Remainder after subtracting the photopeak noise term, Q_all_posf:

  | | 41 V | 40.5 V |
  |---|---|---|
  | observed | 22.2 +- 0.2% | 21.6 +- 0.2% |
  | noise term at the photopeak | 16.5 +- 0.4% | 19.4 +- 0.4% |
  | remainder | 14.9 +- 0.5% | 9.5 +- 1.0% |

  The subtraction closes at both biases and the remainders disagree by
  five points, five sigma, in the wrong direction (41 V has more
  photoelectrons). Same electronics noise in ADC ns, same crystal and
  source; what differs: overvoltage 6.5 against 6.0 V, trigger rate 225
  against 148 kHz, dark rate per channel. Open, with the growth of the
  control scatter with the charge.
- Budget of the 20.7% of Q_0_pos at 41 V, as far as it can be written:
  about two thirds of the variance is electronic noise, the rest is
  the photoelectron statistics (930 photoelectrons, 7.7% by themselves,
  more with the crosstalk and afterpulsing excess) and what the hole
  leaves after the correction, which this chain cannot separate.
- To do: the noise term of the table measured on the photopeak events
  instead of all REAL events; the pre-window check repeated on the
  same events with --adj-before 0.5; the 196K_40V_th20 run as a third
  bias point for the remainder (its 511 keV line is partly sculpted by
  the trigger, ratio 2.2 to 1275); the noise injection as the
  cross-check of the masked estimators.

### The full window at 3 tau, every whole window kept: the sample doubles, the noise term grows (`src/jscan/pass.jl`, `src/py/offline.py`)

- The pass of the afternoon (no common mode, per-channel baseline from
  the 1 us before and the 2 us after the window, window 1 + 4 us) on
  195K_41V_th20, 196K_40d5V_th20 and 184K_40V_th19: SHORT (no room for
  the full pre-window) took 45 to 51% of the events, not a third; and
  the 2 us after-window sat on the tail (23% of the height at 4 us =
  1.5 tau), so Q_0 came out at 487k against 564k with the end-of-record
  baseline: a constant 14% loss, not a resolution term. Results on the
  1903 REAL events without pile-up at 41 V: Q_0_eq 467k, 38.3% FWHM,
  noise term 13.2%, 36.0% after; Q_0_pos 28.5%, 24.9% after (814 of the
  1903 with a seed channel without factor). 40.5 V, 1676 events: Q_0_pos
  24.7%, 18.8% after. 184K_40V_th19: Q_0_pos 29.3%, 23.6% after; both
  lines visible, ratio near 2.5; decay 3.07 us.
- Decided (evening): the window ends at 3 tau after the maximum (8.1 us
  at 196 K), so the after-window sits where the exponential is 5 to 2.4%
  of its height (a constant loss of about 8%, computable per run from
  the fitted decay); the full window with its noise measured on the
  control charges and subtracted in quadrature is the estimator, all the
  light, 95% of the exponential instead of 78%; the pre-window is the
  samples the record has, up to 1 us, with a floor of 0.25 us, below it
  the after-window alone; SHORT is only a window cut by the record start
  or without its 2 us after-window. Statistical error on the remainder
  after subtraction, 1900 events: 16.0 +- 0.7% for the full window
  against 18.0 +- 0.5% for the clipped one, no cost that matters. The
  pass writes `events/short` and `events/nb_used`; offline.py splits the
  photopeak by pre-window group (figure `prewindow`) and falls back to
  one exponential when the two agree within 10% or an error exceeds half
  its value.
- Results, REAL without pile-up, window 1 + 8.1 us (figures
  `src/py/plots/195K_41V_*`, `196K_40d5V_*`):

  | | 41 V | 40.5 V |
  |---|---|---|
  | EMPTY / REAL / TAIL / SHORT | 8.5 / 86.8 / 2.0 / 2.7% | 8.4 / 87.6 / 1.6 / 2.4% |
  | REAL without pile-up | 4115 | 4126 |
  | Q_0_eq peak (ADC ns) | 608 400 +- 4600 | 514 600 +- 3400 |
  | Q_0_eq noise (ADC ns) | 42 800 | 42 400 |
  | Q_0_eq width / noise term / remainder | 38.2 +- 0.4 / 16.6 +- 0.2 / 34.4 +- 0.5% | 35.6 +- 0.4 / 19.4 +- 0.2 / 29.8 +- 0.5% |
  | Q_all_eq width / noise term / remainder | 40.3 +- 0.4 / 18.4 +- 0.2 / 35.9 +- 0.5% | 36.6 +- 0.4 / 21.8 +- 0.2 / 29.3 +- 0.5% |
  | Q_0_pos width / noise term / remainder | 20.7 +- 0.2 / 17.0 +- 0.2 / 11.9 +- 0.5% | 19.5 +- 0.2 / 19.6 +- 0.2 / 0 |
  | decay, one exponential | 2.723 +- 0.008 us | 2.690 +- 0.008 us |

  Q_0 at 41 V 635k against 487k with the 4 us window (the 1.3 of 95%
  over 78% less the 8% of the after-window). The noise of Q_0_eq grew
  from 26 100 to 42 800 ADC ns, a factor 1.64 for a window 1.82 times
  longer, between the square root and the linear scaling. The
  two-exponential fit is degenerate on both runs.
- The noise term of a positive-only sum (Q_0, Q_0_pos) from the control
  charges is an overestimate: in the event the mask keeps the channels
  whose noise came out positive, whose fluctuation is smaller than a
  full noise sample over the same channels, which is what the control
  gives; at 40.5 V the term equals the observed width and the remainder
  is zero. The subtraction is exact for the unmasked sums (Q_all_eq:
  remainder 36 and 29%, with the position spread still inside). The
  estimator with both properties is the response correction applied to
  the full sum, the per-channel correction of the scope (B9), to do.
- The pre-window split tests the trigger, not the baseline: full
  against partial pre-window, Q_0_pos peak 604 400 +- 4600 against
  586 700 +- 4200 at 41 V (2018 and 2054 events), 529 300 +- 4600
  against 524 100 +- 3400 at 40.5 V (1797 and 2262); 43 and 67 events
  with the after-window alone. The 1275 keV line sits almost entirely
  in the full-pre-window group: the maximum of a prompt trigger lands
  about 0.25 us after it, past the 2 us boundary, and the bright events
  trigger promptly; the partial group is the late triggers, the dimmer
  events, with its 511 keV peak 1 to 3% lower. A direct test of the
  baseline is a pass of the same run with `--adj-before 0.5`, the same
  events with two pre-windows.
- Pile-up flagged fell from 16.7% to 5.2%: the search region beyond the
  window shrank, and a second pulse inside the 9 us window goes
  unflagged, about 4% of the events at the 4 kHz real rate.

### The 39 V test: the trigger-rate noise is not the width (`src/jscan/pass.jl`, `src/py/offline.py`)

- Pass of 197K_39V_th20 with the current defaults (dead-channel common
  mode, the 41 V SPE template, the only one that exists), then offline.py
  (figures `src/py/plots/197K_39V_*.png`). REAL 95%, pile-up 13%, 4133
  events used; decay 2.72 +- 0.01 us (2.76 at 41 V).
- Photopeak FWHM/peak, 41 V against 39 V (trigger rate 225 against
  2.7 kHz): Q 49 / 58%, Q_all 47 / 55%, Q_0 46 / 44%, Q_0_eq 45 / 45%,
  Q_kept 62 / 58%. The charge noise of the sum sigma_Q fell from 66 000
  to 36 000 ADC ns, as the baseline scatter predicted; the peak fell to
  437k from 564k (0.77, the overvoltage ratio 4.4 / 6.5 V, PDE
  saturated). The width did not move.
- Why no visible effect: the noise enters in quadrature and was never
  the dominant term. Noise term 2.355 sigma_Q / peak: 28% at 41 V, 19% at
  39 V; subtracted from the observed 46 and 44% it leaves 37 and 40%.
  Quadrature predicts 46 -> 42%, observed 46 -> 44%: the effect is there
  at the 2 to 4 point level and lost in the fit wobble. The hypothesis
  that the trigger disturbance enters the per-channel charges beyond
  what sigma_Q measures is rejected. A quieter trigger would gain a few
  points, not the factor 3.
- The remaining 37 to 40% at both biases is the thing to explain. At
  39 V the multiplicity split repeats the 41 V picture: peak 259k for
  <= 6 lit channels, 378k for 7 to 10, 473k for > 10 (a factor 1.8), 25
  to 34% wide inside each class. Per seed channel: median width 24%,
  36% far from the dead block, 46% next to it.
- The spectra of all estimators are a single broad Gaussian from 150k to
  700k: no Compton continuum below, no 1275 keV line in the plotted
  range. The "photopeak" looks like a mixture of populations, not a
  resolved line.
- Next (the second branch of the handoff): look at the per-channel
  charges of low- and high-multiplicity photopeak events, by eye, to
  decide whether the factor 1.8 is in the light (position, dead block)
  or in our per-channel charge.

### Scope of the full pass, agreed 2026-10-15 (docs/antarctica.tex: "BGO scan: the Antarctica runs", the step-by-step note)

A. The pass (Julia), one file per run; each step with its plot.
- A1 Windows per run from the average pulse: 1 us before the maximum,
  a fixed multiple of the decay time after; a per-run table the pass
  reads; record length and pre-trigger from the file. Plot: the average
  pulse per run, linear and log, with the window shaded and the decay fit.
- A2 Pulse search on the equalized sum of the live channels, no common
  mode: smoothed sum, maximum, height against the smoothed noise, rising
  edge. Plot: per run, the distributions of the height over noise and
  of t_max.
- A3 Classes, in this order: EMPTY (no pulse), TAIL (no rising edge),
  SHORT (less than 1 us of record before the window or 2 us after it),
  REAL (pile-up flag when a second maximum stands beyond the window and
  its veto). Plot: the height spectrum by class and the class fractions
  per run; across runs, the SHORT fraction against the trigger rate.
- A4 Single-photoelectron candidates in the last 10 us of every event
  of every class, skipping events with a second maximum in the search
  region. Plot: candidates per class (SHORT events contribute), the
  amplitude spectrum with the fit per channel and the gain map (spefit).
- A5 Per channel, REAL events only: baseline as the mean of the 1 us
  before and the 2 us after the window; charge over the window; time of
  the channel maximum; a control charge of the same geometry on an
  empty stretch of the tail. Plot: one event, one channel, the waveform
  with the window, the two baseline stretches and the level; per run,
  the control-charge width per channel on the map.
- A6 Output per event: times, window, height, class, flags; for REAL
  events the 64 charges, 64 control charges, 64 channel times; the
  summed waveform; metadata with run conditions, parameters, commit.

B. Offline (Python), per run.
- B7 Sample: REAL without pile-up; the fraction removed as SHORT per
  run; the photopeak position against the SHORT fraction across runs
  as the check that the selection does not favour bright events. Plot:
  that scatter across runs.
- B8 Equalization: the SPE gain table from 41 V; channel 48 and the
  seven without a factor from their light response. Plot: the gain map
  and the light-response map, before and after.
- B9 Response correction: per seed channel first, then per channel from
  the data, applied before the sum. Plot: the per-seed-channel peak map
  before and after, and the photopeak before and after.
- B10 Estimators: equalized sum of all live channels (reference),
  positive-only and zero-suppressed for comparison; the same sums on
  the control charges are the noise of each. Plot: the spectra with
  fits and, beside each, the control-charge spectrum with its width.
- B11 Photopeaks: 511 keV and 1275 keV where present, Gaussian on a
  Compton step; light output in ADC ns and photoelectrons; FWHM over
  peak. Plot: the spectrum with both fits and the step.
- B12 Decay time from the average pulse, one or two exponentials. Plot:
  the average pulse with the fit, linear and log.

C. Resolution budget per run, each term measured: electronic noise
  (control charges), photoelectron statistics (peak in photoelectrons,
  with the crosstalk and afterpulsing excess from the SPE charge against
  its amplitude), response non-uniformity (spread of the per-channel
  peaks before and after correction), equalization error (gain table
  errors), Compton under the peak (the fit's step), the remainder.
  Plot: stacked quadrature terms per run against the measured width.

D. Deliverables (scan.py): per run one table line (T3, bias,
  overvoltage, trigger rate, class fractions, 511 keV peak in ADC ns
  and photoelectrons, resolution and budget, decay time, 511 to 1275
  ratio); light output against temperature raw and corrected for the
  overvoltage (gain linear in overvoltage, breakdown from the datasheet
  slope anchored at the 196 K scan and the 253 K run, PDE saturated,
  the crosstalk and afterpulsing excess the stated systematic); decay
  time against temperature; resolution and budget against temperature
  and bias. Plots: those four, plus the 511 to 1275 ratio per run as the
  trigger check.

Runs: 195K_41V_th20 and 196K_40d5V_th20 only (T3 195.8 and 196.2 K,
-77 C, the Antarctica runs): both hold the 511 keV line intact (ratio
2.5 and 2.4 to 1275 keV), 41 V carries the gain table, 40.5 V has a fifth
of the pile-up and is the control. The 253 K run and the colder runs get
their own notes, each with its window and its trigger check.
Documents kept in step: the pass, offline.py, scan.py, the window table.
Expectation figures for white, independent noise from src/py/noisetoy.py
(noisetoy_doc.tex). The note absorbs pre_analysis.tex (retired
2026-10-15): the detector with the array figure of chanmap.py, the scan,
the acquisition, the single photoelectron; dropped from it the filter
study and the refined template, which gained nothing, and the dark-pulse
search at 39 V. The dead channels' "charge noise of white noise" is not
carried over (in doubt, see above).

### The noise, settled on empty windows; the baseline adopted (`src/py/wavefit.py`, `docs/manual/wavefit_doc.tex`)

- Method: the charge of a 5 us window with nothing in it, in the tail of
  the record (12 to 24 us), per channel and for the sum, over 236 clean
  events of 195K_41V_th20 (245 of 197K_39V_th20, same numbers within
  10%): the scatter of that charge is the noise a baseline method leaves
  on a real charge, with no signal at risk. Scatter in ADC ns, 41 V:

  | Baseline per channel | dead | live | sum of 52 live |
  |---|---|---|---|
  | end of record, last 2 us (the pass) | 3473 | 3581 | 29 000 |
  | mean of the 2 us right after the window | 3044 | 3130 | 26 600 |
  | line on 1 us before and 6 us after | 2540 | 2670 | 26 600 |
  | line on the after side only | 5114 | 5238 | 45 900 |
  | line + sinusoid | 2595 | 2770 | 31 000 |
  | end of record after the dead-channel common mode | 3350 | 3719 | 60 000 |
  | white noise, 4 counts | 756 | 756 | 5 400 |

- The dead-channel common mode of the pass is the noise of the sum: the
  raw sum scatters by 29 000, close to the 25 800 of independent
  channels (the common part is small); subtracting the running mean of
  the 12 dead channels injects their own slow noise, 3473 / sqrt(12) x
  52 = 52 000, and the sum goes to 60 000 (the 66 000 of offline.py).
  Within a record the dead channels share nothing with the live ones
  (correlation 0.00): the subtraction removes nothing. Dropped.
- The slow noise per channel is broadband, 1/f-like, the same in dead
  and live channels (3473 against 3581), independent of bias and trigger
  rate; the single-channel spectra show no line, the autocorrelation one
  lobe to 1.5 us and nothing after. It is in the front-end or the
  digitizer. Samples more than about 1.5 us from the window know
  nothing of the noise inside it: a baseline from adjacent samples on
  both sides gains 1.37, from the after side alone 1.14 (a slope
  extrapolated from one side hurts: 0.68), and no model fitted outside
  can reach the white floor, 3.4 times below.
- Withdrawn from the entry below and from cmnoise_doc.tex: the "slow
  oscillation near 0.33 MHz" (a noisy spectrum of the mean at 0.08 MHz
  resolution; no line in any single channel) and the "wave" wording.
  In doubt: chanstat's "charge noise of pure white noise" for the dead
  channels, contradicted by the 3473 above; to be rechecked there.
- Between events, the per-channel offsets from the run baseline (1.3
  counts rms per channel) are correlated by digitizer group: 0.5 to 0.6
  within each block of 16 channels (0-15, 16-31, 32-47, 48-63), 0.2 to
  0.3 across blocks, dead and live alike (figure
  docs/figs/cmnoise_195K_offsets.png). An offset of the digitizer's
  groups, record to record; any baseline taken inside the record
  removes it.
- Adopted (2026-10-15, evening): no common mode; per channel, the
  baseline from the mean of the 1 us before and the 2 us after the
  window (or the line through them) for the events whose window has
  room before it (maximum at least 2 us into the record, two thirds of
  the events at 41 V; these are also the events whose window is not
  cut by the record start, which cost 8% of the peak); the others kept
  in the output with a flag and excluded from the spectra; then the gain
  table, then the sum. Expected: the sum's noise from 66 000 to about
  27 000, the per-channel noise from 3700 to 2600. For Laura: the
  front-end's 1/f noise, 3.4 times white per channel, is now the floor.

### The noise is not the trigger: a slow oscillation in every channel, a weak common part in the biased ones (`src/py/cmnoise.py`, `docs/manual/cmnoise_doc.tex`)

- cmnoise.py isolates the noise in the tail of the record (last 12 us,
  pulse gone) on 300 events of 195K_41V_th20 and 197K_39V_th20 (236 and
  245 clean events): the mean of the 52 live channels and the mean of
  the 12 dead channels sample by sample, their rms raw and after a 1 us
  running mean against the white-noise expectation, the 64 x 64
  correlation of the running means, the spectra, the average around the
  trigger, and a search for excursions of the live median (figures
  docs/figs/cmnoise_195K_*.png).
- No kicks of the trigger logic. Any excursion of the live median above
  0.94 counts per channel (49 on the sum) lasting 264 ns is found; 0 in
  236 events at 41 V against 637 expected if each of the 225 kHz firings
  kicked the channels. The average around the trigger shows only the
  pulse of the median channel, a step of +3.3 counts per channel in the
  live median, +0.05 in the dead mean: the "excursion of 200 counts on
  the sum at the trigger" of selection_doc.tex is this step times 52,
  light and not pickup.
- What the noise is: a slow oscillation near 0.33 MHz (a 3 us period,
  seen in the traces and as a peak in the spectra, with smaller ones at
  0.6 and 0.9 MHz) present in every channel, dead or live, with its own
  phase per channel: the running-mean rms is 0.117 for the live mean
  against 0.047 white, 0.181 for the dead mean against 0.098, the same
  excess; correlation between channels 0.03 live-live, 0.00 dead-dead
  and live-dead, no block structure (SiPM matrices or digitizer groups
  of 16). On a single channel it is about 0.7 counts rms, the 2400 ADC
  ns of per-channel charge noise against 700 white. It is in the
  front-end or the digitizer, not in the bias or the SiPMs.
- The weak common part among the biased channels (0.03 per pair, uniform)
  still dominates the sum: 1 + 51 x 0.03 = 2.5, so 60% of the slow
  variance of the sum is common, 40% per channel. An 18 MHz line in the
  live mean and not in the dead mean is a coherent pickup on the bias
  side, harmless for the charge.
- Nothing depends on the bias or the trigger rate: every number is the
  same at 41 V (225 kHz) and 39 V (2.7 kHz). The morning's clue, a
  baseline scatter falling as the square root of the trigger rate (47
  counts at 41 V to 5.6 at 39 V), compared a pass with the dead-channel
  common mode to passes with the interpolated median: with the same
  method the pass files give 47 and 29 (clipped 32 and 28), and
  cmnoise.py 44 and 48 on clean events. Withdrawn.
- Between events every channel's offset from the run baseline scatters
  by 1.3 counts, dead or live, with a correlation of 0.4 to 0.5 between
  any two channels: a board-wide offset of about 0.9 counts plus a
  per-channel part; it cancels in the charge (same-record baseline).
- Consequence: the charge noise is the slow oscillation within the
  record, per channel from the front-end and weakly common among the
  biased channels, not the trigger. Software handle: per channel, a
  baseline model over the window and its margins, and fewer channels in
  the sum. For Laura: the low-frequency noise of the front-end (0.33
  MHz), the 18 MHz pickup, and what the dead channels are connected to.

### The 39 V peak is the 1275 keV line: the trigger loses the 511 keV events below 40 V (`src/py/offline.py`, bias scan at 196 K)

- Correction by seed channel (offline.py, new charge Q_0_pos = Q_0_eq
  divided by the photopeak position of the seed channel of the event,
  relative to the median, measured on the run itself; documented in
  offline_doc.tex). Photopeak FWHM Q_0_eq -> Q_0_pos: 41 V 46 -> 22%,
  40.5 V 39 -> 20%, 40 V 38 -> 19%, 39 V 45 -> 41% (39.7% with channel
  48 dropped: dropping it changes little). At 39.5 V the correction is
  meaningless (two lines of similar weight, the per-channel peak is
  either). The 41 V corrected peak is visibly narrower with a shoulder
  on its low side; the fit sits on the core.
- The correction worked at 41 V and not at 39 V because the two runs
  do not show the same line. Q_0 spectra of the 196 K bias scan (REAL
  without pile-up; figure docs/figs/offline_196K_39d5V_spectra.png):

  | Bias | Trigger | 511 peak | 1275 peak | ratio |
  |---|---|---|---|---|
  | 41 V | 225 kHz | 540k | 1350k | 2.5 |
  | 40.5 V | 148 kHz | 420k | 1000k | 2.4 |
  | 40 V | 97 kHz | 350k | 780k | 2.2 |
  | 39.5 V | 29 kHz | 280k | 530k | 1.9 |
  | 39 V | 2.7 kHz | absent (would be near 220k) | 430k | |

  At 39.5 V two clean peaks, the upper with more events (2050 against
  1460); at 39 V the lower is gone (a few dozen events near 220k to
  260k, the "multiplicity <= 6" class), and only 22 of 4133 events lie
  above 1.8 times the peak against 899 of 3863 at 41 V. Everything
  measured on 197K_39V_th20 (the 45% width, the seed-channel spread, the
  multiplicity split) was the 1275 keV line.
- Mechanism: the single-photoelectron pulse falls from 20.7 counts at
  41 V to about 14 at 39 V through the 20-count threshold. Below
  40.5 V a channel fires only on two photoelectrons piled within the
  24 ns pulse; a 511 keV event gives about one photoelectron per 20 ns
  in its brightest channel over the 2.7 us decay and crosses only on a
  fluctuation (hence the late triggers and the TAIL events); a 1275 keV
  event has 2.5 times the rate and always triggers. The ratio of the
  peaks falling from 2.5 to 1.9 is the same effect: the trigger keeps
  the bright side of the 511 peak and pushes its recorded position up.
  The 39 V test of the morning was not the test intended; its
  conclusion on the noise (the width does not follow the noise) stands.
- Resolution of the 511 keV line after the seed-channel correction:
  20% at 40.5 V and 19% at 40 V with a noise term (2.355 sigma_Q / peak)
  of 13%, leaving about 15% in quadrature, inside the expected 10 to
  15%. To be confirmed by eye on the fits (the shoulder, the uncorrected
  events) before quoting. At 41 V the 28% noise term of the full sum
  overestimates the noise of Q_0 (22% observed).
- Consequence for the whole scan: the nominal-bias runs have
  overvoltages of 4.8 to 5.8 V, single-photoelectron pulses of roughly
  15 to 18 counts against a 19-count threshold, the regime of 39.5 to
  40 V here; their 511 keV peaks may be sculpted the same way. The runs
  at 39 V (197 K) and 38.5 V (187 K) are 1275 keV runs for the 511 keV
  yield. NEXT: the Q_0 spectra of every run with the 511 to 1275 ratio
  and the two-peak structure as the check (a pass per run, then
  offline.py), before any light yield versus temperature.
- Figures kept: docs/figs/offline_196K_39d5V_spectra.png (the two
  lines), docs/figs/offline_197K_39V_seedmap.png (the response map at
  39 V, on the 1275 keV line).

### The photopeak position depends on the seed channel: light collection, not gain (`src/py/offline.py`)

- offline.py now locates the photopeak per seed channel (the live channel
  with the largest equalized charge) on Q_0 and on Q_0_eq, and prints one
  line per channel (events, both positions relative to the median, gain
  factor, width). 197K_39V_th20, 44 channels with 40+ events:
  spread of the positions 0.11 rms on both, range 0.75 to 1.22 on Q_0
  and 0.76 to 1.21 on Q_0_eq; median width of one channel 28% on Q_0,
  24% on Q_0_eq. The gain table narrows the peak of each channel a little
  and leaves the spread between channels untouched. The ringers 11, 12,
  49 (factors 1.36 to 1.52) sit at 1.00, 0.88, 1.03 after equalization.
  An imprecise equalization is not the cause of the spread.
- The spread follows the map (figure `src/py/plots/197K_39V_seedmap.png`):
  the lowest peaks ring the dead block, 17 and 18 at 0.76 and 0.79, then
  28, 34, 36, 37, 31, 24, 29 at 0.82 to 0.88 (also 4 and 5 at the top
  middle, 0.82 and 0.85); the highest are far from it, 52 at 1.21, 2 at
  1.18, 8, 10, 16, 55 at 1.12. An interaction next to the block loses 15
  to 25% of its light. The 2% quoted earlier for "next to the block"
  was an average dominated by channel 48 (752 of the 2779 events).
- No top-to-bottom gradient: rows 0 and 1 average 1.03, rows 6 and 7
  1.04. If the array shares the frame gradient, the SPE table absorbed it.
- The width of a single seed channel far from the block is 12 to 20%
  (52, 55, 2, 59, 13), with the 19% noise term still inside.
- Channel 48, a ringer without a gain factor, is the seed of 752 events
  (it wins the equalized comparison whenever it has light); its 41% width
  is a mixture. To be given a factor from its light response, or dropped.
- Conclusion: an uncorrected spread of 11% rms in the response costs 26%
  FWHM on the overall peak, the largest single term. Next: one factor per
  seed channel, measured from the peaks of the run itself, applied to the
  charge of the event; then the fit of the overall peak again.

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
- NEXT (done, see the entry above): pass of 197K_39V_th20 with the current
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
