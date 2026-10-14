"""Make the test fixture test/data/185K_39d5V_th19.h5: the layout of a WebDump
run (datasets, attributes, LZF compression) with 3 events, 64 channels and 300
samples, and waveform values following a known pattern so reads can be checked
exactly: w[e, c, s] = 32768 + 100 e + 3 c + (s mod 11), indices 0-based.
Run once: python3 make_fixture.py (writes next to this script, in data/)."""
from pathlib import Path
import numpy as np
import h5py

NEV, NCH, NS, PRE = 3, 64, 300, 20
out = Path(__file__).parent / "data" / "185K_39d5V_th19.h5"
out.parent.mkdir(exist_ok=True)

e, c, s = np.ogrid[:NEV, :NCH, :NS]
w = (32768 + 100 * e + 3 * c + (s % 11)).astype(np.uint16)
ch = np.arange(NCH, dtype=np.int16)
base = 32768.0 + 0.5 * np.arange(NCH)
thr = np.full(NCH, 19, dtype=np.int64)

with h5py.File(out, "w") as f:
    f.create_dataset("waveforms", data=w, chunks=(2, NCH, NS), compression="lzf")
    f.attrs["waveforms_axes"] = "event, channel (see /channels), sample"
    f.create_dataset("channels", data=ch)
    f.create_dataset("timestamp_ns", data=np.array([0.0, 1.0e6, 2.5e6]))
    f.create_dataset("trigger_id", data=np.array([10, 11, 13], dtype=np.uint32))
    a = f.attrs
    a["model"], a["serial"] = "DT2745", "52694"
    a["fw_type"], a["fpga_fw"], a["software"] = "Scope", "1.0.257", "WebDump"
    a["url"] = "dig2://caendgtz-usb-52694"
    a["start_time"] = "2026-10-02T14:47:45"
    a["events"], a["stop_reason"], a["dropped_events"] = NEV, "fixture", 0
    a["sampling_mhz"], a["tsample_ns"] = 125.0, 8.0
    a["record_length"], a["pretrigger"], a["adc_nbit"] = NS, PRE, 16
    a["adc_to_volts"], a["zin_ohm"] = 6.1455e-05, 100.0
    a["dc_offsets_percent"] = np.full(NCH, 50.001)
    a["vga_gains_db"] = np.ones(4)
    a["baselines"] = base
    a["trigger_mode"], a["trigger_source"] = "self", "ITLA"
    a["trigger_channels"] = str(list(range(NCH)))
    a["itl_logic"], a["itl_majority"], a["edge"] = "Majority", 2, "RISE"
    a["threshold_mode"] = "Relative"
    a["thresholds"] = thr
    a["thresholds_board_absolute"] = (base + thr).astype(np.int64)
print(f"wrote {out} ({out.stat().st_size} bytes)")
