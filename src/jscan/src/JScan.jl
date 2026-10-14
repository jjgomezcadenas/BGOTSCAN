"""
    JScan

Julia side of the BGOTSCAN analysis. So far: the reader of the run files
written by WebDump (CAEN DT2745). See docs/jscan_doc.tex.
"""
module JScan

using HDF5

export register_lzf, parse_runname, openrun, waveforms, time_axis

"Directory of the LZF plugin (src/h5lzf), built by its build.sh."
const LZF_PLUGIN_DIR = normpath(joinpath(@__DIR__, "..", "..", "h5lzf"))

# ---------------------------------------------------------------- plugin ----

"Current plugin search path of the HDF5 library."
function plugin_paths()
    n = Ref{Cuint}(0)
    HDF5.API.h5pl_size(n)
    paths = String[]
    for i in 0:Int(n[]) - 1
        len = HDF5.API.h5pl_get(Cuint(i), C_NULL, Csize_t(0))
        buf = Vector{Cchar}(undef, len + 1)
        GC.@preserve buf begin
            HDF5.API.h5pl_get(Cuint(i), pointer(buf), Csize_t(len + 1))
            push!(paths, unsafe_string(pointer(buf)))
        end
    end
    return paths
end

"""
    register_lzf(dir = LZF_PLUGIN_DIR) -> dir

Add `dir` to the plugin search path of the HDF5 library (once), so that the
LZF-compressed waveforms can be read. `dir` must contain the built plugin
`libh5lzf.dylib` (or `.so`).
"""
function register_lzf(dir::AbstractString = LZF_PLUGIN_DIR)
    isdir(dir) || throw(ArgumentError("plugin directory not found: $dir"))
    any(isfile(joinpath(dir, "libh5lzf.$ext")) for ext in ("dylib", "so")) ||
        throw(ArgumentError("no libh5lzf plugin in $dir; build it with src/h5lzf/build.sh"))
    dir = String(dir)
    dir in plugin_paths() || HDF5.API.h5pl_append(dir)
    return dir
end

# -------------------------------------------------------------- run name ----

"""
    parse_runname(name) -> (T_K, V_bias, thr_adc)

Parse a run name such as "185K_39d5V_th19": temperature in K (Int), SiPM
bias in V (Float64, `d` is the decimal point) and threshold in ADC (Int).
"""
function parse_runname(name::AbstractString)
    m = match(r"^(\d+)K_(\d+)(?:d(\d+))?V_th(\d+)$", name)
    m === nothing && throw(ArgumentError("run name not of the form <T>K_<V>[d<dec>]V_th<thr>: $name"))
    V = parse(Float64, m[3] === nothing ? m[2] : m[2] * "." * m[3])
    return (T_K = parse(Int, m[1]), V_bias = V, thr_adc = parse(Int, m[4]))
end

# ------------------------------------------------------------------ open ----

"""
    openrun(path) -> NamedTuple
    openrun(f, path)

Open a run file read-only and return its parameters together with the open
file (`fid`) and the waveform dataset (`wf`), not read. Close it with
`close(run.fid)`. The do-block form closes it for you.
"""
function openrun(path::AbstractString)
    fid = h5open(path, "r")
    wf = fid["waveforms"]
    nsamp, nch, nev = size(wf)
    attr(k) = read_attribute(fid, k)
    name = first(splitext(basename(path)))
    rn = parse_runname(name)
    return (
        name = name, path = String(path),
        T_K = rn.T_K, V_bias = rn.V_bias, thr_adc = rn.thr_adc,
        nev = nev, nch = nch, nsamp = nsamp,
        channels = Int.(read(fid["channels"])),
        tsample_ns = Float64(attr("tsample_ns")),
        sampling_mhz = Float64(attr("sampling_mhz")),
        pretrigger = Int(attr("pretrigger")),
        adc_to_volts = Float64(attr("adc_to_volts")),
        zin_ohm = Float64(attr("zin_ohm")),
        baselines = Float64.(attr("baselines")),
        thresholds = Int.(attr("thresholds")),
        threshold_mode = String(attr("threshold_mode")),
        thresholds_abs = Int.(attr("thresholds_board_absolute")),
        edge = String(attr("edge")),
        itl_logic = String(attr("itl_logic")),
        itl_majority = Int(attr("itl_majority")),
        start_time = String(attr("start_time")),
        timestamp_ns = Float64.(read(fid["timestamp_ns"])),
        trigger_id = Int.(read(fid["trigger_id"])),
        fid = fid, wf = wf,
    )
end

function openrun(f::Function, path::AbstractString)
    run = openrun(path)
    try
        return f(run)
    finally
        close(run.fid)
    end
end

# ------------------------------------------------------------- waveforms ----

"""
    waveforms(run, events::AbstractUnitRange) -> Array{UInt16,3}
    waveforms(run, event::Integer)            -> Matrix{UInt16}

Raw ADC counts of the given events (1-based): (nsamp, nch, nevents), or
(nsamp, nch) for one event.
"""
function waveforms(run, events::AbstractUnitRange{<:Integer})
    checkbounds(Bool, 1:run.nev, events) || throw(BoundsError(1:run.nev, events))
    return run.wf[:, :, events]::Array{UInt16,3}
end

function waveforms(run, event::Integer)
    checkbounds(Bool, 1:run.nev, event) || throw(BoundsError(1:run.nev, event))
    return run.wf[:, :, event]::Matrix{UInt16}
end

# ------------------------------------------------------------- time axis ----

"""
    time_axis(run) -> Vector{Float64}

Time of each sample in ns from the trigger: t[pretrigger + 1] == 0.
"""
time_axis(run) = [(i - 1 - run.pretrigger) * run.tsample_ns for i in 1:run.nsamp]

end # module
