# The pass: one read of a run that writes, per event, the corrected summed
# waveform, the charge and time of every channel in the scintillation window,
# and the single-photoelectron candidates of every channel.
# Definitions: docs/selection_doc.tex (steps 1 to 9), docs/spe_doc.tex (steps 1
# to 3); conventions and output layout: docs/jscan_doc.tex, section "The pass".

using Statistics
using Dates

export moving_average, common_mode, interpolate_cm, pulse_search, windows,
       channel_charges, matched_filter, spe_candidates, read_template,
       process_event, pass

# ------------------------------------------------------------ pure functions ----

"Population rms about the mean (numpy.std)."
rms(v) = (m = mean(v); sqrt(mean((x - m)^2 for x in v)))

"rms after removing the samples beyond `nsig` times a first rms (pulses excluded)."
function clipped_rms(v, nsig = 3.0)
    s = rms(v)
    s == 0 && return 0.0
    m = mean(v)
    return rms(x for x in v if abs(x - m) < nsig * s)
end

"""
    moving_average(v, n) -> Vector{Float64}

Mean of `v` over a window of `n` samples centred on each sample, the edges using
the samples available; the window of sample `i` is `i - n÷2 : i + (n-1)÷2`, as
numpy.convolve(v, ones(n), "same") / counts in wfview.smooth.
"""
function moving_average(v::AbstractVector, n::Integer)
    n <= 1 && return Float64.(v)
    L = length(v)
    c = cumsum(Float64.(v))
    out = Vector{Float64}(undef, L)
    lo_off, hi_off = n ÷ 2, (n - 1) ÷ 2
    @inbounds for i in 1:L
        lo = max(1, i - lo_off)
        hi = min(L, i + hi_off)
        s = c[hi] - (lo > 1 ? c[lo - 1] : 0.0)
        out[i] = s / (hi - lo + 1)
    end
    return out
end

"Median across channels at each sample: x is (samples, channels)."
common_mode(x::AbstractMatrix) = vec(median(x, dims = 2))

"""
    interpolate_cm(m, lo, hi, edge) -> (mt, L, R)

`m` with the samples `lo:hi` (inclusive) replaced by a straight line from `L`, the
mean of `m` over the `edge` samples before `lo`, to `R`, the mean over the `edge`
samples after `hi`. When the window touches an edge of the record the missing end
is replaced by the other.
"""
function interpolate_cm(m::AbstractVector, lo::Integer, hi::Integer, edge::Integer)
    L_ = length(m)
    left = lo > 1 ? mean(@view m[max(1, lo - edge):lo - 1]) : nothing
    right = hi < L_ ? mean(@view m[hi + 1:min(L_, hi + edge)]) : nothing
    left === nothing && right === nothing && return (copy(m), 0.0, 0.0)
    L = left === nothing ? right : left
    R = right === nothing ? L : right
    mt = copy(m)
    n = hi - lo + 1
    mt[lo:hi] .= n == 1 ? L : range(L, R, length = n)
    return mt, L, R
end

"""
    pulse_search(S, nsm) -> (Sh, imax, H, irise)

Moving average of width `nsm`, its maximum (index, height) and the rising edge:
the last sample before `imax` where `Sh < H/2`, or 0 when there is none.
"""
function pulse_search(S::AbstractVector, nsm::Integer)
    Sh = moving_average(S, nsm)
    imax = argmax(Sh)
    H = Sh[imax]
    irise = 0
    @inbounds for i in imax - 1:-1:1
        if Sh[i] < H / 2
            irise = i
            break
        end
    end
    return Sh, imax, H, irise
end

"""
    windows(imax, a, b, ns) -> (lo, hi, clo, chi, short)

Charge window `lo:hi` = `imax - a : imax + b` clipped to `1:ns`, and a control
window of the same length in the pulse-free part of the record: the longer of
`1:lo-1` and `hi+1:ns`, and within it the samples farthest from the pulse;
`short` is true when that part is shorter than the window.
"""
function windows(imax::Integer, a::Integer, b::Integer, ns::Integer)
    lo, hi = max(1, imax - a), min(ns, imax + b)
    n = hi - lo + 1
    nbefore, nafter = lo - 1, ns - hi
    if nafter >= nbefore
        clo, chi = max(hi + 1, ns - n + 1), ns
        short = nafter < n
    else
        clo, chi = 1, min(lo - 1, n)
        short = nbefore < n
    end
    return lo, hi, clo, chi, short
end

"""
    channel_charges(y, nbase, lo, hi, clo, chi, nsm, dt) -> (q, qctrl, imaxc)

Per channel (columns of `y`): baseline = median of the last `nbase` samples,
charge `dt * sum(y - baseline)` over `lo:hi` and over the control window
`clo:chi`, and the sample of the maximum of the moving average of `y - baseline`
inside the window.
"""
function channel_charges(y::AbstractMatrix, nbase::Integer, lo, hi, clo, chi, nsm, dt)
    ns, nch = size(y)
    q = zeros(Float64, nch)
    qc = zeros(Float64, nch)
    imx = zeros(Int, nch)
    for c in 1:nch
        yc = @view y[:, c]
        Bc = median(@view yc[ns - nbase + 1:ns])
        q[c] = dt * (sum(@view yc[lo:hi]) - Bc * (hi - lo + 1))
        qc[c] = dt * (sum(@view yc[clo:chi]) - Bc * (chi - clo + 1))
        sm = moving_average(yc[lo:hi] .- Bc, nsm)
        imx[c] = lo - 1 + argmax(sm)
    end
    return q, qc, imx
end

"""
    matched_filter(x, T) -> z

Least-squares height of a pulse of shape `T` (odd length, peak 1 at the centre)
at every sample: `z[i] = sum_j T[j] x[i+j] / sum T^2`, zero within `H` of the
edges.
"""
function matched_filter(x::AbstractVector, T::AbstractVector)
    H = (length(T) - 1) ÷ 2
    L = length(x)
    z = zeros(Float64, L)
    ssq = sum(abs2, T)
    @inbounds for i in H + 1:L - H
        s = 0.0
        for j in -H:H
            s += T[j + H + 1] * x[i + j]
        end
        z[i] = s / ssq
    end
    return z
end

"""
    spe_candidates(z, lo, hi, thr, sep) -> Vector{Int}

Local maxima of `z` strictly inside `lo:hi` (neighbours inside the range) above
`thr`, taken in decreasing order of height, at least `sep` samples from every one
already taken; returned sorted by index.
"""
function spe_candidates(z::AbstractVector, lo::Integer, hi::Integer, thr::Real, sep::Integer)
    loc = Int[]
    @inbounds for i in lo + 1:hi - 1
        if z[i] > thr && z[i] >= z[i - 1] && z[i] > z[i + 1]
            push!(loc, i)
        end
    end
    order = sortperm(loc, by = i -> z[i], rev = true)
    kept = Int[]
    for i in loc[order]
        all(abs(i - k) >= sep for k in kept) && push!(kept, i)
    end
    return sort(kept)
end

"Read a template CSV written by template.py (offset_samples,value; comment lines with #)."
function read_template(path::AbstractString)
    vals = Float64[]
    for line in eachline(path)
        startswith(line, "#") && continue
        startswith(line, "offset") && continue
        isempty(strip(line)) && continue
        push!(vals, parse(Float64, split(line, ",")[2]))
    end
    isodd(length(vals)) || throw(ArgumentError("template $path: even number of values"))
    abs(vals[(length(vals) + 1) ÷ 2] - 1) < 1e-6 ||
        throw(ArgumentError("template $path: centre value is not 1"))
    return vals
end

# ------------------------------------------------------------- one event ----

"""
    process_event(w, p) -> NamedTuple

`w`: (samples, channels) raw counts of one event; `p`: the parameters (a
NamedTuple, see `pass`). Returns the per-event quantities, the summed waveform,
the per-channel charges and times, and the SPE candidates.
"""
function process_event(w::AbstractMatrix, p)
    ns, nch = size(w)
    x = Float64.(w) .- p.baselines'                         # (ns, nch)
    # 2. common mode, provisional pulse, W
    m = common_mode(x)
    S0 = vec(sum(x .- m, dims = 2))
    S0 .-= median(@view S0[ns - p.nbase + 1:ns])
    _, imax0, _, _ = pulse_search(S0, p.nsm)
    lo0, hi0, _, _, _ = windows(imax0, p.a, p.b, ns)
    # 3, 4. estimate, subtraction, sum, baseline
    mt = p.cm == "interp" ? interpolate_cm(m, lo0, hi0, p.edge)[1] :
         p.cm == "median" ? m : zeros(ns)
    y = x .- mt
    S = vec(sum(y, dims = 2))
    B = median(@view S[ns - p.nbase + 1:ns])
    S .-= B
    s1 = rms(@view S[ns - p.nbase + 1:ns])
    # 5. pulse search, windows
    Sh, imax, H, irise = pulse_search(S, p.nsm)
    lo, hi, clo, chi, short = windows(imax, p.a, p.b, ns)
    # 6. per-channel charges and times
    q, qctrl, imx = channel_charges(y, p.nbase, lo, hi, clo, chi, p.nsm, p.dt)
    # 7. SPE candidates per channel, on x minus its own median
    Ht = (length(p.T) - 1) ÷ 2
    rlo = max(ns - p.ntail + 1, Ht + 1)
    rhi = ns - Ht
    spe_ev = Int[]; spe_ch = Int[]; spe_i = Int[]; spe_A = Float64[]; spe_res = Float64[]
    sigz = zeros(Float64, nch)
    for c in 1:nch
        xc = x[:, c] .- median(@view x[:, c])
        z = matched_filter(xc, p.T)
        sz = clipped_rms(@view z[rlo:rhi])
        sigz[c] = sz
        for i in spe_candidates(z, rlo, rhi, p.nsig * sz, 2 * Ht)
            A = z[i]
            res = rms(xc[i - Ht:i + Ht] .- A .* p.T)
            push!(spe_ch, c); push!(spe_i, i); push!(spe_A, A); push!(spe_res, res)
        end
    end
    return (t_max = p.t[imax], H = H, t_rise = irise == 0 ? NaN : p.t[irise],
            win = (lo, hi), ctrl = (clo, chi), short = short, baseline = B, rms = s1,
            S = S, q = q, qctrl = qctrl, tc = [p.t[i] for i in imx],
            spe_ch = spe_ch, spe_t = [p.t[i] for i in spe_i], spe_A = spe_A, spe_res = spe_res,
            spe_sigma = sigz)
end

# ------------------------------------------------------------------ pass ----

"Thermometers of a run from temperatures.txt next to the data directory, or nothing."
function thermometers(runpath::AbstractString, name::AbstractString)
    f = joinpath(dirname(dirname(abspath(runpath))), "temperatures.txt")
    isfile(f) || return nothing
    vals = Dict{String,Vector{String}}()
    for line in eachline(f)
        occursin("=", line) || continue
        k, v = strip.(split(line, "=", limit = 2))
        vals[k] = [strip(s, [' ', '\'', '"']) for s in split(strip(v, ['[', ']', ' ']), ",")]
    end
    all(haskey(vals, k) for k in ("names", "T1", "T2", "T3")) || return nothing
    i = findfirst(==(name), vals["names"])
    i === nothing && return nothing
    return (T1 = parse(Float64, vals["T1"][i]), T2 = parse(Float64, vals["T2"][i]),
            T3 = parse(Float64, vals["T3"][i]))
end

"""
    pass(runpath, out; template, events=nothing, cm="interp", edge=32, nbase=250,
         nsm=64, window=(1000.0, 4000.0), tail_us=10.0, nsig=4.5, block=200,
         scale=4, verbose=true) -> out

Run the pass over `runpath` and write `out` (HDF5). `template` is the path of the
SPE template CSV. `events` is a 1-based range or `nothing` for all.
"""
function pass(runpath::AbstractString, out::AbstractString; template::AbstractString,
              events = nothing, cm::AbstractString = "interp", edge::Integer = 32,
              nbase::Integer = 250, nsm::Integer = 64, window = (1000.0, 4000.0),
              tail_us::Real = 10.0, nsig::Real = 4.5, block::Integer = 200,
              scale::Integer = 4, verbose::Bool = true)
    cm in ("interp", "median", "none") || throw(ArgumentError("cm must be interp, median or none"))
    T = read_template(template)
    run = openrun(runpath)
    try
        ev = events === nothing ? (1:run.nev) : events
        checkbounds(Bool, 1:run.nev, ev) || throw(BoundsError(1:run.nev, ev))
        ns, nch, nev = run.nsamp, run.nch, length(ev)
        dt = run.tsample_ns
        p = (baselines = run.baselines, nbase = nbase, nsm = nsm,
             a = round(Int, window[1] / dt), b = round(Int, window[2] / dt),
             cm = cm, edge = edge, dt = dt, t = time_axis(run), T = T,
             ntail = round(Int, tail_us * 1000 / dt), nsig = nsig)
        # per-event arrays
        t_max = zeros(Float32, nev); Hs = zeros(Float32, nev); t_rise = zeros(Float32, nev)
        win_lo = zeros(Int32, nev); win_hi = zeros(Int32, nev)
        ctrl_lo = zeros(Int32, nev); ctrl_hi = zeros(Int32, nev)
        base = zeros(Float32, nev); rmss = zeros(Float32, nev)
        sums = zeros(Int16, ns, nev)
        qc = zeros(Float32, nch, nev); qcc = zeros(Float32, nch, nev); tc = zeros(Float32, nch, nev)
        spes = zeros(Float32, nch, nev)
        spe_ev = Int32[]; spe_ch = Int16[]; spe_t = Float32[]; spe_A = Float32[]; spe_r = Float32[]
        nclip = 0
        lim = 32767 ÷ scale
        k = 0
        for a in first(ev):block:last(ev)
            b = min(last(ev), a + block - 1)
            W = waveforms(run, a:b)
            for e in a:b
                k += 1
                r = process_event(@view(W[:, :, e - a + 1]), p)
                t_max[k] = r.t_max; Hs[k] = r.H; t_rise[k] = r.t_rise
                win_lo[k], win_hi[k] = r.win .- 1          # 0-based in the file
                ctrl_lo[k], ctrl_hi[k] = r.ctrl .- 1
                base[k] = r.baseline; rmss[k] = r.rms
                for i in 1:ns
                    v = round(r.S[i] * scale)
                    if abs(v) > lim
                        nclip += 1
                        v = clamp(v, -lim, lim)
                    end
                    sums[i, k] = Int16(v)
                end
                qc[:, k] .= r.q; qcc[:, k] .= r.qctrl; tc[:, k] .= r.tc
                spes[:, k] .= r.spe_sigma
                for j in eachindex(r.spe_ch)
                    push!(spe_ev, e - 1); push!(spe_ch, run.channels[r.spe_ch[j]])
                    push!(spe_t, r.spe_t[j]); push!(spe_A, r.spe_A[j]); push!(spe_r, r.spe_res[j])
                end
            end
            verbose && println("  events $(a - 1)..$(b - 1) done, $(length(spe_ev)) SPE candidates")
        end
        # write
        mkpath(dirname(abspath(out)))
        h5open(out, "w") do f
            g = create_group(f, "events")
            g["index"] = Int32.(collect(ev) .- 1)
            g["trigger_id"] = Int64.(run.trigger_id[ev])
            g["timestamp_ns"] = run.timestamp_ns[ev]
            g["t_max"] = t_max; g["H"] = Hs; g["t_rise"] = t_rise
            g["win_lo"] = win_lo; g["win_hi"] = win_hi; g["ctrl_lo"] = ctrl_lo; g["ctrl_hi"] = ctrl_hi
            g["baseline"] = base; g["rms"] = rmss
            ds = create_dataset(f, "sum", Int16, (ns, nev), chunk = (ns, min(nev, 50)), deflate = 4)
            ds[:, :] = sums
            attrs(ds)["scale"] = scale
            attrs(ds)["clipped"] = nclip
            for (nm, arr) in (("qc", qc), ("qc_ctrl", qcc), ("tc", tc), ("spe_sigma", spes))
                d = create_dataset(f, nm, Float32, (nch, nev), chunk = (nch, min(nev, 500)), deflate = 4)
                d[:, :] = arr
            end
            s = create_group(f, "spe")
            s["event"] = spe_ev; s["channel"] = spe_ch; s["time"] = spe_t; s["A"] = spe_A
            s["residual"] = spe_r
            A = attrs(f)
            A["run"] = run.name; A["path"] = abspath(runpath)
            A["T_K"] = run.T_K; A["V_bias"] = run.V_bias; A["thr_adc"] = run.thr_adc
            A["start_time"] = run.start_time; A["pretrigger"] = run.pretrigger
            A["tsample_ns"] = dt; A["baselines"] = run.baselines; A["channels"] = run.channels
            th = thermometers(runpath, run.name)
            if th !== nothing
                A["T1"] = th.T1; A["T2"] = th.T2; A["T3"] = th.T3
            end
            tid = run.trigger_id
            dur = (run.timestamp_ns[end] - run.timestamp_ns[1]) * 1e-9
            A["trigger_rate_hz"] = dur > 0 ? mod(tid[end] - tid[1], 2^24) / dur : NaN
            A["cm"] = cm; A["edge"] = edge; A["baseline_window"] = nbase; A["smooth"] = nsm
            A["window_ns"] = collect(Float64, window); A["tail_us"] = Float64(tail_us); A["nsigma"] = Float64(nsig)
            A["template_file"] = abspath(template); A["template"] = T
            A["events_first"] = first(ev) - 1; A["events_last"] = last(ev) - 1
            A["jscan_version"] = string(pkgversion(JScan))
            A["date"] = string(now())
            commit = try
                readchomp(`git -C $(dirname(abspath(runpath))) rev-parse --short HEAD`)
            catch
                try readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`) catch; "" end
            end
            A["git_commit"] = commit
        end
    finally
        close(run.fid)
    end
    return out
end
