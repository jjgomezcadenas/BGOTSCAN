# Tests of the pass (docs/manual/jscan_doc.tex, section "The pass"): the pure functions
# on synthetic arrays, and the end-to-end comparison with the Python prototypes
# on the real-event fixture when it and the dumps are present.

using Statistics
using JScan: rms, clipped_rms

const REAL_FIXTURE = joinpath(@__DIR__, "data", "real", "197K_39V_th20.h5")
const DUMP_SEL = joinpath(@__DIR__, "data", "real", "197K_39V_th20_selection.json")
const DUMP_SPE = joinpath(@__DIR__, "data", "real", "197K_39V_th20_spe.json")
const TEMPLATE = joinpath(@__DIR__, "..", "..", "spe_template_195K_41V_th20.csv")

@testset "moving_average" begin
    v = Float64.(1:10)
    @test moving_average(v, 1) == v
    # window i - n÷2 : i + (n-1)÷2, edges with the samples available
    m3 = moving_average(v, 3)
    @test m3[1] ≈ mean(v[1:2])
    @test m3[5] ≈ mean(v[4:6])
    @test m3[10] ≈ mean(v[9:10])
    m4 = moving_average(v, 4)
    @test m4[5] ≈ mean(v[3:6])                 # 5-2 : 5+1
    @test m4[1] ≈ mean(v[1:2])
    @test moving_average(zeros(20), 64) == zeros(20)
end

@testset "common_mode" begin
    x = [1.0 5.0 3.0; 2.0 2.0 9.0; 0.0 0.0 0.0]
    @test common_mode(x) == [3.0, 2.0, 0.0]
end

@testset "dead_cm" begin
    ns, nch = 400, 6
    wander = 2.0 .* sin.(2π .* (1:ns) ./ 400)          # slow common mode
    x = zeros(ns, nch)
    for c in 1:nch
        x[:, c] .= wander
    end
    x[150:160, 1] .+= 50.0                               # a pulse in a live channel
    x[:, 5] .+= 0.3 .* randn(ns); x[:, 6] .+= 0.3 .* randn(ns)
    cm = dead_cm(x, [5, 6], 25)
    @test maximum(abs.(cm .- wander)) < 0.3              # follows the wander, not the pulse
    @test dead_cm(x, Int[], 25) == zeros(ns)
    S = vec(sum(x .- cm, dims = 2))
    @test S[155] > 45 && abs(S[300]) < 2                 # the pulse survives, the wander is gone
end

@testset "channel_baseline_line!" begin
    ns = 1000
    y = zeros(ns, 2)
    y[:, 1] .= 1.0 .+ 0.01 .* (1:ns)            # a linear drift
    y[:, 2] .= 5.0
    y[300:320, 1] .+= 40.0                        # a pulse, between the two regions
    @test channel_baseline_line!(y, 150, 600, 700) == 0
    @test maximum(abs.(y[1:150, 1])) < 1e-9 && maximum(abs.(y[600:700, 1])) < 1e-9
    @test abs(y[310, 1] - 40.0) < 1e-9             # the pulse survives
    @test maximum(abs.(y[:, 2])) < 1e-9
    z = ones(ns, 1)
    @test channel_baseline_line!(z, 150, 950, 1100) == 1   # post region beyond the record: untouched
    @test all(z .== 1.0)
end

@testset "interpolate_cm" begin
    m = collect(0.0:0.1:9.9)                    # 100 samples, slope 0.1
    mt, L, R = interpolate_cm(m, 41, 60, 4)
    @test L ≈ mean(m[37:40]) && R ≈ mean(m[61:64])
    @test mt[41] ≈ L && mt[60] ≈ R
    @test mt[1:40] == m[1:40] && mt[61:end] == m[61:end]
    @test all(abs.(mt[41:60] .- m[41:60]) .< 0.5)   # a line through a line
    # window at the start: the right end serves for both
    mt2, L2, R2 = interpolate_cm(m, 1, 10, 4)
    @test L2 == R2 ≈ mean(m[11:14])
    @test all(mt2[1:10] .≈ R2)
    # window at the end
    mt3, L3, R3 = interpolate_cm(m, 91, 100, 4)
    @test R3 == L3 ≈ mean(m[87:90])
end

@testset "pulse_search" begin
    S = zeros(500)
    S[200:260] .= [100.0 * exp(-(i - 200) / 20) for i in 200:260]  # a decaying pulse from 200
    Sh, imax, H, irise = pulse_search(S, 1)
    @test imax == 200 && H ≈ 100.0 && irise == 199
    # above half maximum from the first sample: no rising edge
    S2 = [100.0 * exp(-i / 200) for i in 1:500]
    _, imax2, _, irise2 = pulse_search(S2, 1)
    @test imax2 == 1 && irise2 == 0
end

@testset "windows" begin
    lo, hi, clo, chi, short = windows(1000, 125, 500, 3000)
    @test (lo, hi) == (875, 1500) && (clo, chi) == (2375, 3000) && !short
    lo, hi, clo, chi, short = windows(1, 125, 500, 3000)         # at the start
    @test (lo, hi) == (1, 501) && (clo, chi) == (2500, 3000) && !short
    lo, hi, clo, chi, short = windows(2990, 125, 500, 3000)      # at the end
    @test (lo, hi) == (2865, 3000) && clo == 1 && chi == 136 && !short
    lo, hi, clo, chi, short = windows(1500, 125, 1400, 3000)     # no room: short
    @test short && chi - clo + 1 < hi - lo + 1
end

@testset "channel_charges" begin
    ns, nch = 1000, 3
    y = zeros(ns, nch)
    y[:, 2] .= 2.0                      # constant: baseline 2, charge 0
    y[300:400, 3] .= 5.0                # 101 samples of 5 in channel 3
    q, qctrl, imx = channel_charges(y, 100, 250, 450, 700, 900, 1, 8.0)
    @test q[1] == 0 && q[2] == 0 && q[3] ≈ 8.0 * 5 * 101
    @test qctrl[3] == 0
    @test 300 <= imx[3] <= 400
end

@testset "matched_filter and spe_candidates" begin
    T = [0.0, 0.2, 0.6, 1.0, 0.5, 0.2, 0.0]
    x = zeros(200)
    x[100 - 3:100 + 3] .= 3.0 .* T             # a pulse of height 3 at 100
    z = matched_filter(x, T)
    @test z[100] ≈ 3.0
    @test argmax(z) == 100
    @test z[1:3] == zeros(3) && z[end - 2:end] == zeros(3)
    c = spe_candidates(z, 4, 197, 1.0, 7)
    @test c == [100]
    # two pulses closer than sep: only the higher survives; farther: both
    x2 = copy(x); x2[104 - 3:104 + 3] .+= 2.0 .* T
    z2 = matched_filter(x2, T)
    @test length(spe_candidates(z2, 4, 197, 1.0, 7)) == 1
    x3 = copy(x); x3[150 - 3:150 + 3] .= 2.0 .* T
    @test spe_candidates(matched_filter(x3, T), 4, 197, 1.0, 7) == [100, 150]
    # below threshold: nothing
    @test spe_candidates(z, 4, 197, 5.0, 7) == Int[]
end

@testset "rms and clipped_rms" begin
    v = [1.0, 2.0, 3.0, 4.0]
    @test rms(v) ≈ sqrt(mean((v .- 2.5) .^ 2))
    w = vcat(randn(1000), [100.0])
    @test clipped_rms(w) < 1.5 && rms(w) > 2.5
end

@testset "read_template" begin
    tmp = tempname() * ".csv"
    open(tmp, "w") do f
        println(f, "# test"); println(f, "offset_samples,value")
        for (j, v) in zip(-2:2, [0.1, 0.5, 1.0, 0.4, 0.1]); println(f, "$j,$v"); end
    end
    @test read_template(tmp) == [0.1, 0.5, 1.0, 0.4, 0.1]
    open(tmp, "w") do f
        println(f, "offset_samples,value"); println(f, "0,0.5")
    end
    @test_throws ArgumentError read_template(tmp)
    rm(tmp)
end

@testset "process_event on the synthetic fixture" begin
    # the pattern fixture holds no pulse; the pass must run and give finite numbers
    T = read_template(TEMPLATE)
    openrun(FIXTURE) do run
        w = waveforms(run, 1)
        p = (baselines = run.baselines, nbase = 50, nsm = 8, a = 10, b = 40, cm = "interp",
             edge = 8, dead_cols = Int[], cm_smooth = 25, chbase = "end", pre_hi = 15, post_lo = 60, post_hi = 80,
             dt = run.tsample_ns, t = time_axis(run), T = T, ntail = 100, nsig = 4.5)
        r = process_event(w, p)
        @test length(r.S) == run.nsamp && all(isfinite, r.S)
        @test length(r.q) == run.nch && length(r.tc) == run.nch
        @test r.win[1] >= 1 && r.win[2] <= run.nsamp
        @test length(r.spe_sigma) == run.nch
    end
end

@testset "pass end to end against the Python dumps" begin
    if isfile(REAL_FIXTURE) && isfile(DUMP_SEL) && isfile(DUMP_SPE)
        out = tempname() * ".h5"
        pass(REAL_FIXTURE, out; template = TEMPLATE, verbose = false)
        sel = JSON_parse(read(DUMP_SEL, String))
        spe = JSON_parse(read(DUMP_SPE, String))
        h5open(out, "r") do f
            tmax = read(f["events/t_max"]); H = read(f["events/H"])
            B = read(f["events/baseline"]); s1 = read(f["events/rms"])
            wlo = read(f["events/win_lo"]); whi = read(f["events/win_hi"])
            qc = read(f["qc"]); sumw = read(f["sum"]); scale = read_attribute(f["sum"], "scale")
            for (k, ev) in enumerate(sel)
                @test tmax[k] ≈ ev["t_max"] atol = 1e-3
                @test H[k] ≈ ev["H"] rtol = 1e-3
                @test B[k] ≈ ev["baseline"] atol = 1e-2
                @test s1[k] ≈ ev["rms"] rtol = 1e-3
                @test wlo[k] == ev["win"][1] && whi[k] == ev["win"][2]
                @test maximum(abs.(qc[:, k] .- Float32.(ev["q"]))) < 1e-3 * maximum(abs.(ev["q"]))
                @test maximum(abs.(Float64.(sumw[:, k]) ./ scale .- ev["S"])) <= 0.5 / scale + 1e-6
            end
            sev = read(f["spe/event"]); sch = read(f["spe/channel"]); sA = read(f["spe/A"]); st = read(f["spe/time"])
            for (k, ev) in enumerate(spe)
                for (chs, cands) in ev["channels"]
                    c = parse(Int, chs)
                    sel_rows = findall((sev .== ev["event"]) .& (sch .== c))
                    @test length(sel_rows) == length(cands)
                    for (row, cd) in zip(sort(sel_rows, by = i -> st[i]), sort(cands, by = x -> x["time"]))
                        @test st[row] ≈ cd["time"] atol = 1e-3
                        @test sA[row] ≈ cd["A"] rtol = 1e-3
                    end
                end
            end
        end
        rm(out)
    else
        @info "real-event fixture or dumps absent: end-to-end test skipped"
    end
end
