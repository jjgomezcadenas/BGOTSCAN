using Test, HDF5, JScan

const FIXTURE = joinpath(@__DIR__, "data", "185K_39d5V_th19.h5")
const NEV, NCH, NS, PRE = 3, 64, 300, 20
# pattern of make_fixture.py, with 1-based Julia indices (s, c, e)
pattern(s, c, e) = UInt16(32768 + 100 * (e - 1) + 3 * (c - 1) + mod(s - 1, 11))
expected(ss, cs, es) = [pattern(s, c, e) for s in ss, c in cs, e in es]

@testset "parse_runname" begin
    @test parse_runname("128K_37V_th19") == (T_K = 128, V_bias = 37.0, thr_adc = 19)
    @test parse_runname("185K_39d5V_th19") == (T_K = 185, V_bias = 39.5, thr_adc = 19)
    @test parse_runname("196K_40d5V_th20") == (T_K = 196, V_bias = 40.5, thr_adc = 20)
    @test parse_runname("128K_37V_th19").V_bias isa Float64
    @test_throws ArgumentError parse_runname("128K_37V")
    @test_throws ArgumentError parse_runname("128K_37V_th19.h5")
    @test_throws ArgumentError parse_runname("run_20261002_145441")
end

@testset "register_lzf" begin
    # The plugin cannot be unloaded, so this testset must run before any
    # waveform read, and the "fails before" check only makes sense when no
    # plugin path comes from the environment.
    if !haskey(ENV, "HDF5_PLUGIN_PATH")
        h5open(FIXTURE, "r") do f
            @test_throws Exception read(f["waveforms"])
        end
    end
    @test_throws ArgumentError register_lzf(joinpath(@__DIR__, "no_such_dir"))
    @test_throws ArgumentError register_lzf(@__DIR__)      # exists, no plugin
    dir = register_lzf()
    @test dir == JScan.LZF_PLUGIN_DIR
    @test dir in JScan.plugin_paths()
    n = length(JScan.plugin_paths())
    @test register_lzf() == dir                             # second call harmless
    @test length(JScan.plugin_paths()) == n
    h5open(FIXTURE, "r") do f
        @test read(f["waveforms"]) == expected(1:NS, 1:NCH, 1:NEV)
    end
end

@testset "openrun" begin
    run = openrun(FIXTURE)
    @test run.name == "185K_39d5V_th19"
    @test run.path == FIXTURE
    @test (run.T_K, run.V_bias, run.thr_adc) == (185, 39.5, 19)
    @test (run.nev, run.nch, run.nsamp) == (NEV, NCH, NS)
    @test run.channels == collect(0:NCH - 1)
    @test run.channels isa Vector{Int}
    @test run.tsample_ns == 8.0
    @test run.sampling_mhz == 125.0
    @test run.pretrigger == PRE
    @test run.adc_to_volts == 6.1455e-5
    @test run.zin_ohm == 100.0
    @test run.baselines == 32768.0 .+ 0.5 .* (0:NCH - 1)
    @test run.baselines isa Vector{Float64}
    @test run.thresholds == fill(19, NCH)
    @test run.threshold_mode == "Relative"
    @test run.thresholds_abs == floor.(Int, 32768.0 .+ 0.5 .* (0:NCH - 1) .+ 19)
    @test run.edge == "RISE"
    @test run.itl_logic == "Majority"
    @test run.itl_majority == 2
    @test run.start_time == "2026-10-02T14:47:45"
    @test run.timestamp_ns == [0.0, 1.0e6, 2.5e6]
    @test run.trigger_id == [10, 11, 13]
    @test run.trigger_id isa Vector{Int}
    @test run.fid isa HDF5.File && isopen(run.fid)
    @test run.wf isa HDF5.Dataset && size(run.wf) == (NS, NCH, NEV)
    close(run.fid)
    @test !isopen(run.fid)

    # do-block form closes the file, also on error
    fid = nothing
    r = openrun(FIXTURE) do run
        fid = run.fid
        run.nev
    end
    @test r == NEV && !isopen(fid)
    @test_throws ErrorException openrun(FIXTURE) do run
        fid = run.fid
        error("boom")
    end
    @test !isopen(fid)
end

@testset "waveforms" begin
    openrun(FIXTURE) do run
        w = waveforms(run, 1:NEV)
        @test w isa Array{UInt16,3} && size(w) == (NS, NCH, NEV)
        @test w == expected(1:NS, 1:NCH, 1:NEV)
        @test waveforms(run, 2:2) == expected(1:NS, 1:NCH, 2:2)     # inner range
        @test waveforms(run, 2:3) == expected(1:NS, 1:NCH, 2:3)     # last block
        m = waveforms(run, 3)
        @test m isa Matrix{UInt16} && size(m) == (NS, NCH)
        @test m == expected(1:NS, 1:NCH, 3:3)[:, :, 1]
        @test_throws BoundsError waveforms(run, 0:1)
        @test_throws BoundsError waveforms(run, 3:4)
        @test_throws BoundsError waveforms(run, 4)
        @test_throws BoundsError waveforms(run, 0)
    end
end

@testset "time_axis" begin
    openrun(FIXTURE) do run
        t = time_axis(run)
        @test t isa Vector{Float64} && length(t) == NS
        @test t[PRE + 1] == 0.0
        @test t[1] == -PRE * 8.0
        @test all(diff(t) .== 8.0)
    end
end

const REAL = normpath(joinpath(@__DIR__, "..", "..", "..",
    "2026_10_02_bgo_cryo_tempscan", "data_cryo", "128K_37V_th19.h5"))
@testset "real run" begin
    if isfile(REAL)
        openrun(REAL) do run
            @test (run.T_K, run.V_bias, run.thr_adc) == (128, 37.0, 19)
            @test (run.nch, run.nsamp) == (64, 3000) && run.nev >= 10
            @test run.pretrigger == 200 && run.tsample_ns == 8.0
            w = waveforms(run, 1:10)
            @test size(w) == (3000, 64, 10)
            @test all(20000 .< w .< 45000)      # 16-bit ADC near mid-range
        end
    else
        @info "real run not found, skipping" REAL
    end
end

# ---- the pass ----------------------------------------------------------------
"Minimal JSON reader for the dumps of the Python prototypes (objects, arrays, numbers, strings)."
function JSON_parse(s::AbstractString)
    pos = Ref(1)
    skip() = (while pos[] <= lastindex(s) && isspace(s[pos[]]); pos[] = nextind(s, pos[]); end)
    function value()
        skip()
        c = s[pos[]]
        if c == '{'
            pos[] = nextind(s, pos[]); d = Dict{String,Any}()
            skip()
            if s[pos[]] == '}'; pos[] = nextind(s, pos[]); return d; end
            while true
                skip(); k = value(); skip(); @assert s[pos[]] == ':'; pos[] = nextind(s, pos[])
                d[k] = value(); skip()
                if s[pos[]] == ','; pos[] = nextind(s, pos[]); continue; end
                @assert s[pos[]] == '}'; pos[] = nextind(s, pos[]); return d
            end
        elseif c == '['
            pos[] = nextind(s, pos[]); a = Any[]
            skip()
            if s[pos[]] == ']'; pos[] = nextind(s, pos[]); return a; end
            while true
                push!(a, value()); skip()
                if s[pos[]] == ','; pos[] = nextind(s, pos[]); continue; end
                @assert s[pos[]] == ']'; pos[] = nextind(s, pos[]); return a
            end
        elseif c == '"'
            i = nextind(s, pos[]); j = i
            while s[j] != '"'; j = nextind(s, j); end
            pos[] = nextind(s, j); return String(s[i:prevind(s, j)])
        else
            j = pos[]
            while j <= lastindex(s) && (isdigit(s[j]) || s[j] in "+-.eE"); j = nextind(s, j); end
            tok = s[pos[]:prevind(s, j)]; pos[] = j
            tok == "NaN" && return NaN
            return occursin(r"[.eE]", tok) ? parse(Float64, tok) : parse(Int, tok)
        end
    end
    return value()
end

include("test_pass.jl")
