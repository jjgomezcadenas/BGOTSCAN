#!/usr/bin/env julia
# The pass over one run. See docs/manual/jscan_doc.tex, section "The pass".
#
#   julia --project=src/jscan src/jscan/pass.jl RUN [--out DIR] [--events a:b]
#         [--template FILE] [--cm dead|interp|median|none] [--dead LIST] [--cm-smooth 125]
#         [--chbase end|line] [--guard 0.4] [--post-window 10,12] [--edge 32] [--baseline-window 250]
#         [--smooth 64] [--window 1000,4000] [--tail 10] [--nsigma 4.5] [--block 200]

using JScan

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const DATA = joinpath(REPO, "2026_10_02_bgo_cryo_tempscan", "data_cryo")

function usage()
    println(strip(read(@__FILE__, String)[1:findfirst("using JScan", read(@__FILE__, String))[1] - 1]))
    exit(1)
end

function parse_args(args)
    opts = Dict{String,Any}("out" => joinpath(REPO, "pass"), "events" => nothing, "template" => nothing,
                            "cm" => "dead", "dead" => JScan.DEAD, "cm-smooth" => 125,
                            "chbase" => "end", "guard" => 0.4, "post-window" => (10.0, 12.0),
                            "edge" => 32, "baseline-window" => 250, "smooth" => 64,
                            "window" => (1000.0, 4000.0), "tail" => 10.0, "nsigma" => 4.5, "block" => 200)
    pos = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if startswith(a, "--")
            key = a[3:end]
            haskey(opts, key) || (println("unknown option $a"); usage())
            i += 1
            i <= length(args) || (println("missing value for $a"); usage())
            v = args[i]
            opts[key] = key == "events" ? (let p = split(v, ":"); parse(Int, p[1]):parse(Int, p[2]) end) :
                        key in ("window", "post-window") ? Tuple(parse.(Float64, split(v, ","))) :
                        key == "dead" ? parse.(Int, split(v, ",")) :
                        key in ("edge", "baseline-window", "smooth", "block", "cm-smooth") ? parse(Int, v) :
                        key in ("tail", "nsigma", "guard") ? parse(Float64, v) : v
        else
            push!(pos, a)
        end
        i += 1
    end
    length(pos) == 1 || usage()
    return pos[1], opts
end

function main()
    runarg, o = parse_args(ARGS)
    path = isfile(runarg) ? runarg :
           isfile(joinpath(DATA, runarg)) ? joinpath(DATA, runarg) :
           isfile(joinpath(DATA, runarg * ".h5")) ? joinpath(DATA, runarg * ".h5") :
           (println("run not found: $runarg"); exit(1))
    name = first(splitext(basename(path)))
    template = o["template"] === nothing ? joinpath(REPO, "src", "spe_template_$name.csv") : o["template"]
    isfile(template) || (println("template not found: $template (make it with template.py --write, or give --template)"); exit(1))
    register_lzf()
    out = joinpath(o["out"], name * "_pass.h5")
    println("pass over $path -> $out")
    t0 = time()
    pass(path, out; template = template, events = o["events"], cm = o["cm"], dead = o["dead"],
         cm_smooth = o["cm-smooth"], edge = o["edge"], chbase = o["chbase"], guard_us = o["guard"],
         post_window = o["post-window"],
         nbase = o["baseline-window"], nsm = o["smooth"], window = o["window"], tail_us = o["tail"],
         nsig = o["nsigma"], block = o["block"])
    println("done in $(round(time() - t0, digits = 1)) s")
end

main()
