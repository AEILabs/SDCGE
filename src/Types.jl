# Usage:
#   data = init_data()
#   default_sets!(data)              # loads N activities/products and SAM account groups
#   setup_sam_accounts!(data)        # builds named SAM account lists
#
# N defaults to 100 (P001..P100).  For any other sector list, fill `data.sets[:i]`
# and the memberships :cr :lv :e :ft :fd first (see `read_sets_csv!` in SAM.jl);
# `default_sets!` uses `get!` throughout, so pre-populated sets always win.
#
# This file contains only data containers and set/account defaults.

mutable struct LinkageData
    sets::Dict{Symbol,Vector}
    par::Dict{Symbol,Any}
    sam_accounts::Dict{Symbol,Vector{String}}
    sam::Matrix{Float64}
    balanced_sam::Matrix{Float64}
    sam_index::Dict{String,Int}
    metadata::Dict{Symbol,Any}
end

function LinkageData()
    return LinkageData(
        Dict{Symbol,Vector}(),
        Dict{Symbol,Any}(),
        Dict{Symbol,Vector{String}}(),
        zeros(0,0),
        zeros(0,0),
        Dict{String,Int}(),
        Dict{Symbol,Any}(),
    )
end

init_data() = LinkageData()

"""Fill in the LINKAGE sets that have not been supplied yet.

`S[:i]` (activities = products) defaults to the 100 codes P001..P100; any other
sector list must already be in `data.sets[:i]` — every `get!` below then leaves the
caller's sets untouched.  The positional memberships (crops = first 10, livestock =
next 10, energy = 71:75, fertiliser = 76:78) only make sense for that default list,
so with `|S[:i]| != 100` the sets :cr :lv :e :ft :fd must be supplied explicitly.
"""
function default_sets!(data::LinkageData)
    S = data.sets
    products = get!(S, :i, ["P" * lpad(string(n), 3, "0") for n in 1:100])
    get!(S, :j, S[:i])
    get!(S, :k, S[:i])
    get!(S, :r, ["R1", "R2", "R3", "R4"])
    get!(S, :rp, S[:r])
    get!(S, :v, ["Old", "New"])
    get!(S, :l, ["UnSkLab", "SkLab"])
    get!(S, :ul, ["UnSkLab"])
    get!(S, :sl, ["SkLab"])
    get!(S, :h, ["HH"])
    get!(S, :f, ["Gov", "Inv"])
    get!(S, :in, ["HH", "Gov", "Inv"])
    get!(S, :t, [1, 2, 3])

    n = length(products)
    if n != 100 && !all(haskey(S, s) for s in (:cr, :lv, :e, :ft, :fd))
        missing_sets = [s for s in (:cr, :lv, :e, :ft, :fd) if !haskey(S, s)]
        error("default_sets!: with |S[:i]| = $(n) != 100 the sets " *
              join(missing_sets, ", ") * " must be supplied (e.g. via " *
              "read_sets_csv!); the positional defaults apply only to P001..P100.")
    end
    # Lazy `get!(f, dict, key)` so the positional slices are never evaluated when the
    # caller supplied the set — `products[71:75]` would otherwise throw for N < 75.
    get!(() -> products[1:10],  S, :cr)
    get!(() -> products[11:20], S, :lv)
    get!(S, :ag, vcat(S[:cr], S[:lv]))
    get!(S, :ip, [x for x in S[:i] if !(x in S[:ag])])
    get!(() -> products[71:75], S, :e)
    get!(() -> products[76:78], S, :ft)
    get!(() -> products[1:10],  S, :fd)
    get!(S, :nf, [x for x in S[:i] if !(x in S[:ag])])   # == products[21:100] at N = 100
    get!(S, :nnft, [x for x in S[:i] if !(x in S[:ft])])
    get!(S, :nnfd, [x for x in S[:i] if !(x in S[:fd])])
    get!(S, :gz, ["national", "urban", "rural"])
    get!(S, :gs, ["urban", "rural"])
    return data
end
