# Results export utilities for solved LCGE/LINKAGE JuMP models.
#
# Typical use after solving:
#
#     solve_model!(m)
#     export_results!(m, data; outdir="results")
#
# Or:
#
#     solve_and_save!(m, data; outdir="results")
#
# The code avoids a hard dependency on CSV.jl. CSV files are written directly
# using Base I/O so the existing Project.toml remains lightweight.

using Dates

function _lcge_csv_escape(x)
    s = string(x)
    if occursin(",", s) || occursin("\"", s) || occursin("\n", s) || occursin("\r", s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function _lcge_write_csv(path::AbstractString, header, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(_lcge_csv_escape.(header), ","))
        for row in rows
            println(io, join(_lcge_csv_escape.(row), ","))
        end
    end
    return path
end

function _lcge_value(v::JuMP.VariableRef)
    try
        return JuMP.value(v)
    catch
        return missing
    end
end

function _lcge_start_value(v::JuMP.VariableRef)
    try
        sv = start_value(v)
        return sv === nothing ? missing : sv
    catch
        return missing
    end
end

function _lcge_index_string(idx)
    if idx isa Tuple
        return join(string.(idx), "|")
    elseif idx isa CartesianIndex
        return join(string.(Tuple(idx)), "|")
    else
        return string(idx)
    end
end

# Map a container position (as returned by `eachindex`) to its axis labels.
# DenseAxisArray positions are CartesianIndex/Int; SparseAxisArray keys are
# already the label tuples.
function _lcge_index_label(obj, idx)
    if hasproperty(obj, :axes)
        ax = obj.axes
        try
            if idx isa CartesianIndex
                return join((string(ax[d][idx[d]]) for d in 1:length(idx)), "|")
            elseif idx isa Integer && length(ax) == 1
                return string(ax[1][idx])
            end
        catch
        end
    end
    return _lcge_index_string(idx)
end

# True once optimize! has been called; querying values before that emits one
# JuMP warning per variable, so callers skip value lookups instead.
function _lcge_has_solution(m)
    try
        return JuMP.termination_status(m) != JuMP.MOI.OPTIMIZE_NOT_CALLED
    catch
        return false
    end
end

function _lcge_collect_variable_rows!(rows::Vector{Vector{Any}}, name::Symbol, obj; solved::Bool=true)
    if obj isa JuMP.VariableRef
        val = solved ? _lcge_value(obj) : missing
        start = _lcge_start_value(obj)
        diff = (val === missing || start === missing) ? missing : val - start
        pct = (val === missing || start === missing || abs(start) <= 1.0e-12) ? missing : 100.0 * (val - start) / start
        push!(rows, Any[string(name), "", val, start, diff, pct])
        return rows
    end

    # Iterate over JuMP containers (DenseAxisArray, SparseAxisArray, Array).
    # `pairs()` does not work on DenseAxisArray, so use eachindex + indexing.
    try
        for idx in eachindex(obj)
            var = obj[idx]
            if var isa JuMP.VariableRef
                val = solved ? _lcge_value(var) : missing
                start = _lcge_start_value(var)
                diff = (val === missing || start === missing) ? missing : val - start
                pct = (val === missing || start === missing || abs(start) <= 1.0e-12) ? missing : 100.0 * (val - start) / start
                push!(rows, Any[string(name), _lcge_index_label(obj, idx), val, start, diff, pct])
            end
        end
    catch
        # Non-variable JuMP container (e.g., constraints) — skip silently.
    end
    return rows
end

"""Return a flat DataFrame with one row per model variable element.

Columns:
- `variable`: JuMP container name
- `index`: pipe-separated index label
- `value`: solution value
- `start_value`: initialization value
- `change_from_start`: solution minus start value
- `pct_change_from_start`: percentage change from start value
"""
function results_dataframe(m)
    rows = Vector{Vector{Any}}()
    solved = _lcge_has_solution(m)
    for (name, obj) in object_dictionary(m)
        _lcge_collect_variable_rows!(rows, name, obj; solved=solved)
    end
    return DataFrame(
        variable = [r[1] for r in rows],
        index = [r[2] for r in rows],
        value = [r[3] for r in rows],
        start_value = [r[4] for r in rows],
        change_from_start = [r[5] for r in rows],
        pct_change_from_start = [r[6] for r in rows],
    )
end

function _lcge_rows_from_dataframe(df::DataFrame)
    rows = Vector{Vector{Any}}()
    for r in eachrow(df)
        push!(rows, Any[r.variable, r.index, r.value, r.start_value, r.change_from_start, r.pct_change_from_start])
    end
    return rows
end

function _lcge_variable_rows(df::DataFrame, variable_name::AbstractString)
    rows = Vector{Vector{Any}}()
    for r in eachrow(df)
        if string(r.variable) == variable_name
            push!(rows, Any[r.index, r.value, r.start_value, r.change_from_start, r.pct_change_from_start])
        end
    end
    return rows
end

function _lcge_scalar_rows(m)
    rows = Vector{Vector{Any}}()
    solved = _lcge_has_solution(m)
    for (name, obj) in object_dictionary(m)
        if obj isa JuMP.VariableRef
            push!(rows, Any[string(name), solved ? _lcge_value(obj) : missing, _lcge_start_value(obj)])
        end
    end
    return rows
end

function _lcge_metadata_rows(m)
    rows = Vector{Vector{Any}}()
    push!(rows, Any["created_at", string(now())])
    try
        push!(rows, Any["termination_status", string(termination_status(m))])
    catch
        push!(rows, Any["termination_status", "not_available"])
    end
    try
        push!(rows, Any["primal_status", string(primal_status(m))])
    catch
        push!(rows, Any["primal_status", "not_available"])
    end
    try
        push!(rows, Any["dual_status", string(dual_status(m))])
    catch
        push!(rows, Any["dual_status", "not_available"])
    end
    return rows
end

"""Export solution results to CSV files.

Files written:
- `results_all_variables.csv`
- `results_scalars.csv`
- `results_metadata.csv`
- one file per common variable, for example `results_XP.csv`, `results_GDP.csv`

Returns the full results DataFrame.
"""
function export_results!(m, data::LinkageData; outdir::AbstractString="results")
    mkpath(outdir)

    # Always export and verify the SAM balance report alongside model results.
    export_sam_balance_report!(data; outdir=outdir)

    df = results_dataframe(m)

    _lcge_write_csv(
        joinpath(outdir, "results_all_variables.csv"),
        ["variable", "index", "value", "start_value", "change_from_start", "pct_change_from_start"],
        _lcge_rows_from_dataframe(df),
    )

    _lcge_write_csv(
        joinpath(outdir, "results_scalars.csv"),
        ["variable", "value", "start_value"],
        _lcge_scalar_rows(m),
    )

    _lcge_write_csv(
        joinpath(outdir, "results_metadata.csv"),
        ["item", "value"],
        _lcge_metadata_rows(m),
    )

    # Per-variable CSVs. Includes both legacy LINKAGE names (XD, XM, E, INVEST, ...)
    # and the actual JuMP variables in the current model (XDs, XDd, XMT, ES, FDInv, ...).
    # The `if !isempty(rows)` guard skips legacy names that aren't in the model.
    common_variables = [
        # Output & supply / demand
        "XP", "XA", "XDs", "XDd", "XDc", "XDf",
        "XMT", "XM1", "XM2", "XMc", "XMf", "XMgr",
        "ES", "WTFd", "WTFs", "WTFin", "WTFout",
        # Production / intermediate
        "ND", "VA", "XPv", "Kvd", "LV", "Td", "Fd",
        "XAp", "XAc", "XAf", "fert", "feed", "XEp",
        # Prices
        "PP", "PX", "PA", "PD", "PMT", "PE", "PM",
        "PC", "PAc", "PFD", "PVA", "PND", "PT", "PF", "NPT", "PTLnd",
        "W", "NW", "AVGW", "TW", "WMIN",
        "R", "NR", "TR", "CHIv",
        # Macro aggregates
        "GDP", "RGDP", "CPI", "PGDP", "GDPMPr",
        # Income / expenditure / savings
        "YH", "YD", "YC", "YSTAR", "SAV", "DeprY", "CPIH",
        "TY", "FY", "KY", "LY",
        # Government, investment, foreign
        "YG", "Sg", "RSg", "TarY", "RTarY",
        "FD", "FDInv", "GOVDEM", "INVDEM", "InvSh",
        "Sf", "ER",
        # Labor markets
        "LS", "UE", "MIGR", "PS",
        # Capital
        "KS", "KSs", "KActual", "KNorm", "K0", "RR",
        "WRR", "WPMg", "WXMg",
        # Legacy LINKAGE names (kept for backward compatibility; empty if absent)
        "XD", "XQ", "XM", "E", "D",
        "INVEST", "SAVE", "GOVREV", "GEXP", "TAXREV",
    ]

    for vname in common_variables
        rows = _lcge_variable_rows(df, vname)
        if !isempty(rows)
            _lcge_write_csv(
                joinpath(outdir, "results_$(vname).csv"),
                ["index", "value", "start_value", "change_from_start", "pct_change_from_start"],
                rows,
            )
        end
    end

    # Key economic indicators summary (one row per indicator).
    _lcge_write_csv(
        joinpath(outdir, "results_summary.csv"),
        ["indicator", "description", "value", "start_value"],
        _lcge_macro_summary_rows(m, data),
    )

    return df
end

"""Sum the product of two indexed JuMP variables over `idxset`.

Returns `(solution_total, start_total)`, or `missing` if either container is
absent (e.g. `ER` under `trade_closure = :balanced`).
"""
function _sum2(m, sym_a::Symbol, sym_b::Symbol, idxset)
    (haskey(m, sym_a) && haskey(m, sym_b)) || return missing
    total = 0.0; start = 0.0
    try
        for k in idxset
            a  = _lcge_value(m[sym_a][k...]);       b  = _lcge_value(m[sym_b][k...])
            sa = _lcge_start_value(m[sym_a][k...]); sb = _lcge_start_value(m[sym_b][k...])
            (a isa Real && b isa Real)   && (total += float(a) * float(b))
            (sa isa Real && sb isa Real) && (start += float(sa) * float(sb))
        end
        return (total, start)
    catch
        return missing
    end
end

"""Build a named summary of key macro / closure / factor-income indicators.

Each row maps a familiar economic concept to one or more JuMP variables (summing
where needed). Useful for quickly inspecting whether the model reproduces the
SAM benchmark and for comparing pre/post-shock outcomes.
"""
function _lcge_macro_summary_rows(m, data::LinkageData)
    rows = Vector{Vector{Any}}()
    S = data.sets
    i = S[:i]; r = S[:r]; h = S[:h]; ag = S[:ag]
    PAR = try parameters(data) catch; Dict{Symbol,Any}() end

    function _v(sym, idx=())
        haskey(m, sym) || return missing
        try
            return isempty(idx) ? _lcge_value(m[sym]) : _lcge_value(m[sym][idx...])
        catch
            return missing
        end
    end
    function _vs(sym, idx=())
        haskey(m, sym) || return missing
        try
            return isempty(idx) ? _lcge_start_value(m[sym]) : _lcge_start_value(m[sym][idx...])
        catch
            return missing
        end
    end
    function _sum(sym, idxset)
        haskey(m, sym) || return missing
        total = 0.0; start = 0.0
        try
            for k in idxset
                v  = _lcge_value(m[sym][k...])
                sv = _lcge_start_value(m[sym][k...])
                v  isa Real && (total += float(v))
                sv isa Real && (start += float(sv))
            end
            return (total, start)
        catch
            return missing
        end
    end
    function _push(name, desc, value, start)
        push!(rows, Any[name, desc, value, start])
    end

    # ── Output and trade ──────────────────────────────────────────────────────
    let pair = _sum(:XP, [(ii,) for ii in i])
        pair isa Tuple && _push("XP_total",   "Total gross output Σ XP[i]",                       pair[1], pair[2])
    end
    let pair = _sum(:XA, [(ii,) for ii in i])
        pair isa Tuple && _push("XA_total",   "Total Armington demand Σ XA[i]",                   pair[1], pair[2])
    end
    let pair = _sum(:XDs, [(ii,) for ii in i])
        pair isa Tuple && _push("XD_total",   "Total domestic supply Σ XDs[i]  (legacy XD)",      pair[1], pair[2])
    end
    let pair = _sum(:XMT, [(ii,) for ii in i])
        pair isa Tuple && _push("XM_total",   "Total imports Σ XMT[i]  (legacy XM)",              pair[1], pair[2])
    end
    let pair = _sum(:ES, [(ii,) for ii in i])
        pair isa Tuple && _push("E_total",    "Total exports Σ ES[i]  (legacy E)",                pair[1], pair[2])
    end
    let pair = _sum(:ND, [(ii,) for ii in i])
        pair isa Tuple && _push("ND_total",   "Total intermediate demand Σ ND[i]",                pair[1], pair[2])
    end
    let pair = _sum(:VA, [(ii,vv) for ii in i, vv in S[:v]])
        pair isa Tuple && _push("VA_total",   "Total value added Σ VA[i,v]",                      pair[1], pair[2])
    end

    # ── Macro aggregates per region ───────────────────────────────────────────
    for rr in r
        _push("GDP_$(rr)",  "Nominal GDP at producer prices, region $(rr)",                 _v(:GDP,(rr,)),  _vs(:GDP,(rr,)))
        _push("RGDP_$(rr)", "Real GDP (Σ XP), region $(rr)",                                _v(:RGDP,(rr,)), _vs(:RGDP,(rr,)))
        _push("PGDP_$(rr)", "GDP deflator GDP/RGDP, region $(rr)",                          _v(:PGDP,(rr,)), _vs(:PGDP,(rr,)))
        _push("CPI_$(rr)",  "Consumer price index (avg PC), region $(rr)",                  _v(:CPI,(rr,)),  _vs(:CPI,(rr,)))
    end

    # ── Household income / savings ────────────────────────────────────────────
    for hh in h
        _push("YH_$(hh)",   "Household income (Y_5), $(hh)",                                _v(:YH,(hh,)),   _vs(:YH,(hh,)))
        _push("YD_$(hh)",   "Disposable income after tax (Y_7), $(hh)",                     _v(:YD,(hh,)),   _vs(:YD,(hh,)))
        _push("YC_$(hh)",   "Income for consumption (Y_8), $(hh)",                          _v(:YC,(hh,)),   _vs(:YC,(hh,)))
        _push("SAV_$(hh)",  "Household savings (D_3), $(hh)",                               _v(:SAV,(hh,)),  _vs(:SAV,(hh,)))
        _push("DeprY_$(hh)","Depreciation (Y_6), $(hh)",                                    _v(:DeprY,(hh,)),_vs(:DeprY,(hh,)))
    end

    # ── Government, investment, savings totals ────────────────────────────────
    _push("GOVREV",   "Government revenue YG  (legacy GOVREV)",                             _v(:YG),         _vs(:YG))
    _push("GEXP",     "Government expenditure PFD[Gov]*FD[Gov]  (legacy GEXP)",
          (haskey(m,:PFD) && haskey(m,:FD)) ? (_v(:PFD,("Gov",))*_v(:FD,("Gov",))) : missing,
          (haskey(m,:PFD) && haskey(m,:FD)) ? (_vs(:PFD,("Gov",))*_vs(:FD,("Gov",))) : missing)
    _push("Sg",       "Government saving (C_4)",                                            _v(:Sg),         _vs(:Sg))
    _push("INVEST",   "Total investment PFD[Inv]*FD[Inv]  (legacy INVEST)",
          (haskey(m,:PFD) && haskey(m,:FD)) ? (_v(:PFD,("Inv",))*_v(:FD,("Inv",))) : missing,
          (haskey(m,:PFD) && haskey(m,:FD)) ? (_vs(:PFD,("Inv",))*_vs(:FD,("Inv",))) : missing)
    _push("FDInv",    "Investment quantity FDInv (F_31)",                                   _v(:FDInv),      _vs(:FDInv))
    _push("InvSh",    "Investment / GDP ratio (C_10)",                                      _v(:InvSh),      _vs(:InvSh))

    # ── Tax revenue components (TAXREV ≈ YG = sum of all tax flows) ──────────
    _push("TAXREV",   "Total tax revenue ≈ YG (legacy TAXREV)",                             _v(:YG),         _vs(:YG))
    _push("TarY",     "Tariff revenue (C_1)",                                               _v(:TarY),       _vs(:TarY))

    # ── Aggregate savings (SAVE = household + government + foreign) ───────────
    let pair = _sum(:Sf, [(rr,) for rr in r])
        sav_h = sum((v for v in (_v(:SAV,(hh,)) for hh in h) if v isa Real); init=0.0)
        sav_h_start = sum((v for v in (_vs(:SAV,(hh,)) for hh in h) if v isa Real); init=0.0)
        sg = _v(:Sg);  sg_s = _vs(:Sg)
        if pair isa Tuple
            sf, sf_s = pair
            tot = sav_h + (sg isa Real ? sg : 0.0) + sf
            ts  = sav_h_start + (sg_s isa Real ? sg_s : 0.0) + sf_s
            _push("SAVE", "Aggregate savings = HH + Gov + Foreign  (legacy SAVE)", tot, ts)
        end
    end

    # ── Balance of payments ───────────────────────────────────────────────────
    # World-price import and export values, their difference (the trade balance,
    # positive = deficit) and the real exchange rate.  Under trade_closure = :bop
    # C_BOP forces `MW_total - EW_total = Σ Sf`; under :balanced the two are equal
    # good by good by construction.  `SI_gap` is the savings-investment residual
    # of C-9, which is dropped under :bop and must come out at ~0 by Walras' law.
    let rp_ = S[:rp]
        mw = _sum2(m, :WPM, :WTFd, [(rrp,rr,ii) for ii in i, rr in r, rrp in rp_])
        ew = _sum2(m, :WPE, :WTFs, [(rr,rrp,ii) for ii in i, rr in r, rrp in rp_])
        if mw isa Tuple && ew isa Tuple
            _push("MW_total", "Imports at world prices Σ WPM·WTFd",  mw[1], mw[2])
            _push("EW_total", "Exports at world prices Σ WPE·WTFs",  ew[1], ew[2])
            _push("TradeBalance", "World-price trade deficit Σ WPM·WTFd − Σ WPE·WTFs",
                  mw[1]-ew[1], mw[2]-ew[2])
        end
    end
    _push("ER", "Real exchange rate (C_BOP / C_ER; absent under :balanced)",                _v(:ER),         _vs(:ER))
    let sav_h = sum((v for v in (_v(:SAV,(hh,)) for hh in h) if v isa Real); init=0.0),
        dep_h = sum((v for v in (_v(:DeprY,(hh,)) for hh in h) if v isa Real); init=0.0),
        sf    = sum((v for v in (_v(:Sf,(rr,)) for rr in r) if v isa Real); init=0.0),
        sg    = _v(:Sg),
        inv   = (haskey(m,:PFD) && haskey(m,:FD)) ? _v(:PFD,("Inv",))*_v(:FD,("Inv",)) : missing
        if inv isa Real && sg isa Real
            _push("SI_gap", "Savings-investment residual of C-9 (0 by Walras' law under :bop)",
                  inv - (sav_h + dep_h + sg + sf), 0.0)
        end
    end

    # ── Factor incomes ────────────────────────────────────────────────────────
    _push("TY", "Land income (Y_1)",                                                        _v(:TY),         _vs(:TY))
    _push("FY", "Natural-resource income (Y_2)",                                            _v(:FY),         _vs(:FY))
    _push("KY", "Capital income (Y_4)",                                                     _v(:KY),         _vs(:KY))
    for ll in S[:l]
        _push("LY_$(ll)", "Labor income by skill, $(ll)",                                   _v(:LY,(ll,)),   _vs(:LY,(ll,)))
    end

    # ── Land / NR markets ─────────────────────────────────────────────────────
    _push("TLnd",  "Aggregate land supply",                                                 _v(:TLnd),       _vs(:TLnd))
    _push("PTLnd", "Aggregate land price",                                                  _v(:PTLnd),      _vs(:PTLnd))
    _push("KS",    "Aggregate capital supply",                                              _v(:KS),         _vs(:KS))
    _push("TR",    "Aggregate capital return",                                              _v(:TR),         _vs(:TR))
    _push("PNUM",  "Numeraire",                                                             _v(:PNUM),       _vs(:PNUM))

    return rows
end

"""Solve the model and immediately export result files."""
function solve_and_save!(m, data::LinkageData; outdir::AbstractString="results")
    solve_model!(m)
    export_results!(m, data; outdir=outdir)
    return m
end
