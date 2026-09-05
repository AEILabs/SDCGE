# Recursive-dynamic LINKAGE solver.
#
# Each period is solved as a STATIC equilibrium using the cleaned static model.
# Between periods we update the exogenous state (capital stock, labor supply,
# productivity, etc.) using the solved period-t values, then re-build and re-
# solve for period t+1.
#
# Typical usage:
#
#     data, history = run_recursive_dynamic!(periods=10)
#
# Or specify growth/depreciation rates:
#
#     data, history = run_recursive_dynamic!(
#         periods=10,
#         delta=0.05,            # capital depreciation rate
#         g_labor=0.02,          # labor force growth
#         g_tfp=0.015,           # total factor productivity growth
#         g_land=0.0,            # land supply growth
#         outdir="results/dynamic")
#
# Outputs:
# - results/dynamic/period_<t>/...   (full static-model outputs per period)
# - results/dynamic/time_series.csv  (key macro indicators stacked across t)
#
# ─────────────────────────────────────────────────────────────────────────────
# CAPITAL STOCK vs FLOW UNITS  (fixed 2026-09-03)
# ─────────────────────────────────────────────────────────────────────────────
# `PAR[:KSupply][(i,v)]` and `PAR[:K0][i]` are capital RENTAL PAYMENTS: the
# calibration reads them straight off the SAM's CAP row (CAP_i/|v| per vintage),
# so they are a price × quantity flow per period, not a stock.  `FDInv` (the
# solved investment quantity, = FD["Inv"]) is a COMMODITY flow.  Adding one to
# the other, as the previous `K[t+1] = (1−δ)K[t] + FDInv` update did, mixes
# units: at the benchmark CAP = 17,701 while FDInv = 11,519, so `K0` jumped by a
# factor of ~2 per period and PATH failed from period 3 on.
#
# The state variable is now an explicit capital STOCK in commodity units:
#
#   Kstock0      = FDInv0 / delta                     (aggregate; stationary at
#                                                      the benchmark, I0 = δ·K0)
#   share_i      = CAP_i / Σ_j CAP_j                  (benchmark rental shares)
#   Kstock0_i    = Kstock0 · share_i
#   kappa_i      = Σ_v KSupply0[i,v] / Kstock0_i      (rental per unit of stock,
#                                                      constant = the benchmark
#                                                      gross rate of return)
#
# Period update, with I_t the solved `FDInv` (commodity units, the same units
# Kstock0 was defined in):
#
#   Kstock_{t+1,i} = (1 − delta)·Kstock_{t,i} + I_t·share_i
#   KSupply_{t+1}[i,v] = kappa_i · Kstock_{t+1,i} · vshare_{i,v}
#   K0_{t+1}[i]        = KSupply_{t+1}[i,"Old"]        (F-30 / F_K0 anchor)
#
# `vshare_{i,v}` is the benchmark vintage split of sector i's capital (0.5/0.5
# with the current calibration, which has no vintage information in the SAM).
# `vintage_rule = :flow` instead uses the textbook vintage reading
# (Old = kappa·(1−δ)·Kstock_t, New = kappa·I_t·share_i, same sector total).
# `:benchmark_shares` is the default because the two vintages are technologically
# identical in this port (Calibration.jl gives them the same alpha_*), so the
# textbook split adds no economics, while F-24 pins RR = 1 and therefore makes
# F-30 fix XPv[i,"Old"] rigidly: at the (1−δ)/δ = 95/5 split XPv[i,"New"] has
# only 5 % headroom before hitting its zero bound.  See the report for measured
# behaviour of both rules.
#
# With g_labor = g_tfp = 0 the stock is stationary (Kstock_{t+1} = Kstock_t), so
# every period reproduces the benchmark exactly.

"""
`true` when the last solve of `m` converged.  A period that hit its iteration
limit returns whatever point PATH stopped at, which is not an equilibrium; the
state updates below refuse to propagate such a point, so one failed period makes
the path plateau instead of corrupting every period after it.
"""
function _period_converged(m)
    try
        return termination_status(m) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)
    catch
        return false
    end
end

"""
Capital-stock state for the recursive dynamics, created on first use and cached
in `data.metadata[:capital_state]`.

Returns a named tuple `(stock, kappa, share, vshare, stock0, delta)` where
`stock[i]` is sector i's capital stock in commodity units, `kappa[i]` converts
stock to the rental payment used by `PAR[:KSupply]`/`PAR[:K0]`, `share[i]` is
the benchmark capital share used to allocate investment, and `vshare[(i,v)]` is
the benchmark vintage split.  See the units note at the top of this file.
"""
function _capital_state!(data::LinkageData; delta::Float64=0.05)
    st = get(data.metadata, :capital_state, nothing)
    st === nothing || return st

    PAR = data.metadata[:PAR]
    S = data.sets; i = S[:i]; v = S[:v]
    nv = max(length(v), 1)

    # Benchmark rental payments (SAM CAP row) per sector and per vintage.
    rent = Dict{String,Float64}(ii => sum(get(PAR[:KSupply], (ii,vv), 0.0) for vv in v)
                                for ii in i)
    rent_tot = sum(values(rent))

    # Benchmark investment in commodity units: FDInv0 = FD["Inv"] at the benchmark.
    B = get(PAR, :bench, nothing)
    inv0 = B isa Dict ? float(get(B, :INV, 0.0)) : 0.0
    inv0 > 0.0 || (inv0 = float(get(PAR, :INVEST0, 0.0)))
    (inv0 > 0.0 && delta > 0.0) ||
        error("Cannot initialise the capital stock: need a positive benchmark investment " *
              "(got $inv0) and a positive depreciation rate (got $delta).")

    stock_tot = inv0 / delta
    share  = Dict{String,Float64}(ii => rent_tot > 1.0e-9 ? rent[ii]/rent_tot : 1.0/length(i)
                                  for ii in i)
    stock  = Dict{String,Float64}(ii => stock_tot * share[ii] for ii in i)
    kappa  = Dict{String,Float64}(ii => stock[ii] > 1.0e-12 ? rent[ii]/stock[ii] : 0.0
                                  for ii in i)
    vshare = Dict{Any,Float64}((ii,vv) => rent[ii] > 1.0e-12 ?
                get(PAR[:KSupply], (ii,vv), 0.0)/rent[ii] : 1.0/nv
                for ii in i for vv in v)

    st = (stock=stock, kappa=kappa, share=share, vshare=vshare,
          stock0=stock_tot, delta=delta)
    data.metadata[:capital_state] = st
    return st
end

"""
Advance the capital stock one period from the solved model `m` and write the
resulting rental-unit tables (`PAR[:KSupply]`, `PAR[:K0]`, `PAR[:gamma_K]`,
`PAR[:KY0]`).  Shared by `update_period_data!` and
`update_period_data_scenario!`; see the units note at the top of this file.
"""
function _update_capital_stock!(data::LinkageData, m;
        delta::Float64=0.05, vintage_rule::Symbol=:benchmark_shares)

    PAR = data.metadata[:PAR]
    S = data.sets; i = S[:i]; v = S[:v]
    st = _capital_state!(data; delta=delta)

    inv_t = _period_converged(m) ? (try JuMP.value(m[:FDInv]) catch; missing end) : missing
    inv_t === missing && @warn "Period did not converge; capital stock carried forward unchanged."
    if inv_t isa Real && isfinite(inv_t)
        use_flow = vintage_rule === :flow && ("Old" in v) && ("New" in v)
        for ii in i
            old_stock = st.stock[ii]
            new_inv   = inv_t * st.share[ii]
            new_stock = (1 - delta) * old_stock + new_inv
            k = st.kappa[ii]
            if use_flow
                # Textbook vintage reading: Old = surviving capital, New = this
                # period's investment.  Same sector total as the branch below.
                PAR[:KSupply][(ii,"Old")] = max(k * (1 - delta) * old_stock, 1.0e-9)
                PAR[:KSupply][(ii,"New")] = max(k * new_inv,                 1.0e-9)
            else
                for vv in v
                    PAR[:KSupply][(ii,vv)] = max(k * new_stock * st.vshare[(ii,vv)], 1.0e-9)
                end
            end
            # K0 is the F-30 / F_K0 anchor and is in the same rental units.
            PAR[:K0][ii] = PAR[:KSupply][(ii, ("Old" in v) ? "Old" : first(v))]
            st.stock[ii] = new_stock
        end
    end

    # Aggregate rental anchor and the F-21 capital CET shares (must sum to 1).
    PAR[:KY0] = sum(get(PAR[:KSupply], (ii,vv), 0.0) for ii in i for vv in v)
    oldv = ("Old" in v) ? "Old" : first(v)
    tot_k_old = sum(get(PAR[:KSupply], (ii,oldv), 0.0) for ii in i)
    if tot_k_old > 1.0e-9
        for ii in i
            PAR[:gamma_K][ii] = PAR[:KSupply][(ii,oldv)] / tot_k_old
        end
    end
    return data
end

"""
Record realised sectoral labour demand in `PAR[:LV0]`.

`PAR[:LV0]` is a start-value table only — no equation reads it — but
`Initialization.jl` builds the `LV`, `LY` and (through F-10) `UE` start values
from it.  Left at its base-year value it drifts further from the solution every
period, and the `UE` start in particular then sits far from where F-10 puts it,
which is the one place the dynamic path reliably trips PATH up.
"""
function _track_labour_demand!(data::LinkageData, m)
    _period_converged(m) || return data
    PAR = data.metadata[:PAR]
    S = data.sets; i = S[:i]; l = S[:l]
    try
        LVv = m[:LV]
        for ll in l, ii in i
            val = JuMP.value(LVv[ll,ii])
            (val isa Real && isfinite(val)) || continue
            PAR[:LV0][(ll,ii)] = max(float(val), 0.0)
        end
    catch err
        @warn "Could not record realised labour demand" exception=(err, catch_backtrace())
    end
    return data
end

"""
Update exogenous data fields (PAR + supplies) using period-t solution.
This is the recursive-dynamic update step:

* Capital:    Kstock[t+1] = (1−δ)·Kstock[t] + I[t] in commodity units, converted
              back to the rental units of `KSupply`/`K0` by `kappa` (see the
              units note at the top of this file)
* Labor:      LS[t+1] = (1+g_labor)·LS[t]
* TFP:        AT[t+1] = (1+g_tfp)·AT[t]
* Land/NR:    grow at g_land / g_nres (default 0)

`vintage_rule` selects how the sector's capital is split across vintages:
`:benchmark_shares` (default) keeps the calibrated split, `:flow` uses
Old = surviving capital, New = this period's investment.

Returns the modified `data`.
"""
function update_period_data!(data::LinkageData, m;
        delta::Float64=0.05,
        g_labor::Float64=0.02,
        g_tfp::Float64=0.0,
        g_land::Float64=0.0,
        g_nres::Float64=0.0,
        vintage_rule::Symbol=:benchmark_shares)

    PAR = data.metadata[:PAR]
    S = data.sets
    i = S[:i]; l = S[:l]; v = S[:v]; gz = S[:gz]; ag = S[:ag]; h = S[:h]

    # ── 1. CAPITAL ACCUMULATION (explicit stock in commodity units) ──────────
    _update_capital_stock!(data, m; delta=delta, vintage_rule=vintage_rule)

    # ── 2. LABOR SUPPLY GROWTH ───────────────────────────────────────────────
    for ll in l
        PAR[:LSupply][ll] *= (1 + g_labor)
        for gg in gz
            PAR[:LS0][(ll, gg)] *= (1 + g_labor)
        end
        PAR[:LY0][ll] *= (1 + g_labor)
    end
    _track_labour_demand!(data, m)

    # ── 3. TOTAL FACTOR PRODUCTIVITY ─────────────────────────────────────────
    if abs(g_tfp) > 1.0e-12
        for ii in i
            PAR[:AT][ii] *= (1 + g_tfp)
        end
    end

    # ── 4. LAND AND NATURAL-RESOURCE SUPPLIES ────────────────────────────────
    if abs(g_land) > 1.0e-12
        for ii in i
            PAR[:TSupply][ii] *= (1 + g_land)
        end
        PAR[:TY0] *= (1 + g_land)
        if haskey(PAR[:chi_T], :land)
            PAR[:chi_T][:land] *= (1 + g_land)
        end
        # Re-normalize gamma_T to preserve sum to 1
        tot_t = sum(get(PAR[:TSupply], ii, 0.0) for ii in ag)
        if tot_t > 1.0e-9
            for ii in i
                PAR[:gamma_T][ii] = (ii in ag) ? PAR[:TSupply][ii] / tot_t : 0.0
            end
        end
    end
    if abs(g_nres) > 1.0e-12
        for ii in i
            PAR[:FSupply][ii] *= (1 + g_nres)
            PAR[:chi_F][ii] *= (1 + g_nres)
        end
        PAR[:FY0] *= (1 + g_nres)
    end

    # ── 4b. EXOGENOUS REAL INVESTMENT (trade_closure = :bop) ─────────────────
    # C-INV pins FD[Inv] to PAR[:FDInv0] inside a period (see Closure.jl for why
    # a share-of-GDP rule cannot be used).  Re-base it on this period's real GDP
    # so it still tracks a growing economy; with zero growth RGDP is unchanged
    # and FDInv0 stays exactly at its benchmark value.  `PAR[:Sfbar]` is
    # deliberately NOT grown: foreign saving is an exogenous level throughout.
    if Symbol(get(PAR, :trade_closure, :bop)) !== :balanced
        rgdp = try JuMP.value(m[:RGDP][first(S[:r])]) catch; nothing end
        if rgdp isa Real && isfinite(rgdp) && rgdp > 0
            PAR[:FDInv0] = get(PAR, :chi_inv, 0.2) * rgdp
        end
    end

    # ── 5. UPDATE AGGREGATE ANCHORS USED BY INITIALIZATION ────────────────────
    # KY0 and gamma_K are set by _update_capital_stock! above.
    # YH benchmark: factor incomes (assuming PF/W/R = 1 at next benchmark).
    PAR[:YH0] = Dict(hh => max(
        PAR[:TY0] + PAR[:FY0] + sum(values(PAR[:LY0])) + PAR[:KY0], 1.0e-9)
        for hh in h)

    return data
end

"""
Extract a compact named-tuple of key macro indicators for one period.
"""
function _period_summary(m, data::LinkageData, period::Int)
    S = data.sets
    function _v(sym, idx=())
        haskey(m, sym) || return missing
        try
            return isempty(idx) ? JuMP.value(m[sym]) : JuMP.value(m[sym][idx...])
        catch
            return missing
        end
    end
    function _sum(sym, keys)
        haskey(m, sym) || return missing
        total = 0.0
        try
            for k in keys
                v = JuMP.value(m[sym][k...])
                v isa Real && (total += float(v))
            end
            return total
        catch
            return missing
        end
    end

    return (
        period         = period,
        GDP_R1         = _v(:GDP, ("R1",)),
        RGDP_R1        = _v(:RGDP, ("R1",)),
        PGDP_R1        = _v(:PGDP, ("R1",)),
        CPI_R1         = _v(:CPI, ("R1",)),
        XP_total       = _sum(:XP, [(ii,) for ii in S[:i]]),
        YH_HH          = _v(:YH, ("HH",)),
        YC_HH          = _v(:YC, ("HH",)),
        SAV_HH         = _v(:SAV, ("HH",)),
        GOVREV         = _v(:YG),
        GEXP           = (_v(:PFD, ("Gov",)) isa Real && _v(:FD, ("Gov",)) isa Real) ?
                          _v(:PFD, ("Gov",)) * _v(:FD, ("Gov",)) : missing,
        INVEST         = (_v(:PFD, ("Inv",)) isa Real && _v(:FD, ("Inv",)) isa Real) ?
                          _v(:PFD, ("Inv",)) * _v(:FD, ("Inv",)) : missing,
        FDInv          = _v(:FDInv),
        Sg             = _v(:Sg),
        InvSh          = _v(:InvSh),
        KS             = _v(:KS),
        KActual        = _v(:KActual),
        TR             = _v(:TR),
        TLnd           = _v(:TLnd),
        PTLnd          = _v(:PTLnd),
        PNUM           = _v(:PNUM),
        TY             = _v(:TY),
        FY             = _v(:FY),
        KY             = _v(:KY),
        LY_UnSkLab     = _v(:LY, ("UnSkLab",)),
        LY_SkLab       = _v(:LY, ("SkLab",)),
    )
end

"""
Write the per-period summaries as a single time-series CSV (one row per period).
"""
function _write_time_series(history::Vector, outdir::AbstractString)
    isempty(history) && return nothing
    cols = collect(keys(history[1]))
    path = joinpath(outdir, "time_series.csv")
    mkpath(outdir)
    open(path, "w") do io
        println(io, join(string.(cols), ","))
        for row in history
            println(io, join((string(getfield(row, c)) for c in cols), ","))
        end
    end
    return path
end

"""
Collect ALL variable values for one solved period into a flat dict
`(variable_symbol, index_tuple) => value`. Used to build the consolidated
Excel workbook in `_write_dynamic_xlsx`.

For JuMP `DenseAxisArray`s (the common case), `eachindex` returns
`CartesianIndex` values that don't preserve the original axis labels (e.g.
"R1", "P001"). We translate to axis labels via `obj.axes` so the workbook
shows the readable sector / region / vintage names.
"""
function _collect_period_values(m)
    snap = Dict{Tuple{Symbol,Any},Float64}()
    for (name, obj) in object_dictionary(m)
        if obj isa JuMP.VariableRef
            try
                v = JuMP.value(obj)
                v isa Real && (snap[(name, ())] = float(v))
            catch
            end
            continue
        end
        # Try to get axis labels for DenseAxisArray; fall back to raw index.
        labelfn = if hasproperty(obj, :axes)
            ax = obj.axes
            (idx) -> begin
                if idx isa CartesianIndex
                    Tuple(ax[d][idx[d]] for d in 1:length(idx))
                elseif idx isa Tuple
                    Tuple(ax[d][idx[d]] for d in 1:length(idx))
                else
                    (ax[1][idx],)
                end
            end
        else
            (idx) -> idx
        end
        try
            for idx in eachindex(obj)
                var = obj[idx]
                var isa JuMP.VariableRef || continue
                try
                    v = JuMP.value(var)
                    v isa Real || continue
                    lbl = try labelfn(idx) catch; idx end
                    snap[(name, lbl)] = float(v)
                catch
                end
            end
        catch
        end
    end
    return snap
end

"""
Write the consolidated dynamic workbook with one sheet per variable family.

Sheet layout (per variable):
- Row 1: ["index", "period_1", "period_2", ..., "period_T"]
- Subsequent rows: [index_label, value_t1, value_t2, ..., value_tT]

Plus two summary sheets:
- "macro_summary": named indicators from `history`
- "scalars": all scalar variables (no index) tracked across periods
"""
function _write_dynamic_xlsx(snapshots::Vector, history::Vector, outdir::AbstractString)
    isempty(snapshots) && return nothing
    mkpath(outdir)
    path = joinpath(outdir, "dynamic_results.xlsx")
    T = length(snapshots)

    # Discover all variable names across all snapshots.
    var_indices = Dict{Symbol, Set{Any}}()
    for snap in snapshots
        for ((name, idx), _) in snap
            push!(get!(var_indices, name, Set{Any}()), idx)
        end
    end

    period_cols = ["period_$(t)" for t in 1:T]

    XLSX.openxlsx(path, mode="w") do xf
        # ── Sheet: macro_summary ──────────────────────────────────────────────
        if !isempty(history)
            cols = collect(keys(history[1]))
            sheet_name = "macro_summary"
            sheet = if XLSX.sheetnames(xf)[1] == "Sheet1"
                XLSX.rename!(xf[1], sheet_name); xf[1]
            else
                XLSX.addsheet!(xf, sheet_name)
            end
            # Header row: ["indicator", "period_1", ..., "period_T"]
            sheet["A1"] = "indicator"
            for t in 1:T
                sheet[XLSX.CellRef(1, t+1)] = "period_$(t)"
            end
            for (i, ind) in enumerate(cols)
                sheet[XLSX.CellRef(i+1, 1)] = string(ind)
                for t in 1:T
                    v = getfield(history[t], ind)
                    if v isa Real && isfinite(v)
                        sheet[XLSX.CellRef(i+1, t+1)] = float(v)
                    end
                end
            end
        end

        # ── One sheet per variable ────────────────────────────────────────────
        # Sort variable names so sheet ordering is deterministic.
        for vname in sort(collect(keys(var_indices)); by=string)
            sheet_name = string(vname)
            # Excel sheet names max 31 chars; truncate if needed.
            length(sheet_name) > 31 && (sheet_name = sheet_name[1:31])
            sheet = XLSX.addsheet!(xf, sheet_name)

            indices = sort(collect(var_indices[vname]); by=_index_sort_key)

            # Header row.
            sheet["A1"] = "index"
            for t in 1:T
                sheet[XLSX.CellRef(1, t+1)] = "period_$(t)"
            end

            # Data rows.
            for (i, idx) in enumerate(indices)
                sheet[XLSX.CellRef(i+1, 1)] = _index_label(idx)
                for t in 1:T
                    v = get(snapshots[t], (vname, idx), nothing)
                    if v isa Real && isfinite(v)
                        sheet[XLSX.CellRef(i+1, t+1)] = float(v)
                    end
                end
            end
        end
    end
    return path
end

_index_label(idx) = idx isa Tuple ? join(string.(idx), "|") : (idx === () ? "" : string(idx))
_index_sort_key(idx) = idx isa Tuple ? join(string.(idx), "|") : string(idx)

"""
    run_recursive_dynamic!(; periods=10, delta=0.05, g_labor=0.02, g_tfp=0.015,
                           g_land=0.0, g_nres=0.0,
                           outdir="results/dynamic", show_solver_output=false)

Solve the LINKAGE model recursively over `periods` periods. Each period is a
static-equilibrium solve with capital, labor and productivity updated between
periods from the previous solution.

Returns `(data, history)` where `history` is a vector of named tuples — one per
period — containing key macro indicators. The full per-period CSVs are written
to `outdir/period_<t>/`, and a consolidated `time_series.csv` is written to
`outdir/`.
"""
function run_recursive_dynamic!(;
        periods::Int=10,
        delta::Float64=0.05,
        g_labor::Float64=0.02,
        g_tfp::Float64=0.015,
        g_land::Float64=0.0,
        g_nres::Float64=0.0,
        outdir::AbstractString="results/dynamic",
        show_solver_output::Bool=false,
        vintage_rule::Symbol=:benchmark_shares)

    mkpath(outdir)
    data = init_data()
    prepare_data!(data)

    history    = Vector{NamedTuple}()
    snapshots  = Vector{Dict{Tuple{Symbol,Any},Float64}}()

    for t in 1:periods
        println("\n========================================")
        println("  Recursive-Dynamic Period t = $t / $periods")
        println("========================================")

        m = model(data; show_solver_output=show_solver_output)
        solve_model!(m; show_diagnostics=false)

        # Capture per-period values for the consolidated workbook.
        push!(snapshots, _collect_period_values(m))
        push!(history,   _period_summary(m, data, t))

        # Update exogenous state for next period (skip after the last period).
        if t < periods
            update_period_data!(data, m;
                delta=delta, g_labor=g_labor, g_tfp=g_tfp,
                g_land=g_land, g_nres=g_nres, vintage_rule=vintage_rule)
        end
    end

    # ── Consolidated outputs ──────────────────────────────────────────────────
    # One Excel workbook: one sheet per variable family, columns = periods.
    xlsx_path = _write_dynamic_xlsx(snapshots, history, outdir)
    println("\nDynamic results written to ", xlsx_path)

    # Time-series CSV summary (handy for quick inspection / downstream tools).
    _write_time_series(history, outdir)
    println("Time-series summary written to ", joinpath(outdir, "time_series.csv"))

    return data, history, snapshots
end
