# Policy-experiment driver.
#
# An Excel workbook acts as the experiment design.  Each sheet stores one
# exogenous block over the simulation horizon.  The Julia code reads the
# workbook, builds a `Scenario` per row in the `scenarios` sheet and runs the
# recursive-dynamic solver for each active scenario.
#
# Workflow:
#
#     write_policy_template("policy_experiments.xlsx";
#                           periods=10, n_scenarios=10)
#     # edit the file by hand to design experiments, then:
#     run_policy_experiments!("policy_experiments.xlsx";
#                             outdir="results/scenarios")
#
# Each scenario writes its results to `outdir/<scenario_name>/` with the
# same consolidated XLSX + plot files as the standalone dynamic run.
#
# The workbook is shaped so every exogenous variable that varies across
# scenarios and periods is fully visible — no aggregation that hides
# heterogeneity (productivity AT and labor growth are by activity / skill).

# ─── Scenario container ──────────────────────────────────────────────────────

struct Scenario
    sim_id::Int
    name::String
    description::String
    periods::Int
    delta::Float64                      # capital depreciation rate (scalar)
    AT::Dict{Tuple{Any,Int},Float64}     # AT[(activity, period)] – productivity level
    g_labor::Dict{Tuple{Any,Int},Float64}# g_labor[(skill, period)] – growth rate
    g_land::Vector{Float64}              # length = periods
    g_nres::Vector{Float64}              # length = periods
end

# ─── Sheet names used by the template ────────────────────────────────────────
const _POL_SHEETS = (
    scenarios = "scenarios",
    AT        = "AT_by_activity",
    g_labor   = "g_labor_by_skill",
    g_land    = "g_land",
    g_nres    = "g_nres",
)

# ─── Template writer ─────────────────────────────────────────────────────────

"""
    write_policy_template(path::AbstractString; periods=10, n_scenarios=10)

Create a starter Excel workbook for policy experiments.  All values default to
the benchmark (`AT = 1`, `g_labor = 0.02`, `g_land = g_nres = 0`) so that
running the file unchanged reproduces 10 identical baseline runs.  Open the
file in Excel and overwrite cells to design experiments.

Sheets:
- `scenarios`        – one row per simulation (sim_id, name, description,
                        periods, delta, active)
- `AT_by_activity`   – productivity level by (sim_id, activity, period)
- `g_labor_by_skill` – labor-supply growth by (sim_id, skill, period)
- `g_land`           – land-supply growth by (sim_id, period)
- `g_nres`           – natural-resource growth by (sim_id, period)
"""
function write_policy_template(path::AbstractString;
        periods::Int=10, n_scenarios::Int=10)

    # We need the sector / skill names to size the sheets; the SAM defaults are
    # the canonical 100 products and 2 skill types.
    data = init_data(); default_sets!(data)
    activities = data.sets[:i]   # 100 sectors
    skills     = data.sets[:l]   # ["UnSkLab", "SkLab"]

    period_cols = ["period_$(t)" for t in 1:periods]
    mkpath(dirname(path))

    XLSX.openxlsx(path, mode="w") do xf
        # ── scenarios ────────────────────────────────────────────────────────
        scen = if XLSX.sheetnames(xf)[1] == "Sheet1"
            XLSX.rename!(xf[1], _POL_SHEETS.scenarios); xf[1]
        else
            XLSX.addsheet!(xf, _POL_SHEETS.scenarios)
        end
        header = ["sim_id","name","description","periods","delta","active"]
        for (j,h) in enumerate(header); scen[XLSX.CellRef(1,j)] = h; end
        defaults = [
            (1, "baseline",       "benchmark; no growth shocks"),
            (2, "high_TFP",       "1.5% TFP growth in all activities"),
            (3, "low_TFP",        "0.5% TFP growth in all activities"),
            (4, "fast_labor",     "labor force grows 3%/period"),
            (5, "slow_labor",     "labor force grows 1%/period"),
            (6, "land_growth",    "land supply grows 1%/period"),
            (7, "TFP_in_crops",   "2% TFP growth, crops only (P001-P010)"),
            (8, "TFP_in_industry","2% TFP growth, industrial sectors P021-P100"),
            (9, "high_depreciation","capital depreciation 10% (vs 5% baseline)"),
            (10,"all_growth",     "TFP 2% + labor 3% + land 0.5%"),
        ]
        for r in 1:n_scenarios
            id, nm, desc = r <= length(defaults) ? defaults[r] :
                           (r, "scenario_$(r)", "user-defined scenario $(r)")
            scen[XLSX.CellRef(r+1, 1)] = id
            scen[XLSX.CellRef(r+1, 2)] = nm
            scen[XLSX.CellRef(r+1, 3)] = desc
            scen[XLSX.CellRef(r+1, 4)] = periods
            scen[XLSX.CellRef(r+1, 5)] = (nm == "high_depreciation") ? 0.10 : 0.05
            scen[XLSX.CellRef(r+1, 6)] = "yes"
        end

        # ── AT_by_activity ───────────────────────────────────────────────────
        at = XLSX.addsheet!(xf, _POL_SHEETS.AT)
        at["A1"] = "sim_id"; at["B1"] = "activity"
        for t in 1:periods; at[XLSX.CellRef(1, t+2)] = period_cols[t]; end
        row = 2
        for sid in 1:n_scenarios
            for ii in activities
                at[XLSX.CellRef(row, 1)] = sid
                at[XLSX.CellRef(row, 2)] = ii
                vals = _default_AT_path(sid, ii, periods)
                for t in 1:periods
                    at[XLSX.CellRef(row, t+2)] = vals[t]
                end
                row += 1
            end
        end

        # ── g_labor_by_skill ─────────────────────────────────────────────────
        gl = XLSX.addsheet!(xf, _POL_SHEETS.g_labor)
        gl["A1"] = "sim_id"; gl["B1"] = "skill"
        for t in 1:periods; gl[XLSX.CellRef(1, t+2)] = period_cols[t]; end
        row = 2
        for sid in 1:n_scenarios
            for ll in skills
                gl[XLSX.CellRef(row, 1)] = sid
                gl[XLSX.CellRef(row, 2)] = ll
                vals = _default_g_labor_path(sid, ll, periods)
                for t in 1:periods
                    gl[XLSX.CellRef(row, t+2)] = vals[t]
                end
                row += 1
            end
        end

        # ── g_land ───────────────────────────────────────────────────────────
        gland = XLSX.addsheet!(xf, _POL_SHEETS.g_land)
        gland["A1"] = "sim_id"
        for t in 1:periods; gland[XLSX.CellRef(1, t+1)] = period_cols[t]; end
        for sid in 1:n_scenarios
            gland[XLSX.CellRef(sid+1, 1)] = sid
            vals = _default_g_land_path(sid, periods)
            for t in 1:periods
                gland[XLSX.CellRef(sid+1, t+1)] = vals[t]
            end
        end

        # ── g_nres ───────────────────────────────────────────────────────────
        gnres = XLSX.addsheet!(xf, _POL_SHEETS.g_nres)
        gnres["A1"] = "sim_id"
        for t in 1:periods; gnres[XLSX.CellRef(1, t+1)] = period_cols[t]; end
        for sid in 1:n_scenarios
            gnres[XLSX.CellRef(sid+1, 1)] = sid
            for t in 1:periods
                gnres[XLSX.CellRef(sid+1, t+1)] = 0.0
            end
        end
    end
    return path
end

# ─── Defaults for each scenario ──────────────────────────────────────────────

function _default_AT_path(sid::Int, activity, periods::Int)
    # AT level path; cumulative growth from period 1.
    if sid == 1                               # baseline
        return ones(Float64, periods)
    elseif sid == 2                           # high TFP (all activities, 1.5%/period)
        return [(1.015)^(t-1) for t in 1:periods]
    elseif sid == 3                           # low TFP (all activities, 0.5%)
        return [(1.005)^(t-1) for t in 1:periods]
    elseif sid == 7                           # crops-only TFP (2%)
        is_crop = startswith(string(activity), "P0") &&
                  parse(Int, string(activity)[2:end]) <= 10
        return is_crop ? [(1.02)^(t-1) for t in 1:periods] :
                          ones(Float64, periods)
    elseif sid == 8                           # industrial sectors TFP (2%)
        ai = parse(Int, string(activity)[2:end])
        return (ai >= 21) ? [(1.02)^(t-1) for t in 1:periods] :
                            ones(Float64, periods)
    elseif sid == 10                          # all growth: 2% TFP everywhere
        return [(1.02)^(t-1) for t in 1:periods]
    else
        return ones(Float64, periods)
    end
end

function _default_g_labor_path(sid::Int, skill, periods::Int)
    if sid == 4                                # fast labor
        return fill(0.03, periods)
    elseif sid == 5                            # slow labor
        return fill(0.01, periods)
    elseif sid == 10                           # all growth
        return fill(0.03, periods)
    else
        return fill(0.02, periods)             # baseline labor growth
    end
end

function _default_g_land_path(sid::Int, periods::Int)
    if sid == 6
        return fill(0.01, periods)
    elseif sid == 10
        return fill(0.005, periods)
    else
        return zeros(Float64, periods)
    end
end

# ─── Reader ──────────────────────────────────────────────────────────────────

"""
    read_policy_scenarios(path::AbstractString) -> Vector{Scenario}

Parse the Excel workbook at `path` and return only the scenarios whose
`active` cell is `"yes"` (case-insensitive).
"""
function read_policy_scenarios(path::AbstractString)
    isfile(path) || error("read_policy_scenarios: file not found: $path")
    scenarios = Scenario[]

    XLSX.openxlsx(path, mode="r") do xf
        scen_data = xf[_POL_SHEETS.scenarios][:]
        at_data    = xf[_POL_SHEETS.AT][:]
        gl_data    = xf[_POL_SHEETS.g_labor][:]
        gland_data = xf[_POL_SHEETS.g_land][:]
        gnres_data = xf[_POL_SHEETS.g_nres][:]

        # Each row of `scenarios` after the header is one scenario.
        for r in 2:size(scen_data, 1)
            sid_raw = scen_data[r, 1]
            sid_raw === missing && continue
            sid = sid_raw isa Real ? Int(sid_raw) : parse(Int, string(sid_raw))
            name = string(coalesce(scen_data[r, 2], "scenario_$(sid)"))
            desc = string(coalesce(scen_data[r, 3], ""))
            T_cell = scen_data[r, 4]
            T = T_cell isa Real ? Int(T_cell) : parse(Int, string(T_cell))
            delta_cell = scen_data[r, 5]
            delta = delta_cell isa Real ? float(delta_cell) : parse(Float64, string(delta_cell))
            active = lowercase(string(coalesce(scen_data[r, 6], "no"))) in ("yes","y","true","1")
            active || continue

            AT       = _read_indexed_block(at_data,    T)
            g_labor  = _read_indexed_block(gl_data,    T)
            g_land   = _read_per_sim_row(gland_data, sid, T)
            g_nres   = _read_per_sim_row(gnres_data, sid, T)

            # Filter AT / g_labor to this scenario only.
            AT_filtered = Dict{Tuple{Any,Int},Float64}()
            for ((s, key, t), val) in AT
                s == sid && (AT_filtered[(key, t)] = val)
            end
            gl_filtered = Dict{Tuple{Any,Int},Float64}()
            for ((s, key, t), val) in g_labor
                s == sid && (gl_filtered[(key, t)] = val)
            end

            push!(scenarios, Scenario(
                sid, name, desc, T, delta,
                AT_filtered, gl_filtered, g_land, g_nres))
        end
    end
    return scenarios
end

# Two-key indexed block: rows = [sim_id, key, period_1..period_T].
# Returns Dict{(sim_id, key, t) => value}.
function _read_indexed_block(data::Matrix, T::Int)
    out = Dict{Tuple{Int,Any,Int},Float64}()
    nrows, ncols = size(data)
    nrows < 2 && return out
    for r in 2:nrows
        sid_cell = data[r, 1]
        sid_cell === missing && continue
        sid = sid_cell isa Real ? Int(sid_cell) : parse(Int, string(sid_cell))
        key = data[r, 2]
        key === missing && continue
        for t in 1:T
            c = data[r, 2+t]
            c === missing && continue
            v = c isa Real ? float(c) : (try parse(Float64, string(c)) catch; NaN end)
            isfinite(v) || continue
            out[(sid, key, t)] = v
        end
    end
    return out
end

# Per-sim single row: rows = [sim_id, period_1..period_T]. Returns Vector.
function _read_per_sim_row(data::Matrix, sid::Int, T::Int)
    out = zeros(Float64, T)
    for r in 2:size(data, 1)
        sid_cell = data[r, 1]
        sid_cell === missing && continue
        rid = sid_cell isa Real ? Int(sid_cell) : parse(Int, string(sid_cell))
        rid == sid || continue
        for t in 1:T
            c = data[r, 1+t]
            c === missing && continue
            out[t] = c isa Real ? float(c) : (try parse(Float64, string(c)) catch; 0.0 end)
        end
        break
    end
    return out
end

# ─── Update step using Scenario values ───────────────────────────────────────

"""
    update_period_data_scenario!(data, m, scen::Scenario, t_next::Int)

Update `data.metadata[:PAR]` from the period-`t` solve using the scenario's
period-`t_next` exogenous values.  This is the scenario-aware version of
`update_period_data!`.
"""
function update_period_data_scenario!(data::LinkageData, m, scen::Scenario, t_next::Int;
        vintage_rule::Symbol=:benchmark_shares,
        delta::Float64=scen.delta)
    PAR = data.metadata[:PAR]
    S   = data.sets
    i = S[:i]; l = S[:l]; v = S[:v]; gz = S[:gz]; ag = S[:ag]; h = S[:h]

    # Capital accumulation: explicit stock in commodity units, converted back to
    # the rental units of KSupply/K0.  Shared with update_period_data!; see the
    # units note at the top of RecursiveDynamic.jl.
    _update_capital_stock!(data, m; delta=delta, vintage_rule=vintage_rule)

    # Labor supply per skill: use scenario growth rates.
    for ll in l
        g = get(scen.g_labor, (ll, t_next), 0.02)
        PAR[:LSupply][ll] *= (1 + g)
        for gg in gz
            PAR[:LS0][(ll, gg)] *= (1 + g)
        end
        PAR[:LY0][ll] *= (1 + g)
    end
    _track_labour_demand!(data, m)

    # Productivity AT per activity: scenario provides absolute LEVELS, not growth.
    # PAR[:AT][i] is replaced with the period-t_next level (1.0 = baseline).
    for ii in i
        at_level = get(scen.AT, (ii, t_next), 1.0)
        # AT[i] in the static model defaults to 1.0; we set it relative to
        # the period-1 level (== AT[1] = 1.0 normalises everything).
        PAR[:AT][ii] = at_level
    end

    # Land supply growth.
    g_t = scen.g_land[t_next]
    if abs(g_t) > 1.0e-12
        for ii in i; PAR[:TSupply][ii] *= (1 + g_t); end
        PAR[:TY0] *= (1 + g_t)
        haskey(PAR[:chi_T], :land) && (PAR[:chi_T][:land] *= (1 + g_t))
        tot_t = sum(get(PAR[:TSupply], ii, 0.0) for ii in ag)
        if tot_t > 1.0e-9
            for ii in i
                PAR[:gamma_T][ii] = (ii in ag) ? PAR[:TSupply][ii] / tot_t : 0.0
            end
        end
    end

    # Natural-resource supply growth.
    n_t = scen.g_nres[t_next]
    if abs(n_t) > 1.0e-12
        for ii in i
            PAR[:FSupply][ii] *= (1 + n_t)
            PAR[:chi_F][ii]   *= (1 + n_t)
        end
        PAR[:FY0] *= (1 + n_t)
    end

    # Update aggregate anchors (KY0 and gamma_K are set by _update_capital_stock!).
    PAR[:YH0] = Dict(hh => max(
        PAR[:TY0] + PAR[:FY0] + sum(values(PAR[:LY0])) + PAR[:KY0], 1.0e-9)
        for hh in h)
    return data
end

# ─── Single-scenario, in-memory runner ───────────────────────────────────────

"""
    run_scenario!(data, scen::Scenario; periods = scen.periods, delta = scen.delta,
                  on_period = nothing, abort = Ref(false), period_modifier = nothing,
                  vintage_rule = :benchmark_shares, show_solver_output = false,
                  show_diagnostics = false, convergence_tolerance = 1.0e-8,
                  time_limit = 3600, verbose = true)
        -> (history, snapshots)

Run ONE `Scenario` over `periods` periods against already-prepared `data`, entirely in
memory: no Excel round trip, no workbook, no plots, nothing written to disk.  This is the
loop body of [`run_policy_experiments!`](@ref) — which now calls it — exposed for callers
that build a `Scenario` in Julia (a web backend, a batch driver, a test).

`data` must already have been through `prepare_data!`; the scenario's period-1 `AT` levels
are applied here, exactly as `run_policy_experiments!` did before every scenario.

Hooks:
- `period_modifier(t, data)` runs BEFORE period `t` is built, so a caller can edit
  `parameters(data)` for levers the `Scenario` struct does not carry (tariffs `:tau_m`,
  taxes `:tau_p`/`:kappa_h`, the government-spending share `:chi_gov`, …).  This is the
  only supported way to give those levers a time path.
- `on_period(t, status, elapsed, history_t)` runs after every period with the PATH
  termination status as a `String`, the wall-clock seconds the build+solve took, and that
  period's `_period_summary` named tuple.
- `abort[]` is checked before each period; when it is `true` the loop stops and returns
  the periods solved so far (possibly none).

`delta` overrides the scenario's depreciation rate in the capital update; `periods`
overrides its horizon.  Non-converged periods behave exactly as in the other drivers: the
state is not propagated, so the trajectory plateaus instead of being corrupted.
"""
function run_scenario!(data::LinkageData, scen::Scenario;
        periods::Int=scen.periods,
        delta::Float64=scen.delta,
        on_period=nothing,
        abort::Ref{Bool}=Ref(false),
        period_modifier=nothing,
        vintage_rule::Symbol=:benchmark_shares,
        show_solver_output::Bool=false,
        show_diagnostics::Bool=false,
        convergence_tolerance::Float64=1.0e-8,
        time_limit::Real=3600,
        verbose::Bool=true)

    PAR = parameters(data)
    # Initial AT (period 1) — apply the scenario's period-1 levels right away.
    for ii in data.sets[:i]
        PAR[:AT][ii] = get(scen.AT, (ii, 1), 1.0)
    end

    history   = Vector{NamedTuple}()
    snapshots = Vector{Dict{Tuple{Symbol,Any},Float64}}()

    for t in 1:periods
        abort[] && break
        period_modifier === nothing || period_modifier(t, data)
        verbose && println("\n--- Scenario $(scen.name) — Period t = $(t) / $(periods) ---")
        elapsed = @elapsed begin
            m = model(data; show_solver_output=show_solver_output)
            solve_model!(m; show_diagnostics=show_diagnostics,
                         convergence_tolerance=convergence_tolerance,
                         time_limit=time_limit)
        end
        push!(snapshots, _collect_period_values(m))
        summary = _period_summary(m, data, t)
        push!(history, summary)
        if on_period !== nothing
            status = try string(termination_status(m)) catch; "UNKNOWN" end
            on_period(t, status, elapsed, summary)
        end
        if t < periods
            update_period_data_scenario!(data, m, scen, t+1;
                vintage_rule=vintage_rule, delta=delta)
        end
    end

    return history, snapshots
end

# ─── Main scenario-runner ────────────────────────────────────────────────────

"""
    run_policy_experiments!(xlsx_path; outdir="results/scenarios", make_plots=true)

Read the policy-experiment workbook and run every active scenario.  Each
scenario writes its own `dynamic_results.xlsx`, `time_series.csv`, and
optional plot folder under `outdir/<scenario_name>/`.

Returns `Dict(name => (data, history, snapshots))` for all scenarios.
"""
function run_policy_experiments!(xlsx_path::AbstractString;
        outdir::AbstractString="results/scenarios",
        make_plots::Bool=true,
        show_solver_output::Bool=false)

    scenarios = read_policy_scenarios(xlsx_path)
    isempty(scenarios) && error("No active scenarios in $xlsx_path")
    mkpath(outdir)
    results = Dict{String, Tuple}()

    for scen in scenarios
        println("\n################################################################")
        println("  SCENARIO ", scen.sim_id, " : ", scen.name)
        println("  ", scen.description)
        println("  periods=", scen.periods, "  delta=", scen.delta)
        println("################################################################")

        data = init_data()
        prepare_data!(data)
        history, snapshots = run_scenario!(data, scen;
            show_solver_output=show_solver_output)

        results[scen.name] = (data, history, snapshots)
    end

    # ── Single consolidated workbook holding ALL scenarios ───────────────────
    xlsx_path = _write_all_scenarios_xlsx(scenarios, results, outdir)
    println("\nConsolidated results: ", xlsx_path)

    # ── Cross-scenario comparison plots ───────────────────────────────────────
    if make_plots
        try
            files = plot_all_scenarios(xlsx_path; outdir=joinpath(outdir, "plots"))
            println("Comparison plots:     ", joinpath(outdir, "plots"), "  (", length(files), " files)")
        catch err
            @warn "Cross-scenario plotting failed" exception=(err, catch_backtrace())
        end
    end

    return results
end

"""
Write ALL scenarios into one Excel workbook.

Sheet layout (per variable family):
- Row 1: ["scenario", "index", "period_1", ..., "period_T"]
- Subsequent rows: [scenario_name, index_label, val_t1, val_t2, ..., val_tT]

Sheets included:
- `scenarios`       — copy of the scenario design table
- `macro_summary`   — named macro indicators (one row per scenario × indicator)
- One sheet per JuMP variable family (XP, GDP, YH, ...)
"""
function _write_all_scenarios_xlsx(scenarios::Vector{Scenario}, results::Dict, outdir::AbstractString)
    mkpath(outdir)
    path = joinpath(outdir, "all_scenarios_results.xlsx")

    # Discover every variable name across every scenario × period snapshot,
    # and every index seen for that variable.
    var_indices = Dict{Symbol, Set{Any}}()
    max_T = 0
    for scen in scenarios
        haskey(results, scen.name) || continue
        _, _, snaps = results[scen.name]
        max_T = max(max_T, length(snaps))
        for snap in snaps
            for ((name, idx), _) in snap
                push!(get!(var_indices, name, Set{Any}()), idx)
            end
        end
    end

    XLSX.openxlsx(path, mode="w") do xf
        # ── Sheet: scenarios design ──────────────────────────────────────────
        sd = if XLSX.sheetnames(xf)[1] == "Sheet1"
            XLSX.rename!(xf[1], "scenarios"); xf[1]
        else
            XLSX.addsheet!(xf, "scenarios")
        end
        sd["A1"] = "sim_id"; sd["B1"] = "name"; sd["C1"] = "description"
        sd["D1"] = "periods"; sd["E1"] = "delta"
        for (i, s) in enumerate(scenarios)
            sd[XLSX.CellRef(i+1, 1)] = s.sim_id
            sd[XLSX.CellRef(i+1, 2)] = s.name
            sd[XLSX.CellRef(i+1, 3)] = s.description
            sd[XLSX.CellRef(i+1, 4)] = s.periods
            sd[XLSX.CellRef(i+1, 5)] = s.delta
        end

        # ── Sheet: macro_summary (named indicators, one row per scenario) ────
        ms = XLSX.addsheet!(xf, "macro_summary")
        ms["A1"] = "scenario"; ms["B1"] = "indicator"
        for t in 1:max_T; ms[XLSX.CellRef(1, t+2)] = "period_$(t)"; end
        row = 2
        for scen in scenarios
            haskey(results, scen.name) || continue
            history = results[scen.name][2]
            isempty(history) && continue
            indicators = collect(keys(history[1]))
            for ind in indicators
                ms[XLSX.CellRef(row, 1)] = scen.name
                ms[XLSX.CellRef(row, 2)] = string(ind)
                for t in 1:length(history)
                    v = getfield(history[t], ind)
                    if v isa Real && isfinite(v)
                        ms[XLSX.CellRef(row, t+2)] = float(v)
                    end
                end
                row += 1
            end
        end

        # ── One sheet per variable family ────────────────────────────────────
        for vname in sort(collect(keys(var_indices)); by=string)
            sheet_name = string(vname)
            length(sheet_name) > 31 && (sheet_name = sheet_name[1:31])
            sh = XLSX.addsheet!(xf, sheet_name)

            sh["A1"] = "scenario"
            sh["B1"] = "index"
            for t in 1:max_T; sh[XLSX.CellRef(1, t+2)] = "period_$(t)"; end

            indices = sort(collect(var_indices[vname]); by=_index_sort_key)
            row = 2
            for scen in scenarios
                haskey(results, scen.name) || continue
                _, _, snaps = results[scen.name]
                for idx in indices
                    sh[XLSX.CellRef(row, 1)] = scen.name
                    sh[XLSX.CellRef(row, 2)] = _index_label(idx)
                    for t in 1:length(snaps)
                        v = get(snaps[t], (vname, idx), nothing)
                        if v isa Real && isfinite(v)
                            sh[XLSX.CellRef(row, t+2)] = float(v)
                        end
                    end
                    row += 1
                end
            end
        end
    end
    return path
end

# Older per-scenario column-layout writer removed — replaced by
# `_write_all_scenarios_xlsx` which keeps one sheet per variable family with
# scenarios as rows (cleaner for filtering and plotting).
