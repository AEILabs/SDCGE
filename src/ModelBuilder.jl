# Usage:
#   data = init_data()
#   prepare_data!(data)                         # default 100-sector SAM → balance → calibrate
#   prepare_data!(data; sets_path="data/csv/sets.csv", source=:csv, sam_path="data/csv/sam.csv")
#                                               # N-sector run: sets from file, then the SAM
#   m = model(data)                             # creates Model(PATHSolver.Optimizer) and builds LINKAGE
#   solve_model!(m)                             # calls optimize!(m)
#
# One-call alternatives:
#   m, data = run_linkage!()                    # full default pipeline and solve
#   m = build_model(data; optimizer=PATHSolver.Optimizer)
#
# This file is the clean public API. User scripts should call these functions instead
# of manually sequencing every low-level SAM/calibration/build routine.

"""Prepare LINKAGE data for model construction.

Keyword arguments:
- `sets_path = nothing`: read the sets from a `set,item` CSV with `read_sets_csv!`
  before any default is applied — this is how the model runs with N != 100 sectors.
  The file must define `i` (activities = products) and, unless N == 100, also the
  memberships `cr`, `lv`, `e`, `ft`, `fd` (with `ft`/`fd` disjoint from `e`); `ag`,
  `ip`, `nf`, `nnft`, `nnfd` are derived, and anything the file omits (`r`, `v`, `l`,
  ...) keeps its `default_sets!` value.  The SAM then needs 2N+16 accounts.
- `source = :default`: use the internally generated N-sector SAM.
- `source = :csv`: read SAM from `sam_path` using `read_sam_csv!`.
- `source = :excel`: read SAM from `sam_path` using `read_sam_excel!`.
- `balance = :ras`: apply RAS balancing. Use `balance = :none` to skip.
- `calibrate = true`: run `calibrate_from_sam!`.
- `precompute = true`: store `data.metadata[:PAR] = precompute_parameters(data)`.
- `outdir = "results"`: directory for the SAM balance report and the balanced SAM.
  Pass `outdir = nothing` to skip those writes entirely — needed by embedding callers
  (a web backend, a batch driver) that must not write into the process's working
  directory, and by anything that calls `prepare_data!` concurrently, since the two
  files would otherwise be overwritten by every call.
"""
function prepare_data!(data::LinkageData=init_data();
        source::Symbol=:default,
        sam_path::Union{Nothing,String}=nothing,
        sets_path::Union{Nothing,AbstractString}=nothing,
        balance::Symbol=:ras,
        calibrate::Bool=true,
        precompute::Bool=true,
        outdir::Union{Nothing,AbstractString}="results")

    sets_path === nothing || read_sets_csv!(data, sets_path)
    default_sets!(data)
    setup_sam_accounts!(data)

    if source == :default
        build_default_large_sam!(data)
    elseif source == :csv
        sam_path === nothing && error("prepare_data!(; source=:csv) requires sam_path.")
        read_sam_csv!(data, sam_path)
    elseif source == :excel
        sam_path === nothing && error("prepare_data!(; source=:excel) requires sam_path.")
        read_sam_excel!(data, sam_path)
    else
        error("Unknown SAM source: $(source). Use :default, :csv, or :excel.")
    end

    # The raw SAM read from file or generated internally may be unbalanced.
    # Check only the accounting structure here; do not require balance until
    # after the balancing step below.
    validate_sam!(data; require_balanced=false)

    if balance == :ras
        balance_sam_ras!(data)
    elseif balance == :none
        data.balanced_sam = copy(data.sam)
    else
        error("Unknown balance method: $(balance). Use :ras or :none.")
    end

    # The model is calibrated and initialized only from the balanced SAM.
    # `balance_sam_ras!` also replaces `data.sam` with the balanced matrix,
    # so downstream code cannot accidentally use the unbalanced raw SAM.
    assert_balanced_sam!(data; tol=1.0e-6)

    # Write diagnostics and the balanced SAM for inspection/re-use (skipped when
    # `outdir === nothing`).
    if outdir !== nothing
        export_sam_balance_report!(data; outdir=outdir)
        export_balanced_sam!(data; path=joinpath(outdir, "balanced_sam.csv"))
    end

    if calibrate
        calibrate_from_sam!(data)
    end

    if precompute
        data.metadata[:PAR] = precompute_parameters(data)
    end

    return data
end

"""Return the precomputed parameter table, computing it if needed."""
function parameters(data::LinkageData)
    if !haskey(data.metadata, :PAR)
        data.metadata[:PAR] = precompute_parameters(data)
    end
    return data.metadata[:PAR]
end

"""Low-level in-place model builder.

This function assumes a JuMP model has already been initialized. It is kept for
advanced users who want to control the solver object manually.
"""
function build_linkage_model!(m, data::LinkageData)
    default_sets!(data)
    if length(data.sam_accounts) == 0
        setup_sam_accounts!(data)
    end
    if size(data.balanced_sam, 1) == 0
        build_default_large_sam!(data)
        balance_sam_ras!(data)
    end
    assert_balanced_sam!(data; tol=1.0e-6)
    if !haskey(data.par, :intermediate_share)
        calibrate_from_sam!(data)
    end

    PAR = parameters(data)

    add_variables!(m, data)
    initialize_from_sam!(m, data)
    check_initialization!(m; error_on_bad=true)
    add_production_equations!(m, data, PAR)
    add_income_equations!(m, data, PAR)
    add_demand_equations!(m, data, PAR)
    add_trade_equations!(m, data, PAR)
    add_equilibrium_equations!(m, data, PAR)
    add_closure_equations!(m, data, PAR)
    add_factor_equations!(m, data, PAR)
    add_other_equations!(m, data, PAR)
    # Dynamics: not part of the static MCP. Recursive dynamics is implemented
    # externally by RecursiveDynamic.jl (single trajectory) and
    # PolicyScenarios.jl (multi-scenario experiments).
    return m
end


"""Print a compact model-size summary before solving."""
function print_model_diagnostics(m)
    println("\n================ LINKAGE MCP MODEL ================")
    try
        println("Variables:   ", num_variables(m))
    catch err
        println("Variables:   unavailable (", typeof(err), ")")
    end
    try
        println("Constraints: ", num_constraints(m; count_variable_in_set_constraints=false))
    catch err
        println("Constraints: unavailable (", typeof(err), ")")
    end
    try
        println("Solver:      ", solver_name(m))
    catch
        println("Solver:      PATHSolver.Optimizer")
    end
    println("===================================================\n")
    return nothing
end

"""Print the solve status returned by PATH."""
function print_solver_status(m, status=nothing)
    println("\n================ PATH SOLVER STATUS ================")
    if status !== nothing
        println("optimize! returned: ", status)
    end
    try
        println("Solver:             ", solver_name(m))
    catch
        println("Solver:             PATHSolver.Optimizer")
    end
    try
        println("Termination status: ", termination_status(m))
    catch err
        println("Termination status unavailable (", typeof(err), ")")
    end
    try
        println("Primal status:      ", primal_status(m))
    catch err
        println("Primal status unavailable (", typeof(err), ")")
    end
    println("====================================================\n")
    return nothing
end

"""Create and build a JuMP/MCP model from calibrated LINKAGE data.

`show_solver_output = true` is stored on the model and used by `solve_model!`.
The LINKAGE-PATH version has no objective function.
"""
function build_model(data::LinkageData=prepare_data!();
        optimizer=PATHSolver.Optimizer,
        optimizer_attributes=Dict{String,Any}(),
        show_solver_output::Bool=true)

    # Build the JuMP model with PATH attached immediately.
    # This avoids the confusing display `solver: none` and ensures
    # complementarity constraints are sent directly to PATH.
    opt = optimizer === nothing ? PATHSolver.Optimizer : optimizer
    m = Model(opt)

    # Store and apply PATH options.
    path_output = show_solver_output ? "yes" : "no"
    m.ext[:path_output] = path_output
    m.ext[:path_optimizer_attributes] = optimizer_attributes

    # PATHSolver versions differ in accepted option names, so ignore only
    # unsupported attributes while keeping the model solver attached.
    try
        set_optimizer_attribute(m, "output", path_output)
    catch err
        @warn "Could not set PATH output option" exception=(err, catch_backtrace())
    end
    for (k, v) in optimizer_attributes
        try
            set_optimizer_attribute(m, k, v)
        catch err
            @warn "Could not set optimizer attribute" attribute=k value=v exception=(err, catch_backtrace())
        end
    end

    build_linkage_model!(m, data)
    return m
end

"""Convenience alias requested by the user: `model(data)` returns a built JuMP/MCP model."""
function model(data::LinkageData=prepare_data!();
        optimizer=nothing,
        optimizer_attributes=Dict{String,Any}(),
        show_solver_output::Bool=true)
    return build_model(data; optimizer=optimizer,
        optimizer_attributes=optimizer_attributes,
        show_solver_output=show_solver_output)
end

"""Optimize an already-built MCP model and show PATH progress on screen.

Examples:
```julia
m, data = run_linkage!(solve=false)
solve_model!(m; output="yes")
solve_model!(m; output="verbose")
```
"""
function solve_model!(m;
        solver::Symbol=:PATH,
        convergence_tolerance::Float64=1.0e-8,
        output::Union{Nothing,AbstractString}=nothing,
        time_limit::Real=3600,
        show_diagnostics::Bool=true)

    check_initialization!(m; error_on_bad=true)

    path_output = output === nothing ? get(m.ext, :path_output, "yes") : output

    # Make sure PATH output is visible at solve time.
    try
        set_optimizer_attribute(m, "output", path_output)
    catch err
        @warn "Could not set PATH output option at solve time" exception=(err, catch_backtrace())
    end
    try
        set_optimizer_attribute(m, "convergence_tolerance", convergence_tolerance)
    catch
    end
    try
        set_optimizer_attribute(m, "time_limit", time_limit)
    catch
    end

    if show_diagnostics
        print_model_diagnostics(m)
        println("PATH options:")
        println("  output                = ", path_output)
        println("  convergence_tolerance = ", convergence_tolerance)
        println("  time_limit            = ", time_limit)
        println("\nStarting PATH solve...\n")
    end

    optimize!(m)

    if show_diagnostics
        print_solver_status(m, nothing)
    end

    return m
end

"""Build and solve LINKAGE from prepared data."""
function solve_linkage!(data::LinkageData=prepare_data!();
        optimizer=nothing,
        optimizer_attributes=Dict{String,Any}(),
        show_solver_output::Bool=true)
    m = model(data; optimizer=optimizer,
        optimizer_attributes=optimizer_attributes,
        show_solver_output=show_solver_output)
    solve_model!(m; output=show_solver_output ? "yes" : "no")
    return m
end

"""Full one-call pipeline: prepare data, build model, solve, and return `(model, data)`."""
function run_linkage!(; source::Symbol=:default,
        sam_path::Union{Nothing,String}=nothing,
        balance::Symbol=:ras,
        optimizer=nothing,
        optimizer_attributes=Dict{String,Any}("tol" => 1.0e-6),
        solve::Bool=true,
        write_results::Bool=false,
        make_plots::Bool=false,
        results_dir::AbstractString="results",
        show_solver_output::Bool=true,
        path_output::AbstractString="yes")

    data = prepare_data!(init_data(); source=source, sam_path=sam_path, balance=balance)
    m = model(data; optimizer=optimizer,
        optimizer_attributes=optimizer_attributes,
        show_solver_output=show_solver_output)
    if solve
        solve_model!(m; output=show_solver_output ? path_output : "no")
        if write_results || make_plots
            export_results!(m, data; outdir=results_dir)
        end
        if make_plots
            plot_results(results_dir; outdir=joinpath(results_dir, "plots"))
        end
    end
    return m, data
end

# Backward-compatible alias.
run_model! = run_linkage!
