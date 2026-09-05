# Smoke test: SAM preparation (default, CSV, sets-from-file and a 12-sector
# GTAP-coded run), model build, squareness and results export. Does not call
# PATH unless LCGE_TEST_SOLVE=true is set.
# Run with:  julia --project=. test/runtests.jl   (or Pkg.test())

using Test
include(joinpath(@__DIR__, "..", "src", "LinkageModel.jl"))
using .LinkageModel
using JuMP
using DataFrames

@testset "LCGE-V4 smoke test" begin
    data = prepare_data!()
    # 2N + 16 accounts: N activities, N commodities, 5 factors, 6 taxes,
    # 4 institutions, 1 margin account.
    @test length(data.sam_accounts[:all]) == 2 * length(data.sets[:i]) + 16
    @test sam_balance_summary(data)[:max_abs_gap] < 1e-6

    # Bring-your-own SAM path (CSV).
    d2 = init_data()
    default_sets!(d2)
    setup_sam_accounts!(d2)
    read_sam_csv!(d2, joinpath(@__DIR__, "..", "data", "csv", "sam.csv"))
    n2 = 2 * length(d2.sets[:i]) + 16
    @test size(d2.sam) == (n2, n2)

    # Sets read from file: the shipped sets.csv must reproduce the default sets
    # and, through them, exactly the same balanced SAM and benchmark.
    sam_csv  = joinpath(@__DIR__, "..", "data", "csv", "sam.csv")
    sets_csv = joinpath(@__DIR__, "..", "data", "csv", "sets.csv")
    d3 = prepare_data!(init_data(); sets_path=sets_csv, source=:csv, sam_path=sam_csv, outdir=nothing)
    d4 = prepare_data!(init_data(); source=:csv, sam_path=sam_csv, outdir=nothing)
    @test Set(keys(d3.sets)) == Set(keys(d4.sets))
    @test all(d3.sets[k] == d4.sets[k] for k in keys(d4.sets))
    @test d3.balanced_sam == d4.balanced_sam
    @test d3.par[:bench][:kappa] == d4.par[:bench][:kappa]

    # N != 100: 12 GTAP-coded sectors with the memberships supplied explicitly
    # (as read_sets_csv! supplies them).  |e| = 3, |ft| = 1, |fd| = 2, and a
    # single region — the other end of the |r| range the trade block must handle.
    d5 = init_data()
    d5.sets[:i]  = ["pdr","wht","gro","ctl","oap","coa","oil","ely","chm","tex","trd","osg"]
    d5.sets[:cr] = ["pdr","wht","gro"]
    d5.sets[:lv] = ["ctl","oap"]
    d5.sets[:e]  = ["coa","oil","ely"]
    d5.sets[:ft] = ["chm"]
    d5.sets[:fd] = ["wht","gro"]
    d5.sets[:r]  = ["R1"]
    prepare_data!(d5; outdir=nothing)
    @test length(d5.sam_accounts[:all]) == 2 * 12 + 16
    @test d5.sets[:ag] == vcat(d5.sets[:cr], d5.sets[:lv])
    @test d5.sets[:ip] == d5.sets[:nf] == ["coa","oil","ely","chm","tex","trd","osg"]
    @test d5.sets[:rp] == ["R1"]          # rp follows r, so |r| = 1 stays square
    m12 = model(d5; show_solver_output=false)
    @test num_variables(m12) == num_constraints(m12; count_variable_in_set_constraints=false)

    # A SAM whose labels are not the 2N+16 the sets imply must be rejected by the
    # reader (naming the mismatches), not deep inside calibrate_from_sam!.
    @test_throws ErrorException read_sam_csv!(d5, sam_csv)

    # A non-default sector list without :cr/:lv/:e/:ft/:fd must fail loudly:
    # the positional defaults would otherwise be disjoint from :i.
    d6 = init_data(); d6.sets[:i] = ["a", "b", "c"]
    @test_throws ErrorException default_sets!(d6)

    m = model(data; show_solver_output=false)
    nv = num_variables(m)
    nc = num_constraints(m; count_variable_in_set_constraints=false)
    println("variables = $nv, constraints = $nc")
    @test nv == nc

    df = results_dataframe(m)
    @test nrow(df) == nv
    @test !any(occursin("CartesianIndex", string(x)) for x in df.index)

    if get(ENV, "LCGE_TEST_SOLVE", "false") == "true"
        # The 12-sector run must replicate its own benchmark just like the 100-sector one.
        solve_model!(m12; output="no", show_diagnostics=false)
        println("12-sector termination_status = ", termination_status(m12))
        @test termination_status(m12) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)
        @test maximum(abs(value(m12[:XP][p]) / d5.par[:bench][:XP][p] - 1) for p in d5.sets[:i]) < 1e-4

        solve_model!(m; output="no", show_diagnostics=false)
        println("termination_status = ", termination_status(m))
        @test termination_status(m) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)

        # Recursive dynamics with zero growth must reproduce the benchmark every
        # period: the capital stock is stationary when I = delta*K, which only
        # holds if update_period_data! keeps stock and rental units apart.
        # Driven directly (not through run_recursive_dynamic!) to skip the XLSX
        # writer and keep the test fast.
        dd = prepare_data!(init_data())
        rgdp = Float64[]; ks = Float64[]; yg = Float64[]
        for t in 1:2
            md = model(dd; show_solver_output=false)
            solve_model!(md; output="no", show_diagnostics=false)
            @test termination_status(md) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)
            push!(rgdp, value(md[:RGDP]["R1"]))
            push!(ks,   value(md[:KS]))
            push!(yg,   value(md[:YG]))
            t < 2 && update_period_data!(dd, md; delta=0.05, g_labor=0.0, g_tfp=0.0)
        end
        println("zero-growth dynamics: RGDP = ", rgdp, "  KS = ", ks)
        @test abs(rgdp[2]/rgdp[1] - 1) < 1e-6
        @test abs(ks[2]/ks[1]     - 1) < 1e-6
        @test abs(yg[2]/yg[1]     - 1) < 1e-6
    end
end
