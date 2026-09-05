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

    # Trade closure.  The default :bop keeps the SAM's own imports, exports and
    # final demand (no rescale), books the trade deficit as exogenous foreign
    # saving and adds the real exchange rate ER; :balanced is the legacy
    # convention (imports = (1+tau_m)(1+tau_e) * exports good by good, final
    # demand rescaled, no ER).  Both must stay square.
    db = prepare_data!(init_data(); outdir=nothing, trade_closure=:balanced)
    @test data.par[:trade_closure] == :bop && db.par[:trade_closure] == :balanced
    @test haskey(m, :ER) && !haskey(model(db; show_solver_output=false), :ER)
    mb = model(db; show_solver_output=false)
    @test num_variables(mb) == num_constraints(mb; count_variable_in_set_constraints=false)
    @test all(v == 0.0 for v in values(db.par[:Sfbar]))
    let M = data.balanced_sam, idx = data.sam_index, iset = data.sets[:i]
        hh_sam = sum(max(M[idx["COM_"*p], idx["HH"]], 0.0) for p in iset)
        @test isapprox(data.par[:bench][:HH], hh_sam; rtol=1e-10)      # :bop: no final-demand rescale
        @test db.par[:bench][:HH] != data.par[:bench][:HH]              # :balanced rescales it
        # CIF imports (incl. margins on imports) - exports - margin sales - export tax
        deficit = sum(M[idx["ROW"], idx["COM_"*p]] + M[idx["TRD_MRG"], idx["COM_"*p]] -
                      M[idx["COM_"*p], idx["ROW"]] - M[idx["COM_"*p], idx["TRD_MRG"]] for p in iset) -
                  M[idx["TAX_EXP"], idx["ROW"]]
        @test isapprox(sum(values(data.par[:Sfbar])), deficit; rtol=1e-8)
    end
    @test_throws ErrorException prepare_data!(init_data(); outdir=nothing, trade_closure=:foo)
    # :fixed_er swaps the complementary variable of the balance of payments.
    df_ = prepare_data!(init_data(); outdir=nothing, bop_closure=:fixed_er)
    mf = model(df_; show_solver_output=false)
    @test haskey(mf, :C_ER) && num_variables(mf) == num_constraints(mf; count_variable_in_set_constraints=false)

    # A real country SAM (Kenya, EMERGING/GTAP hybrid, 65 goods) with its own
    # trade deficit: loads, calibrates without a final-demand rescale and builds square.
    ken_sam  = joinpath(@__DIR__, "data", "KEN_2018_hybrid_sam.csv")
    ken_sets = joinpath(@__DIR__, "data", "KEN_2018_hybrid_sets.csv")
    dk = prepare_data!(init_data(); sets_path=ken_sets, source=:csv, sam_path=ken_sam, balance=:none, outdir=nothing)
    @test sum(values(dk.par[:Sfbar])) > 0                                # Kenya runs a trade deficit
    mk = model(dk; show_solver_output=false)
    @test num_variables(mk) == num_constraints(mk; count_variable_in_set_constraints=false)

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
        @test abs(value(m[:ER]) - 1) < 1e-6                              # ER = 1 at the benchmark

        # Both closures and the Kenya SAM replicate their benchmarks (the real SAM to
        # 0.5 %: its calibrated start point has a 1e-3 residual in P-5).
        for (mm, dd_, tol) in ((mb, db, 1e-4), (mf, df_, 1e-4), (mk, dk, 5e-3))
            solve_model!(mm; output="no", show_diagnostics=false)
            @test termination_status(mm) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)
            @test maximum(abs(value(mm[:XP][p]) / dd_.par[:bench][:XP][p] - 1) for p in dd_.sets[:i]) < tol
        end

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
