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

    # An economy without livestock (or fertiliser) production: the set is empty, the
    # model still builds square (ag = cr).
    d5b = init_data()
    d5b.sets[:i]  = d5.sets[:i]; d5b.sets[:cr] = d5.sets[:cr]; d5b.sets[:lv] = String[]
    d5b.sets[:e]  = d5.sets[:e]; d5b.sets[:ft] = String[]; d5b.sets[:fd] = d5.sets[:fd]; d5b.sets[:r] = ["R1"]
    prepare_data!(d5b; outdir=nothing)
    @test d5b.sets[:ag] == d5b.sets[:cr]
    m12b = model(d5b; show_solver_output=false)
    @test num_variables(m12b) == num_constraints(m12b; count_variable_in_set_constraints=false)

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

    # :wage_floor's benchmark unemployment (no PATH needed): the CSV reader, the default read
    # next to a CSV SAM, and set_benchmark_unemployment! (idempotent: employment stays fixed).
    @testset "benchmark unemployment" begin
        tmp = mktempdir()
        for f in ("sam.csv", "sets.csv"); cp(joinpath(@__DIR__, "..", "data", "csv", f), joinpath(tmp, f)); end
        open(joinpath(tmp, "unemployment.csv"), "w") do io
            println(io, "labour,rate"); println(io, "UnSkLab,0.07"); println(io, "SkLab,0.03")
        end
        @test read_unemployment_csv(joinpath(tmp, "unemployment.csv")) == Dict("UnSkLab" => 0.07, "SkLab" => 0.03)
        du = prepare_data!(init_data(); sets_path=joinpath(tmp, "sets.csv"), source=:csv,
                           sam_path=joinpath(tmp, "sam.csv"), outdir=nothing)
        @test du.par[:ue_data] == Dict("UnSkLab" => 0.07, "SkLab" => 0.03)
        Pu = parameters(du)
        @test all(Pu[:UE0][ll] == 0.0 for ll in du.sets[:l])        # nothing applied before a :wage_floor build
        emp(P, ll) = P[:LSupply][ll] * (1 - P[:UE0][ll])
        e0 = Dict(ll => emp(Pu, ll) for ll in du.sets[:l])
        set_benchmark_unemployment!(du, 0.10); set_benchmark_unemployment!(du, du.par[:ue_data])
        @test all(isapprox(Pu[:UE0][ll], du.par[:ue_data][ll]) for ll in du.sets[:l])
        @test all(isapprox(emp(Pu, ll), e0[ll]; rtol=1e-12) for ll in du.sets[:l])
        @test all(isapprox(Pu[:LS0][(ll,"national")], e0[ll] / (1 - Pu[:UE0][ll]); rtol=1e-12) for ll in du.sets[:l])
        @test_throws ErrorException set_benchmark_unemployment!(du, 0.97)
    end
end

# ── Numeric regression: the default closure on a real database ───────────────
# Solves the committed Kenya database (test/data/KEN_2018_hybrid_*, 65 goods)
# under the shipped defaults (:fixed_wage, bop_closure = :flex_er, trade_closure
# = :bop), at the benchmark and under a uniform -10 % TFP shock, and compares
# the count, sum and sum of squares of every variable family (XP, VA, PABS, ...)
# against test/reference/KEN_2018_hybrid_default.tsv at rtol 1e-8.  Any change
# to a solved value moves its family's sums, so this guards the promise that a
# change to the model leaves the default closure's solutions untouched.  It
# needs PATH, so it runs only with LCGE_TEST_SOLVE=true.  After a change that is
# MEANT to move the default solution, regenerate the reference deliberately:
#     LCGE_TEST_SOLVE=true LCGE_WRITE_REFERENCE=true julia --project=. test/runtests.jl
if get(ENV, "LCGE_TEST_SOLVE", "false") == "true"
    @testset "default-closure regression (KEN_2018_hybrid)" begin
        ref_path = joinpath(@__DIR__, "reference", "KEN_2018_hybrid_default.tsv")
        ken_sam  = joinpath(@__DIR__, "data", "KEN_2018_hybrid_sam.csv")
        ken_sets = joinpath(@__DIR__, "data", "KEN_2018_hybrid_sets.csv")
        write_ref = get(ENV, "LCGE_WRITE_REFERENCE", "false") == "true"

        family_stats(m) = begin
            acc = Dict{String,Vector{Float64}}()
            for v in all_variables(m)
                fam = String(first(split(name(v), '[')))
                x = value(v)
                s = get!(acc, fam, [0.0, 0.0, 0.0]); s[1] += 1; s[2] += x; s[3] += x * x
            end
            acc
        end

        got = Dict{Tuple{String,String},Vector{Float64}}()
        for (tag, at) in (("bench", 1.0), ("tfp090", 0.90))
            d = prepare_data!(init_data(); sets_path=ken_sets, source=:csv, sam_path=ken_sam,
                              balance=:none, outdir=nothing)
            PAR = parameters(d)
            @test PAR[:labour_closure] == :fixed_wage           # the default this test protects
            for k in keys(PAR[:AT]); PAR[:AT][k] = at; end
            mr = model(d; show_solver_output=false)
            solve_model!(mr; output="no", show_diagnostics=false)
            println("regression ", tag, ": termination_status = ", termination_status(mr))
            @test termination_status(mr) == JuMP.LOCALLY_SOLVED
            for (fam, s) in family_stats(mr); got[(tag, fam)] = s; end
        end

        if write_ref
            mkpath(dirname(ref_path))
            open(ref_path, "w") do io
                println(io, "scenario\tfamily\tn\tsum\tsumsq")
                for ((tag, fam), s) in sort!(collect(got); by=first)
                    println(io, tag, '\t', fam, '\t', Int(s[1]), '\t', repr(s[2]), '\t', repr(s[3]))
                end
            end
            println("wrote ", ref_path)
        else
            @test isfile(ref_path)
            ref = Dict{Tuple{String,String},Vector{Float64}}()
            for line in Iterators.drop(eachline(ref_path), 1)
                f = split(line, '\t')
                ref[(String(f[1]), String(f[2]))] = [parse(Float64, f[3]), parse(Float64, f[4]), parse(Float64, f[5])]
            end
            @test Set(keys(ref)) == Set(keys(got))               # same variable families, both scenarios
            nbad = 0
            for (k, r) in ref
                haskey(got, k) || continue
                g = got[k]
                ok = g[1] == r[1] && isapprox(g[2], r[2]; rtol=1e-8, atol=1e-12) &&
                     isapprox(g[3], r[3]; rtol=1e-8, atol=1e-12)
                ok || (nbad += 1; println("  regression mismatch ", k, ": got ", g, " ref ", r))
            end
            @test nbad == 0
        end
    end
end

# ── :wage_floor on a real database (needs PATH) ──────────────────────────────
# KEN_2018_hybrid at the benchmark scale the simulator calibrates at (largest flow 1e5),
# trade_closure = :bop, bop_closure = :fixed_er.  A null run reproduces the benchmark with 5 %
# unemployment; a 0 % start is :full_employment exactly; the wage stays on its floor while a
# skill has unemployment and rises once it has none; the :flex_er pairing is refused.
if get(ENV, "LCGE_TEST_SOLVE", "false") == "true"
    @testset "wage_floor closure (KEN_2018_hybrid)" begin
        ken_sam  = joinpath(@__DIR__, "data", "KEN_2018_hybrid_sam.csv")
        ken_sets = joinpath(@__DIR__, "data", "KEN_2018_hybrid_sets.csv")
        function prep(bop)
            d = prepare_data!(init_data(); sets_path=ken_sets, source=:csv, sam_path=ken_sam, balance=:none,
                              outdir=nothing, bop_closure=bop, calibrate=false, precompute=false)
            sc = 1e5 / maximum(abs, d.balanced_sam); d.balanced_sam .*= sc; d.sam === d.balanced_sam || (d.sam .*= sc)
            LinkageModel.calibrate_from_sam!(d); d.metadata[:PAR] = LinkageModel.precompute_parameters(d)
            return d
        end
        function solve_closure(closure; u=nothing, tm=1.0, bop=:fixed_er)
            d = prep(bop); P = parameters(d); P[:labour_closure] = closure
            u === nothing || set_benchmark_unemployment!(d, u)
            for k in keys(P[:tau_m]); P[:tau_m][k] *= tm; end
            m = model(d; show_solver_output=false)
            solve_model!(m; output="no", show_diagnostics=false)
            @test termination_status(m) == JuMP.LOCALLY_SOLVED
            return m, d
        end
        L = ["UnSkLab", "SkLab"]
        ue(m, ll) = value(m[:UE][ll, "national"])
        floorgap(m, ll) = value(m[:TW][ll, "national"]) / value(m[:WMIN][ll, "national"]) - 1

        # refused without a fixed exchange rate
        dflex = prep(:flex_er); parameters(dflex)[:labour_closure] = :wage_floor
        @test_throws ErrorException model(dflex; show_solver_output=false)

        # null run: the benchmark, with 5 % unemployment at the floor wage
        m, d = solve_closure(:wage_floor; u=0.05)
        xp0 = parameters(d)[:bench][:XP]
        @test maximum(abs(value(m[:XP][ii]) / xp0[ii] - 1) for ii in d.sets[:i] if xp0[ii] > 1e-9) < 1e-5
        @test all(isapprox(ue(m, ll), 0.05; atol=1e-8) for ll in L)
        @test all(abs(floorgap(m, ll)) < 1e-8 for ll in L)

        # 0 % start = :full_employment (20 % tariff cut)
        mf, df = solve_closure(:full_employment; tm=0.8)
        mw, dw = solve_closure(:wage_floor; u=0.0, tm=0.8)
        @test isapprox(value(mw[:RGDP]["R1"]), value(mf[:RGDP]["R1"]); rtol=1e-8)
        @test maximum(abs(value(mw[:XP][ii]) / value(mf[:XP][ii]) - 1) for ii in dw.sets[:i] if value(mf[:XP][ii]) > 1e-6) < 1e-7
        @test all(abs(floorgap(mw, ll) - floorgap(mf, ll)) < 1e-7 for ll in L)

        # the regime switch: 5 % start, 20 % tariff cut -> unemployment falls, wage at the floor;
        # 2 % start, zero tariffs -> unskilled unemployment exhausted and its wage above the floor,
        # skilled still unemployed at the floor
        m5, _ = solve_closure(:wage_floor; u=0.05, tm=0.8)
        @test all(1e-4 < ue(m5, ll) < 0.05 for ll in L)
        @test all(abs(floorgap(m5, ll)) < 1e-8 for ll in L)
        m2, _ = solve_closure(:wage_floor; u=0.02, tm=0.0)
        @test ue(m2, "UnSkLab") < 1e-7 && floorgap(m2, "UnSkLab") > 1e-4
        @test ue(m2, "SkLab") > 1e-4 && abs(floorgap(m2, "SkLab")) < 1e-8
        println("wage_floor: 5 % start, -20 % tariffs: UE ", round.([ue(m5, ll) for ll in L]; digits=4),
                "; 2 % start, zero tariffs: UE ", round.([ue(m2, ll) for ll in L]; digits=4),
                ", TW/WMIN - 1 ", round.([floorgap(m2, ll) for ll in L]; sigdigits=3))
    end
end

# ── Negative government revenue: households that consume more than they earn ─────────────────
# Calibration.jl convention (6).  The 12-sector synthetic economy plus a remittance R from
# abroad (an HH x ROW transfer) that households spend on imports, sized so that the calibrated
# revenue YG0 = Tother + kappa_h·YH0 is -0.5 x investment: kappa_h < 0 is a net transfer larger
# than every other tax, the case of Lebanon, Syria, Kyrgyzstan, ... (2023 databases).  The
# transfer is not read on its own; it is part of the trade deficit and of household spending,
# and the benchmark must still reproduce the SAM, with C-9 closing.  The solves fail on the
# YG >= 0 bound that preceded 2026-10-06: C-3 then held only at YG = 1e-8, C-9 missed by
# 0.5 x investment under :bop, and under inv_closure = :savings investment came out 50 % high.
@testset "negative government revenue (transfer-financed consumption)" begin
    function remittance_economy(; kw...)
        d = init_data()
        d.sets[:i]  = ["pdr","wht","gro","ctl","oap","coa","oil","ely","chm","tex","trd","osg"]
        d.sets[:cr] = ["pdr","wht","gro"]; d.sets[:lv] = ["ctl","oap"]; d.sets[:e] = ["coa","oil","ely"]
        d.sets[:ft] = ["chm"]; d.sets[:fd] = ["wht","gro"]; d.sets[:r] = ["R1"]
        prepare_data!(d; outdir=nothing, calibrate=false, precompute=false, kw...)
        M = d.balanced_sam; ix = d.sam_index; hh = ix["HH"]; row = ix["ROW"]
        coms = [ix["COM_"*p] for p in d.sets[:i]]
        calibrate_from_sam!(d)
        R = d.par[:bench][:YG] + 0.5 * d.par[:bench][:INV]    # YG0 falls one for one with C
        w = [M[row, c] for c in coms]; w ./= sum(w)
        Δ = zeros(size(M))
        Δ[hh, row] += R                                         # the transfer
        for (c, wc) in zip(coms, w)
            Δ[c, hh] += R * wc; Δ[row, c] += R * wc             # spent on imported consumption
        end
        M .+= Δ; d.sam === M || (d.sam .+= Δ)
        calibrate_from_sam!(d); d.metadata[:PAR] = precompute_parameters(d)
        return d, R
    end
    c9(m) = (I = value(m[:PFD]["Inv"]) * value(m[:FD]["Inv"]);
             (I - (sum(value(m[:SAV][h]) + value(m[:DeprY][h]) for h in ("HH",)) + value(m[:Sg]) + value(m[:Sf]["R1"]))) / I)

    d, R = remittance_economy()
    M = d.balanced_sam; ix = d.sam_index; B = d.par[:bench]; iset = d.sets[:i]
    @test R > 0 && isapprox(M[ix["HH"], ix["ROW"]], R; rtol=1e-12)
    @test maximum(abs.(vec(sum(M; dims=2)) .- vec(sum(M; dims=1)))) < 1e-8 * sum(M)   # still balanced
    @test isapprox(B[:YG], -0.5 * B[:INV]; rtol=1e-8) && d.par[:kappa_h]["HH"] < 0
    # the SAM's household consumption and its trade deficit (which includes the transfer)
    @test isapprox(B[:HH], sum(M[ix["COM_"*p], ix["HH"]] for p in iset); rtol=1e-12)
    deficit = sum(M[ix["ROW"], ix["COM_"*p]] - M[ix["COM_"*p], ix["ROW"]] for p in iset) - M[ix["TAX_EXP"], ix["ROW"]]
    @test isapprox(B[:Sf], deficit; rtol=1e-8)
    # C-9 closes in the calibrated table (to the 1e-6 floor on SAV0)
    @test abs(B[:INV] - (B[:SAV] + B[:DeprY] + B[:Sg] + B[:Sf])) < 1e-5
    m = model(d; show_solver_output=false)
    @test !has_lower_bound(m[:YG]) && start_value(m[:YG]) == B[:YG]
    @test num_variables(m) == num_constraints(m; count_variable_in_set_constraints=false)

    if get(ENV, "LCGE_TEST_SOLVE", "false") == "true"
        xp0 = B[:XP]
        solve_model!(m; output="no", show_diagnostics=false)
        @test termination_status(m) == JuMP.LOCALLY_SOLVED
        @test maximum(abs(value(m[:XP][p]) / xp0[p] - 1) for p in iset) < 1e-4
        @test isapprox(value(m[:YG]), B[:YG]; rtol=1e-6)
        @test abs(c9(m)) < 1e-6
        # a 20 % tariff cut: C-9 still implied (Walras' law) to the ~1e-5 a shocked solve leaves on the
        # unmodified economy too (-8e-6 there; the bounded YG left -0.5 here), and revenue still < 0
        PAR = parameters(d); tm0 = copy(PAR[:tau_m])
        for (k, v) in tm0; PAR[:tau_m][k] = 0.8 * v; end
        mt = model(d; show_solver_output=false)
        solve_model!(mt; output="no", show_diagnostics=false)
        @test termination_status(mt) == JuMP.LOCALLY_SOLVED
        @test abs(c9(mt)) < 1e-4 && value(mt[:YG]) < 0
        for (k, v) in tm0; PAR[:tau_m][k] = v; end
        # savings-driven investment (C-9 imposed): the benchmark investment replicates
        ds, _ = remittance_economy(; bop_closure=:fixed_er, inv_closure=:savings)
        ms = model(ds; show_solver_output=false)
        solve_model!(ms; output="no", show_diagnostics=false)
        @test termination_status(ms) == JuMP.LOCALLY_SOLVED
        @test isapprox(value(ms[:FD]["Inv"]), ds.par[:bench][:INV]; rtol=1e-6)
        println("negative revenue: YG0/I = ", round(B[:YG] / B[:INV]; digits=3), ", kappa_h = ",
                round(d.par[:kappa_h]["HH"]; digits=3), "; C-9/I benchmark ", round(c9(m); sigdigits=2),
                ", -20 % tariffs ", round(c9(mt); sigdigits=2))
    end
end

# ── An intermediate-input tax on a sector without intermediate inputs ────────────────────────
# Four 133-sector hybrids (Kyrgyzstan, Laos, Nepal, Pakistan 2023) book TAX_INT on activities
# whose COM x ACT column is empty; tau_Ap = txi/EPS (~1e11) left C-3 ~1e6 off at the start point
# and none of their benchmarks solved.  The tax is now an output tax of that sector.  Built from
# the 12-sector synthetic economy: sector "trd" buys its intermediates no more (their cost goes
# to capital, the goods to household consumption), but keeps its input tax.
@testset "input tax without intermediate inputs" begin
    d = init_data()
    d.sets[:i]  = ["pdr","wht","gro","ctl","oap","coa","oil","ely","chm","tex","trd","osg"]
    d.sets[:cr] = ["pdr","wht","gro"]; d.sets[:lv] = ["ctl","oap"]; d.sets[:e] = ["coa","oil","ely"]
    d.sets[:ft] = ["chm"]; d.sets[:fd] = ["wht","gro"]; d.sets[:r] = ["R1"]
    prepare_data!(d; outdir=nothing, calibrate=false, precompute=false)
    M = d.balanced_sam; ix = d.sam_index; a = ix["ACT_trd"]; iset = d.sets[:i]
    Δ = zeros(size(M))
    for p in iset
        c = ix["COM_"*p]; x = M[c, a]
        Δ[c, a] -= x; Δ[ix["CAP"], a] += x; Δ[ix["HH"], ix["CAP"]] += x; Δ[c, ix["HH"]] += x
    end
    M .+= Δ; d.sam === M || (d.sam .+= Δ)
    txi = M[ix["TAX_INT"], a]; txo = M[ix["TAX_OUT"], a]
    @test txi > 0 && sum(M[ix["COM_"*p], a] for p in iset) == 0
    calibrate_from_sam!(d); d.metadata[:PAR] = precompute_parameters(d)
    P = parameters(d); X = P[:bench][:XP]["trd"]
    @test all(P[:tau_Ap][(p, "trd")] == 0 for p in iset)
    @test isapprox(P[:tau_p]["trd"] * X / (1 + P[:tau_p]["trd"]), txo + txi; rtol=1e-12)
    if get(ENV, "LCGE_TEST_SOLVE", "false") == "true"
        m = model(d; show_solver_output=false)
        solve_model!(m; output="no", show_diagnostics=false)
        @test termination_status(m) == JuMP.LOCALLY_SOLVED
        @test maximum(abs(value(m[:XP][p]) / P[:bench][:XP][p] - 1) for p in iset) < 1e-4
        @test isapprox(value(m[:YG]), P[:bench][:YG]; rtol=1e-6)
    end
end
