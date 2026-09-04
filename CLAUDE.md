# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

LCGE-V4: a LINKAGE-style CGE model (N sectors, 100 by default, × 4 regions × 2 labour skills × 2 capital vintages) written as a square mixed complementarity problem in JuMP and solved with PATH. Module `LinkageModel` in `src/LinkageModel.jl`. There is no package name/UUID in `Project.toml`, so the module is always loaded from source:

```julia
include("src/LinkageModel.jl"); using .LinkageModel   # never `using LinkageModel`
```

The `envcge-v9` branch of this repo is an unrelated model (EnvCGE) with its own CLAUDE.md. Never merge the branches.

## Commands

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'        # once; 5–15 min (JuMP, Plots precompile)
julia --project=. test/runtests.jl                          # smoke test, ~30 s, no PATH call
LCGE_TEST_SOLVE=true julia --project=. test/runtests.jl     # + benchmark solve, must be LOCALLY_SOLVED
julia examples/04_build_model.jl                            # build only, prints counts
julia examples/05_solve_and_export.jl                       # build + solve + results/ + plots (~1 min)
```

Timings: `include` ≈7 s, `prepare_data!` ≈6 s, `model()` ≈14 s, `solve_model!` ≈1–4 s. Do not run `Pkg.resolve()`/`Pkg.update()` casually — the committed `Manifest.toml` (resolved on Julia 1.11, runs on 1.12) is the known-good environment.

## Pipeline and architecture

`prepare_data!` → `model` → `solve_model!` → `export_results!` (`run_linkage!` does all four; note its defaults `write_results=false, make_plots=false`).

- `Types.jl` — `LinkageData` (sets, SAM, calibrated tables, `metadata`). `default_sets!` fills in whatever is missing (every assignment is a `get!`, so pre-populated sets win): 100 products `P001`–`P100` (also used as activities), regions `R1`–`R4`, vintages `Old`/`New`, labour `UnSkLab`/`SkLab`, one household `HH`, final-demand agents `Gov`/`Inv`. Energy = P071–P075, fertiliser = P076–P078, crops P001–P010, livestock P011–P020, `S[:ag]` = crops ∪ livestock. Those positional groups only apply to `P001`–`P100`: with `|S[:i]| ≠ 100` the sets `:cr :lv :e :ft :fd` must be supplied or `default_sets!` errors.
- `SAM.jl` — 2N+16-account SAM (216 at N=100; not 406, whatever old text says), synthetic default (`build_default_large_sam!`), CSV/Excel readers, RAS balancing, and `read_sets_csv!` (two-column `set,item`, e.g. `data/csv/sets.csv`) which `prepare_data!(...; sets_path=...)` calls before `default_sets!` — the only way to run N ≠ 100. `ft`/`fd` must be disjoint from `e` or the model stops being square; `|r| = 1` is supported. The SAM has **no regional dimension**; regional/bilateral structure is synthesised with uniform shares.
- **`Calibration.jl::calibrate_from_sam!` is the single source of truth for the benchmark.** It derives every share/scale parameter from SAM value shares (CES `α_j = s_j (P/P_j)^(1−σ)`, CET `β_j = s_j (P/P_j)^(1+σ)`) *and* the complete start point `PAR[:bench]`, which `Initialization.jl::initialize_from_sam!` reads verbatim (no hard-coded fractions anywhere). The elasticities used in the α formulas are written into `par` so equations and calibration cannot drift. Its header documents the conventions: benchmark prices 1 (`PX = 1/(1+tau_p)`, `PND = 1+tau_Ap`); `tau_p = TAX_OUT/(output−TAX_OUT)`, `tau_Ap = TAX_INT/intermediate`; balanced trade (imports = `(1+tau_m)(1+tau_e)`·exports per good, final demand scaled ≈12 % to close the SAM's deficit, `lambda_w = (1+tau_m)(1+tau_e)`); land ag-only (non-ag land payments → capital); `theta = 0`; `kappa_h` solved so C-9 delivers benchmark investment; all elasticities 0.5; `lambda_* = 1`; `GammaC` diagonal.
- `ParameterTables.jl::precompute_parameters` builds `PAR::Dict{Symbol,Any}` (defaults, overwritten by the calibrated tables; stored as `data.metadata[:PAR]`, accessed via `parameters(data)`). Tables are `Dict` keyed by product or tuples such as `(r, rp, product)` for `:tau_m`, `(product, vintage)` for `:KSupply`. Policy shocks are applied by editing `PAR` between `prepare_data!` and `model`.
- `Variables.jl` declares every JuMP container; `Initialization.jl` sets start values/bounds (`enforce_nlp_safe_bounds_and_starts!` moves zero lower bounds to 1e-8 — except `TauPR` and `UE`, which must be able to reach exactly 0 for T-12/F-10 to hold); the equation files `Production/Income/Demand/Trade/Equilibrium/Closure/Factors/Other.jl` each add one paper-numbered block (`P-`, `Y-`, `D-`, `T-`, `E-`, `C-`, `F-`) as `@constraint(m, F ⟂ x)` pairs. Constraint names follow `<Block>_<n>[indices]`, e.g. `T_5[R1,P075]`.
- CET/CES convention: demand/supply equations are in share form `X_k = β_k (P_k/P)^σ X` with Σβ = 1; primal aggregators use `b_k = β_k^(-1/σ)` (CET, T-16) so they hold at the benchmark. T-7 aggregates over source regions; T-18 is the CET FOC `WTFs = β_z (PE/PET)^σ ES` (the ratio was inverted before 2026-09-03 and made tariff shocks perverse).
- `ModelBuilder.jl` — `build_linkage_model!` calls the blocks in order; `solve_model!` runs `check_initialization!` first and errors on bad start values.
- `Results.jl` — `results_dataframe(m)` (one arg) and `export_results!`; per-variable CSVs are written only for containers that exist (legacy names in `common_variables` such as `GOVREV` silently produce nothing — the real variable is `YG`).
- `RecursiveDynamic.jl` re-solves the static model period by period, updating K/L/A between periods; `PolicyScenarios.jl` drives many trajectories from an Excel workbook (`write_policy_template` → edit → `run_policy_experiments!`) and, for one scenario in memory, `run_scenario!(data, scen; periods, delta, on_period, abort, period_modifier)` — the same loop body with a per-period callback, a cooperative abort `Ref`, and a `period_modifier(t, data)` hook that runs before each period's `model(data)` so levers outside the `Scenario` struct (`:tau_m`, `:tau_p`, `:kappa_h`, `:chi_gov`) can be given a time path. `prepare_data!(...; outdir=nothing)` skips the SAM-report writes; `LCGE_NO_PLOTS=1` skips `Plotting.jl` (and Plots/GR) at load time. `Dynamics.jl` (paper G-equations) was removed; recursive dynamics live only in these two files.
- `Diagnostics.jl` — `diagnose_model` pairs each complementarity constraint with its variable and evaluates `residual_at_start`. JuMP stores `F ⟂ x` as a positional `Vector` of scalar expressions (`func[1]` = F, `func[2]` = x), not MOI vector functions.

## Known state of the model (verified 2026-09-03)

- Build is square: 48,099 variables = 48,099 constraints. **Benchmark replicates**: max `residual_at_start` 8.6e-6 (all in equations summing thousands of 1e-8-bounded zeros), PATH `LOCALLY_SOLVED` in 1 major iteration, ≈4 s, max |pct change| 2.7e-4 %.
- Tariff shock `PAR[:tau_m][("R1","R2","P001")] = 0.20` (from 0.113): `LOCALLY_SOLVED` in ≈1 s; taxed flow −1.9 %, `PM` +3.8 %, exporter FOB −3.7 %, `RGDP[R1]` −1e-4 %.
- Judge any calibration/initialisation change by the residual-by-family table (`diagnostic_equation_matches(m)` → aggregate `residual_at_start` by constraint family) — max must stay ≈1e-5 — and by a benchmark solve with variables unchanged.
- Recursive dynamics (`RecursiveDynamic.jl`, shared by `PolicyScenarios.jl`): explicit capital stock in `data.metadata[:capital_state]` — `Kstock0 = FDInv0/δ` allocated by benchmark rental shares, `κ_i = ΣKSupply0/Kstock0_i` (≈0.0768), `K_{t+1} = (1−δ)K_t + I_t·share_i`, `KSupply = κ·K·vshare`, `K0 = KSupply[Old]`; `vintage_rule=:benchmark_shares` (default; `:flow` breaks zero-growth replication by +3.3 % because vintages carry no technology). Non-converged periods do not update the state (warning). `PAR[:LV0]` tracks solved labour demand as next period's start. Zero growth reproduces the benchmark (test `LCGE_TEST_SOLVE=true`); TFP growth alone solves each period (≈5 s).
- Labour closure switch `PAR[:labour_closure]` (ParameterTables.jl default `:fixed_wage`; `Factors.jl::add_factor_equations!` branches on it; header documents both regimes with labels):
  - `:fixed_wage` (default, byte-identical to the original equations): `W = 1` (F-12), F-10 makes `UE` absorb supply−demand, F-7 `(TW−WMIN)·UE = 0 ⟂ UE` is degenerate with `TW ≡ WMIN ≡ 1`, F-4/F-6/F-11 circular (`AVGW·1e-9 = 0`); `W` and `NW` are unlinked; F-21 collapses to `TR = TR` (F-24 `RR = 1` ⇒ `R[i,Old] ≡ TR`), so `TR` has no determining equation — hidden because the benchmark start is already the solution. Labour supply is not binding: `g_labor` only creates unemployment; default 10-period run fails 7/10 periods.
  - `:full_employment` (experimental): F-6 `Σ LV + N·LF_d = LS·(1−UE0) ⟂ TW[national]`, F-7 `TW[z] = TW[national]`, F-9 `WMIN = χ·PS^ω·PABS^ω`, F-10 `UE = UE0` (=0, from the SAM), F-11 `NW = φ·TW`, F-12 `W = (1+τ_l)·NW`, F-25 `KS = Σ KSupply` exogenous and F-21 capital-market clearing `Σ Kvd + Σ N·KF_d = KS ⟂ TR`; numeraire `PAR[:numeraire]` `:pabs` (`PABS = 1`, default) or `:cpi` (`CPI[R1] = 1 ⟂ PABS`; `PNUM` cannot anchor — its Jacobian column is empty). Square (48,099), start residual 8.6e-6, but PATH stalls (`SLOW_PROGRESS`, ~110 s, drift 0.004 %): every step is cut to ~1e-3 — a near-singular Jacobian in the capital/vintage block (`F_29`, `P_1`, `P_2`, `F_30`, `F_23`), where F-24 `RR = 1` leaves the vintage decomposition no slack. **Next lead: relax F-24 / make old-vintage output adjustable before flipping the default.** Tried and reverted (comments in `Initialization.jl`): warm starts from the previous period (worse), `LCGE_START_EPS = 1e-10` (worsens default drift), dropping the redundant `RR` upper bound (no effect).
- Other unresolved: `AT` enters P-1/P-2 as `AT·…` but P-3 as `…/AT` (CES demand should carry `AT^(σ−1)`; only matters away from `AT = 1`); trade margins are zero at the benchmark (`zeta_t = 0`).
- PATH licence is set in `LinkageModel.jl` from `PATH_LICENSE_STRING` if present, else a built-in courtesy licence.

## Conventions

- `results/` is git-ignored generated output; never commit it.
- Example scripts activate the project themselves (`Pkg.activate(joinpath(@__DIR__, ".."))`) so they run with plain `julia examples/xx.jl` from the repo root.
- Equation comments cite the LINKAGE paper numbering; keep the labels when editing an equation, and record any convention decision in the `Calibration.jl` header.
