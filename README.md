# LCGE-V4 — LINKAGE-style economy-wide policy simulation model

LCGE-V4 is a Julia implementation of a LINKAGE-style computable general
equilibrium (CGE) model: N sectors (100 by default), 4 regions, two labour
skills with rural–urban migration, old/new capital vintages, nested Armington/CET trade
with tariff-rate quotas, and a recursive-dynamic extension. It simulates how a
policy change (a tariff, a productivity shock, labour growth, …) ripples
through prices, production, employment, trade and household income.

The model is written as a square **mixed complementarity problem (MCP)** in
[JuMP](https://jump.dev) and solved with the **PATH** solver.

> **Sibling branch.** The `envcge-v9` branch of this repository holds an
> unrelated model (EnvCGE V9, an ENVISAGE-style climate CGE). The two share no
> code or history; do not merge them.

---

## Current status (2026-09-03)

- The default synthetic economy prepares, balances and builds correctly:
  216-account SAM (2N+16 for N = 100 sectors), 48,099 variables = 48,099
  complementarity constraints.
- **The benchmark replicates**: every equation holds at the calibrated start
  values (max residual 9e-6), PATH reports `LOCALLY_SOLVED` after one major
  iteration (≈4 s) and all variables stay within 3e-4 % of their start
  values.
- A 20 % import tariff on one bilateral flow (`examples` / README below)
  solves in ≈1 s with textbook incidence: the taxed import flow falls 1.9 %,
  its tariff-inclusive price rises 3.8 %, the exporter's FOB price falls
  3.7 %, competing sources are barely affected, real GDP falls by 1e-4 %.
- Recursive dynamics keep an explicit capital stock (`Kstock0 = I0/δ`,
  `K_{t+1} = (1−δ)K_t + I_t`, rental supply `KSupply = κ·K`). With zero
  growth every period reproduces the benchmark (max change 8e-8); with TFP
  growth alone each period solves. **Runs with labour growth still fail in
  most periods** because of the labour closure — see "Known limitations".
- Always check `termination_status(m)` after `solve_model!`. Now that the
  benchmark converges, an `ITERATION_LIMIT` after a shock means the shock is
  too large for a single step — apply it in smaller increments, re-solving
  from the previous solution.

---

## Requirements

| Requirement | Notes |
|---|---|
| Julia ≥ 1.10 | `Manifest.toml` was resolved with 1.11; tested on 1.12.4 |
| PATH licence | A public courtesy licence is built in (valid to 31 Dec 2035). Set `PATH_LICENSE_STRING` to override it. |

Julia packages (installed automatically by `Pkg.instantiate()`):

| Package | Purpose |
|---|---|
| `JuMP`, `PATHSolver`, `Complementarity` | MCP formulation and PATH solver interface |
| `DataFrames` | Results tables |
| `XLSX` | Excel SAM input, dynamic/scenario workbooks |
| `Plots` | Optional charts of results and trajectories — set `LCGE_NO_PLOTS=1` before loading `src/LinkageModel.jl` to skip `Plotting.jl` entirely (the `plot_*` functions are then undefined; nothing else changes) |

---

## Installation

```bash
git clone https://github.com/SebKrantz/SDCGE
cd SDCGE
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The first `instantiate` precompiles JuMP and Plots and can take 5–15 minutes.

---

## Quick start

```julia
include("src/LinkageModel.jl")
using .LinkageModel
using JuMP

data = init_data()        # empty data container
prepare_data!(data)       # built-in synthetic 100-sector SAM → balance → calibrate
m = model(data)           # build the JuMP/PATH model (~15 s)
solve_model!(m)           # solve with PATH (~4 s at the benchmark)
termination_status(m)     # LOCALLY_SOLVED
export_results!(m, data)  # write results/*.csv
```

One-call equivalent (results are only written if you ask for them):

```julia
m, data = run_linkage!(write_results=true)
```

Rough timings on a laptop: `include` 7 s, `prepare_data!` 6 s, `model` 14 s,
`solve_model!` 1–4 s.

---

## Output

`export_results!(m, data)` writes to `results/` (git-ignored):

| File | Contents |
|---|---|
| `results_all_variables.csv` | Every variable element: value, start value, change, % change |
| `results_summary.csv` | Named macro indicators (GDP, CPI, household income, government revenue `YG`, investment, savings, …) |
| `results_scalars.csv` | Scalar variables |
| `results_<VAR>.csv` | One file per major variable container, e.g. `results_XP.csv`, `results_GDP.csv`, `results_YH.csv`, `results_YG.csv` |
| `results_metadata.csv` | Solver termination/primal status and timestamp |
| `balanced_sam.csv`, `sam_balance_table.csv`, `sam_balance_summary.csv` | The RAS-balanced SAM and its row/column gaps |

In Julia, `results_dataframe(m)` returns the same table as
`results_all_variables.csv`. `pct_change_from_start` is the change relative
to the calibrated benchmark (≈ 0 for a benchmark run).

---

## Using your own SAM

The economy is defined by a square Social Accounting Matrix with **2N + 16
accounts** for N sectors, in the standard order (N activities, N commodities,
5 factors, 6 taxes, 4 institutions, 1 margin account) — 216 accounts for the
default N = 100; see `data/csv/sam_accounts.csv`. The SAM is RAS-balanced
automatically and the balance report is written next to the results.

```julia
# CSV: first row / first column hold the account labels
data = prepare_data!(init_data(); source=:csv, sam_path="data/csv/sam.csv")

# Excel: sheet "SAM" of the workbook
data = prepare_data!(init_data(); source=:excel, sam_path="data/linkage_100sector_data.xlsx")

m = model(data)
```

### Using your own sectors (N ≠ 100)

Pass `sets_path` to read the sector list and its groupings from a two-column
`set,item` CSV (`data/csv/sets.csv` is the shipped 100-sector file):

```julia
data = prepare_data!(init_data(); sets_path="data/csv/sets.csv",
                     source=:csv, sam_path="data/csv/sam.csv")
```

The file must define `i` (activities = products). Unless N = 100 it must also
define the memberships `cr` (crops), `lv` (livestock), `e` (energy), `ft`
(fertiliser) and `fd` (feed) — the built-in positional defaults (crops = first
10, energy = 71:75, …) only make sense for `P001`–`P100`, and `default_sets!`
now errors instead of silently applying them. Constraints:

- every item of a sector set must be a member of `i`;
- `ft ∩ e = ∅` and `fd ∩ e = ∅` (the `XAp` columns are partitioned into
  fertiliser/feed, energy and "other" blocks; an overlap breaks squareness);
- `|e| ≥ 1`, and `|ft| ≥ 1` whenever `|cr| ≥ 1`, `|fd| ≥ 1` whenever
  `|lv| ≥ 1`, `|ag| ≥ 1`;
- `ag = cr ∪ lv`, `ip = i ∖ ag`, `nf = i ∖ ag` are derived when the file omits
  them, and `nnft = i ∖ ft`, `nnfd = i ∖ fd` are always recomputed.

Sets the file does not mention (`r`, `v`, `l`, `h`, `f`, `t`, …) keep their
defaults, so the file can also shrink the pseudo-regions to a single `R1`
(`r,R1`) — the SAM has no regional dimension, the four regions only spread the
national trade totals over 16 uniform `(r, rp)` cells, and `|r| = 1` builds
square, replicates the benchmark and cuts the bilateral blocks by `|r|²`.

The SAM must then carry exactly the 2N + 16 accounts derived from `i`
(`ACT_<code>`, `COM_<code>`, …); `read_sam_csv!`/`read_sam_excel!` check the
labels and name the mismatches.

`examples/02_read_csv_sam.jl` and `03_read_excel_sam.jl` show the explicit
step-by-step version (`read_sam_csv!` → `validate_sam!` → `balance_sam_ras!`
→ `calibrate_from_sam!`). Read "Calibration conventions" below before using
a real SAM: the model cannot represent a trade deficit, and the calibration
adjusts final demand to close it.

---

## Applying a policy shock

All calibrated parameters live in `parameters(data)` (a `Dict{Symbol,Any}`,
stored in `data.metadata[:PAR]`). Modify them after `prepare_data!` and before
`model`:

```julia
data = prepare_data!(init_data())
PAR = parameters(data)
PAR[:tau_m][("R1", "R2", "P001")] = 0.20   # import tariff, index (r, rp, product)
m = model(data)
solve_model!(m)
```

Commonly used entries:

| Key | Description | Index |
|---|---|---|
| `:tau_m`, `:tau_e` | Import tariff / export tax rate | `(r, rp, product)` |
| `:tau_p` | Output tax rate | `product` |
| `:AT` | Total factor productivity | `product` |
| `:LSupply` | Labour supply | `"UnSkLab"` or `"SkLab"` |
| `:KSupply` | Capital supply (rental units) | `(product, vintage)` |

---

## Recursive dynamics and policy experiments

```julia
# 10 periods; each period is one static solve, then K, L and A are updated
data, history, snapshots = run_recursive_dynamic!(periods=10, delta=0.05,
                                                  g_labor=0.02, g_tfp=0.015,
                                                  outdir="results/dynamic")
plot_dynamic_results("results/dynamic/dynamic_results.xlsx")

# Excel-driven batch of scenarios (one sheet per exogenous block)
write_policy_template("data/policy_experiments.xlsx"; periods=10, n_scenarios=10)
results = run_policy_experiments!("data/policy_experiments.xlsx";
                                  outdir="results/scenarios", make_plots=true)
```

See `examples/07_recursive_dynamics.jl` and `08_policy_experiments.jl`.

### One scenario, in memory

`run_scenario!` runs a single `Scenario` without the Excel round trip and without
writing anything, and gives an embedding caller (a web backend, a batch driver) the
three hooks it needs: a per-period callback, a cooperative abort flag, and a
`period_modifier` for the levers the `Scenario` struct does not carry — tariffs
(`:tau_m`), the output tax (`:tau_p`), the direct tax (`:kappa_h`), the government
spending share (`:chi_gov`).

```julia
data = init_data()
prepare_data!(data; source=:csv, sam_path="data/csv/sam.csv",
              outdir=nothing)               # nothing ⇒ write no SAM report

scen = Scenario(1, "tariff_cut", "", 10, 0.05,
                Dict((i, t) => 1.0 for i in data.sets[:i], t in 1:10),   # AT levels
                Dict((l, t) => 0.0 for l in data.sets[:l], t in 1:10),   # g_labor rates
                zeros(10), zeros(10))                                    # g_land, g_nres

stop = Ref(false)
history, snapshots = run_scenario!(data, scen;
    on_period = (t, status, secs, h) -> println("period $t: $status in $secs s"),
    abort = stop,
    period_modifier = (t, d) -> (parameters(d)[:tau_m][("R1","R2","P001")] = 0.0))
```

`history` is one named tuple of macro indicators per period; `snapshots` is every
solved variable of every period, keyed `(variable, index_labels)`.

---

## Examples and tests

`examples/` contains nine numbered scripts from data preparation to the full
dynamic/policy pipeline; see `examples/README.md` for the list and runtimes.

```bash
julia examples/04_build_model.jl          # build only, prints variable/constraint counts
julia --project=. test/runtests.jl        # smoke test: SAM, CSV import, build, results labels (~30 s)
LCGE_TEST_SOLVE=true julia --project=. test/runtests.jl   # additionally solves the benchmark (LOCALLY_SOLVED)
```

---

## Diagnostics

```julia
print_model_diagnostics(m)   # variable/constraint counts and solver
diagnose_model(m)            # equation ↔ variable matching and residual_at_start per equation
```

`diagnose_model` writes `results/diagnostics/*.csv`; the `residual_at_start`
column of the equation table is the quantity to watch when changing
calibration or initialisation — it must be ≈ 0 for every equation at the
benchmark.

---

## Calibration conventions

`calibrate_from_sam!` (src/Calibration.jl) derives every share parameter and
the complete benchmark start point (`parameters(data)[:bench]`) from the
balanced SAM, so that each equation holds exactly at the start values. The
conventions, and the three places where the SAM cannot supply what the
equations need, are documented in the header of `Calibration.jl`:

- Benchmark prices are 1 (net-of-tax producer price `PX = 1/(1+tau_p)`);
  CES shares `α_j = s_j (P/P_j)^(1−σ)`, CET shares `β_j = s_j (P/P_j)^(1+σ)`
  from SAM value shares; all elasticities default to 0.5 (no elasticity data
  in the SAM); technical-change indices `λ = 1`.
- **Balanced trade.** The trade block (E-2 with T-21) makes the CIF value of
  imports of each good identically equal its FOB export value, and there is no
  balance-of-payments equation, so a trade deficit cannot be represented. The
  synthetic SAM has a 7,439 deficit at border prices. Exports are kept at
  their SAM values, imports are set to `(1+tau_m)(1+tau_e)` × exports, and
  household, government and investment demand are scaled down (≈12 %) so that
  absorption equals output minus exports plus imports. Intermediate demand,
  production and factor payments keep their SAM values. Consequently the
  direct-tax rate is solved so that investment is financed (`kappa_h ≈ 0.24`),
  household saving is ≈ 0 and government saving funds investment.
- **Land is agricultural only** (the factor equations force zero land outside
  `S[:ag]`); land payments the synthetic SAM assigns to other sectors are
  reassigned to capital.
- **Subsistence quantities `theta = 0`** (LES collapses to proportional
  budget shares); trade margins are zero at the benchmark (`zeta_t = 0`).

---

## Known limitations

- **Labour closure.** The default regime (`parameters(data)[:labour_closure]
  = :fixed_wage`) fixes the wage `W` at 1 and lets unemployment `UE` absorb
  any gap between labour supply and demand (F-10); F-4, F-6, F-7 and F-11 are
  circular/degenerate in this regime, so once `UE` leaves its bound PATH
  loses its footing. Consequences: `g_labor` in the dynamic runs raises
  unemployment instead of output, and 7 of 10 periods of the default
  10-period run end in `ITERATION_LIMIT` (non-converged periods do not update
  the state, so the path plateaus; a warning is printed). TFP growth alone
  solves every period.
  A market-clearing alternative is implemented and selectable
  (`PAR[:labour_closure] = :full_employment`, with `PAR[:numeraire] ∈
  {:pabs (default), :cpi}`): F-6 becomes labour-market clearing ⟂ `TW`,
  `UE` is fixed at its benchmark, `NW = φ·TW`, `W = (1+τ_l)·NW`, F-21 becomes
  capital-market clearing ⟂ `TR` (in the default regime `TR` has no
  determining equation — F-21 collapses to `TR = TR`, hidden because the
  benchmark start is already the solution). It is square and replicates the
  benchmark to 0.004 %, **but PATH stalls (`SLOW_PROGRESS`)** because the
  vintage block is rigid: F-24 pins `RR = 1`, so old-vintage output cannot
  adjust to a wage change. That is the next thing to fix before switching
  the default; see the `Factors.jl` header.
- **Vintages carry no technology**: `Calibration.jl` gives Old and New capital
  identical shares and prices, so the dynamic update keeps the benchmark
  Old/New split (`vintage_rule=:benchmark_shares`); the flow-based split
  (`:flow`) is available but shifts output between vintages without economic
  content.
- `AT` (productivity) enters P-1/P-2 as `AT·…` but P-3 as `…/AT`; the CES
  demand form should carry `AT^(σ−1)`. Harmless at the benchmark (`AT = 1`),
  but check before relying on TFP shocks quantitatively.
- Zero-valued variables sit at a 1e-8 safety bound, which leaves residual
  floors of ~1e-8 on ~180 equation families; harmless for PATH.

---

## Troubleshooting

**`termination_status(m)` is `ITERATION_LIMIT` after a shock** — the shock is
too large for one Newton path from the benchmark. Apply it in steps (e.g.
tariff 0.11 → 0.15 → 0.20), re-solving each time; the previous solution is
kept as the start point when you modify `PAR` and rebuild.

**PATH licence error** — set `PATH_LICENSE_STRING` before starting Julia
(see [PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl)); by default the
built-in courtesy licence is used.

**`SAM is not balanced`** — RAS balancing is applied automatically; large
residual gaps in `results/sam_balance_table.csv` mean the input SAM is
inconsistent.

**`SAM file ... does not carry the expected accounts`** — the file's account
labels must be exactly the 2N + 16 the model derives from `data.sets[:i]`
(the 216 in `data/csv/sam_accounts.csv` for the default sectors). The error
lists what is missing and what is unexpected.

**`default_sets!: with |S[:i]| = N != 100 ...`** — a non-default sector list
needs the memberships `cr`, `lv`, `e`, `ft`, `fd`; supply them via
`sets_path` (see "Using your own sectors").

---

## Project layout

```
SDCGE/  (branch main)
├── src/
│   ├── LinkageModel.jl     module entry point; include order and exports
│   ├── Types.jl            LinkageData container, default_sets!
│   ├── SAM.jl              SAM accounts, synthetic SAM, CSV/Excel readers, RAS balancing
│   ├── Calibration.jl      calibrate_from_sam!: shares and benchmark start point from the SAM
│   ├── ParameterTables.jl  precompute_parameters → PAR dictionary (defaults + calibrated tables)
│   ├── Functions.jl        CES / CET / Armington helper functions
│   ├── Initialization.jl   start values and bounds for all variables (from PAR[:bench])
│   ├── Variables.jl        @variables declarations
│   ├── Production.jl, Income.jl, Demand.jl, Trade.jl, Equilibrium.jl,
│   │   Closure.jl, Factors.jl, Other.jl   paper-numbered equation blocks (P-, Y-, D-, T-, E-, C-, F-)
│   ├── ModelBuilder.jl     prepare_data!, model/build_model, solve_model!, run_linkage!
│   ├── Results.jl          results_dataframe, export_results!
│   ├── Plotting.jl         plot_results, plot_dynamic_results, plot_all_scenarios
│   ├── Diagnostics.jl      diagnose_model, print_equation_diagnostics
│   ├── RecursiveDynamic.jl run_recursive_dynamic!
│   └── PolicyScenarios.jl  write_policy_template, run_policy_experiments!
├── data/                   synthetic SAM (CSV + Excel), policy_experiments.xlsx template
├── examples/               numbered walkthrough scripts (start here)
├── test/runtests.jl        smoke test
└── results/                generated output (git-ignored)
```

---

## What the model covers

| Feature | Detail |
|---|---|
| Sectors | N, from `data.sets[:i]`; 100 by default (crops P001–P010, livestock P011–P020, energy P071–P075, fertiliser P076–P078, other industry and services). Any other list comes from `sets_path` |
| Regions | 4 (`R1`–`R4`) with bilateral trade; `|r| = 1` is supported |
| Labour | Unskilled and skilled, rural/urban zones, migration |
| Capital | Old (installed) and new (investment) vintages |
| Trade | Nested Armington import demand, CET export supply, tariff-rate quotas |
| Government | Output, intermediate, trade, factor and income taxes; fiscal balance |
| Households | One representative household; ELES/AIDADS-style demand system |
| Dynamics | Recursive: capital accumulation, labour growth, productivity growth |
