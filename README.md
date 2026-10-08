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

## Current status (2026-09-05)

- The default synthetic economy prepares, balances and builds correctly:
  216-account SAM (2N+16 for N = 100 sectors), 48,100 variables = 48,100
  complementarity constraints (the extra one is the real exchange rate `ER`).
- **Any number of sectors** from a `sets.csv` file, and a single region
  (`|r| = 1`), see "Using your own sectors" below.
- **Balance-of-payments trade closure** (`trade_closure = :bop`, default): the
  SAM's own imports, exports and final demand are the benchmark, the trade
  deficit is exogenous foreign saving and the real exchange rate clears the
  current account; `:balanced` keeps the legacy per-good trade balance. See
  "Closures".
- **Real country databases**: 333 SAMs for base year 2023 built by the pipeline in
  `~/Documents/Data/CGE` ([AEILabs/CGE-SAMs](https://github.com/AEILabs/CGE-SAMs)) under
  `data/<ISO3>_<year>_<family>/` (git-ignored, regenerate with its `R/09_export_sdcge.R`;
  `data/registry.csv` lists them): `gtap12` (GTAP 12, 144 countries, 65 sectors), `hybrid`
  (EMERGING 2023 levels × GTAP 12 composite structure, 25 economies GTAP 12 lacks) and
  `hybrid133` (EMERGING's own 89–133 sectors, 164 economies). Every one loads, balances and
  calibrates without rescaling. `gtap12`: PATH solves 142/144, 139 within 1 % of the benchmark,
  2-period run 142 (Mauritius and Hong Kong stall in the land CET: GTAP 12 books no crop land
  there). `hybrid`: 25/25, 23 within 1 %. `hybrid133`: 156/164 solve but only 74 within 1 %
  (median 1.2 %; production-nest residuals of tiny activities). The 2017 GTAP Africa V3 and 2018
  hybrid databases cited in older notes (e.g. `KEN_2017_gtap11afr`) are no longer exported;
  regenerate one with `Rscript R/09_export_sdcge.R gtap11afr 2017 KEN` if needed.
- **Household income, saving, transfers and GDP** (2026-10-07, branch `fix/household-income`; see
  "Calibration conventions", "GDP and gross output" and "What changed for users"). The benchmark
  reproduces the SAM's household and government saving and closes C-9 to ~1e-9 of investment
  (solved, 1e5 and raw scale); transfers from abroad are a lump sum fixed in real domestic terms, so
  `kappa_h` is a tax rate again (it was a net transfer on 38 databases); `GDP`/`RGDP`/`PGDP` are
  GDP at market prices (they were gross output, now `GO`/`RGO`/`PGO`), and C-6 ties government
  demand to real GDP.
- **The savings–investment account closes on every database** (2026-10-06). The 133-sector
  hybrids of Kyrgyzstan, Laos, Nepal and Pakistan, which no closure solved, solve their
  benchmark under all three labour closures, and a −20 % tariff run under all three (Nepal: not
  under `:full_employment`).
- **The benchmark replicates**: every equation holds at the calibrated start
  values, PATH reports `LOCALLY_SOLVED` after one major iteration (≈4 s).
- Recursive dynamics keep an explicit capital stock (`Kstock0 = I0/δ`,
  `K_{t+1} = (1−δ)K_t + I_t`, rental supply `KSupply = κ·K`). With zero
  growth every period reproduces the benchmark; with TFP growth alone each
  period solves. **Runs with labour growth still fail in most periods**
  because of the labour closure — see "Known limitations".
- Always check `termination_status(m)` after `solve_model!`. An
  `ITERATION_LIMIT` after a shock means the shock is too large for a single
  step — apply it in smaller increments, re-solving from the previous solution.

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
- `|e| ≥ 1`; `lv`, `ft` and `fd` may be empty (an economy without livestock
  or fertiliser production — the bundle price then stays at 1), and a set the
  file does not mention is empty;
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
→ `calibrate_from_sam!`). Read "Calibration conventions" and "Closures" below before using
a real SAM: under the default `:bop` closure the SAM's trade deficit becomes
exogenous foreign saving and final demand is used as is; under `:balanced`
the calibration rescales final demand to close it.

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
- **Trade closure** (`prepare_data!(...; trade_closure = :bop | :balanced)`,
  see the Closures section below). Under the default `:bop` the SAM's imports,
  exports and final demand are used verbatim and the trade deficit becomes
  exogenous foreign saving. Under `:balanced` (the legacy convention) the trade
  block (E-2 with T-21) makes the CIF value of imports of each good identically
  equal its FOB export value: exports are kept at their SAM values, imports are
  set to `(1+tau_m)(1+tau_e)` × exports, and household, government and
  investment demand are scaled down (≈12 % for the synthetic SAM, which has a
  7,439 deficit at border prices) so that absorption equals output minus
  exports plus imports. Intermediate demand, production and factor payments
  keep their SAM values in both cases.
- **Land is agricultural only** (the factor equations force zero land outside
  `S[:ag]`); land payments the synthetic SAM assigns to other sectors are
  reassigned to capital.
- **Subsistence quantities `theta = 0`** (LES collapses to proportional
  budget shares); trade margins are zero at the benchmark (`zeta_t = 0`).
- **Saving and taxes** (convention (6) in the header, 2026-10-07). The benchmark reproduces
  the SAM's household saving `S_H` (`INV × HH − HH × INV`) and government saving `S_G`
  (`INV × GOV − GOV × INV`): `SAV0 = S_H − DeprY0` (household gross saving `SAV + DeprY` is
  `S_H`; `DeprY = 0.05·KY`), disposable income `YD0 = C + SAV0`, and the direct-tax rate
  `kappa_h = 1 − YD0/YH0`, which is the households' final-demand tax (`TAX_OUT × HH`; the
  exported SAMs carry no other direct tax) over their income: 0.4 % on Kenya `gtap12`, 5.5 % on
  Lebanon. The final-demand taxes of government and investment (`TAX_OUT × GOV/INV`) are read
  as `tau_Af`, so `PFD = 1 + tau_Af` at the benchmark. With zero subsistence, `1 − Σ mu_c` is
  the household saving rate. C-9 then holds identically at the benchmark (to ~1e-14 of
  investment). Until 2026-10-07 household saving was the macro residual (≈ 0) and `kappa_h`
  took all household non-consumption, negative on 38 of the 333 databases.
- **Transfers from abroad** (convention (7)). Net current transfers to households are a lump
  sum fixed in real domestic terms, `PAR[:WTRbar]`, valued at `PNUM·PABS` (Y-5) and paid
  through the balance of payments at that value (`C_BOP`). It is neutral to the exchange rate;
  under `bop_closure = :fixed_er` it is identical to a foreign-currency lump sum to first
  order. (Valued at the exchange rate, the imputed gap, 54 % of household income on Lebanon,
  made a real appreciation cut household income one for one: +20 % foreign saving under
  `:fixed_wage`/`:flex_er` gave Lebanon −59 % real GDP.) They come from the SAM's `HH × ROW − ROW × HH`. The
  exported SAMs carry none: GTAP books remittances, aid and foreign borrowing in
  `S − I = X − M`. So where households dissave (`S_H < DeprY0`, 40 databases), the gap is
  booked as the transfer and `SAV0 = 0`, with foreign saving lower by the same amount; foreign
  saving `Sf` is then the current-account deficit. A `transfers.csv` next to `sam.csv`
  (`institution,value`, `HH`, SAM units; or `prepare_data!(...; transfers_path=)`) sets the
  amount explicitly; it is folded into the `HH × ROW` cell, and any dissaving it leaves stays
  in `SAV0`. Under `:balanced` (no exchange rate) no transfer is modelled.
- **An input tax on a sector that buys no intermediate inputs** (`TAX_INT` × `ACT` with an
  empty `COM` × `ACT` column) is booked as that sector's output tax: as a rate on inputs it
  would be `txi/1e-9`. Four 133-sector hybrids (Kyrgyzstan, Laos, Nepal, Pakistan 2023) carry
  one; their benchmarks failed under every closure until 2026-10-06.

## GDP and gross output

Since 2026-10-07 `GDP`, `RGDP` and `PGDP` are GDP (`Other.jl` M-1..M-3), on the expenditure
side at market prices:

    GDP  = Σ PAc·XAc + Σ_f PFD·FD + Σ WPE·WTFs − Σ WPM·WTFd            (M-1, nominal)
    RGDP = Σ PAc0·XAc + Σ_f PFD0·FD + Σ WPE0·WTFs − Σ WPM0·WTFd        (M-2, benchmark prices)
    PGDP = GDP / RGDP                                                  (M-3)

i.e. household purchases, government and investment demand, FOB exports less CIF imports; the
benchmark prices `PAc0 = 1 + tau_Ac`, `PFD0 = Σ a_f·(1 + tau_Af)`, `WPE0`, `WPM0` are calibrated
into `PAR`, so a shock to a tax rate does not move them. This is the measure the EPS app has
computed itself since its P29. Until then the three were gross output; they are now `GO = Σ PP·XP`,
`RGO = Σ XP` and `PGO = GO/RGO` (in `results_dataframe`, `export_results!` and the period
summaries as `GO_R1`, `RGO_R1`, `PGO_R1`). C-6 sets real government demand
`FD[Gov] = chi_gov · RGDP`, so `chi_gov` is government demand's share of real GDP (Kenya
`gtap12` 0.124; it was 0.066 of gross output), and `chi_inv`, which `update_period_data!` uses
to re-base investment between periods, is investment's. C-10 `InvSh` is investment's share of
nominal GDP. Every region of a multi-region build carries the national value.

### What changed for users (2026-10-07)

- `GDP_R1`, `RGDP_R1`, `PGDP_R1` (and the `GDP`, `RGDP`, `PGDP` variable families) are GDP now;
  the old gross-output numbers are `GO_R1`, `RGO_R1`, `PGO_R1` (`GO`, `RGO`, `PGO`).
- `PAR[:chi_gov]` is government demand's share of real GDP, about twice its old value (Kenya
  `gtap12` 0.124, was 0.066); a lever that sets it to an absolute level must be re-based.
  `PAR[:chi_inv]` likewise; between periods investment now follows real GDP. `InvSh` is
  investment over nominal GDP.
- `PAR[:kappa_h]` is a tax rate (the households' final-demand tax over their income, ≥ 0; Kenya
  0.4 %, Lebanon 5.5 %); it used to be the residual of household non-consumption, a net transfer
  on 38 databases.
- `SAV` is the SAM's household saving net of depreciation (a free variable), `Sg` the SAM's
  government saving (Kenya −7,416 at the 1e5 scale, was 1,065), `YG` includes the final-demand
  taxes of government and investment (`tau_Af`; `PFD = 1 + tau_Af` at the benchmark) and no
  longer the overstated export tax.
- `Sf` and `PAR[:Sfbar]` are the current-account deficit net of transfers to households
  (`PAR[:WTRbar]`, new). Where households dissave the transfer is imputed and `Sfbar` falls by
  it, below zero on Lebanon and Kyrgyzstan, so a foreign-saving shock stated as a percentage
  of `Sfbar` changes sign there: state it on `Sfbar + Σ WTRbar` (the SAM's `INV × ROW`).
- New: `transfers.csv` / `prepare_data!(...; transfers_path=)`, `read_transfers_csv`,
  `fold_transfers!`.
- Under `:bop` with fixed investment the saving fixes alone leave every real result unchanged
  (consumption is `(C0/YH0)·YH` either way); C-6 on real GDP and the lump-sum transfers do
  change them.

---

## Closures

`parameters(data)[:trade_closure]` selects how the current account is closed
(set it through `prepare_data!(...; trade_closure=, bop_closure=)`; the choice is
stored in `data.par` and read by the equation files, so it must be made before
calibration).

- **`:bop`** (default) — small open economy with a balance-of-payments
  equation. World prices are exogenous in foreign currency (`PWE0`, `PWM0`,
  chosen so every benchmark price is 1) and converted with the real exchange
  rate `ER` (E-2: `WPE = ER·PWE0`; T-21: `WPM = (1+zeta_t)·ER·PWM0`; T-20 then
  gives the producer's export price `PE = WPE/(1+tau_e)`). Imports (T-9) and
  exports (T-18) are independent, the SAM's own trade flows and final demand
  are the benchmark (no rescale), and the current-account deficit is booked as
  exogenous foreign saving `Sfbar` (= CIF imports − FOB exports − export tax − net
  transfers to households, the SAM's `INV × ROW` net of `ROW × INV` less the transfer of
  convention (7)). `C_BOP` (CIF imports − FOB exports = `Sf` + `PNUM·PABS·Σ WTRbar`)
  is imposed and the savings–investment balance C-9 is dropped, since by
  Walras' law it is the same restriction (`export_results!` reports its
  residual as `SI_gap`). Investment is then pinned by `C_INV`
  (`FD[Inv] = PAR[:FDInv0]`, an exogenous real level that
  `update_period_data!` re-bases to `chi_inv·RGDP` between periods; a
  share-of-GDP rule was tried and left the model's real scale nearly
  unanchored under the fixed-wage labour closure). `bop_closure` picks what
  clears the current account: `:flex_er` (default; `Sf = PNUM·ER·Sfbar`, `ER`
  adjusts) or `:fixed_er` (`ER = ER0`, the home region's `Sf` adjusts).
- **`:balanced`** — the legacy convention: E-2 `WTFd = lambda_w·WTFs` and
  T-21 `WPM = WPE/lambda_w` tie the CIF value of imports to the FOB value of
  exports good by good, `Sf = 0`, C-9 is imposed, and there is no `ER`
  variable. The benchmark must be trade-balanced, which `calibrate_from_sam!`
  enforces by rescaling final demand (Calibration conventions above). Kept so
  that old results reproduce; the two regimes have identical equation counts.

**Investment under `:bop`** (`inv_closure`). `:fixed` (default) is the `C_INV` rule
above. `:savings` keeps C-9 instead (investment = household + government +
foreign saving) and leaves `C_BOP` as the equation of `ER`; the system stays
square and replicates its benchmark, but with `bop_closure = :flex_er` the real
scale of the economy is weakly anchored (a uniform +10 pp tariff moved real gross output, M-2
`RGDP` = ΣXP, by −46 % on the synthetic SAM and −26 % for South Africa, against −2 % for
Kenya), so use it only with `bop_closure = :fixed_er`, where it gives the
Keynesian fixed-wage response (tariff revenue is saved and invested: real gross output
+0.4–0.6 %, investment +4–14 %). Under the default the same shock lowers real
gross output by 6–11 % (Kenya, Benin, South Africa, synthetic SAM): the extra government
saving is not invested and, with fixed wages, demand falls. The comparison is in
`~/Documents/Data/CGE/validation/sdcge_branch_check/inv_closure_comparison.txt`.

### Labour closures

`parameters(data)[:labour_closure]`, set after `prepare_data!` and before `model(data)`:

| closure | wage | unemployment | capital | needs |
|---|---|---|---|---|
| `:fixed_wage` (default) | `W = 1` | absorbs supply − demand; labour demand is **uncapped** | perfectly elastic at `TR = 1` | any trade closure |
| `:full_employment` | clears the labour market | fixed at its benchmark rate | clears on `TR` (stock exogenous) | `bop_closure = :fixed_er` (warns otherwise) |
| `:wage_floor` | at its floor `WMIN` (the benchmark wage) while there is unemployment, rises once there is none | `TW ≥ WMIN ⟂ UE ≥ 0`, starting from the database's rate | clears on `TR` | `trade_closure = :bop`, `bop_closure = :fixed_er` (errors otherwise) |

```julia
data = prepare_data!(init_data(); source=:csv, sam_path=…, sets_path=…, bop_closure = :fixed_er)
parameters(data)[:labour_closure] = :wage_floor
# optional override of the database's benchmark unemployment (one rate, or by skill):
set_benchmark_unemployment!(data, Dict("UnSkLab" => 0.08, "SkLab" => 0.04))
```

`:wage_floor` is the LINKAGE minimum-wage regime: labour clears with unemployment
(F-6: `Σ LV = LS·(1 − UE)` ⟂ `TW`) and F-10 is the complementarity `TW − WMIN ≥ 0 ⟂ UE ≥ 0`, so an
expansion first hires the unemployed at the benchmark wage and raises the wage only once a skill's
unemployment is exhausted; a contraction raises unemployment at the floor. A SAM records only the
employed, so the benchmark unemployment rate comes from data: `prepare_data!` reads
`unemployment.csv` (`labour,rate`, UnSkLab/SkLab) next to a CSV SAM — the country databases from
`~/Documents/Data/CGE` carry the ILO modelled 2023 rate, split by skill from labour-force surveys —
and the first `:wage_floor` build applies it (labour force = employment / (1 − u)); without a file
the rate is 0, and `set_benchmark_unemployment!` overrides it. With a 0 % start the closure is
`:full_employment` (tested to 1e-8). Capital, the old vintage, `W = (1+τ_l)·NW`, the
natural-resource routing and the numeraire are `:full_employment`'s, so it needs a fixed real
exchange rate: under `:flex_er` the nominal floor and the CPI numeraire pin the real consumer wage,
and a 20 % tariff cut on Kenya 2023 *lowers* real gross output (ΣXP) 4.3 % with unemployment rising from 5 % to
13 % (zero tariffs does not converge). The first `:wage_floor` build scales the labour force in the
parameter table; build another closure from a fresh `prepare_data!` if it must stay unscaled.

Batch on the 333 country databases (**pre-fix**: measured before the household-income, transfer and GDP
fixes of 2026-10-07, which change these numbers; a re-validation on the merged branch follows) (2026-10-08, main `16fd6c7` = `4de52d8` plus docs, on the
data pipeline's exports with the value-added floor; SAM at the simulator's benchmark scale,
largest flow 1e5; `~/Documents/Data/CGE/validation/sdcge_closure_check/`): real GDP on the
expenditure side at benchmark prices (`real_gdp_pct`, as the simulator computes it since P29;
SDCGE's own M-2 `RGDP = ΣXP` is real GROSS output, in brackets after the median) against each
closure's own benchmark solve, median [90th percentile, maximum] over the solved cases. A case
counts as solved when PATH certifies it AND it passes the simulator's Walras-law gate: C-9,
dropped under `:bop`, within 1e-3 of investment of the benchmark solve's residual
(`si_gap_rel`, `walras_ok`).

| | benchmark solves | −20 % tariffs | zero tariffs |
|---|---|---|---|
| `:fixed_wage` (`:flex_er`) | 328 | 315: +0.46 % (gross +0.53 %) [3.2 %, 70 %] | 287: +2.3 % (gross +2.6 %) [16 %, 131 %]; labour demand above the labour force in 275 (median 2.8 %, up to 153 %) |
| `:full_employment` (`:fixed_er`) | 328 (328) | 323 (323): +0.02 % (gross +0.07 %) [0.16 %, 1.2 %] | 301 (303): +0.08 % (gross +0.32 %) [0.63 %, 4.0 %] |
| `:wage_floor` (`:fixed_er`, database rates) | 328 (328) | 326 (326): +0.10 % (gross +0.17 %) [0.65 %, 3.3 %] | 308 (312): +0.44 % (gross +0.72 %) [2.1 %, 7.1 %] |

In brackets after the counts: solved when a failed solve is retried with the SAM at its own
scale (the 1e5 normalisation that rescues `:fixed_wage` puts the two market-clearing closures'
benchmark residuals, ~5e-6, above PATH's absolute 1e-6, so PATH has to move on a near-singular
Jacobian; at the raw scale the start point already passes; where both solve, real GDP agrees to
5e-6 points). The gate rejects 20 cases PATH certifies at 1e5 and 18 at the raw scale, all under
zero tariffs except Mauritania's −20 % runs: C-9 misses by 0.1–0.4 of investment on Mauritania (`:full_employment`
−20 % tariffs at the raw scale −0.39, zero tariffs −0.42; `:fixed_wage` −20 % at 1e5 −0.07, its
+21 % gross output a false solution), Comoros (zero tariffs, raw scale, −0.27), Burundi,
Madagascar and Slovenia (`:wage_floor` zero tariffs −0.07, the one case that put `:wage_floor`
below `:full_employment`), and by 1e-3–9e-3 on Gabon, Mozambique, Russia and Togo. Paired,
`:wage_floor` real GDP is never below `:full_employment` and is below `:fixed_wage` in 91 %
(−20 %) / 95 % (zero tariffs) of the databases where all three solve. After zero tariffs both
skills are still unemployed at the floor wage in 253 databases, one skill has reached full
employment in 39 and both in 16. Kenya 2023 (database rates 5.2 % / 8.2 %), real GDP:
−20 % tariffs +0.87 % (gross output +1.12 %), zero tariffs +3.0 % (+5.2 %; unemployment 0, wage
+20 %, nominal GDP +26 %); its 133-sector hybrid +2.2 %; South Africa +1.5 %, Cameroon +1.9 %
(fixed wage: +65 %). A three-period run on Kenya keeps unemployment and the wage flat across
periods (no zig-zag). Every closure solves the same 328 benchmarks; the five that fail under all
three have no crop land in GTAP 12 (Hong Kong `gtap12` and `hybrid133`, Mauritius `gtap12` and
`hybrid133`, Iceland `hybrid133`). Against the batch before the value-added floor (2026-10-06,
first run on the `wage-floor` branch plus 14 re-run databases), 32 `:full_employment` and 38
`:wage_floor` cases solve that did not, all but one on the 90 databases whose sectors without
value added (Syria `gtap12` `nfm`/`i_s`, Mongolia ten manufactures, Belgium's 133-sector coal)
the pipeline floors at 0.1 % of output; Mongolia `gtap12` and Syria `hybrid133` solve every case
at 1e5. Tajikistan `gtap12` under `:fixed_wage`, −20 % tariffs, gives −5.1 % (−16.8 % in gross
output before the 2026-10-06 revenue fix): `YG` and `Sg` feed only C-9, which `:bop` drops, so
that was PATH reaching a different solution of the weakly anchored fixed-wage,
flexible-exchange-rate equilibrium (unemployment 17 % against 5 %), not an effect of the fix.

The data pipeline in `~/Documents/Data/CGE` ships each country with `sam.csv`
(the real SAM, for `:bop`) and `sam_balanced_trade.csv` (pre-balanced good by
good, for `:balanced`).

---

## Known limitations

- **Labour closure.** (`:wage_floor`, above, is the supported alternative that respects the
  labour force; the notes below predate it.) The default regime (`parameters(data)[:labour_closure]
  = :fixed_wage`) fixes the wage `W` at 1 and lets unemployment `UE` absorb
  any gap between labour supply and demand (F-10); F-4, F-6, F-7 and F-11 are
  circular/degenerate in this regime, so once `UE` leaves its bound PATH
  loses its footing. Consequences: `g_labor` in the dynamic runs raises
  unemployment instead of output, and 7 of 10 periods of the default
  10-period run end in `ITERATION_LIMIT` (non-converged periods do not update
  the state, so the path plateaus; a warning is printed). TFP growth alone
  solves every period.
  A market-clearing alternative is implemented and selectable — **run it with
  `prepare_data!(...; bop_closure = :fixed_er)`**, which is what anchors the
  nominal price level (see below):

  ```julia
  data = prepare_data!(init_data(); source=:csv, sam_path=…, sets_path=…,
                       bop_closure = :fixed_er)
  parameters(data)[:labour_closure] = :full_employment
  ```

  F-6 becomes labour-market clearing ⟂ `TW`, `UE` is fixed at its benchmark,
  `NW = φ·TW`, `W = (1+τ_l)·NW`, F-21 becomes capital-market clearing ⟂ `TR`
  (in the default regime `TR` has no determining equation — F-21 collapses to
  `TR = TR`, hidden because the benchmark start is already the solution), and
  F-24 fixes the old vintage's capital-output ratio so that `RR` is a genuine
  relative return. `PAR[:numeraire]` defaults to `:pabs` under `:fixed_er` and
  to `:cpi` under `:flex_er`; do not mix them up.

  Measured on the 65-sector GTAP-Africa databases: square, benchmark replicates
  (max |ΔXP| 0.0004 % on KEN/TCD, 0.00008 % on ZMB), and a +10 % TFP shock is
  `LOCALLY_SOLVED` in 1–10 s with **gross output up 13.9–27.2 %, real household
  consumption up 15.8–31.0 %, employment unchanged and `UE` = 0** — where the
  default `:fixed_wage` regime turns the same shock into a 14–20 % fall in
  employment (`UE` reaches 0.95 on Chad) and flat-to-falling output. A
  23-period TFP path at 1.2 %/yr converges **23/23 on Kenya, Chad, Zambia and
  Lesotho** and 21/23 on Cabo Verde (two late periods end in `SLOW_PROGRESS`).

  **It is still not the default** because of that last case, because the
  `gtap11afr`/`hybrid` batch has not been re-run under it, and because it needs
  the non-default `:fixed_er` trade closure, which `PAR[:labour_closure]` cannot
  set (`bop_closure` is a `prepare_data!` argument). With `:flex_er`
  the model is homogeneous of degree one in every nominal price, `PABS = 1` is
  not a numeraire at all, and `PABS` ends up clearing the land market instead
  (land supply becomes demand-determined: +141 % on KEN). See the `Factors.jl`
  header for the measurements and for the closure experiment that was tried and
  rejected.
- **Transfers imputed from dissaving.** No exported SAM carries a transfer cell, so the transfer
  booked where households dissave (convention (7)) is everything that finances the gap:
  remittances, aid and foreign borrowing alike, all held fixed in real domestic terms. GTAP 12's
  balance-of-payments flows (remittances `rmi`/`rmo` in `summary.json`'s `bop_flows`) would
  separate the remittances; the pipeline would export them as `transfers.csv`. Aid to the
  government is not separated (it stays in foreign saving and `S_G`). The households'
  final-demand tax is inside `kappa_h` (`tau_Ac = 0`), so `GDP` is at market prices except for
  that tax.
- **Vintages carry no technology**: `Calibration.jl` gives Old and New capital
  identical shares and prices, so the dynamic update keeps the benchmark
  Old/New split (`vintage_rule=:benchmark_shares`); the flow-based split
  (`:flow`) is available but shifts output between vintages without economic
  content.
- `AT` (productivity) enters P-1/P-2 as `AT·…` but P-3 as `…/AT`; the CES
  demand form should carry `AT^(σ−1)`. Harmless at the benchmark (`AT = 1`),
  but check before relying on TFP shocks quantitatively.
- `solve_model!` uses PATH's default convergence tolerance (1e-6, absolute on the
  residual norm). The earlier 1e-8 made PATH stop with `ITERATION_LIMIT`
  ("cumulative minor iterlim met") or `SLOW_PROGRESS` on solutions it had
  reached to 1e-8 relative, e.g. the second period of a zero-growth run.
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
├── data/                   synthetic SAM (CSV + Excel), policy_experiments.xlsx template;
│                           country databases <ISO3>_<year>_<family>/ (git-ignored, from AEILabs/CGE-SAMs)
├── examples/               numbered walkthrough scripts (start here)
├── test/runtests.jl        smoke test (+ test/data: a real Kenya SAM fixture)
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
