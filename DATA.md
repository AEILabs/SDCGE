# DATA.md — building a real-country SAM for LCGE-V4 (`LinkageModel`)

This document specifies what data this model needs and the exact format/shape it must take:
files, the 2N+16 account list and order (216 for the default 100 sectors), what each block of
the matrix must contain, sign
conventions, the balance requirement, every derived parameter and its formula, the scenario
workbook, and the in-memory `Scenario` struct. It does **not** prescribe where your source data
should live — every loader below takes a path argument, so location is the caller's choice. A
short, non-prescriptive "possible sources" section at the end suggests what national accounts,
input-output and GTAP data *could* supply; it is not part of the format contract.

All citations are `file:line` against branch `n-sector-sets` of this repo as checked out.

---

## 1. What the model needs, conceptually

LCGE-V4 needs one square, balanced **Social Accounting Matrix (SAM)** with:

- **N activities** and **N commodities** (one-to-one, LINKAGE-style — activity `i` produces
  commodity `i`), grouped internally into crops (`:cr`), livestock (`:lv`), energy (`:e`),
  fertiliser (`:ft`), feed (`:fd`) and "rest" (everything else) — this grouping drives which
  production nest a sector uses, not just a label (`src/Calibration.jl:69-71,121-127`; sets
  defined in `src/Types.jl`). N defaults to 100 with the positional groups crops `P001`–`P010`,
  livestock `P011`–`P020`, energy `P071`–`P075`, fertiliser `P076`–`P078`; **any other sector
  list is supplied by file** with `prepare_data!(...; sets_path=...)`, which must then also
  carry `cr`, `lv`, `e`, `ft`, `fd` (`ft` and `fd` disjoint from `e`) — see §4.
- **5 factors**: two labour skills (unskilled, skilled), capital, land, one sector-specific
  "natural resource" (`src/SAM.jl:15`).
- **2 capital vintages** (Old/New) — not a SAM account; a per-sector split applied after loading
  (`src/Types.jl:43`, `src/Calibration.jl:217-240`).
- **6 tax accounts**: output tax, intermediate-input tax, import tariff, export tax, a factor
  tax and an income tax account that exist in the schema but are **not read by the calibrator**
  today (`src/SAM.jl:16`; see §2.3, §4).
- **4 institutions**: one representative household, government, investment/savings, rest of
  world (`src/SAM.jl:17`).
- **1 trade-margin account** (`src/SAM.jl:18`).

That is N+N+5+6+4+1 = **2N + 16 accounts**; for the default N = 100 that is **216**, verified
directly against `data/csv/sam_accounts.csv` (217 lines = 1 header + 216 accounts).

Beyond the SAM, the **dynamic** extension needs: a labour-supply growth rate per skill per
period, a total-factor-productivity path per activity per period, an economy-wide land-supply
growth rate per period, an economy-wide natural-resource-supply growth rate per period, and a
capital depreciation rate (§3). It also needs **elasticities** — the model ships none tied to
data; every one defaults to 0.5 (§2.3).

### What is synthetic in the bundled dataset today

- **The whole SAM is synthetic.** `build_default_large_sam!` (`src/SAM.jl:41-143`) generates
  every cell from a deterministic formula (`output = 900.0 + 3.0*pos`, `src/SAM.jl:52`, plus
  fixed cost-share constants) — no national accounts, survey or GTAP extract underlies it.
- **Regional/bilateral trade shares are uniform, not data-driven.** The SAM itself carries **no
  regional dimension** (it is a single national account matrix); the model's 4 pseudo-regions
  (`R1`–`R4`) and every bilateral trade weight (`beta_1`, `beta_2`, `beta_w`, `beta_z`) default
  to `1/|r|` or `1/(|r|·|rp|)` regardless of what SAM is loaded
  (`src/ParameterTables.jl:227-229,242`), and `calibrate_from_sam!` never touches them.
- **All elasticities default to 0.5.** A single constant `LCGE_SIGMA = 0.5`
  (`src/Calibration.jl:57`) drives every CES/CET share the calibrator derives; every other
  elasticity (Armington tiers, CET export tiers, consumption, migration, factor supply) is a
  flat `0.5` default in `src/ParameterTables.jl` (§2.3 has the full list) that calibration never
  overrides.
- **No base year, no currency.** Nothing in `LinkageData`, `SAM.jl` or `Calibration.jl`
  references a year or a currency/unit; SAM cell values are plain dimensionless magnitudes.
- **The loader requires exactly the 2N + 16 accounts the model derives from `data.sets[:i]`**
  (the 216 of `data/csv/sam_accounts.csv` for the default sectors). `read_sam_csv!`/
  `read_sam_excel!` compare row labels to column labels, then compare the file's labels to
  `data.sam_accounts[:all]` as a set and name the mismatches. The sector list itself comes from
  `sets_path`, not from the SAM file (§4).

---

## 2. The exact data contract

### 2.1 File formats accepted

Two loaders, both invoked through `prepare_data!(data; source=..., sam_path=...)`
(`src/ModelBuilder.jl:29-79`):

| `source` | Function | Format |
|---|---|---|
| `:csv` | `read_sam_csv!(data, sam_path)` | A CSV whose first row and first column hold the 2N+16 account labels; cell `[r+1,c+1]` is the flow from column account `c` to row account `r`. |
| `:excel` | `read_sam_excel!(data, sam_path; sheet="SAM")` | A workbook with a sheet named `SAM`; the reader sizes the range `A1:<(2N+17)th column><2N+17>` from `length(data.sam_accounts[:all])`, i.e. from `data.sets[:i]`, so the sheet must be exactly (2N+17)×(2N+17) — 217×217 at N = 100. |

A third keyword, `sets_path`, reads the sector list and its groupings from a two-column
`set,item` CSV *before* any of this (`read_sets_csv!`); it is what makes N ≠ 100 possible (§4).
The file needs the exact header row `set,item` and one row per (set, item) pair; item order
inside a set is preserved, and the order of `i` fixes the order of the `ACT_`/`COM_` accounts
(§2.2). `data/csv/sets.csv` is the shipped 100-sector example.

| Set | Meaning | Required | If the file omits it |
|---|---|---|---|
| `i` | the N activity/product codes (activity `i` produces commodity `i`); items must be unique | **yes** | error |
| `j`, `k` | equation-index aliases of `i` | no | copied from `i` |
| `cr` | crops — get the crop production nest (land + fertiliser) | yes unless N = 100 | positional default `i[1:10]`; error if N ≠ 100 |
| `lv` | livestock — get the livestock nest (feed + land) | yes unless N = 100 | positional default `i[11:20]`; error if N ≠ 100 |
| `e` | energy goods — the energy bundle inside every sector's intermediate nest | yes unless N = 100 | positional default `i[71:75]`; error if N ≠ 100 |
| `ft` | fertiliser goods — the crop nest's fertiliser input | yes unless N = 100 | positional default `i[76:78]`; error if N ≠ 100 |
| `fd` | feed goods — the livestock nest's feed input | yes unless N = 100 | positional default `i[1:10]`; error if N ≠ 100 |
| `ag` | agriculture — the only sectors that keep their `LAND` payments; land outside `ag` is silently moved to capital (`Calibration.jl:94-100`) | no | `cr ∪ lv` |
| `ip` | non-agricultural ("industry") sectors — `F_PT_nonag`/`F_NPT_nonag`/`F_Td_nonag`/`F_Ts_nonag` fix their land and natural-resource blocks (`Factors.jl:268-272,422`), and policy scenario 8 shocks them | no | `i ∖ ag` |
| `nf` | non-farm sectors (destructured in `Production.jl`/`Other.jl`; no current equation indexes it) | no | `i ∖ ag` |
| `nnft`, `nnfd` | complements `i ∖ ft`, `i ∖ fd` | no | **always** recomputed, even if supplied |
| `r` | regions; `rp` follows `r` | no | `R1`–`R4`; a single line `r,R1` gives `\|r\| = 1`, which builds square and replicates the benchmark |

Any other set (`v`, `l`, `ul`, `sl`, `h`, `f`, `in`, `t`, `gz`, `gs`) is left to `default_sets!`
and should stay at its default: the calibrator writes `LY0`/`LV0` under the literal skill names
`UnSkLab`/`SkLab` (`Calibration.jl:182-183,367-368`), `a_f`/`XAf0` under `Gov`/`Inv`
(`Calibration.jl:417-418,431`), the household is the SAM's single `HH` account
(`Calibration.jl:78`) and P-5 reads the `Old` vintage by name (`Production.jl:32`).

`read_sets_csv!` validates that every item of `j k cr lv ag ip e ft fd nf nnft nnfd` is a member
of `i`, and that `ft ∩ e = ∅` and `fd ∩ e = ∅` (§4); `default_sets!` then errors if `|i| ≠ 100`
and any of `cr lv e ft fd` is still missing, rather than applying the positional defaults to a
sector list they do not describe.

Both readers require `row_accounts == col_accounts` exactly, i.e. the same 2N+16 labels in the
same order on both axes, and both hand off to `set_sam!`
(`src/SAM.jl:31-39`), which re-derives `data.sam_index` from whatever label list was read — so
a caller *can* supply its own account list/order, but every other part of the pipeline
(`calibrate_from_sam!`, all equation files) looks up `ACT_<code>`/`COM_<code>` for the codes in
`data.sets[:i]` plus the fixed `LAB_UNSK`…`TRD_MRG` names, so the labels must match the account
list `setup_sam_accounts!` derives — `read_sam_csv!`/`read_sam_excel!` now check this and name
the missing/unexpected labels instead of failing later inside `calibrate_from_sam!`.

The repository's own bundled example illustrates the shape: `data/csv/sam.csv` (216×216 data +
label row/column) paired with `data/csv/sam_accounts.csv`, or the `SAM` sheet of
`data/linkage_100sector_data.xlsx` (confirmed by direct inspection with `XLSX.jl`: sheet `SAM`
is 217×217, sheet `SAM_Accounts` is 217×3 with columns `group,account,label`, sheet `Sets` is
447×2 with columns `set,item`, sheet `BalanceCheck` is 217×5 with
`account,row_sum,column_sum,gap,abs_gap`). These are illustrations of the format, not a
prescribed location — `sam_path` is a caller-supplied argument.

### 2.2 The 2N+16-account list, in order — shown for the default N = 100 (from `data/csv/sam_accounts.csv`, verified)

For any other sector list the shape is the same: `ACT_<code>` for every code of `data.sets[:i]`
in file order, then `COM_<code>` in the same order, then the fixed 16 non-sector accounts.

| Rows (1-indexed) | Group | Codes | Count |
|---|---|---|---|
| 2–11 | activities | `ACT_P001`–`ACT_P010` (crops) | 10 |
| 12–21 | activities | `ACT_P011`–`ACT_P020` (livestock) | 10 |
| 22–71 | activities | `ACT_P021`–`ACT_P070` (rest) | 50 |
| 72–76 | activities | `ACT_P071`–`ACT_P075` (energy) | 5 |
| 77–79 | activities | `ACT_P076`–`ACT_P078` (fertiliser) | 3 |
| 80–101 | activities | `ACT_P079`–`ACT_P100` (rest) | 22 |
| 102–201 | commodities | `COM_P001`–`COM_P100`, same crop/livestock/rest/energy/fertiliser/rest sub-ranges as activities | 100 |
| 202–206 | factors | `LAB_UNSK`, `LAB_SK`, `CAP`, `LAND`, `NRES` | 5 |
| 207–212 | taxes | `TAX_OUT`, `TAX_INT`, `TAX_IMP`, `TAX_EXP`, `TAX_FACT`, `TAX_INC` | 6 |
| 213–216 | institutions | `HH`, `GOV`, `INV`, `ROW` | 4 |
| 217 | margins | `TRD_MRG` | 1 |

Total 2N+16 accounts (N+N+5+6+4+1) = 216 at N = 100, constructed by `setup_sam_accounts!`
from `data.sets[:i]`, asserted in `test/runtests.jl`
(`@test length(data.sam_accounts[:all]) == 2 * length(data.sets[:i]) + 16`).

### 2.3 What each block must contain, and orientation

**Convention: rows receive payments, columns pay** (`src/SAM.jl:7`, and directly verified in
`calibrate_from_sam!`: `M[idx["COM_"*jj], idx["ACT_"*ii]]` is read as intermediate use — the
commodity **row** receives from the activity **column**, `src/Calibration.jl:83-84`;
`M[idx["COM_"*p], hh_col]` is household consumption — the commodity row receives from the
household column, `src/Calibration.jl:152`).

| Block (row × column) | Content |
|---|---|
| `COM_j` × `ACT_i` | Intermediate use: purchases of commodity `j` by activity `i` (`src/Calibration.jl:83-84`) |
| `LAB_UNSK`/`LAB_SK`/`CAP`/`LAND`/`NRES` × `ACT_i` | Factor payments by activity `i` (`src/Calibration.jl:86-90`) |
| `TAX_OUT`/`TAX_INT` × `ACT_i` | Output tax and intermediate-input tax paid by activity `i` (`src/Calibration.jl:91-92`) |
| `ACT_i` × `COM_i` | Activity `i`'s sales of its own output to the matching commodity account (`src/SAM.jl:79`, illustrative generator) |
| `ROW` × `COM_i` | Imports of commodity `i` (`src/Calibration.jl:135`) |
| `TAX_IMP` × `COM_i` | Import tariff revenue on commodity `i` (`src/Calibration.jl:136`) |
| `TRD_MRG` × `COM_i` | Trade-margin services embodied in commodity `i` (`src/Calibration.jl:137`) |
| `COM_i` × `HH`/`GOV`/`INV` | Final demand for commodity `i` by household / government / investment (`src/Calibration.jl:152-154`) |
| `COM_i` × `ROW` | Exports of commodity `i` (`src/Calibration.jl:133`) |
| `TAX_EXP` × `ROW` | Export tax revenue (`src/Calibration.jl:144`) |
| `HH` × factor rows | Factor income paid out to the household (all factor income routes through `HH`) |
| `GOV` × tax rows | Tax revenue paid out to government |
| institution × institution (`HH`/`GOV`/`INV`/`ROW`) | Transfers, savings, and the current-account closure |

### 2.4 Balance requirement and RAS

- `validate_sam!(data; require_balanced, tol)` (`src/SAM.jl:189-200`) checks: square
  (`SAM.jl:192`), account count matches `sam_accounts.csv` length (`SAM.jl:193`), no negative
  cells (`SAM.jl:194`), and — only if `require_balanced=true` — that
  `max|colsum − rowsum| ≤ tol` (default `tol=1e-6`, `SAM.jl:196-197`).
- `prepare_data!` first validates the raw SAM **without** requiring balance
  (`require_balanced=false`, `src/ModelBuilder.jl:53`), then, unless you pass `balance=:none`,
  calls `balance_sam_ras!(data; maxiter=2000, tol=1e-8)` (`src/SAM.jl:202-229`): a biproportional
  RAS to a symmetric row/column target (`target = 0.5·(rowsum+colsum)`, rescaled to preserve the
  matrix total, floored at `1e-9`, `src/SAM.jl:211-215`), alternating row- and column-scaling
  until the max gap is below `tol` or 2000 iterations elapse (`SAM.jl:217-224`). The raw,
  unbalanced input is kept at `data.metadata[:raw_sam]` for diagnostics (`SAM.jl:207`); the
  RAS output replaces both `data.sam` and `data.balanced_sam` and is re-validated at `tol=1e-5`
  (`SAM.jl:227`).
- `assert_balanced_sam!` (`SAM.jl:279-283`, default `tol=1e-6`) and `sam_balance_summary`
  (`SAM.jl:262-276`: `balanced` boolean at `1e-6`, `max_abs_gap`, `max_gap_account`,
  `total_abs_gap`, `relative_max_gap`) are the checks to run on a new SAM. `prepare_data!` calls
  both automatically and (unless `outdir=nothing`) writes `sam_balance_table.csv`/
  `sam_balance_summary.csv` via `export_sam_balance_report!` (`SAM.jl:291-313`,
  `ModelBuilder.jl:65-67`).
- A correctly calibrated benchmark then replicates to `max residual_at_start ≈ 1e-5` or better
  (measured `8.6e-6`/`9e-6` on the shipped SAM — `CLAUDE.md:44`, `README.md:24`); this is the
  number `diagnose_model` reports (`src/Diagnostics.jl:352-394`) and is the pass/fail bar for
  any new SAM (§5).

### 2.5 Parameters `calibrate_from_sam!` derives from the SAM (with formula)

Everything below is computed inside `calibrate_from_sam!` (`src/Calibration.jl:59-469`) and
stored in `data.par`, then overlaid onto the `PAR` defaults by `precompute_parameters`
(`src/ParameterTables.jl:276-278`, "Overlay SAM-calibrated parameters. Keys in data.par replace
defaults.").

| Parameter(s) | Formula (in words) | Cite |
|---|---|---|
| `tau_p[i]` | Output tax rate = `TAX_OUT_i / (output_i − TAX_OUT_i)` (net-of-tax mark-up) | `Calibration.jl:110` |
| `tau_Ap[·,i]` | Intermediate-input tax = `TAX_INT_i / intermediate_i`, applied uniformly to sector `i`'s whole ND/energy/fertiliser/feed bill | `Calibration.jl:113` |
| `tau_m[·,·,i]` | Import tariff = `TAX_IMP_i / (imports_i + margin-on-imports_i)` | `Calibration.jl:143` |
| `tau_e[·,·,i]` | **One economy-wide** export tax rate = `TAX_EXP_total / Σ_i ES0_i` (not sector-specific, despite the 3-key index) | `Calibration.jl:144` |
| `lambda_w` | `:balanced` only: `(1+tau_m_i)(1+tau_e)` — forces CIF imports = FOB exports good-by-good (convention 2); `1.0` under `:bop`, where E-2/T-21 no longer use it | `Calibration.jl` (trade closure block) |
| `Sfbar[r]`, `PWE0`, `PWM0`, `ER0`, `FDInv0`, `chi_inv` | `:bop` (default): exogenous foreign saving = `Σ XMT0/(1+tau_m) − Σ (1+tau_e)·ES0` (the SAM's trade deficit at world prices, booked on the home region), world prices `1+tau_e` / `1/(1+tau_m)` so every benchmark price is 1, `ER0 = 1`, exogenous real investment `FDInv0 = INVEST0` and its GDP share (used to re-base it between periods). `Sfbar = 0` under `:balanced` | `Calibration.jl` (trade closure and balance-of-payments blocks) |
| `alpha_nd`, `alpha_va`, `alpha_l`, `alpha_hktef`, `alpha_fert`, `alpha_hkte`, `alpha_e`, `alpha_hkt`, `alpha_h`, `alpha_kt`, `alpha_k`, `alpha_t`, `alpha_ff`, `alpha_ktel`, `alpha_tfd`, `alpha_feed`, `alpha_hkte_liv` | Every CES nest share = benchmark cost/value share, scaled `s_j·(P/P_j)^(1−σ)` with `σ=0.5`; collapses to the plain cost share when all nest prices are 1 | `Calibration.jl:259-301` |
| `a_nd`, `alpha_ep`, `alpha_ft`, `alpha_fd` | Composition of the intermediate/energy/fertiliser/feed bundles = `IOc[j,i]/bundle_i`, or uniform `1/n` if the bundle is empty for that sector | `Calibration.jl:303-320` |
| `lambda_k`, `lambda_t`, `lambda_f`, `lambda_ep`, `lambda_ft`, `lambda_fd`, `lambda_l`, `lambda_w` (efficiency, not the tariff `lambda_w` above) | Technical-change/efficiency indices = 1.0 in the benchmark year | `Calibration.jl:340-349` |
| `sigma_p`, `sigma_v`, `sigma_f`, `sigma_e`, `sigma_h`, `sigma_k`, `sigma_feed`, `sigma_ep`, `sigma_ft`, `sigma_fd` | Locked to `LCGE_SIGMA = 0.5` (the constant the shares above were built with) | `Calibration.jl:57,351-357` |
| `LSupply`, `KSupply`, `TSupply`, `FSupply`, `K0`, `LV0` | Benchmark factor supplies = the SAM's `LAB_UNSK`/`LAB_SK`/`CAP`/`LAND`/`NRES` payment rows, divided by the vintage count where relevant | `Calibration.jl:359-369` |
| `UE0[l]` | `0.0` for every skill — "the SAM records labour payments, not an unemployment rate, so full employment is the only rate consistent with it" | `Calibration.jl:370-373` |
| `chi_T[:land]`, `gamma_T[i]` | Ag-only land CET scale/shares; non-agricultural land payments in the SAM are **reassigned to capital** (convention 1) before this is computed | `Calibration.jl:94-100,381-384` |
| `chi_F[i]`, `gamma_K[i]` | Sector-specific-factor scale and capital CET shares from the `NRES`/`CAP` rows | `Calibration.jl:386-388` |
| `kappa_h` | Direct/income tax rate, **solved** (not read from `TAX_INC`) as the residual that makes C-9 (savings=investment) hold exactly: `kappa = 1 − (HH0+2·SAV0)/YH0` | `Calibration.jl:196-206` |
| `beta_m[i]`, `beta_d[i]`, `beta_es[i]`, `beta_xd[i]`, `alpha_dc/mc/df/mf` | Armington import/domestic and CET export/domestic-sales shares = benchmark value shares (nest prices are 1); the household/gov/inv-specific variants are just the aggregate share broadcast to every agent (no agent-specific data in the SAM) | `Calibration.jl:401-408` |
| `theta[k,h]` | `0.0` — ELES subsistence quantities (convention 3, "the SAM carries no information on them") | `Calibration.jl:411` |
| `mu_c[k,h]` | Household consumption budget share = `C0_k/YC0` | `Calibration.jl:412` |
| `GammaC[i,k,h]` | `1` if `i==k` else `0` — one Armington composite per consumption good (diagonal) | `Calibration.jl:413-414` |
| `a_f[i,"Gov"/"Inv"]`, `chi_gov` | Final-demand composition shares and government-spending share of GDP = `GOV0/GDP0` | `Calibration.jl:416-420` |
| `par[:bench]` (the `B` dict, ~35 keys) | The complete benchmark start-value table for every JuMP variable (`XP`, `ND`, `VA`, prices, nest quantities, capital/land/factor demands, trade flows, income/saving/GDP aggregates) — consumed **verbatim** by `initialize_from_sam!`, never re-derived | `Calibration.jl:436-466` |

**Not derived from the SAM at all** (loader limitation, not a data question — see §4): `TAX_FACT`
and `TAX_INC` are declared SAM accounts (`SAM.jl:16`) but `calibrate_from_sam!` never reads
either one (confirmed by grep — the only source references to these two names are the account
list itself and one comment, `Calibration.jl:193`: "TAX_FACT = TAX_INC = 0" ). A real SAM that
populates them today would have that information **silently discarded**; the payroll tax
(`tau_l`) and capital tax (`tau_k`) stay at their `ParameterTables.jl` default of `0.0`
regardless.

### 2.6 Parameters left at `ParameterTables.jl` defaults (elasticities and otherwise)

`precompute_parameters` (`src/ParameterTables.jl:84-281`) fills a generic default for every
parameter first, then the calibration overlay above replaces only the keys `data.par` sets.
Everything in this table is **untouched by calibration** — a real SAM does not affect these
unless the code is extended:

| Parameter(s) | Default | Cite |
|---|---|---|
| `sigma_top_m`, `sigma_w1`, `sigma_w2`, `sigma_w3` (Armington tiers) | `0.5` | `ParameterTables.jl:226,230-232` |
| `sigma_z`, `sigma_z2`, `sigma_TT` (CET export tiers) | `0.5` | `ParameterTables.jl:235,243,251` |
| `sigma_c`, `sigma_mc`, `sigma_mf` (consumption nests) | `0.5` | `ParameterTables.jl:210,216,221` |
| `sigma_ul`, `sigma_sl`, `alpha_ul`, `alpha_sl` (labour-skill disaggregation) | `0.5`; shares uniform `1/|ul|`, `1/|sl|` | `ParameterTables.jl:115-119` |
| `omega_migr` (migration elasticity), `eta_T`, `omega_T` (land supply), `omega_F` (nat.-resource supply), `omega_K`, `eta_k` (capital supply) | `0.5` (or `Inf` regime switches) | `ParameterTables.jl:165,173,177,180` |
| `beta_1`, `beta_2`, `beta_w` (regional/bilateral Armington weights), `beta_z` (CET bilateral export weights) | Uniform `1/|r|` or `1/(|r|·|rp|)` | `ParameterTables.jl:227-229,242` |
| `tau_l` (payroll tax), `tau_t` (land tax), `tau_k` (capital tax) | `0.0` (never read from the SAM — §2.5) | `ParameterTables.jl:167,175,181` |
| `chi_wmin` (minimum wage), `phi_wage` (wage-dispersion shares) | `1.0` | `ParameterTables.jl:168,166` |
| `TRG`, `WTR`, `WTRgov_in/out`, `WTRinv_in/out`, `Sfbar` (transfers, foreign saving) | `0.0` | `ParameterTables.jl:140,144,255-258,260` |
| `tau_pr`, `tau_in`, `tau_out`, `tau_trq_share` (TRQ), `zeta_t` (trade-margin cost) | `0.0` (inactive) | `ParameterTables.jl:244-248` |
| `labour_closure` | `:fixed_wage` | `ParameterTables.jl:157` |
| `numeraire` | `:pabs` | `ParameterTables.jl:162` |
| `sigma_b`, and the vestigial `t`-indexed dynamic defaults (`g_pop`, `g_L`, `g_T`, `pi_dyn`, `gamma_s/t/f/e`, `RGDPMP0`, `Pop0`, `ChiL0/T0/F0`, `lambda_l0/k0/t0/f0/ep0`, `alpha_p_share`) | `0.5` / assorted; **not read by any active equation file** — dead parameters from a removed dynamic block, per `CLAUDE.md:39` | `ParameterTables.jl:183-193,196-201,273` |

**Contradiction found in the code, flagged for the author.** `src/Factors.jl`'s header (lines
20-35) documents `:full_employment` as "(default)", and the closure-selection line itself falls
back to `:full_employment` if the key is absent (`labour_closure = get(PAR, :labour_closure,
:full_employment)`, `Factors.jl:106`). But `precompute_parameters` always sets the key, and it
sets it to **`:fixed_wage`** (`ParameterTables.jl:157`), which is also what `CLAUDE.md:48-50` and
`README.md:277-296` describe as the actual default and the only regime that currently converges
reliably. In practice `:fixed_wage` is the live default; the `Factors.jl` header comment and its
`get(...)` fallback default are stale/misleading. The same pattern repeats one level down:
`Factors.jl:411` falls back to `numeraire = :cpi` if unset, while `ParameterTables.jl:162` always
sets `:pabs` — again `:pabs` is what actually runs.

---

## 3. The scenario workbook and the in-memory `Scenario`

### 3.1 Workbook contract (`write_policy_template`/`read_policy_scenarios`, `PolicyScenarios.jl`)

Five fixed sheet names (`_POL_SHEETS`, `PolicyScenarios.jl:38-44`):

| Sheet | Columns | Grain |
|---|---|---|
| `scenarios` | `sim_id, name, description, periods, delta, active` | one row per scenario |
| `AT_by_activity` | `sim_id, activity, period_1..period_T` | one row per (scenario, activity) — productivity **level** (1.0 = baseline) |
| `g_labor_by_skill` | `sim_id, skill, period_1..period_T` | one row per (scenario, skill) — growth **rate** |
| `g_land` | `sim_id, period_1..period_T` | one row per scenario — economy-wide growth rate |
| `g_nres` | `sim_id, period_1..period_T` | one row per scenario — economy-wide growth rate |

`read_policy_scenarios(path)` (`PolicyScenarios.jl:224-270`) returns only rows whose `active`
cell parses to `"yes"/"y"/"true"/"1"` (case-insensitive, `PolicyScenarios.jl:246-247`); `T`
(periods) and `delta` are read per-scenario from the `scenarios` sheet
(`PolicyScenarios.jl:242-245`). Direct inspection of the shipped `data/policy_experiments.xlsx`
with `XLSX.jl` confirms: `scenarios` sheet is 4×6 (header + 3 rows: `baseline`, `high_TFP`,
`low_TFP`, all `periods=15, delta=0.05, active="yes"`); `AT_by_activity` is 301×17 (header + 300
rows = 3 scenarios × 100 activities, 15 period columns); `g_labor_by_skill` is 7×17 (header + 6
rows = 3 scenarios × 2 skills); `g_land`/`g_nres` are each 4×16 (header + 3 scenario rows, 15
period columns). This is fewer than `write_policy_template`'s own default of 10 scenarios ×
10 periods (`PolicyScenarios.jl:64-65`) — the shipped file is a smaller illustrative example,
not the template's own default shape.

### 3.2 The `Scenario` struct (`PolicyScenarios.jl:25-35`)

```julia
struct Scenario
    sim_id::Int
    name::String
    description::String
    periods::Int
    delta::Float64                        # capital depreciation rate (scalar, per scenario)
    AT::Dict{Tuple{Any,Int},Float64}       # AT[(activity, period)] – productivity LEVEL
    g_labor::Dict{Tuple{Any,Int},Float64}  # g_labor[(skill, period)] – growth RATE
    g_land::Vector{Float64}                # length = periods, economy-wide scalar
    g_nres::Vector{Float64}                # length = periods, economy-wide scalar
end
```

Only these four exogenous channels are first-class: `AT` is a **level replacement** each period
(`PolicyScenarios.jl:347-352`, `update_period_data_scenario!`), `g_labor` **compounds** onto
`LSupply`/`LS0`/`LY0` (`PolicyScenarios.jl:335-342`), `g_land`/`g_nres` compound onto
`TSupply`/`FSupply` with land shares renormalised to sum to 1
(`PolicyScenarios.jl:354-376`). Tariffs, output/direct taxes, government spending share, and
every other `PAR` entry in §2.6 have **no `Scenario` field**.

### 3.3 Building and running a scenario without the Excel round trip

`run_scenario!(data, scen; periods, delta, on_period, abort, period_modifier, ...)`
(`PolicyScenarios.jl:418-464`, documented `README.md:187-215`) runs one `Scenario` entirely in
memory against already-`prepare_data!`d `data` — no workbook, no plots, nothing written to disk.
Three hooks matter for a caller that builds scenarios programmatically (e.g. a web app):

- `period_modifier(t, data)` runs **before** period `t` is built, so a caller can set any `PAR`
  entry outside the `Scenario` struct — tariffs (`:tau_m`), output tax (`:tau_p`), direct tax
  (`:kappa_h`), government spending share (`:chi_gov`), or anything else in §2.6 — giving it a
  time path (`PolicyScenarios.jl:404-407,442`).
- `on_period(t, status, elapsed, history_t)` runs after every period with the PATH termination
  status, wall-clock seconds, and that period's `_period_summary` named tuple
  (`PolicyScenarios.jl:408-410,453-456`).
- `abort::Ref{Bool}` is checked before each period; setting it `true` stops the loop and returns
  whatever periods solved so far (`PolicyScenarios.jl:411-412,441`).

Returns `(history, snapshots)`: `history` is one named-tuple-per-period of ~25 macro indicators
(`RGDP_R1`, `CPI_R1`, `YH_HH`, `KS`, `TR`, `LY_UnSkLab`, `LY_SkLab`, …,
`RecursiveDynamic.jl:292-341`); `snapshots` is every solved variable of every period, keyed
`(variable_symbol, index_tuple)` (`RecursiveDynamic.jl:375-417`). `prepare_data!(...;
outdir=nothing)` (`ModelBuilder.jl:35,68-71`) skips the SAM-report disk writes so a concurrent
caller doesn't race on `"results/"`.

---

## 4. What a real country needs beyond the current loader

| Need | Loader limitation (code) | Data-availability question |
|---|---|---|
| **Sector concordance / N≠100 sectors** | **Implemented.** `prepare_data!(...; sets_path=<set,item CSV>)` reads the sector list into `data.sets[:i]` and derives the 2N+16 SAM accounts from it; the equation layer is generic in N and in the cardinality of `cr`/`lv`/`e`/`ft`/`fd`, so your own codes (GTAP `pdr`, `wht`, …) can be used directly. What you still supply yourself: the memberships `cr`, `lv`, `e`, `ft`, `fd` — the LINKAGE production nests branch on them, not on a label, and the positional defaults (crops = first 10, energy = 71:75, …) only apply to `P001..P100`; `default_sets!` errors rather than applying them to a different N. Constraints the equations impose: `ft ∩ e = ∅` and `fd ∩ e = ∅` (an overlap gives an `XAp` column two equations and breaks squareness), `|e| ≥ 1`, `|ft| ≥ 1` if `|cr| ≥ 1`, `|fd| ≥ 1` if `|lv| ≥ 1`, `|ag| ≥ 1`. `ag`, `ip`, `nf`, `nnft`, `nnfd` are derived. | An ISIC/HS↔your-codes correspondence table is a data-engineering task the caller must do; the model provides no crosswalk. |
| **Real regional/bilateral trade shares** | The SAM has no regional dimension at all — a single national SAM in, and `beta_1`/`beta_2`/`beta_w`/`beta_z` stay at their uniform `ParameterTables.jl` defaults (§2.6) regardless of what is loaded, because `calibrate_from_sam!` never writes them. The four pseudo-regions therefore only spread the national trade totals over 16 uniform `(r, rp)` cells; a `sets.csv` line `r,R1` collapses them to one region, which builds square and replicates the benchmark while cutting the bilateral blocks by `|r|²`. They enter `src/Trade.jl` at T-3/T-5/T-7 (source-region tiering, `Trade.jl:49-75`) and T-18/T-19 (CET bilateral export allocation, `Trade.jl:119-136`). **Not implemented: no loader path writes region-specific shares from data.** | Whether a country-level 4-region or bilateral breakdown even exists (e.g. from a regionalised IO table) is a separate availability question from the missing code path. |
| **Base year and currency** | Nothing in `LinkageData`/`SAM.jl`/`Calibration.jl` carries a year or currency field; SAM values are unitless. | Purely a metadata-tracking task for the caller (the web app's `dataset.toml`/`registry.toml`, §6, carry `base_year`/`currency` fields precisely because the model doesn't). |
| **Real elasticities** | `sigma_p`, `sigma_v`, `sigma_f`, `sigma_e`, `sigma_h`, `sigma_k`, `sigma_feed`, `sigma_ep`, `sigma_ft`, `sigma_fd` are **calibration-locked**: changing them means editing `Calibration.jl`'s single `LCGE_SIGMA` constant (or generalising it to a per-nest value) and re-deriving every `alpha_*`/`beta_*` share, since the share formulas embed `σ` (`Calibration.jl:57,76,259-301`). The remaining elasticities in §2.6 (Armington/CET tiers, migration, factor supply, consumption) are free-standing `PAR` entries a caller can overwrite post-`prepare_data!` with no recalibration needed. | No elasticity file format exists at all (`CLAUDE.md:39`, confirmed — no elasticity data anywhere in `data/`); a real-country elasticity set (e.g. GTAP-derived Armington elasticities) would need a small hand-authored table and a short loader you write yourself. |
| **Labour-force projections for `g_labor`** | `g_labor` is always a **compounding rate** applied uniformly per skill per period, either as a scalar (`update_period_data!`, `RecursiveDynamic.jl:222-238`) or a `Dict{(skill,period)=>rate}` (`Scenario.g_labor`, `PolicyScenarios.jl:32,335-342`). There is no "supply an absolute level path" option (documented gap; still true on this branch). | Converting a population/labour-force projection (e.g. by education level, as a proxy for skill) into per-period compounding rates is a transformation the caller must do before calling `run_recursive_dynamic!`/`run_scenario!`. |

---

## 5. Possible sources (suggestions only — not part of the format contract)

These are conceptual pointers to what kinds of data could inform each SAM block; none of them
prescribe a directory layout or a specific file path, since every loader above takes a path
argument you choose.

- **National-accounts-style data** (GDP by expenditure, value added by industry/ISIC,
  institutional-sector accounts, government spending by function) can inform: value added by
  activity (splitting each activity's `LAB_UNSK`/`LAB_SK`/`CAP`/`LAND`/`NRES` rows), the
  household/government/investment columns of final demand, and aggregate tax revenue (`TAX_OUT`,
  `TAX_INT`, `TAX_IMP`, `TAX_EXP`). It carries no bilateral or regional trade detail and no
  factor-skill split on its own.
- **A multi-regional input-output / ICIO-style table** (intermediate-use matrices by sector,
  bilateral imports/exports by sector and partner) can inform: the `COM_j`×`ACT_i` intermediate-
  use block, the `ROW`×`COM_i` import row and `COM_i`×`ROW` export column, and — only if the
  loader is extended per §4 — genuine bilateral trade shares in place of the uniform `beta_*`
  defaults.
- **A GTAP-style database** (documentation only for this exercise) would additionally offer:
  factor-payment detail split by labour skill/capital/land, tariff and export-tax rate detail by
  bilateral pair (finer than this model's per-sector, non-bilateral `tau_m`/`tau_e`), and
  substitution-elasticity estimates that could replace the flat `0.5` defaults if the
  calibration formulas in §2.5/§4 are generalised to accept them.

### Step-by-step recipe to assemble, balance and validate a new SAM

1. Fix your sector list. Either map your source classification onto `P001..P100`, respecting the
   crop/livestock/energy/fertiliser positional groups (§4), or write a `set,item` CSV listing `i`
   and the memberships `cr`, `lv`, `e`, `ft`, `fd` for your own N codes and pass it as
   `sets_path`. Then lay out the 2N+16 account labels — `ACT_<code>`, `COM_<code>`, then the 16
   fixed accounts — in the order of `data/csv/sam_accounts.csv` (§2.2).
2. Populate each block per §2.3: intermediate use, factor payments, output/intermediate taxes,
   imports/tariffs/margins, exports/export tax, final demand by institution. Leave `TAX_FACT`/
   `TAX_INC` at zero unless you plan to extend `calibrate_from_sam!` to read them (§2.5).
3. Write the matrix as a CSV (labelled first row/column) or the `SAM` sheet of an XLSX workbook
   with a `(2N+17)×(2N+17)` range (§2.1), and call
   `prepare_data!(init_data(); sets_path=<your sets CSV or nothing>, source=:csv_or_:excel,
   sam_path=<your path>)`.
4. Let `balance_sam_ras!` reconcile row/column gaps automatically (default), or pre-balance the
   matrix yourself and pass `balance=:none` to `prepare_data!`.
5. Check `sam_balance_summary(data)` — `balanced == true` and `max_abs_gap` at machine-epsilon
   scale — before trusting anything downstream (§2.4).
6. `prepare_data!` runs `calibrate_from_sam!`/`precompute_parameters` automatically; build and
   solve the benchmark (`model(data)` → `solve_model!(m)`) and check
   `termination_status(m)` is `LOCALLY_SOLVED`/`OPTIMAL` and `diagnose_model(m)`'s
   `residual_at_start` is ≈`1e-5` or smaller (§2.4) before running any shock or dynamic path.

---

## 6. The web-app registry entry

The companion web app (`CGE-intelligence`) reads
`/Users/sebastiankrantz/Documents/web-projects/CGE-intelligence/data/registry.toml` via
`src/Registry.jl`. Each `[[dataset]]` table requires `id, model, label, path, loader, base_year`
(`Registry.jl:15`); every other key (`n_sectors`, `horizon_default`, `resolution`, …) is carried
through unread into `DatasetEntry.options` (`Registry.jl:14,42`). `path` is resolved relative to
the app root and `available = ispath(path)` (`Registry.jl:41,47`) — a dataset need not exist on
disk yet to be listed (disabled) in the UI.

The bundled dynamic dataset's entry (`registry.toml:21-31`):

```toml
[[dataset]]
id          = "synthetic100"
model       = "dynamic"
label       = "Synthetic 100-sector economy (illustrative, no real country)"
country     = ""
base_year   = 2020          # app-chosen mapping for period 1; disclosed in UI
currency    = "SAM units"
path        = "data/dynamic/synthetic100"
loader      = "linkage_csv"
horizon_default = 10
caveats     = ["Synthetic SAM, uniform bilateral shares, all elasticities 0.5", "No BOP closure"]
```

A sidecar `data/dynamic/synthetic100/dataset.toml` documents the same entry plus a longer
`caveats` list and a `source` provenance string; `Registry.jl` itself only reads
`registry.toml`. The `linkage_csv` loader (`src/models/DynamicCGE.jl:452-463`) requires
`entry.path` to be a directory containing `sam.csv` (or a direct path to one), and calls
`SDCGE.prepare_data!(data; source=:csv, sam_path=sam_path, outdir=nothing)` — i.e. exactly the
CSV contract in §2.1, with `outdir=nothing` so the web app never races on `"results/"` writes.

**To register a new (real-country) dataset**: add a `[[dataset]]` table to `registry.toml` with
a unique `id`, `model = "dynamic"`, a `path` pointing at wherever your `sam.csv` lives, `loader =
"linkage_csv"`, a real `base_year`/`currency`, and a `caveats` list disclosing what in your SAM
is still synthetic or assumption-driven (regional shares, elasticities, etc., per §1/§4) — the
same pattern as the `synthetic100` entry above. No code change is needed as long as the SAM
satisfies the §2 contract; `Registry.jl` and `DynamicCGE.jl` are generic over any dataset that
does.
