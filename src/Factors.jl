# Usage: add_factor_equations!(model, data, PAR)
# Complete paper-numbered LINKAGE factor-market equations (F-1)--(F-33).
# MCP form: residual ⟂ left-hand/market-clearing variable.
#
# Changes from the original draft:
#  - F_20_cet removed (duplicate KSs equation; F_22 now the sole KSs definition).
#  - F_24_old_supply removed (duplicate CHIv["Old"] equation; F_23 is canonical).
#  - F_29 repurposed: now defines XPv[i,"New"] = XP[i] - XPv[i,"Old"]
#    (XP itself is defined by T_16 in Trade.jl via the CET aggregation).
#  - F_16 changed to ⟂ NPT (land market clearing determines net land price).
#  - F_17 kept as ⟂ PT (gross land price = net land price × tax wedge).
#  - F_19 split: ⟂ Fs for sectors with exogenous price (omega_F = Inf),
#                ⟂ PF for sectors with elastic supply (omega_F finite) to
#                avoid duplicating F_18_supply_factor on Fs.
#  - K0[i] fixed to calibrated value PAR[:K0][i].
#  - PS[gz] fixed to PABS (price level for minimum-wage determination).
#  - PT and NPT for non-agricultural sectors fixed to 1.0 (no land demand).
#  - Td for non-agricultural sectors fixed to 0.0.
#
# ─────────────────────────────────────────────────────────────────────────────
# LABOUR-MARKET CLOSURE  (PAR[:labour_closure], added 2026-09-03)
# ─────────────────────────────────────────────────────────────────────────────
# :full_employment (default)
#   The wage clears the labour market.  Per skill l:
#     (F-6)  sum_i (LV[l,i] + N_i·LF_d[l,i]) = LS[l,"national"]·(1 − UE0[l])
#                                                        ⟂ TW[l,"national"]
#     (F-7)  TW[l,z]   = TW[l,"national"]     (z in urban/rural)   ⟂ TW[l,z]
#     (F-9)  WMIN[l,z] = chi_wmin·PS^omega_ps·PABS^omega_p·(1−UE)^omega_ue
#                                                                  ⟂ WMIN[l,z]
#     (F-10) UE[l,z]   = UE0[l]               (all zones)          ⟂ UE[l,z]
#     (F-11) NW[l,i]   = phi_wage[l,i]·TW[l,"national"]            ⟂ NW[l,i]
#     (F-12) W[l,i]    = (1 + tau_l[l,i])·NW[l,i]                  ⟂ W[l,i]
#     (F-4)  AVGW[l,z] = employment-weighted mean of NW             ⟂ AVGW[l,z]
#   Labour demand responds to W through P-72/P-73/P-74/P-75 (UW, SW), so F-6 is
#   a genuine downward-sloping-demand market-clearing condition on TW.
#
# :fixed_wage (the previous behaviour, kept for reference/back-compat)
#   F-12 pins W[l,i] = 1 and the wage level is then determined by the chain
#   F-4 (AVGW = mean of NW) -> F-6_segmented (TW = AVGW) -> F-11 (NW = phi·TW),
#   which is CIRCULAR: substituting gives AVGW·1e-9 = 0, so the nominal net wage
#   is pinned only through the 1e-9 guard inside F-4.  On top of that
#   F-7_segmented is (TW − WMIN)·UE = 0 ⟂ UE with TW ≡ WMIN ≡ 1, a complementarity
#   satisfied with zero gradient either way.  Together these leave a
#   near-singular direction that PATH handles badly as soon as the state moves
#   away from the benchmark: labour supply is not binding in this regime (growth
#   in LS only raises UE, which enters no other active equation), and 10-period
#   growth runs failed with ITERATION_LIMIT in 7 of 10 periods.
#
# NUMERAIRE (PAR[:numeraire]).  Under :fixed_wage, W = 1 (F-12) is the de facto
# numeraire — labour is 29 % of output, so it carries the whole cost chain.
# Under :full_employment W is endogenous and one absolute condition has to take
# its place.  M-5 (PNUM = 1) can NOT do it: PNUM only ever multiplies transfer
# terms that are all zero here.  WHICH condition is right depends on the TRADE
# closure, so `precompute_parameters` picks the default from `bop_closure`:
#
#   bop_closure = :fixed_er  ->  :pabs   (THE COMBINATION THAT WORKS)
#       C-ER fixes the real exchange rate at ER0 and world prices WPE/WPM are
#       parameters, so the domestic price level is anchored through the border —
#       the standard small-open-economy closure.  PABS = 1 (F_PABS) is then just
#       a unit normalisation and F-13/F-18 are genuine REAL supply schedules.
#       Measured at AT = 1.10 (65-sector GTAP-Africa databases, :full_employment):
#         KEN  LOCALLY_SOLVED  9.3 s  Sum XP +13.9 %  cons +17.6 %  TR 1.23
#                                     CPI 0.891  TLnd +4.2 %
#         TCD  LOCALLY_SOLVED  1.4 s  Sum XP +15.4 %  cons +15.8 %  TR 1.47
#                                     CPI 1.010  TLnd +1.3 %
#         ZMB  LOCALLY_SOLVED  1.9 s  Sum XP +19.9 %  cons +25.6 %  TR 1.57
#                                     CPI 1.055  TLnd +6.2 %
#       Employment is unchanged and UE stays 0 in all three, and RR now spans
#       0.08-0.66, i.e. the vintage margin is actually doing work.
#
#   bop_closure = :flex_er  ->  :cpi     (usable, not sound)
#       With a flexible real exchange rate and foreign saving exogenous in
#       foreign-currency terms the model is homogeneous of degree one in every
#       nominal price, so PABS = 1 anchors NOTHING: PABS appears only as the
#       deflator of F-13 (TLnd = chi_T·(PTLnd/PABS)^eta_T) and F-18
#       (Fs = chi_F·(PF/PABS)^omega_F), in F_PS and in F-9 (WMIN, inert here).
#       Measured on KEN at AT = 1.10 the :pabs solution is exactly the :cpi
#       solution times lambda = 18.7 (CPI 18.6996 vs 1.0000, PABS 1.0000 vs
#       0.0535, PTLnd/PABS = 5.342 and TLnd = 76.48 identical to five digits);
#       sigma_min for that direction is 2.0e-6 against ||J||_1 = 1.6e3, nonzero
#       only because of the 1e-9 guards in the denominators.  CPI[first(r)] = 1
#       ⟂ PABS restores the rank, but PABS is then the only variable left to
#       clear the land/natural-resource markets and runs to 0.04 (KEN), 14.8
#       (TCD) or 190 (ZMB): land supply becomes demand-determined (KEN TLnd
#       +141 %) and TCD's real household consumption actually FALLS 4.6 % under
#       a +10 % TFP shock.  Raising eta_T 0.5 -> 3.0 changes the real solution by
#       < 0.05 pp and only rescales PABS, confirming eta_T is no longer acting as
#       a supply elasticity.  Prefer :fixed_er for anything quantitative.
#
#   Do NOT combine :fixed_er with :cpi: the border already anchors the level, so
#   CPI = 1 is a second absolute condition and PABS absorbs the difference
#   (ZMB at AT = 1.10: PABS 5.2e-7, TLnd x36, SLOW_PROGRESS in 99 s).
#
# TRIED AND REJECTED (2026-09-14).  Spending the capital-market clearing
#   condition F-21 on the numeraire (CPI = 1 ⟂ TR, on the theory that one
#   market-clearing equation is Walras-redundant) does NOT work: this model is
#   not Walras-closed.  The benchmark still replicates, but at AT = 1.10 total
#   capital demand misses KS by +124 % (KEN), -18.2 % (TCD) and -33.0 % (ZMB).
#   There is therefore no redundant market-clearing equation to spend.

function _lcge_badfinite(x)
    return !isfinite(float(x))
end

# A benchmark factor payment at or below this is "no factor at all".  The
# calibrator writes chi_F[i] = max(nrs[i], 1e-9), so a sector with no
# natural-resource row in the SAM gets 1e-9; SAM totals are 1e4-1e5 and the
# smallest genuine sectoral payment in the shipped databases is ~1e-3.
const LCGE_FACTOR_ZERO = 1.0e-6

function add_factor_equations!(model, data::LinkageData, PAR)
    S = data.sets
    default_sets!(data)
    i  = S[:i];  v  = S[:v];  l  = S[:l]
    r  = S[:r];  rp = S[:rp]; gz = S[:gz]; gs = S[:gs]
    ag = S[:ag]; ip = S[:ip]

    W    = model[:W];    LV   = model[:LV];   LF_d  = model[:LF_d]
    Nfirm= model[:Nfirm]
    LS   = model[:LS];   MIGR = model[:MIGR]; AVGW  = model[:AVGW]
    NW   = model[:NW];   TW   = model[:TW];   WMIN  = model[:WMIN]; UE = model[:UE]
    PS   = model[:PS];   PABS = model[:PABS]
    Td   = model[:Td];   Ts   = model[:Ts];   PT    = model[:PT];   NPT = model[:NPT]
    TLnd = model[:TLnd]; PTLnd= model[:PTLnd]
    Fd   = model[:Fd];   Fs   = model[:Fs];   PF    = model[:PF]
    Kvd  = model[:Kvd];  KF_d = model[:KF_d]; R     = model[:R];    NR  = model[:NR]
    KSs  = model[:KSs];  KS   = model[:KS];   TR    = model[:TR]
    RR   = model[:RR];   CHIv = model[:CHIv]; K0    = model[:K0];   XP  = model[:XP]
    XPv  = model[:XPv];  FDInv= model[:FDInv]; KActual=model[:KActual]
    GammaInv=model[:GammaInv]; KNorm=model[:KNorm]

    # ── Pre-compute regime switches outside JuMP macros ──────────────────────
    migr_integrated = [ll for ll in l if  _lcge_badfinite(PAR[:omega_migr][ll])]
    migr_segmented  = [ll for ll in l if !_lcge_badfinite(PAR[:omega_migr][ll])]
    eta_T_inf  = _lcge_badfinite(PAR[:eta_T])
    lndmax_inf = _lcge_badfinite(PAR[:LndMAX])
    omega_T_inf= _lcge_badfinite(PAR[:omega_T])
    omega_K_inf= _lcge_badfinite(PAR[:omega_K])

    # ── Labour-market closure switch (see the header) ────────────────────────
    labour_closure = get(PAR, :labour_closure, :full_employment)
    labour_closure in (:full_employment, :fixed_wage) ||
        error("Unknown PAR[:labour_closure] = $(labour_closure); " *
              "use :full_employment or :fixed_wage.")
    full_employment = labour_closure === :full_employment
    ue0 = get(PAR, :UE0, Dict{Any,Float64}())

    # :full_employment needs a nominal anchor and only the trade closure can give
    # it one (see the NUMERAIRE section of the header).  Warn rather than error so
    # the pre-2026-09-14 combination can still be reproduced.
    if full_employment && Symbol(get(data.par, :bop_closure, :flex_er)) !== :fixed_er
        @warn string("PAR[:labour_closure] = :full_employment with bop_closure = ",
                     get(data.par, :bop_closure, :flex_er),
                     ": the model is then homogeneous of degree one in every nominal ",
                     "price and PABS ends up clearing the land market instead of ",
                     "being a price level (land supply becomes demand-determined). ",
                     "Rebuild with prepare_data!(...; bop_closure = :fixed_er).") maxlog=1
    end

    # Sector-specific ("natural resource") factor.  A sector with NO such factor
    # in the SAM still got the elastic supply schedule F-18 with the calibrator's
    # floor chi_F = 1e-9, while its demand share alpha_ff is exactly 0: F-18 reads
    # Fs = 1e-9·PF^0.5 and P-23/P-47/P-65 give Fd = 0, so PF[i]'s entire Jacobian
    # column is a single entry of order 5e-10.  On data/KEN_2017_gtap11afr 59 of
    # 65 sectors are in that state and they set the smallest singular value of the
    # Jacobian to 2.9e-10 against ||J||_1 = 1.6e3 — a condition number of 1e13 and
    # the largest single source of PATH's pivoting failures under
    # :full_employment (raising sigma_min to 2.0e-6 turned a +10 % TFP shock on
    # KEN from ITERATION_LIMIT into LOCALLY_SOLVED).  Under :full_employment such
    # sectors therefore take the fixed-real-price branch instead
    # (PF = PABS·PF0, Fs = Fd = 0), which leaves the benchmark solution unchanged.
    # :fixed_wage keeps the original routing so it stays byte-identical.
    _lcge_nofactor(ii) = full_employment && get(PAR[:chi_F], ii, 0.0) <= LCGE_FACTOR_ZERO
    omega_F_inf= [ii for ii in i if  (_lcge_badfinite(PAR[:omega_F][ii]) || _lcge_nofactor(ii))]
    omega_F_fin= [ii for ii in i if !(_lcge_badfinite(PAR[:omega_F][ii]) || _lcge_nofactor(ii))]

    # ── (F-1) Rural labor supply ──────────────────────────────────────────────
    @constraint(model, F_1[ll in l],
        (LS[ll,"rural"]) - ((1 + PAR[:g_l][(ll,"rural")]) * PAR[:LS0][(ll,"rural")] - MIGR[ll]) ⟂ LS[ll,"rural"])

    # ── (F-2) Urban labor supply ──────────────────────────────────────────────
    @constraint(model, F_2[ll in l],
        (LS[ll,"urban"]) - ((1 + PAR[:g_l][(ll,"urban")]) * PAR[:LS0][(ll,"urban")] + MIGR[ll]) ⟂ LS[ll,"urban"])

    # ── (F-3) National labor supply ───────────────────────────────────────────
    @constraint(model, F_3[ll in l],
        (LS[ll,"national"]) - (sum(LS[ll,gg] for gg in gs)) ⟂ LS[ll,"national"])

    # ── (F-4) Average wage per zone ───────────────────────────────────────────
    @constraint(model, F_4[ll in l, gg in gz],
        (AVGW[ll,gg] * (sum(LV[ll,ii] + Nfirm[ii]*LF_d[ll,ii] for ii in i) + 1.0e-9)) -
        (sum(NW[ll,ii] * (LV[ll,ii] + Nfirm[ii]*LF_d[ll,ii]) for ii in i)) ⟂ AVGW[ll,gg])

    # ── (F-5) Migration ───────────────────────────────────────────────────────
    @constraint(model, F_5_fixed_migration[ll in migr_integrated],
        (MIGR[ll]) - 0.0 ⟂ MIGR[ll])
    @constraint(model, F_5_migration[ll in migr_segmented],
        (MIGR[ll]) -
        (PAR[:chi_migr][ll] *
         (((1 - UE[ll,"urban"]) * AVGW[ll,"urban"]) /
          (((1 - UE[ll,"rural"]) * AVGW[ll,"rural"]) + 1.0e-9)) ^ PAR[:omega_migr][ll]) ⟂ MIGR[ll])

    if full_employment
        # ── (F-6) Labour-market clearing: the wage clears each skill market ───
        # Demand responds to W through P-72/P-75 (UW, SW), so this is a genuine
        # downward-sloping-demand condition on the economy-wide skill wage TW.
        @constraint(model, F_6_clearing[ll in l],
            (sum(LV[ll,ii] + Nfirm[ii]*LF_d[ll,ii] for ii in i)) -
            (LS[ll,"national"] * (1 - get(ue0, ll, 0.0))) ⟂ TW[ll,"national"])

        # ── (F-7) Zone wage equals the national wage (one integrated market) ──
        @constraint(model, F_7_zone_wage[ll in l, gg in gs],
            (TW[ll,gg]) - (TW[ll,"national"]) ⟂ TW[ll,gg])

        # ── (F-9) Minimum wage: not binding in this regime, kept for reporting ─
        @constraint(model, F_9_wmin[ll in l, gg in gz],
            (WMIN[ll,gg]) - (PAR[:chi_wmin][(ll,gg)] *
            PS[gg]^PAR[:omega_ps][gg] *
            PABS^PAR[:omega_p][gg]) ⟂ WMIN[ll,gg])

        # ── (F-10) Unemployment exogenous at its benchmark rate ───────────────
        @constraint(model, F_10_ue_fixed[ll in l, gg in gz],
            (UE[ll,gg]) - (get(ue0, ll, 0.0)) ⟂ UE[ll,gg])
    else
        # ── (F-6) National wage condition ─────────────────────────────────────
        @constraint(model, F_6_integrated[ll in migr_integrated],
            ((TW[ll,"national"] - WMIN[ll,"national"]) * UE[ll,"national"]) - 0.0 ⟂ UE[ll,"national"])
        @constraint(model, F_6_segmented[ll in migr_segmented],
            (TW[ll,"national"]) - (AVGW[ll,"national"]) ⟂ TW[ll,"national"])

        # ── (F-7) Zone-specific wage condition ────────────────────────────────
        @constraint(model, F_7_integrated[ll in migr_integrated, gg in gs],
            (TW[ll,gg]) - (TW[ll,"national"]) ⟂ TW[ll,gg])
        @constraint(model, F_7_segmented[ll in migr_segmented, gg in gs],
            ((TW[ll,gg] - WMIN[ll,gg]) * UE[ll,gg]) - 0.0 ⟂ UE[ll,gg])

        # ── (F-8) National minimum wage ───────────────────────────────────────
        @constraint(model, F_8[ll in migr_integrated],
            (WMIN[ll,"national"]) - (PAR[:chi_wmin][(ll,"national")] *
            PS["national"]^PAR[:omega_ps]["national"] *
            PABS^PAR[:omega_p]["national"] *
            (1 - UE[ll,"national"])^PAR[:omega_ue]["national"]) ⟂ WMIN[ll,"national"])

        # ── (F-9) Zone minimum wage ───────────────────────────────────────────
        @constraint(model, F_9_integrated[ll in migr_integrated, gg in gs],
            (WMIN[ll,gg]) - (WMIN[ll,"national"]) ⟂ WMIN[ll,gg])
        @constraint(model, F_9_segmented[ll in migr_segmented, gg in gs],
            (WMIN[ll,gg]) - (PAR[:chi_wmin][(ll,gg)] *
            PS[gg]^PAR[:omega_ps][gg] *
            PABS^PAR[:omega_p][gg] *
            (1 - UE[ll,gg])^PAR[:omega_ue][gg]) ⟂ WMIN[ll,gg])

        # ── (F-10) Unemployment rate ──────────────────────────────────────────
        # Restricted to "national" zone for segmented-migration labor types:
        #   - For segmented l: F_7_segmented already defines UE[l,"urban"] and UE[l,"rural"].
        #     This constraint adds only UE[l,"national"].
        #   - For integrated l: F_6_integrated defines UE[l,"national"];
        #     urban/rural UE is handled by F_10_int_gs below.
        @constraint(model, F_10[ll in migr_segmented],
            (UE[ll,"national"] * (LS[ll,"national"] + 1.0e-9)) -
            (LS[ll,"national"] - sum(LV[ll,ii] + Nfirm[ii]*LF_d[ll,ii] for ii in i)) ⟂ UE[ll,"national"])
        @constraint(model, F_10_int_gs[ll in migr_integrated, gg in gs],
            (UE[ll,gg] * (LS[ll,gg] + 1.0e-9)) -
            (LS[ll,gg] - sum(LV[ll,ii] + Nfirm[ii]*LF_d[ll,ii] for ii in i)) ⟂ UE[ll,gg])
    end

    # ── (F-11) Sectoral net wages ─────────────────────────────────────────────
    @constraint(model, F_11[ll in l, ii in i],
        (NW[ll,ii]) - (PAR[:phi_wage][(ll,ii)] * TW[ll,"national"]) ⟂ NW[ll,ii])

    if full_employment
        # ── (F-12) Gross employer wage = net wage × payroll-tax wedge ─────────
        # W is endogenous: it is what the demand side (P-72/P-75) responds to and
        # what F-6 clears.  With tau_l = 0 (no SAM payroll-tax account) this is
        # W = NW, and both equal phi_wage · TW.  Note this also repairs a latent
        # gap in the fixed-wage regime, where W (labour COST) and NW (labour
        # INCOME) were determined by unrelated equations.
        @constraint(model, F_12[ll in l, ii in i],
            (W[ll,ii]) - ((1 + PAR[:tau_l][(ll,ii)]) * NW[ll,ii]) ⟂ W[ll,ii])
    else
        # ── (F-12) Producer wage anchored to numeraire ────────────────────────
        # All gross employer wages W[ll,ii] are fixed at 1. This pins the absolute
        # wage level across all skills and sectors, breaking the price homogeneity
        # of the CGE system. NW (net wages) is then determined by F_11 from TW.
        # Payroll tax revenue in C_3 still uses tau_l × NW × labor demand.
        @constraint(model, F_12[ll in l, ii in i],
            (W[ll,ii]) - (1.0) ⟂ W[ll,ii])
    end

    # ── (F-13) Aggregate land supply ──────────────────────────────────────────
    if eta_T_inf
        @constraint(model, F_13_inf_eta,
            (PTLnd - PABS * PAR[:PTLnd0]) - 0.0 ⟂ PTLnd)
    elseif lndmax_inf
        @constraint(model, F_13_unbounded_land,
            (TLnd - PAR[:chi_T][:land] * (PTLnd / (PABS + 1.0e-9))^PAR[:eta_T]) - 0.0 ⟂ TLnd)
    else
        @constraint(model, F_13_bounded_land,
            (TLnd - PAR[:LndMAX] / (1 + PAR[:chi_T][:land] * exp(-PAR[:gamma_ts] * (PTLnd / (PABS + 1.0e-9))))) - 0.0 ⟂ TLnd)
    end

    # ── (F-14) Aggregate land price / closure ─────────────────────────────────
    # CET case: PTLnd is the CES dual of sectoral PT[ag]; with gamma_T calibrated
    # to land shares (Σ_ag gamma_T = 1), this gives PTLnd = 1 at benchmark.
    if omega_T_inf
        @constraint(model, F_14_mobile_land,
            (TLnd) - (sum(Td[ii] for ii in ag)) ⟂ TLnd)
    else
        @constraint(model, F_14_cet_land,
            (PTLnd - (sum(PAR[:gamma_T][ii] * PT[ii]^(1 + PAR[:omega_T]) for ii in ag))^(1/(1 + PAR[:omega_T]))) - 0.0 ⟂ PTLnd)
    end

    # ── (F-15) Sectoral land allocation ───────────────────────────────────────
    if omega_T_inf
        @constraint(model, F_15_mobile_land[ii in ag], (PT[ii]) - (PTLnd) ⟂ PT[ii])
    else
        @constraint(model, F_15_cet_land[ii in ag],
            (Ts[ii] - PAR[:gamma_T][ii] * (PT[ii] / (PTLnd + 1.0e-9))^PAR[:omega_T] * TLnd) - 0.0 ⟂ Ts[ii])
    end

    # ── (F-16) Land-market equilibrium → determines net land price NPT ────────
    # In MCP form: excess supply (Ts - Td) ⟂ NPT ≥ 0.
    # Either the market clears (Ts = Td) or the land price is zero.
    @constraint(model, F_16[ii in ag],
        (Ts[ii]) - (Td[ii]) ⟂ NPT[ii])

    # ── (F-17) Gross land price = net price × (1 + land tax) ─────────────────
    @constraint(model, F_17[ii in ag],
        (PT[ii]) - ((1 + PAR[:tau_t][ii]) * NPT[ii]) ⟂ PT[ii])

    # ── PT and NPT fixed at 1 for non-agricultural sectors ────────────────────
    @constraint(model, F_PT_nonag[ii in ip],  (PT[ii])  - (1.0) ⟂ PT[ii])
    @constraint(model, F_NPT_nonag[ii in ip], (NPT[ii]) - (1.0) ⟂ NPT[ii])

    # ── Td fixed at 0 for non-agricultural sectors ────────────────────────────
    @constraint(model, F_Td_nonag[ii in ip],  (Td[ii])  - (0.0) ⟂ Td[ii])

    # ── (F-18) Sector-specific factor supply / fixed real factor price ────────
    @constraint(model, F_18_fixed_factor[ii in omega_F_inf],
        (PF[ii] - PABS * PAR[:PF0][ii]) - 0.0 ⟂ PF[ii])
    @constraint(model, F_18_supply_factor[ii in omega_F_fin],
        (Fs[ii] - PAR[:chi_F][ii] * (PF[ii] / (PABS + 1.0e-9))^PAR[:omega_F][ii]) - 0.0 ⟂ Fs[ii])

    # ── (F-19) Sector-specific factor market equilibrium ─────────────────────
    # For sectors with exogenous price (omega_F = Inf): Fs adjusts to equal Fd.
    # For sectors with elastic supply (omega_F finite): PF adjusts to clear
    #   the market (avoids duplicating F_18_supply_factor on Fs).
    @constraint(model, F_19_inf[ii in omega_F_inf],
        (Fs[ii]) - (Fd[ii]) ⟂ Fs[ii])
    @constraint(model, F_19_fin[ii in omega_F_fin],
        (Fs[ii]) - (Fd[ii]) ⟂ PF[ii])

    # ── (F-20) Single-vintage capital supply allocation ────────────────────────
    # CET branch removed: F_22 is the sole equation for KSs.
    # Mobile-capital branch kept when omega_K = Inf.
    if omega_K_inf
        @constraint(model, F_20_mobile_capital[ii in i],
            (R[ii,"Old"]) - (TR) ⟂ R[ii,"Old"])
    end

    # ── (F-21) Economy-wide capital return / aggregate capital-market closure ──
    # NOTE (2026-09-04): F_21_cet_capital is a DEGENERATE identity in this port.
    # F-24 forces RR[i] = 1, F-26 gives NR[i,"Old"] = RR[i]·TR and F-28 gives
    # R = (1+tau_k)·NR, so R[i,"Old"] == TR for every i; with sum_i gamma_K = 1
    # the CET dual collapses to TR = TR·1 = TR, i.e. a zero row in the Jacobian,
    # and TR has no determining equation.  Under :fixed_wage that is invisible —
    # W = 1 makes the benchmark start already the solution, so PATH terminates
    # before the singular row matters — but with an endogenous wage PATH must
    # actually solve and stalls (443 major iterations of backtracking steps,
    # SLOW_PROGRESS at a residual of 1.6e-6).
    # Under :full_employment, capital is therefore given a real market: the
    # aggregate stock is exogenous (F-25 below, from the calibrated/dynamically
    # updated PAR[:KSupply]) and TR clears total capital demand against it.
    # That also makes investment bite in the recursive dynamics, where only
    # K0 = KSupply["Old"] used to reach the model at all.
    if omega_K_inf
        @constraint(model, F_21_mobile_capital,
            (sum(sum(Kvd[ii,vv] for vv in v) + Nfirm[ii]*KF_d[ii] for ii in i)) - (KS) ⟂ KS)
    elseif full_employment
        @constraint(model, F_21_capital_clearing,
            (sum(Kvd[ii,vv] for ii in i for vv in v) +
             sum(Nfirm[ii]*KF_d[ii] for ii in i)) - (KS) ⟂ TR)
    else
        @constraint(model, F_21_cet_capital,
            (TR - (sum(PAR[:gamma_K][ii] * R[ii,"Old"]^(1 + PAR[:omega_K]) for ii in i))^(1/(1 + PAR[:omega_K]))) - 0.0 ⟂ TR)
    end

    # ── (F-22) Sectoral capital-market equilibrium ────────────────────────────
    @constraint(model, F_22[ii in i],
        (sum(Kvd[ii,vv] for vv in v) + Nfirm[ii]*KF_d[ii]) - (KSs[ii]) ⟂ KSs[ii])

    # ── (F-23) Capital-output ratio by vintage ────────────────────────────────
    @constraint(model, F_23[ii in i, vv in v],
        (CHIv[ii,vv] * XPv[ii,vv]) - (Kvd[ii,vv]) ⟂ CHIv[ii,vv])

    # ── (F-24) Old-vintage technology / relative return on old capital ────────
    # F_24_old_supply removed (would duplicate F_23 for CHIv["Old"]).
    #
    # :fixed_wage — the original pin RR[i] = 1.  Combined with F-23 and F-30 it
    #   forces Kvd[i,"Old"] = K0[i] exactly: half the capital stock is frozen in
    #   place AND earns the economy-wide return, so the putty-clay margin the
    #   vintage structure exists to represent does not exist.
    #
    # :full_employment — putty-clay proper.  The OLD vintage's technology is the
    #   one installed when it was new, so its capital-output ratio is fixed at the
    #   calibrated benchmark value; RR[i] = R[i,"Old"]/TR is the relative return
    #   that sustains it, and F-30 then makes old-vintage OUTPUT adjust through
    #   the disinvestment schedule XPv·CHIv = K0·RR^eta_k.  Sectors that want to
    #   contract shed old capital (RR < 1) instead of being forced to keep using
    #   K0; sectors that want to expand bid RR above 1.
    #   Why this matters: with RR pinned, a uniform +10 % TFP shock made the old
    #   vintage's forced output RISE (it has to absorb K0) while F-29
    #   (XPv[New] = XP - XPv[Old]) drove new-vintage output to its lower bound in
    #   4-5 sectors on data/TCD_2017_gtap11afr and data/ZMB_2017_gtap11afr, which
    #   capped gross output at the old vintage's forced level and destroyed those
    #   sectors (ZMB "ctl": XPv[Old] 32.5 -> 52.6, XPv[New] 32.5 -> 0, XP 65 -> 15).
    #   At the benchmark CHIv[i,"Old"] is already chi0 and RR[i] = 1 solves this
    #   equation exactly, so the benchmark is untouched.
    if full_employment
        chi0 = get(PAR[:bench], :CHIv, Dict{Any,Float64}())
        @constraint(model, F_24_old_technology[ii in i],
            (CHIv[ii,"Old"]) - (get(chi0, (ii,"Old"), 1.0)) ⟂ RR[ii])
    else
        @constraint(model, F_24_rr_bound[ii in i],
            (1.0) - (RR[ii]) ⟂ RR[ii])
    end

    # ── (F-25) Aggregate capital supply ───────────────────────────────────────
    # :full_employment — the aggregate stock is exogenous, taken from the
    #   calibrated PAR[:KSupply] (which RecursiveDynamic.jl advances between
    #   periods).  Market clearing is F-21 above, complementary to TR.
    # :fixed_wage — the previous behaviour: KS is *defined* as total capital
    #   demand, so capital is in perfectly elastic supply and TR is undetermined.
    if full_employment
        ks_total = sum(get(PAR[:KSupply], (ii,vv), 0.0) for ii in i for vv in v)
        @constraint(model, F_25, (KS) - (ks_total) ⟂ KS)
    else
        # ⟂ variable changed from KActual to KS: defines the aggregate capital supply KS
        # as the total installed capital demand. F_32 already defines KActual via dynamics.
        @constraint(model, F_25,
            (sum(Kvd[ii,vv] for ii in i for vv in v) + sum(Nfirm[ii]*KF_d[ii] for ii in i)) - (KS) ⟂ KS)
    end

    # ── (F-26) Net return on Old capital ──────────────────────────────────────
    @constraint(model, F_26[ii in i],
        (NR[ii,"Old"]) - (RR[ii] * TR) ⟂ NR[ii,"Old"])

    # ── (F-27) Net return on New capital ──────────────────────────────────────
    @constraint(model, F_27[ii in i],
        (NR[ii,"New"]) - (TR) ⟂ NR[ii,"New"])

    # ── (F-28) Gross capital return inclusive of capital tax ──────────────────
    @constraint(model, F_28[ii in i, vv in v],
        (R[ii,vv]) - ((1 + PAR[:tau_k][(ii,vv)]) * NR[ii,vv]) ⟂ R[ii,vv])

    # ── (F-29) New-vintage output = residual of total and old-vintage output ───
    # XP[ii] is determined by T_16 (CET aggregation in Trade.jl).
    # XPv[ii,"Old"] is determined by F_30 (capital supply schedule).
    # This equation closes the vintage decomposition.
    @constraint(model, F_29[ii in i],
        (XPv[ii,"New"]) - (XP[ii] - XPv[ii,"Old"]) ⟂ XPv[ii,"New"])

    # ── (F-30) Old-vintage output from installed capital ──────────────────────
    @constraint(model, F_30[ii in i],
        (XPv[ii,"Old"] * CHIv[ii,"Old"]) - (K0[ii] * RR[ii]^PAR[:eta_k]) ⟂ XPv[ii,"Old"])

    # ── (F-31) Investment quantity ────────────────────────────────────────────
    # Static-model identity: investment quantity equals the savings-determined value FD[Inv].
    # Avoids the dynamic (1+GammaInv)^nstep formulation that conflicts with C_9.
    @constraint(model, F_31, (FDInv) - (model[:FD]["Inv"]) ⟂ FDInv)

    # ── (F-32) Capital stock identity (static model) ──────────────────────────
    # KActual equals the economy-wide capital supply KS. The dynamic accumulation
    # formula is reserved for the recursive-dynamic extension.
    @constraint(model, F_32, (KActual) - (KS) ⟂ KActual)

    # ── (F-33) Capital index anchor (static model) ────────────────────────────
    # In the benchmark (static) period the normalised capital index equals 1.
    @constraint(model, F_33, (KNorm) - (1.0) ⟂ KNorm)

    # ── Fixing equations for variables without paper-numbered MCP conditions ──

    # K0[i]: baseline old-capital stock fixed at SAM-calibrated value.
    @constraint(model, F_K0[ii in i],
        (K0[ii]) - (get(PAR[:K0], ii, 100.0)) ⟂ K0[ii])

    # PS[gz]: zone price level used in minimum-wage determination; equals PABS.
    @constraint(model, F_PS[gg in gz],
        (PS[gg]) - (PABS) ⟂ PS[gg])

    # Numeraire / absolute price level.  Under :fixed_wage, W = 1 (F-12) is the
    # numeraire and PABS is simply fixed.  Under :full_employment W is endogenous
    # and one absolute condition has to take its place; PAR[:numeraire] selects
    # which (see the header and the report for the measured comparison):
    #   :pabs — keep PABS = 1.  PABS reaches the price system only through the
    #           land and natural-resource supply schedules F-13/F-18 (~3.8 % of
    #           value added, elasticity 0.5), so this anchor is weak.
    #   :cpi  — fix the consumer price index instead (CPI is the mean of the PC
    #           bundle prices and so carries the whole cost chain) and let PABS
    #           be the endogenous absolute price level it is named for, which
    #           also turns F-13/F-18 into genuine REAL supply schedules.
    numeraire = get(PAR, :numeraire, :cpi)
    if full_employment && numeraire === :cpi
        @constraint(model, F_PABS_numeraire, (model[:CPI][first(r)]) - (1.0) ⟂ PABS)
    else
        @constraint(model, F_PABS, (PABS) - (1.0) ⟂ PABS)
    end

    # GammaInv: investment growth rate; zero in the static (one-period) model.
    @constraint(model, F_GammaInv, (GammaInv) - (0.0) ⟂ GammaInv)

    # Ts[ip]: land supply for non-agricultural sectors is zero (no land demand).
    @constraint(model, F_Ts_nonag[ii in ip], (Ts[ii]) - (0.0) ⟂ Ts[ii])

    # Under :full_employment the TW[l,gs] and WMIN[l,"national"] equations are
    # already supplied by F_7_zone_wage and F_9_wmin above.
    if !full_employment
        # TW[l,gs] for segmented migration: zone threshold wage equals average zone wage.
        # For integrated migration this is handled by F_7_integrated (TW[l,gs] = TW[l,"national"]).
        @constraint(model, F_7_tw_segmented[ll in migr_segmented, gg in gs],
            (TW[ll,gg]) - (AVGW[ll,gg]) ⟂ TW[ll,gg])

        # WMIN[l,"national"] for segmented migration: national minimum wage equals average national wage.
        # For integrated migration this is handled by F_8 (WMIN[l,"national"] via chi_wmin formula).
        @constraint(model, F_8_segmented[ll in migr_segmented],
            (WMIN[ll,"national"]) - (AVGW[ll,"national"]) ⟂ WMIN[ll,"national"])
    end

    return model
end
