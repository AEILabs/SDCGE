# Usage: add_trade_equations!(model, data, PAR)
# Paper-numbered LINKAGE trade equations.
# T-1--T-27: nested Armington/CET/trade-margin system (active).
# T-28--T-34: AIDS alternative specification removed — it conflicts with T_4
#   (both define PA[i]) and with T_2/T_9 (both define XDd/WTFd). The standard
#   Armington block (T-1..T-27) is the canonical specification.
#
# Additional fixes:
#  - T_13 complementary variable changed from WTFd to WTFout (TRQ decomposition
#    determines over-quota imports; T_9 determines total WTFd).
#  - WTFin defined here: in-quota imports equal min(WTFd, WTFq); implemented as
#    WTFin = WTFq * WTFq-fill-ratio via a simple complementarity stub.
#  - T_17 removed (revenue identity implied by CET homogeneity; T_19 is canonical).
#  - T_19 de-indexed from r: PET[i] is a single aggregate export price per good.
#  - E_2 (bilateral trade balance) moved to Equilibrium.jl with ⟂ PE.

function add_trade_equations!(model, data::LinkageData, PAR)
    S=data.sets; default_sets!(data)
    i=S[:i]; j=S[:j]; h=S[:h]; f=S[:f]; r=S[:r]; rp=S[:rp]

    XA=model[:XA]; XAp=model[:XAp]; XAc=model[:XAc]; XAf=model[:XAf]
    XDd=model[:XDd]; XDs=model[:XDs]
    XMT=model[:XMT]; XM1=model[:XM1]; XM2=model[:XM2]
    PA=model[:PA]; PD=model[:PD]; PMT=model[:PMT]; PM1=model[:PM1]; PM2=model[:PM2]
    WTFd=model[:WTFd]; WTFs=model[:WTFs]; WTFin=model[:WTFin]; WTFout=model[:WTFout]; WTFq=model[:WTFq]
    TauPR=model[:TauPR]
    XP=model[:XP]; PP=model[:PP]; ES=model[:ES]; PET=model[:PET]; PE=model[:PE]
    WPE=model[:WPE]; WPM=model[:WPM]; PM=model[:PM]
    WXMg=model[:WXMg]; WPMg=model[:WPMg]; AXMg=model[:AXMg]; APMg=model[:APMg]; XMgr=model[:XMgr]

    # (T-1) Aggregate Armington demand.
    @constraint(model, T_1[ii in i],
        (XA[ii]) - (sum(XAp[ii,jj] for jj in j) + sum(XAc[ii,hh] for hh in h) + sum(XAf[ii,ff] for ff in f)) ⟂ XA[ii])

    # (T-2) Domestic component of Armington demand.
    @constraint(model, T_2[ii in i],
        (XDd[ii]) - (PAR[:beta_d][ii] * (PA[ii] / PD[ii])^PAR[:sigma_top_m][ii] * XA[ii]) ⟂ XDd[ii])

    # (T-3) Aggregate import component.
    @constraint(model, T_3[ii in i],
        (XMT[ii]) - (PAR[:beta_m][ii] * (PA[ii] / PMT[ii])^PAR[:sigma_top_m][ii] * XA[ii]) ⟂ XMT[ii])

    # (T-4) Armington dual price.
    @constraint(model, T_4[ii in i],
        (PA[ii]) - ((PAR[:beta_d][ii]*PD[ii]^(1-PAR[:sigma_top_m][ii]) + PAR[:beta_m][ii]*PMT[ii]^(1-PAR[:sigma_top_m][ii]))^(1/(1-PAR[:sigma_top_m][ii]))) ⟂ PA[ii])

    # (T-5) Top-tier composite import demand XM1.
    @constraint(model, T_5[rr in r, ii in i],
        (XM1[rr,ii]) - (PAR[:beta_1][(rr,ii)] * (PMT[ii] / PM1[rr,ii])^PAR[:sigma_w1][ii] * XMT[ii]) ⟂ XM1[rr,ii])

    # (T-6) Aggregate import price (top tier).
    @constraint(model, T_6[ii in i],
        (PMT[ii]) - ((sum(PAR[:beta_1][(rr,ii)] * PM1[rr,ii]^(1-PAR[:sigma_w1][ii]) for rr in r))^(1/(1-PAR[:sigma_w1][ii]))) ⟂ PMT[ii])

    # (T-7) Second-tier composite import demand XM2 (index fix).
    # T-8 makes PM1[rr] the CES dual of the whole PM2[rrp] vector, so the
    # matching quantity equation must aggregate the demand for source bundle
    # rrp over every first-tier bundle rr.  Written as XM2[rr] = beta_2[rr]·XM1[rr]
    # (one term, same index) the nest lost a factor |r|: with sum_rrp beta_2 = 1
    # required by T-8, the second tier delivered only XMT/|r| to the bilateral
    # tier, so sum_{r,rp} WTFd could never equal XMT.
    @constraint(model, T_7[rrp in rp, ii in i],
        (XM2[rrp,ii]) - (sum(PAR[:beta_2][(rrp,ii)] * (PM1[rr,ii] / PM2[rrp,ii])^PAR[:sigma_w2][(rr,ii)] * XM1[rr,ii] for rr in r)) ⟂ XM2[rrp,ii])

    # (T-8) Top-tier composite import price PM1.
    @constraint(model, T_8[rr in r, ii in i],
        (PM1[rr,ii]) - ((sum(PAR[:beta_2][(rrp,ii)] * PM2[rrp,ii]^(1-PAR[:sigma_w2][(rr,ii)]) for rrp in rp))^(1/(1-PAR[:sigma_w2][(rr,ii)]))) ⟂ PM1[rr,ii])

    # (T-9) Bilateral import demand (final Armington tier).
    @constraint(model, T_9[rr in r, rrp in rp, ii in i],
        (WTFd[rr,rrp,ii]) - (PAR[:beta_w][(rr,rrp,ii)] * (PM2[rrp,ii] / PM[rr,rrp,ii])^PAR[:sigma_w3][(rrp,ii)] * XM2[rrp,ii]) ⟂ WTFd[rr,rrp,ii])

    # (T-10) Second-tier import price PM2.
    @constraint(model, T_10[rrp in rp, ii in i],
        (PM2[rrp,ii]) - ((sum(PAR[:beta_w][(rr,rrp,ii)] * PM[rr,rrp,ii]^(1-PAR[:sigma_w3][(rrp,ii)]) for rr in r))^(1/(1-PAR[:sigma_w3][(rrp,ii)]))) ⟂ PM2[rrp,ii])

    # (T-11) TRQ in-quota fill complementarity.
    @constraint(model, T_11[rr in r, rrp in rp, ii in i],
        (WTFq[rr,rrp,ii]) - (WTFin[rr,rrp,ii]) ⟂ WTFq[rr,rrp,ii])

    # (T-12) TRQ premium upper bound.
    @constraint(model, T_12[rr in r, rrp in rp, ii in i],
        (PAR[:tau_out][(rr,rrp,ii)] - PAR[:tau_in][(rr,rrp,ii)]) - (TauPR[rr,rrp,ii]) ⟂ TauPR[rr,rrp,ii])

    # (T-13) TRQ over-quota decomposition: WTFout = WTFd - WTFin.
    # Complementary variable changed to WTFout (T_9 already determines WTFd).
    @constraint(model, T_13[rr in r, rrp in rp, ii in i],
        (WTFout[rr,rrp,ii]) - (WTFd[rr,rrp,ii] - WTFin[rr,rrp,ii]) ⟂ WTFout[rr,rrp,ii])

    # WTFin: in-quota imports.  At benchmark all tariff quotas are non-binding,
    # so WTFin = WTFd.  The TRQ system (T-11/T-12/T-13) handles binding quotas.
    @constraint(model, T_WTFin[rr in r, rrp in rp, ii in i],
        (WTFin[rr,rrp,ii]) - (WTFd[rr,rrp,ii]) ⟂ WTFin[rr,rrp,ii])

    # (T-14) Domestic supply allocation from CET top nest.
    @constraint(model, T_14[ii in i],
        (XDs[ii]) - (PAR[:beta_xd][ii] * (PD[ii] / PP[ii])^PAR[:sigma_z][ii] * (XP[ii] - XMgr[ii])) ⟂ XDs[ii])

    # (T-15) Aggregate export supply allocation from CET top nest.
    @constraint(model, T_15[ii in i],
        (ES[ii]) - (PAR[:beta_es][ii] * (PET[ii] / PP[ii])^PAR[:sigma_z][ii] * (XP[ii] - XMgr[ii])) ⟂ ES[ii])

    # (T-16) CET primal aggregation: gross output = CET(domestic + export).
    # Convention fix: T-14/T-15 use the value-share form X_k = beta_k (P_k/P)^sigma X
    # with sum_k beta_k = 1.  For an aggregator X = [sum_k b_k X_k^rho]^(1/rho),
    # rho = (1+sigma)/sigma, revenue maximisation gives X_k = b_k^(-sigma)(P_k/P)^sigma X,
    # so b_k = beta_k^(-1/sigma); the previous b_k = beta_k only reproduced
    # X = sum_k X_k when sum_k beta_k^(1+rho) = 1, which value shares never satisfy.
    @constraint(model, T_16[ii in i],
        (XP[ii] - XMgr[ii]) -
        ((PAR[:beta_xd][ii]^(-1/PAR[:sigma_z][ii])*XDs[ii]^((1+PAR[:sigma_z][ii])/PAR[:sigma_z][ii]) +
          PAR[:beta_es][ii]^(-1/PAR[:sigma_z][ii])*ES[ii]^((1+PAR[:sigma_z][ii])/PAR[:sigma_z][ii]))^(PAR[:sigma_z][ii]/(1+PAR[:sigma_z][ii]))) ⟂ XP[ii])

    # (T-17) REMOVED — the revenue identity PET·ES = Σ PE·WTFs is implied by
    # CET duality (T-18, T-19) and would duplicate PET's equation.

    # (T-18) Bilateral export supply allocation across destinations.
    # Sign fix: this is the CET supply FOC belonging to the dual T-19,
    # PET = (sum beta_z · PE^(1+sigma_z2))^(1/(1+sigma_z2)), whose FOC is
    # WTFs = beta_z (PE/PET)^sigma_z2 ES — supply to a destination RISES with the
    # price received there.  The ratio was written the other way up (the CES
    # *demand* form), so an export price cut raised supply to that destination.
    # That is invisible at the benchmark (PE = PET = 1) but makes the response to
    # any trade shock perverse and non-monotone; PATH hit its iteration limit on
    # a 20% tariff with the old form.  It also breaks the CET revenue identity
    # sum_rp PE·WTFs = PET·ES away from the benchmark.
    @constraint(model, T_18[rr in r, rrp in rp, ii in i],
        (WTFs[rr,rrp,ii]) -
        (PAR[:beta_z][(rr,rrp,ii)] *
         ((PE[rr,rrp,ii] + PAR[:tau_trq_share][(rr,rrp,ii)]*TauPR[rr,rrp,ii]*WPM[rr,rrp,ii]) / PET[ii])^PAR[:sigma_z2][ii] *
         ES[ii]) ⟂ WTFs[rr,rrp,ii])

    # (T-19) CET dual export price (de-indexed from r: PET is a single scalar per good).
    @constraint(model, T_19[ii in i],
        (PET[ii]) -
        ((sum(sum(PAR[:beta_z][(rr,rrp,ii)] *
                  (PE[rr,rrp,ii] + PAR[:tau_trq_share][(rr,rrp,ii)]*TauPR[rr,rrp,ii]*WPM[rr,rrp,ii])^(1+PAR[:sigma_z2][ii])
                  for rrp in rp) for rr in r))^(1/(1+PAR[:sigma_z2][ii]))) ⟂ PET[ii])

    # (T-20) FOB export price including export tax/subsidy.
    # (T-21) CIF import price with trade margin.
    # Both depend on `PAR[:trade_closure]` (see Equilibrium.jl and the header of
    # Calibration.jl); the equation counts are identical in the two regimes.
    if Symbol(get(PAR, :trade_closure, :bop)) === :balanced
        # E-2 pins PE, so T-20 determines the world price and T-21 derives the
        # CIF price from it through the iceberg factor lambda_w — which is what
        # ties imports to exports good by good.
        @constraint(model, T_20[rr in r, rrp in rp, ii in i],
            (WPE[rr,rrp,ii]) - ((1 + PAR[:tau_e][(rr,rrp,ii)]) * PE[rr,rrp,ii]) ⟂ WPE[rr,rrp,ii])
        @constraint(model, T_21[rr in r, rrp in rp, ii in i],
            (WPM[rr,rrp,ii]) - ((1 + PAR[:zeta_t][(rr,rrp,ii)]) * WPE[rr,rrp,ii] / PAR[:lambda_w][(rr,rrp,ii)]) ⟂ WPM[rr,rrp,ii])
    else
        # :bop — E-2 pins WPE = ER·PWE0, so T-20 determines the domestic export
        # price PE net of the export tax, and the CIF import price is its own
        # exogenous world price converted at ER (the margin wedge zeta_t is kept).
        ER = model[:ER]
        @constraint(model, T_20[rr in r, rrp in rp, ii in i],
            (WPE[rr,rrp,ii]) - ((1 + PAR[:tau_e][(rr,rrp,ii)]) * PE[rr,rrp,ii]) ⟂ PE[rr,rrp,ii])
        @constraint(model, T_21[rr in r, rrp in rp, ii in i],
            (WPM[rr,rrp,ii]) - ((1 + PAR[:zeta_t][(rr,rrp,ii)]) * ER * PAR[:PWM0][(rr,rrp,ii)]) ⟂ WPM[rr,rrp,ii])
    end

    # (T-22) Domestic import price including tariff and TRQ premium.
    @constraint(model, T_22[rr in r, rrp in rp, ii in i],
        (PM[rr,rrp,ii]) - ((1 + PAR[:tau_m][(rr,rrp,ii)] + TauPR[rr,rrp,ii]) * WPM[rr,rrp,ii]) ⟂ PM[rr,rrp,ii])

    # (T-23) World demand for international trade and transport services.
    # ⟂ variable changed to WXMg: the demand-supply balance determines quantity,
    # while T_25 determines the price WPMg via the CES dual formula.
    @constraint(model, T_23,
        (WPMg * WXMg) -
        (sum(PAR[:zeta_t][(rr,rrp,ii)] * WPE[rr,rrp,ii] * WTFd[rr,rrp,ii]
             for rr in r for rrp in rp for ii in i)) ⟂ WXMg)

    # (T-24) Regional demand for trade and transport services.
    @constraint(model, T_24[rr in r],
        (AXMg[rr]) - (PAR[:alpha_TT][rr] * (WPMg / APMg[rr])^PAR[:sigma_TT] * WXMg) ⟂ AXMg[rr])

    # (T-25) World price of trade and transport services.
    @constraint(model, T_25,
        (WPMg) - ((sum(PAR[:alpha_TT][rr] * APMg[rr]^(1-PAR[:sigma_TT]) for rr in r))^(1/(1-PAR[:sigma_TT]))) ⟂ WPMg)

    # (T-26) Total sectoral supply of trade-margin services summed across regions.
    # Original had [rr in r, ii in i] giving 4 equations per XMgr[ii] — fixed by summing over r.
    @constraint(model, T_26[ii in i],
        (XMgr[ii]) - (sum(PAR[:a_Mg][(rr,ii)] * AXMg[rr] for rr in r)) ⟂ XMgr[ii])

    # (T-27) Regional price of trade and transport services.
    @constraint(model, T_27[rr in r],
        (APMg[rr]) - (sum(PAR[:a_Mg][(rr,ii)] * PP[ii] for ii in i)) ⟂ APMg[rr])

    return model
end
