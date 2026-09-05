# Usage: add_closure_equations!(model, data, PAR)
# Paper-numbered LINKAGE domestic closure equations C-1--C-12.

function add_closure_equations!(model, data::LinkageData, PAR)
    S=data.sets; default_sets!(data)
    i=S[:i]; r=S[:r]; rp=S[:rp]; h=S[:h]; f=S[:f]; l=S[:l]; v=S[:v]; ins=S[:in]

    TarY=model[:TarY]; RTarY=model[:RTarY]; YG=model[:YG]; Sg=model[:Sg]; RSg=model[:RSg]
    PGDP=model[:PGDP]; PFD=model[:PFD]; FD=model[:FD]; Sf=model[:Sf]; SAV=model[:SAV]; DeprY=model[:DeprY]
    InvSh=model[:InvSh]; WRR=model[:WRR]; PNUM=model[:PNUM]; FDInv=model[:FDInv]; GDPMPr=model[:GDPMPr]
    PP=model[:PP]; PX=model[:PX]; XP=model[:XP]; PA=model[:PA]; XAp=model[:XAp]; XAc=model[:XAc]; XAf=model[:XAf]
    WPE=model[:WPE]; WPM=model[:WPM]; WTFd=model[:WTFd]; WTFout=model[:WTFout]; WTFs=model[:WTFs]; TauPR=model[:TauPR]
    NW=model[:NW]; LV=model[:LV]; Nfirm=model[:Nfirm]; LF_d=model[:LF_d]; NPT=model[:NPT]; Td=model[:Td]
    NR=model[:NR]; Kvd=model[:Kvd]; KF_d=model[:KF_d]; YH=model[:YH]; FDInvVar=model[:FDInv]
    RGDP=model[:RGDP]

    oldv = ("Old" in v) ? "Old" : first(v)
    gov = ("Gov" in f) ? "Gov" : first(f)
    inv = ("Inv" in f) ? "Inv" : last(f)
    rr0 = first(r)

    # (C-1) Nominal tariff revenue, including in- and over-quota imports.
    @constraint(model, C_1, (TarY) - (sum(WPM[rrp,rr,ii] * (PAR[:tau_m][(rrp,rr,ii)]*WTFd[rrp,rr,ii] + PAR[:tau_out][(rrp,rr,ii)]*WTFout[rrp,rr,ii]) for ii in i for rr in r for rrp in rp)) ⟂ TarY)

    # (C-2) Real tariff revenue.
    @constraint(model, C_2, (RTarY) - (TarY / PGDP[rr0]) ⟂ RTarY)

    # (C-3) Gross government revenues.
    @constraint(model, C_3, (YG) - (sum(PAR[:tau_p][ii] * (1 + PAR[:pi][ii]) * PX[ii] * XP[ii] for ii in i)
          + sum(PAR[:chi_kappa] * PAR[:kappa_h][hh] * YH[hh] for hh in h)
          + sum(PA[ii] * (sum(PAR[:tau_Ap][(ii,jj)]*XAp[ii,jj] for jj in i) + sum(PAR[:tau_Ac][(ii,hh)]*XAc[ii,hh] for hh in h) + sum(PAR[:tau_Af][(ii,ff)]*XAf[ii,ff] for ff in f)) for ii in i)
          + TarY
          + sum(PAR[:tau_e][(rr,rrp,ii)] * WPE[rr,rrp,ii] * WTFs[rr,rrp,ii] for ii in i for rr in r for rrp in rp)
          + sum(PAR[:tau_trq_share][(rrp,rr,ii)] * TauPR[rrp,rr,ii] * WPM[rrp,rr,ii] * WTFd[rrp,rr,ii] for ii in i for rr in r for rrp in rp)
          + sum(PAR[:tau_l][(ll,ii)] * NW[ll,ii] * (LV[ll,ii] + Nfirm[ii]*LF_d[ll,ii]) for ii in i for ll in l)
          + sum(PAR[:tau_t][ii] * NPT[ii] * Td[ii] for ii in i)
          + sum(PAR[:tau_k][(ii,vv)] * NR[ii,vv] * Kvd[ii,vv] for ii in i for vv in v)
          + sum(PAR[:tau_k][(ii,oldv)] * NR[ii,oldv] * Nfirm[ii] * KF_d[ii] for ii in i)) ⟂ YG)

    # (C-4) Government saving / net fiscal position.
    @constraint(model, C_4, (Sg) - (YG - PFD[gov]*FD[gov] - sum(PGDP[rr0]*PAR[:TRG][hh] for hh in h)
             + PNUM*sum(PAR[:WTRgov_in][(rrp,inn)] for rrp in rp for inn in ins)
             - PNUM*sum(PAR[:WTRgov_out][(rrp,inn)] for rrp in rp for inn in ins)) ⟂ Sg)

    # (C-5) Real government saving.
    @constraint(model, C_5, (RSg) - (Sg / PGDP[rr0]) ⟂ RSg)

    # (C-6) Government expenditure volume as share of real GDP at market prices.
    @constraint(model, C_6, (FD[gov]) - (PAR[:chi_gov] * GDPMPr) ⟂ FD[gov])

    # ── C-7 / C-9 / C-BOP: the macro closure ─────────────────────────────────
    # Walras' law makes exactly one of {savings-investment balance, balance of
    # payments, government balance} redundant: summing the household, government
    # and investment budgets with goods-market clearing (E-1) and the Armington /
    # CET duality identities gives
    #     Σ WPM·WTFd − Σ WPE·WTFs  =  PFD[Inv]·FD[Inv] − (Σ SAV + Σ DeprY + Sg)
    # identically.  C-9 and C-BOP are therefore the same restriction, and only one
    # of them may be imposed.
    #
    #   :balanced — E-2 forces world-price trade balance, Sf = PNUM·Sfbar (= 0),
    #               and C-9 is the imposed balance (savings-driven investment).
    #               No BoP equation; nothing changes from the original model.
    #
    #   :bop      — C-BOP is imposed and C-9 is dropped (it holds by Walras' law
    #               and is reported by `export_results!` as a diagnostic).
    #               Investment then needs its own rule, C-INV, which mirrors C-6
    #               for government: a fixed volume share of real GDP.
    #               `PAR[:bop_closure]` picks what clears the current account:
    #                 :flex_er  — Sf = PNUM·ER·Sfbar for every region and C-BOP ⟂ ER
    #                             (foreign saving exogenous in foreign-currency
    #                             terms, real exchange rate adjusts).
    #                 :fixed_er — ER = ER0 and C-BOP ⟂ Sf[rr0] (the home region's
    #                             current account is endogenous).
    if Symbol(get(PAR, :trade_closure, :bop)) === :balanced
        # (C-7) Foreign saving value at world numeraire price.
        @constraint(model, C_7[rr in r], (Sf[rr]) - (PNUM * PAR[:Sfbar][rr]) ⟂ Sf[rr])

        # C-8 REMOVED: all Sf[rr] are already pinned by C_7 (exogenous foreign saving per region).
        # C_8 used Sf[first(r)] as its ⟂ variable, duplicating C_7 for that region.

        # (C-9) Savings-investment balance; one region normally dropped by Walras law.
        @constraint(model, C_9, (PFD[inv]*FD[inv]) - (sum(SAV[hh] + DeprY[hh] for hh in h) + Sg + Sf[rr0]
            + PNUM*sum(PAR[:WTRinv_in][(rrp,inn)] for rrp in rp for inn in ins)
            - PNUM*sum(PAR[:WTRinv_out][(rrp,inn)] for rrp in rp for inn in ins)) ⟂ FD[inv])
    else
        ER = model[:ER]
        flex_er = Symbol(get(PAR, :bop_closure, :flex_er)) === :flex_er

        # (C-7) Domestic-currency value of exogenous (real) foreign saving.
        # Under :fixed_er the home region is skipped: C-BOP determines Sf[rr0].
        @constraint(model, C_7[rr in r; flex_er || rr != rr0],
            (Sf[rr]) - (PNUM * ER * PAR[:Sfbar][rr]) ⟂ Sf[rr])

        # (C-BOP) Balance of payments of the home region, in world prices:
        # CIF import value minus FOB export value equals foreign saving.
        bop_gap = @expression(model,
            sum(WPM[rrp,rr,ii] * WTFd[rrp,rr,ii] for ii in i for rr in r for rrp in rp)
            - sum(WPE[rr,rrp,ii] * WTFs[rr,rrp,ii] for ii in i for rr in r for rrp in rp)
            - Sf[rr0])
        if flex_er
            @constraint(model, C_BOP, (bop_gap) - (0.0) ⟂ ER)
        else
            @constraint(model, C_BOP, (bop_gap) - (0.0) ⟂ Sf[rr0])
            # (C-ER) Fixed real exchange rate.
            @constraint(model, C_ER, (ER) - (PAR[:ER0]) ⟂ ER)
        end

        # (C-9, inv_closure = :savings) Savings-driven investment: C-9 is kept and
        # determines FD[Inv]; C_BOP above is then implied by Walras' law but is kept
        # as the equation of ER (an experiment: see README "Closures").
        if Symbol(get(PAR, :inv_closure, :fixed)) === :savings
            @constraint(model, C_9, (PFD[inv]*FD[inv]) - (sum(SAV[hh] + DeprY[hh] for hh in h) + Sg + Sf[rr0]
                + PNUM*sum(PAR[:WTRinv_in][(rrp,inn)] for rrp in rp for inn in ins)
                - PNUM*sum(PAR[:WTRinv_out][(rrp,inn)] for rrp in rp for inn in ins)) ⟂ FD[inv])
        else
        # (C-INV) Exogenous real investment, replacing C-9 (see the Walras note).
        # It must be an exogenous LEVEL, not a share of real GDP: `chi_inv·GDPMPr`
        # is homogeneous of degree one in quantities, and with C-6 doing the same
        # for government it leaves the model's real scale almost unanchored — the
        # capital market is perfectly elastic under :fixed_wage (F-25) and labour
        # is not binding, so only land and natural resources pin the level.  With
        # the share form a uniform +10 pp tariff on the shipped SAM moved real GDP
        # by -49 % and the response was not monotone in the shock; with the level
        # form it is -5.8 % and monotone.  `PAR[:FDInv0]` is calibrated to the
        # SAM's investment and re-based to `chi_inv · RGDP` between periods by
        # `update_period_data!`, so it still tracks a growing economy.
        @constraint(model, C_INV, (FD[inv]) - (PAR[:FDInv0]) ⟂ FD[inv])
        end
    end

    # (C-10) Investment share of GDP at market prices.
    @constraint(model, C_10, (InvSh) - (PFD[inv] * FD[inv] / GDPMPr) ⟂ InvSh)

    # C-11 REMOVED: PNUM is already pinned to 1 by M_5 in Other.jl (numeraire).
    # Including C_11 here would give PNUM two equations simultaneously.

    # (C-12) World average rate of return to capital.
    @constraint(model, C_12, (WRR) - (sum(PAR[:TR_region][rr] * PAR[:K_region][rr] for rr in r) / sum(PAR[:K_region][rr] for rr in r)) ⟂ WRR)

    # GDPMPr: real GDP at market prices used in C_6 for government expenditure.
    # Defined as real GDP of the first region (RGDP is the regional real GDP
    # index computed in M_2 of Other.jl).
    @constraint(model, C_GDPMPr, (GDPMPr) - (RGDP[rr0]) ⟂ GDPMPr)

    return model
end
