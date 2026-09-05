# Usage:
#   add_variables!(model, data)
# All endogenous variables for the static LINKAGE MCP model.
# Dynamic (G-1..G-30), accounting (Annex H), and AIDS trade (T-28..T-34)
# variables have been removed because their equation blocks are commented out
# in ModelBuilder.jl. Unused legacy variables (XD, XM, E, D, PQ, etc.) have
# also been removed. The remaining set is exactly the set that receives
# complementarity equations in the active equation files.

function add_variables!(model, data::LinkageData)
    S=data.sets; default_sets!(data)
    i=S[:i]; j=S[:j]; k=S[:k]; r=S[:r]; rp=S[:rp]; v=S[:v]; l=S[:l]; h=S[:h]; f=S[:f]; t=S[:t]; gz=S[:gz]; e=S[:e]

    @variables(model, begin
        # ── Production / pricing ─────────────────────────────────────────────
        XPv[i,v] >= 0
        XP[i]    >= 0
        ND[i]    >= 0
        VA[i,v]  >= 0
        PND[i]   >= 0
        PVA[i,v] >= 0
        UVCv[i,v] >= 0
        UVC[i]   >= 0
        AC[i]    >= 0
        Nfirm[i] >= 0
        PX[i]    >= 0
        PP[i]    >= 0
        PROFIT[i] >= 0

        # ── Crop-sector value-added nest (P-9..P-30) ─────────────────────────
        ULD[i]      >= 0
        SLD[i]      >= 0
        HKTEF[i,v]  >= 0
        PHKTEF[i,v] >= 0
        fert[i]     >= 0
        Pfert[i]    >= 0
        HKTE[i,v]   >= 0
        PHKTE[i,v]  >= 0
        XEp[i,v]    >= 0
        PEp[i,v]    >= 0
        HKT[i,v]    >= 0
        PHKT[i,v]   >= 0
        KT[i,v]     >= 0
        PKT[i,v]    >= 0
        Kvd[i,v]    >= 0
        Td[i]       >= 0
        Fd[i]       >= 0

        # ── Intermediate demand ───────────────────────────────────────────────
        XAp[j,i]  >= 0
        PA[i]     >= 0
        PAp[j,i]  >= 0      # tax-inclusive intermediate price (P_aux_PAp_wedge)

        # ── Livestock nest extras (P-31..P-54) ───────────────────────────────
        Pfeed[i]  >= 0
        feed[i]   >= 0
        KTEL[i,v] >= 0
        PKTEL[i,v] >= 0
        TFD[i,v]  >= 0
        PTFD[i,v] >= 0

        # ── Labor market ──────────────────────────────────────────────────────
        LV[l,i]    >= 0
        W[l,i]     >= 0
        UW[i]      >= 0
        SW[i]      >= 0
        LF_d[l,i]  >= 0
        KF_d[i]    >= 0

        # ── Land / sector-specific-factor prices ─────────────────────────────
        PT[i]    >= 0
        NPT[i]   >= 0
        PF[i]    >= 0

        # ── Capital returns ───────────────────────────────────────────────────
        R[i,v]   >= 0
        NR[i,v]  >= 0

        # ── Factor-income aggregates ──────────────────────────────────────────
        TY  >= 0
        FY  >= 0
        LY[l] >= 0
        KY  >= 0

        # ── Household income / demand ─────────────────────────────────────────
        YH[h]     >= 0
        DeprY[h]  >= 0
        YD[h]     >= 0
        YC[h]     >= 0
        XH[k,h]   >= 0
        PC[k]     >= 0
        SAV[h]    >= 0
        YSTAR[h]  >= 0
        CPIH[h]   >= 0

        # ── Government / investment demand ────────────────────────────────────
        GOVDEM[i] >= 0      # government demand by good  (stub equation in Other.jl)
        INVDEM[i] >= 0      # investment demand by good  (stub equation in Other.jl)
        FD[f]     >= 0
        PFD[f]    >= 0
        XAc[i,h]  >= 0
        PAc[i,h]  >= 0
        XAf[i,f]  >= 0

        # ── Armington / domestic-import split ────────────────────────────────
        XA[i]    >= 0
        XDd[i]   >= 0
        XDs[i]   >= 0
        XMT[i]   >= 0
        PMT[i]   >= 0
        PD[i]    >= 0
        XDc[i]   >= 0
        XMc[i]   >= 0
        XDf[i,f] >= 0
        XMf[i,f] >= 0

        # ── Bilateral trade ───────────────────────────────────────────────────
        PM[r,rp,i]    >= 0
        XM1[r,i]      >= 0
        PM1[r,i]      >= 0
        XM2[r,i]      >= 0
        PM2[r,i]      >= 0
        WTFd[r,rp,i]  >= 0
        WTFs[r,rp,i]  >= 0
        WTFin[r,rp,i] >= 0
        WTFout[r,rp,i] >= 0
        WTFq[r,rp,i]  >= 0
        TauPR[r,rp,i] >= 0
        WPE[r,rp,i]   >= 0
        WPM[r,rp,i]   >= 0
        PE[r,rp,i]    >= 0      # bilateral export price; determined by E_2 (equilibrium)

        # ── Export-side CET ───────────────────────────────────────────────────
        ES[i]    >= 0
        PET[i]   >= 0

        # ── Trade-margin services ─────────────────────────────────────────────
        WXMg     >= 0
        WPMg     >= 0
        XMgr[i]  >= 0
        AXMg[r]  >= 0
        APMg[r]  >= 0

        # ── Macro aggregates ──────────────────────────────────────────────────
        GDP[r]   >= 0
        RGDP[r]  >= 0
        CPI[r]   >= 0
        PGDP[r]  >= 0
        GDPMPr   >= 0       # real GDP at market prices; defined in Closure.jl

        # ── Closure / fiscal ──────────────────────────────────────────────────
        TarY   >= 0
        RTarY  >= 0
        YG     >= 0
        Sg
        RSg
        Sf[r]
        InvSh  >= 0
        WRR    >= 0
        PNUM   >= 0

        # ── Labor-market stocks / wages ───────────────────────────────────────
        LS[l,gz]   >= 0
        MIGR[l]
        AVGW[l,gz] >= 0
        NW[l,i]    >= 0
        TW[l,gz]   >= 0
        WMIN[l,gz] >= 0
        UE[l,gz]   >= 0
        PS[gz]     >= 0
        PABS       >= 0

        # ── Land market ───────────────────────────────────────────────────────
        TLnd   >= 0
        PTLnd  >= 0
        Ts[i]  >= 0

        # ── Sector-specific-factor supply ─────────────────────────────────────
        Fs[i]  >= 0

        # ── Capital stocks / dynamics ─────────────────────────────────────────
        KSs[i]     >= 0
        KS         >= 0
        TR         >= 0
        RR[i]      >= 0
        CHIv[i,v]  >= 0
        K0[i]      >= 0
        FDInv      >= 0
        KActual    >= 0
        GammaInv
        KNorm      >= 0
    end)

    # Real exchange rate: only exists under the balance-of-payments closure.
    # Under :balanced the trade block carries no world-price variable, and
    # declaring ER would leave it without an equation (the model must stay square).
    if Symbol(get(parameters(data), :trade_closure, :bop)) !== :balanced
        @variable(model, ER >= 0)
    end

    initialize_from_sam!(model, data)
    return model
end
