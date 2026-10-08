# Usage: add_income_equations!(model, data, PAR)
# Paper-numbered LINKAGE income-distribution equations.
# IMPORTANT: These equations follow the paper equations Y-1--Y-8 exactly in structure.

function add_income_equations!(model, data::LinkageData, PAR)
    S=data.sets; default_sets!(data)
    i=S[:i]; v=S[:v]; l=S[:l]; h=S[:h]

    # JuMP variables declared in Variables.jl
    LV      = model[:LV]
    Kvd     = model[:Kvd]
    Td      = model[:Td]
    Fd      = model[:Fd]
    LF_d    = model[:LF_d]
    KF_d    = model[:KF_d]
    Nfirm   = model[:Nfirm]
    PROFIT  = model[:PROFIT]

    NPT     = model[:NPT]
    PF      = model[:PF]
    NW      = model[:NW]
    NR      = model[:NR]

    TY      = model[:TY]
    FY      = model[:FY]
    LY      = model[:LY]
    KY      = model[:KY]
    YH      = model[:YH]
    DeprY   = model[:DeprY]
    YD      = model[:YD]
    YC      = model[:YC]
    SAV     = model[:SAV]
    PNUM    = model[:PNUM]

    # Paper uses the old capital vintage in the fixed capital cost term of Y-4.
    oldv = ("Old" in v) ? "Old" : first(v)

    # (Y-1) Land remuneration: TY = sum_i NPT_i * T^d_i
    @constraint(model, Y_1, (TY) - (sum(NPT[ii] * Td[ii] for ii in i)) ⟂ TY)

    # (Y-2) Sector-specific factor remuneration: FY = sum_i PF_i * F^d_i
    @constraint(model, Y_2, (FY) - (sum(PF[ii] * Fd[ii] for ii in i)) ⟂ FY)

    # (Y-3) Labor remuneration by skill:
    # LY_l = sum_i NW_{l,i} * (LV^d_{l,i} + N_i * LF^d_{l,i})
    # Includes labor payments to the fixed-cost component under increasing returns.
    @constraint(model, Y_3[ll in l], (LY[ll]) - (sum(NW[ll,ii] * (LV[ll,ii] + Nfirm[ii] * LF_d[ll,ii]) for ii in i)) ⟂ LY[ll])

    # (Y-4) Capital remuneration:
    # KY = sum_i [ sum_v NR_{i,v} * Kv^d_{i,v} + NR_{i,Old} * N_i * KF^d_i + Π_i ]
    # Includes fixed capital cost payments and profits/markups.
    @constraint(model, Y_4, (KY) - (sum(
            sum(NR[ii,vv] * Kvd[ii,vv] for vv in v)
            + NR[ii,oldv] * Nfirm[ii] * KF_d[ii]
            + PROFIT[ii]
            for ii in i
        )) ⟂ KY)

    # (Y-5) Household income allocation across factor incomes, fiscal depreciation,
    # government transfers, and net foreign transfers.
    # Note: phi_* and TRG/WTRbar are calibrated/precomputed tables in PAR.  WTRbar is a lump
    # sum in foreign currency valued at the exchange rate (Calibration.jl convention (7); no ER
    # and WTRbar = 0 under :balanced).
    ERv = haskey(model, :ER) ? model[:ER] : 1.0
    @constraint(model, Y_5[hh in h], (YH[hh]) - (PAR[:phi_T][hh] * TY
          + PAR[:phi_F][hh] * FY
          + sum(PAR[:phi_L][(hh,ll)] * LY[ll] for ll in l)
          + PAR[:phi_K][hh] * (KY - sum(DeprY[hhh] for hhh in h))
          + PAR[:TRG][hh]
          + PNUM * ERv * PAR[:WTRbar][hh]) ⟂ YH[hh])

    # (Y-6) Fiscal depreciation allocated using capital-income shares.
    @constraint(model, Y_6[hh in h], (DeprY[hh]) - (PAR[:phi_K][hh] *
            sum(PAR[:delta_f][(ii,vv)] * NR[ii,vv] * Kvd[ii,vv] for ii in i for vv in v)) ⟂ DeprY[hh])

    # (Y-7) Disposable income after direct tax.
    @constraint(model, Y_7[hh in h], (YD[hh]) - ((1 - PAR[:chi_kappa] * PAR[:kappa_h][hh]) * YH[hh]) ⟂ YD[hh])

    # (Y-8) Income allocated by the ELES: all of disposable income.  The ELES (D-1..D-3) splits it
    # into consumption and saving, D-3 being the saving equation, so YD = Σ PC·XH + SAV.  Until
    # 2026-10-07 this was YC = YD − SAV, which with D-3 subtracted saving twice
    # (YD = Σ PC·XH + 2·SAV; Calibration.jl convention (6)).
    @constraint(model, Y_8[hh in h], (YC[hh]) - (YD[hh]) ⟂ YC[hh])

    return model
end
