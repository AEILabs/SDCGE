# Usage: add_equilibrium_equations!(model, data, PAR)
# Paper-numbered LINKAGE goods-market equilibrium equations.
#
# E-2 has two forms, selected by `PAR[:trade_closure]` (see the header of
# Calibration.jl and the "Closures" section of README.md).  Both occupy the same
# |r|·|rp|·|i| slots, so the model stays square either way:
#
#   :balanced — the legacy bilateral trade-balance condition
#               `WTFd = lambda_w · WTFs`, complementary to the bilateral export
#               price PE.  Together with T-21 (`WPM = WPE/lambda_w`) it makes the
#               CIF value of imports identically equal to the FOB value of
#               exports for every good, so the benchmark must be trade-balanced.
#
#   :bop      — small open economy.  The world FOB export price is exogenous in
#               foreign currency and converted with the real exchange rate:
#               `WPE = ER · PWE0`, complementary to WPE.  T-20 then determines
#               the domestic export price PE = WPE/(1+tau_e) (its complementary
#               variable moves from WPE to PE in Trade.jl), so an export tax
#               lowers what the exporter receives instead of raising the world
#               price.  Import demand (T-9) and export supply (T-18) become
#               independent, and the current account is closed by C-BOP.

function add_equilibrium_equations!(model, data::LinkageData, PAR)
    S=data.sets; default_sets!(data)
    i=S[:i]; r=S[:r]; rp=S[:rp]
    XDs=model[:XDs]; XDd=model[:XDd]; PD=model[:PD]
    WTFd=model[:WTFd]; WTFs=model[:WTFs]; PE=model[:PE]; WPE=model[:WPE]

    # (E-1) Domestic market equilibrium: supply equals demand; domestic price PD clears market.
    # Changed ⟂ variable from XDs to PD: T_14 already defines XDs via CET supply allocation.
    @constraint(model, E_1[ii in i], (XDs[ii]) - (XDd[ii]) ⟂ PD[ii])

    if Symbol(get(PAR, :trade_closure, :bop)) === :balanced
        # (E-2, :balanced) Bilateral trade-flow balance: imports (demand side)
        # equal exports adjusted for iceberg transport loss.  The equilibrium
        # bilateral export price PE[rr,rrp,ii] adjusts to clear each bilateral
        # market.
        @constraint(model, E_2[rr in r, rrp in rp, ii in i],
            (WTFd[rr,rrp,ii]) - (PAR[:lambda_w][(rr,rrp,ii)] * WTFs[rr,rrp,ii]) ⟂ PE[rr,rrp,ii])
    else
        # (E-2, :bop) Exogenous world FOB export price in domestic currency.
        ER = model[:ER]
        @constraint(model, E_2[rr in r, rrp in rp, ii in i],
            (WPE[rr,rrp,ii]) - (ER * PAR[:PWE0][(rr,rrp,ii)]) ⟂ WPE[rr,rrp,ii])
    end

    return model
end
