# Initialization helpers for the static LINKAGE MCP model.
#
# Every start value comes from `parameters(data)[:bench]`, the benchmark table
# built by `calibrate_from_sam!` from the balanced SAM.  There are no hard-coded
# fractions here: a start value and the parameter that determines it always come
# from the same calibration, so each equation holds exactly at the start point.
# See the header of `Calibration.jl` for the price normalisation and for the
# three benchmark conventions (ag-only land, balanced trade, zero subsistence).

function _safe_start_value!(model, name::Symbol, index_tuple::Tuple, value)
    haskey(model, name) || return nothing
    container = model[name]
    try
        var = isempty(index_tuple) ? container : container[index_tuple...]
        set_start_value(var, max(float(value), 1.0e-9))
    catch
    end
    return nothing
end

function _safe_start_value_raw!(model, name::Symbol, index_tuple::Tuple, value)
    haskey(model, name) || return nothing
    container = model[name]
    try
        var = isempty(index_tuple) ? container : container[index_tuple...]
        set_start_value(var, float(value))
    catch
    end
    return nothing
end

# Lower bound substituted for a zero lower bound, and the offset used when a
# start value would otherwise sit on that bound.  Equations that sum thousands of
# structurally-zero quantities (C-3 over ~10^4 XAp cells, Y-3, Y-4, T-1, F-25)
# inherit `count · eps` as their start residual, which is why the max start
# residual is ~8.6e-6 rather than ~1e-10.  Tightening this to 1e-10 was tried:
# it cuts the max start residual to 5.5e-7 but makes the DEFAULT (fixed-wage)
# benchmark drift worse (max |pct change| 2.7e-4 % -> 0.23 %, because the
# near-singular wage direction is excited differently), and it does not fix the
# full-employment stall.  Kept at 1e-8.
const LCGE_START_EPS = 1.0e-8
# Start value for variables the benchmark forces to exactly zero: the smallest
# value that still passes `check_initialization!` (start strictly above the
# lower bound).
const LCGE_ZERO_START = 1.0e-12
# UE start: F-10 says UE·LS = LS - LV, and at the benchmark labour demand equals
# labour supply, so UE = 0.  Starting at 1e-6 left an F-10 residual of ~2e-2.
const LCGE_UE_START  = LCGE_ZERO_START
const LCGE_UE_MAX    = 0.95
# F-24 forces RR = 1 (full capacity utilisation) at the benchmark; RR's upper
# bound is 1, and a variable sitting on its bound is a valid MCP solution.
const LCGE_RR_START  = 1.0

function _is_finite_number(x)
    return x isa Real && isfinite(float(x))
end

function enforce_nlp_safe_bounds_and_starts!(model; eps::Float64=LCGE_START_EPS)
    # Move nonneg lower bounds slightly inside the feasible region and repair
    # missing/non-finite start values.
    for var in all_variables(model)
        try
            if has_lower_bound(var) && lower_bound(var) == 0.0
                set_lower_bound(var, eps)
            end
        catch
        end
        sv = try start_value(var) catch; nothing end
        if sv === nothing || !(sv isa Real) || !isfinite(float(sv))
            try
                if has_lower_bound(var)
                    lb = lower_bound(var)
                    set_start_value(var, max(float(lb) + 1.0, eps))
                else
                    set_start_value(var, 0.0)
                end
            catch
            end
        elseif has_lower_bound(var) && float(sv) <= lower_bound(var)
            try
                set_start_value(var, lower_bound(var) + eps)
            catch
            end
        end
    end

    # Unemployment rates: UE = 0 exactly at the benchmark, and UE only ever
    # appears as (1-UE) or (1-UE)^0, so a zero lower bound is numerically safe.
    # A start set earlier by `initialize_from_sam!` (the dynamic path, where
    # labour supply has grown away from labour demand) is kept.
    if haskey(model, :UE)
        UE = model[:UE]
        try
            for key in eachindex(UE)
                set_lower_bound(UE[key], 0.0)
                set_upper_bound(UE[key], LCGE_UE_MAX)
                sv = try start_value(UE[key]) catch; nothing end
                if !(sv isa Real) || !isfinite(float(sv)) || float(sv) <= 0.0
                    set_start_value(UE[key], LCGE_UE_START)
                end
            end
        catch
        end
    end

    # TauPR (tariff-rate-quota premium): T-12 is `(tau_out - tau_in) - TauPR ⟂ TauPR`,
    # so with no TRQ in the data the MCP requires TauPR = 0 exactly.  A strictly
    # positive lower bound makes that impossible — at the bound the residual is
    # -TauPR < 0, which violates the complementarity condition — and PATH stalls
    # on T_12 with SLOW_PROGRESS at a residual of ~1e-6.  TauPR only ever enters
    # additively (T-18/T-19 base `PE + ...`, T-22 `1 + tau_m + TauPR`, C-1, C-3),
    # so a zero lower bound is numerically safe.
    if haskey(model, :TauPR)
        TP = model[:TauPR]
        try
            for key in eachindex(TP)
                set_lower_bound(TP[key], 0.0)
                set_start_value(TP[key], LCGE_ZERO_START)
            end
        catch
        end
    end

    # RR is a utilisation ratio in (0, 1]. F_24 forces RR = 1 at the benchmark.
    # (Dropping the redundant upper bound was tried as a fix for the
    # full-employment stall — F_24 `1 - RR ⟂ RR` has bound and equation active
    # together — and made no measurable difference; the bound is kept because it
    # is economically meaningful.)
    if haskey(model, :RR)
        RR = model[:RR]
        try
            for key in eachindex(RR)
                set_lower_bound(RR[key], eps)
                set_upper_bound(RR[key], 1.0)
                set_start_value(RR[key], LCGE_RR_START)
            end
        catch
        end
    end

    # GammaInv appears inside (1+GammaInv)^nstep; keep the base positive.
    # F_GammaInv fixes it to 0 in the static model, so that is also its start.
    if haskey(model, :GammaInv)
        try
            set_lower_bound(model[:GammaInv], -0.95)
            set_start_value(model[:GammaInv], 0.0)
        catch
        end
    end

    return model
end

function _maybe_lower_bound(var)
    try; return has_lower_bound(var) ? lower_bound(var) : missing; catch; return missing; end
end
function _maybe_upper_bound(var)
    try; return has_upper_bound(var) ? upper_bound(var) : missing; catch; return missing; end
end
function _maybe_start_value(var)
    try; return start_value(var); catch; return nothing; end
end

function initialization_diagnostics(model; max_items::Int=25)
    rows = NamedTuple[]
    n_missing = 0; n_bad = 0; n_at_or_below_lb = 0
    for var in all_variables(model)
        nm  = name(var)
        sv  = _maybe_start_value(var)
        bad = sv === nothing || !(sv isa Real) || !isfinite(float(sv))
        lb  = _maybe_lower_bound(var)
        ub  = _maybe_upper_bound(var)
        atlb= lb !== missing && sv !== nothing && sv isa Real && float(sv) <= float(lb)
        n_missing += sv === nothing ? 1 : 0
        n_bad     += bad ? 1 : 0
        n_at_or_below_lb += atlb ? 1 : 0
        if (bad || atlb) && length(rows) < max_items
            push!(rows, (variable=nm, start=sv, lower=lb, upper=ub,
                         issue=bad ? "missing_or_nonfinite_start" : "start_at_or_below_lower_bound"))
        end
    end
    return (missing_starts=n_missing, bad_starts=n_bad,
            starts_at_or_below_lower_bound=n_at_or_below_lb, examples=rows)
end

function check_initialization!(model; error_on_bad::Bool=true)
    d = initialization_diagnostics(model)
    if error_on_bad && (d.bad_starts > 0 || d.starts_at_or_below_lower_bound > 0)
        error("Bad model initialization: $(d). Call initialization_diagnostics(model) for details.")
    end
    return d
end

function initialize_from_sam!(model, data::LinkageData)
    default_sets!(data)
    S  = data.sets
    i  = S[:i]; j = S[:j]; k = S[:k]; r = S[:r]; rp = S[:rp]
    v  = S[:v]; l = S[:l]; h = S[:h]; f = S[:f]
    gz = S[:gz]

    PAR = parameters(data)
    haskey(PAR, :bench) || error("Benchmark table missing: run calibrate_from_sam! before initializing.")
    B = PAR[:bench]
    b1(key, idx, dflt=1.0) = get(B[key], idx, dflt)

    nr = length(r)

    # ── Labour-market closure (see the header of Factors.jl) ─────────────────
    # Under :full_employment F-6 forces sum_i LV[l,i] = LS[l,"national"]·(1-UE0),
    # so labour demand tracks labour SUPPLY and the LV/ULD/SLD/LY starts have to
    # be rescaled to it whenever LS0 has been grown by the dynamic drivers.
    # The scale is exactly 1 at the benchmark (LV0 sums to LSupply, UE0 = 0), so
    # the benchmark start point is unchanged.  Under :fixed_wage labour demand is
    # unconstrained and the starts stay at the last realised demand.
    full_employment = get(PAR, :labour_closure, :full_employment) === :full_employment
    ue0    = get(PAR, :UE0, Dict{Any,Float64}())
    lscale = Dict{Any,Float64}()
    for ll in l
        ld = sum(get(PAR[:LV0], (ll,ii), 0.0) for ii in i)
        ls = get(PAR[:LS0], (ll,"national"), ld) * (1 - get(ue0, ll, 0.0))
        lscale[ll] = (full_employment && ld > 1.0e-9 && ls > 0.0) ? ls/ld : 1.0
    end

    # ── Prices at the benchmark normalisation ────────────────────────────────
    # Unit prices everywhere except the tax wedges PX = 1/(1+tau_p) and
    # PND = PEp = Pfert = Pfeed = PAp = 1 + tau_Ap (see Calibration.jl header).
    for nm in [:PP, :PA, :PD, :PMT, :PET, :PT, :NPT, :PF, :PC, :CPI, :PGDP,
               :PVA, :PHKTEF, :PHKTE, :PHKT, :PKT, :PKTEL, :PTFD, :R, :NR]
        haskey(model, nm) || continue
        obj = model[nm]
        try; for key in eachindex(obj); set_start_value(obj[key], 1.0); end; catch; end
    end
    for nm in [:PNUM, :PABS, :PTLnd, :WPMg, :WRR, :TR, :KNorm]
        _safe_start_value!(model, nm, (), 1.0)
    end

    # ── Sector-level quantities and prices ───────────────────────────────────
    for ii in i
        _safe_start_value!(model, :XP,   (ii,), B[:XP][ii])
        _safe_start_value!(model, :XA,   (ii,), B[:XA][ii])
        _safe_start_value!(model, :XDs,  (ii,), B[:XDs][ii])
        _safe_start_value!(model, :XDd,  (ii,), B[:XDd][ii])
        _safe_start_value!(model, :XMT,  (ii,), B[:XMT][ii])
        _safe_start_value!(model, :ES,   (ii,), B[:ES][ii])
        _safe_start_value!(model, :ND,   (ii,), B[:ND][ii])
        _safe_start_value!(model, :PX,   (ii,), B[:PX][ii])
        _safe_start_value!(model, :UVC,  (ii,), B[:PX][ii])
        _safe_start_value!(model, :AC,   (ii,), B[:PX][ii])
        _safe_start_value!(model, :PND,  (ii,), B[:PND][ii])
        _safe_start_value!(model, :Pfert,(ii,), B[:Pfert][ii])
        _safe_start_value!(model, :Pfeed,(ii,), B[:Pfeed][ii])
        _safe_start_value!(model, :Nfirm,(ii,), 1.0)
        # P-8: PROFIT = XP·(PX - AC) = 0 at the competitive benchmark.
        _safe_start_value!(model, :PROFIT, (ii,), 0.0)
        _safe_start_value!(model, :ULD,  (ii,), B[:ULD][ii] * get(lscale, "UnSkLab", 1.0))
        _safe_start_value!(model, :SLD,  (ii,), B[:SLD][ii] * get(lscale, "SkLab",   1.0))
        _safe_start_value!(model, :fert, (ii,), B[:fert][ii])
        _safe_start_value!(model, :feed, (ii,), B[:feed][ii])
        # Factor quantities are read from the live PAR supply tables rather than
        # from :bench, so that RecursiveDynamic/PolicyScenarios (which rescale
        # KSupply/TSupply/FSupply/LSupply between periods) still get start values
        # consistent with the period they are solving.  At the benchmark the two
        # are identical by construction.
        _safe_start_value!(model, :Td,   (ii,), get(PAR[:TSupply], ii, B[:Td][ii]))
        _safe_start_value!(model, :Ts,   (ii,), get(PAR[:TSupply], ii, B[:Ts][ii]))
        _safe_start_value!(model, :Fd,   (ii,), get(PAR[:FSupply], ii, B[:Fd][ii]))
        _safe_start_value!(model, :Fs,   (ii,), get(PAR[:FSupply], ii, B[:Fs][ii]))
        _safe_start_value!(model, :K0,   (ii,), get(PAR[:K0], ii, 1.0))
        _safe_start_value!(model, :KSs,  (ii,),
            sum(get(PAR[:KSupply], (ii,vv), 0.0) for vv in v))
        _safe_start_value!(model, :KF_d, (ii,), 0.0)
        _safe_start_value!(model, :GOVDEM, (ii,), B[:G][ii])
        _safe_start_value!(model, :INVDEM, (ii,), B[:I][ii])
        # T-23/T-24/T-26: zeta_t = 0, so international transport margins and the
        # sectoral margin supply XMgr are zero at the benchmark.
        _safe_start_value!(model, :XMgr, (ii,), 0.0)

        for vv in v
            _safe_start_value!(model, :XPv,  (ii,vv), B[:XPv][(ii,vv)])
            _safe_start_value!(model, :VA,   (ii,vv), B[:VA][(ii,vv)])
            _safe_start_value!(model, :UVCv, (ii,vv), B[:PX][ii])
            _safe_start_value!(model, :Kvd,  (ii,vv), get(PAR[:KSupply], (ii,vv), B[:Kvd][(ii,vv)]))
            _safe_start_value!(model, :KT,   (ii,vv), B[:KT][(ii,vv)])
            _safe_start_value!(model, :HKT,  (ii,vv), B[:HKT][(ii,vv)])
            _safe_start_value!(model, :HKTE, (ii,vv), B[:HKTE][(ii,vv)])
            _safe_start_value!(model, :HKTEF,(ii,vv), B[:HKTEF][(ii,vv)])
            _safe_start_value!(model, :KTEL, (ii,vv), B[:KTEL][(ii,vv)])
            _safe_start_value!(model, :TFD,  (ii,vv), B[:TFD][(ii,vv)])
            _safe_start_value!(model, :XEp,  (ii,vv), B[:XEp][(ii,vv)])
            _safe_start_value!(model, :PEp,  (ii,vv), B[:PEp][ii])
            _safe_start_value!(model, :CHIv, (ii,vv), B[:CHIv][(ii,vv)])
        end

        # P-72/P-74: with one labour type per bundle, LV = ULD / SLD.
        for ll in l
            _safe_start_value!(model, :LF_d, (ll,ii), 0.0)
            _safe_start_value!(model, :LV,   (ll,ii),
                get(PAR[:LV0], (ll,ii), 1.0) * get(lscale, ll, 1.0))
            _safe_start_value!(model, :W,    (ll,ii), 1.0)
            _safe_start_value!(model, :NW,   (ll,ii), 1.0)
        end
        _safe_start_value!(model, :UW, (ii,), 1.0)
        _safe_start_value!(model, :SW, (ii,), 1.0)
    end

    # ── Intermediate demand: XAp = SAM input-output cell, PAp = (1+tau_Ap)·PA ──
    for jj in j, ii in i
        _safe_start_value!(model, :XAp, (jj,ii), get(B[:XAp], (jj,ii), 0.0))
        _safe_start_value!(model, :PAp, (jj,ii), B[:PND][ii])
    end

    # ── Household income and demand (all from the calibrated benchmark) ──────
    for hh in h
        _safe_start_value!(model, :YH,    (hh,), B[:YH])
        _safe_start_value!(model, :DeprY, (hh,), B[:DeprY])
        _safe_start_value!(model, :YD,    (hh,), B[:YD])
        _safe_start_value!(model, :YC,    (hh,), B[:YC])
        _safe_start_value!(model, :SAV,   (hh,), B[:SAV])
        _safe_start_value!(model, :YSTAR, (hh,), B[:YSTAR])
        _safe_start_value!(model, :CPIH,  (hh,), 1.0)
        for ii in i
            xac = B[:C][ii]
            _safe_start_value!(model, :XAc, (ii,hh), xac)
            _safe_start_value!(model, :PAc, (ii,hh), 1.0)
            # D-10/D-11: household Armington split at the calibrated shares.
            _safe_start_value!(model, :XDc, (ii,), get(PAR[:beta_d], ii, 0.9) * xac)
            _safe_start_value!(model, :XMc, (ii,), get(PAR[:beta_m], ii, 0.1) * xac)
        end
    end
    # D-2 with theta = 0: XH[k,h] = mu_c[k,h]·YC = benchmark consumption of k.
    for kk in k, hh in h
        _safe_start_value!(model, :XH, (kk,hh), B[:C][kk])
    end
    _safe_start_value!(model, :TY, (), B[:TY])
    _safe_start_value!(model, :FY, (), B[:FY])
    _safe_start_value!(model, :KY, (), B[:KY])
    # Y-3: LY = sum_i NW·LV, so labour income follows labour DEMAND, not supply.
    # Identical at the benchmark (LV0 sums to LSupply there), but in a dynamic
    # period labour supply has been grown while demand has not.
    for ll in l
        ld = sum(get(PAR[:LV0], (ll,ii), 0.0) for ii in i) * get(lscale, ll, 1.0)
        _safe_start_value!(model, :LY, (ll,), ld > 0 ? ld : B[:LY][ll])
    end

    # ── Other final demand (Gov, Inv) ────────────────────────────────────────
    fd_start = Dict("Gov" => B[:GOV], "Inv" => B[:INV])
    fd_qty   = Dict("Gov" => B[:G],   "Inv" => B[:I])
    for ff in f
        _safe_start_value!(model, :FD,  (ff,), get(fd_start, ff, B[:INV]))
        _safe_start_value!(model, :PFD, (ff,), 1.0)
        qty = get(fd_qty, ff, B[:I])
        for ii in i
            xafv = qty[ii]
            _safe_start_value!(model, :XAf, (ii,ff), xafv)
            _safe_start_value!(model, :XDf, (ii,ff), get(PAR[:beta_d], ii, 0.9) * xafv)
            _safe_start_value!(model, :XMf, (ii,ff), get(PAR[:beta_m], ii, 0.1) * xafv)
        end
    end

    # ── Labour market ────────────────────────────────────────────────────────
    # F-10 gives UE[l,"national"]·LS = LS − Σ_i LV, so the start follows from
    # the (exogenous) labour supply and the benchmark labour demand.  At the
    # benchmark the two are equal and UE starts at ~0; in a dynamic period where
    # LS0 has been grown by RecursiveDynamic/PolicyScenarios this puts UE close
    # to its solution instead of leaving it pinned to its lower bound, which
    # otherwise forces PATH through a basis change it handles badly.
    for ll in l
        _safe_start_value_raw!(model, :MIGR, (ll,), 0.0)
        if full_employment
            # F-10 fixes UE at its benchmark rate in this regime.
            ue_nat = max(get(ue0, ll, 0.0), LCGE_ZERO_START)
            ue_gs  = ue_nat
        else
            # F-10 makes UE absorb supply minus demand; start it there so it is
            # not pinned to its lower bound when LS0 has been grown.
            ls_nat = max(get(PAR[:LS0], (ll,"national"), 1.0), 1.0e-9)
            ld_nat = sum(get(PAR[:LV0], (ll,ii), 0.0) for ii in i)
            ue_nat = clamp(1.0 - ld_nat/ls_nat, LCGE_UE_START, LCGE_UE_MAX)
            ue_gs  = LCGE_UE_START
        end
        for gg in gz
            _safe_start_value!(model, :LS,   (ll,gg), get(PAR[:LS0], (ll,gg), 1.0))
            _safe_start_value!(model, :AVGW, (ll,gg), 1.0)
            _safe_start_value!(model, :TW,   (ll,gg), 1.0)
            _safe_start_value!(model, :WMIN, (ll,gg), 1.0)
            _safe_start_value_raw!(model, :UE, (ll,gg),
                gg == "national" ? ue_nat : ue_gs)
        end
    end

    # ── Bilateral trade ──────────────────────────────────────────────────────
    # Every bilateral quantity is the top-down split of the aggregate benchmark
    # (XMT on the import side, ES on the export side) by the calibrated shares,
    # so T-5..T-10, T-18/T-19 and E-2 all hold exactly.  Border prices carry the
    # export tax and the iceberg factor lambda_w = (1+tau_m)(1+tau_e), which
    # makes WPM = 1/(1+tau_m) and therefore PM = PM2 = PM1 = PMT = 1 (T-22).
    for ii in i
        xmt0 = B[:XMT][ii]
        es0  = B[:ES][ii]
        # T-5: XM1[r] = beta_1[r]·XMT (sum_r beta_1 = 1).
        xm1 = Dict(rr => get(PAR[:beta_1], (rr,ii), 1.0/nr) * xmt0 for rr in r)
        xm1tot = sum(values(xm1))
        # T-7: XM2[rp] = beta_2[rp]·sum_r XM1[r]  (sum_rp beta_2 = 1).
        xm2 = Dict(rr => get(PAR[:beta_2], (rr,ii), 1.0/nr) * xm1tot for rr in r)
        for rr in r
            _safe_start_value!(model, :XM1, (rr,ii), xm1[rr])
            _safe_start_value!(model, :PM1, (rr,ii), 1.0)
            _safe_start_value!(model, :XM2, (rr,ii), xm2[rr])
            _safe_start_value!(model, :PM2, (rr,ii), 1.0)
        end
        for rr in r, rrp in rp
            # T-9: WTFd[r,rp] = beta_w[r,rp]·XM2[rp]  (sum_r beta_w = 1).
            wtfd_val = get(PAR[:beta_w], (rr,rrp,ii), 1.0/nr) * xm2[rrp]
            _safe_start_value!(model, :WTFd,   (rr,rrp,ii), wtfd_val)
            _safe_start_value!(model, :WTFin,  (rr,rrp,ii), wtfd_val)
            _safe_start_value!(model, :WTFq,   (rr,rrp,ii), wtfd_val)
            _safe_start_value!(model, :WTFout, (rr,rrp,ii), 0.0)
            _safe_start_value!(model, :TauPR,  (rr,rrp,ii), 0.0)
            # T-18: WTFs[r,rp] = beta_z[r,rp]·ES  (sum over the whole grid = 1).
            _safe_start_value!(model, :WTFs, (rr,rrp,ii),
                get(PAR[:beta_z], (rr,rrp,ii), 1.0/(nr*nr)) * es0)
            _safe_start_value!(model, :PM,  (rr,rrp,ii), 1.0)
            _safe_start_value!(model, :PE,  (rr,rrp,ii), 1.0)
            _safe_start_value!(model, :WPE, (rr,rrp,ii), B[:WPE][ii])
            _safe_start_value!(model, :WPM, (rr,rrp,ii), B[:WPM][ii])
        end
    end

    # ── International trade and transport margins (zero at the benchmark) ────
    for rr in r
        _safe_start_value!(model, :AXMg, (rr,), 0.0)
        _safe_start_value!(model, :APMg, (rr,), 1.0)
    end
    _safe_start_value!(model, :WXMg, (), 0.0)
    _safe_start_value!(model, :WPMg, (), 1.0)

    # ── Land market ──────────────────────────────────────────────────────────
    # F-13: TLnd = chi_T[:land], calibrated as total agricultural land (land
    # outside S[:ag] is re-assigned to capital, see Calibration.jl).
    _safe_start_value!(model, :TLnd,  (), get(PAR[:chi_T], :land, B[:TLnd]))
    _safe_start_value!(model, :PTLnd, (), 1.0)

    # ── Capital aggregates ───────────────────────────────────────────────────
    ks_start = sum(get(PAR[:KSupply], (ii,vv), 0.0) for ii in i for vv in v)
    ks_start = ks_start > 0 ? ks_start : B[:KS]
    _safe_start_value!(model, :KS,      (), ks_start)
    _safe_start_value!(model, :KActual, (), ks_start)
    _safe_start_value!(model, :FDInv,   (), B[:INV])

    # ── Macro and fiscal closure ─────────────────────────────────────────────
    _safe_start_value!(model, :GDPMPr, (), B[:GDP])
    _safe_start_value!(model, :YG,     (), B[:YG])
    _safe_start_value!(model, :TarY,   (), B[:TarY])
    _safe_start_value!(model, :RTarY,  (), B[:TarY])
    _safe_start_value!(model, :InvSh,  (), B[:InvSh])
    _safe_start_value_raw!(model, :GammaInv, (), 0.0)
    _safe_start_value_raw!(model, :Sg,  (), B[:Sg])
    _safe_start_value_raw!(model, :RSg, (), B[:Sg])
    # Real exchange rate (only present under the :bop trade closure) and foreign
    # saving.  Sfbar is 0 under :balanced, so this is the old `Sf = 0` there.
    _safe_start_value!(model, :ER, (), get(B, :ER, get(PAR, :ER0, 1.0)))
    for rr in r
        _safe_start_value_raw!(model, :Sf, (rr,), get(PAR[:Sfbar], rr, 0.0))
        _safe_start_value!(model, :GDP,  (rr,), B[:GDP])
        _safe_start_value!(model, :RGDP, (rr,), B[:GDP])
        _safe_start_value!(model, :CPI,  (rr,), 1.0)
        _safe_start_value!(model, :PGDP, (rr,), 1.0)
    end

    for gg in gz
        _safe_start_value!(model, :PS, (gg,), 1.0)
    end
    _safe_start_value!(model, :PABS, (), 1.0)

    # NOTE: warm-starting the dynamic periods from the previous period's solved
    # values was tried and removed.  It makes things worse, not better: with
    # zero growth it took a period from 0.8 s to 30 s and ended in
    # SLOW_PROGRESS, because the model has near-singular directions (the
    # F-4/F-6/F-11 wage closure, see CLAUDE.md) and starting exactly on that
    # flat manifold leaves PATH with a badly conditioned basis.  Every period
    # therefore starts from the calibrated benchmark, rescaled through the PAR
    # supply tables above.

    enforce_nlp_safe_bounds_and_starts!(model)
    return model
end
