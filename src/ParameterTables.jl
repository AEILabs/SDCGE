# Usage:
#   PAR = precompute_parameters(data)
# This file precomputes all scalar/table parameters used by JuMP equations.
# JuMP equation files index `PAR` directly; they do not call helper functions inside constraints.

function _fill!(PAR, key::Symbol, combos, val)
    d = Dict{Any,Float64}()
    for c in combos
        d[c] = val
    end
    PAR[key] = d
end



function sanitize_parameters!(PAR; eps::Float64=1.0e-8)
    # Avoid numerical singularities in nonlinear CES/CET expressions.
    # Elasticities exactly equal to 1 make expressions like 1/(1-sigma) blow up.
    for (key, val) in PAR
        if val isa Dict
            for (idx, x) in collect(val)
                if x isa Real
                    xf = float(x)
                    if !isfinite(xf)
                        # Keep only deliberate regime switches as Inf; other Inf/NaN values are unsafe.
                        if !(key in [:omega_migr, :omega_F])
                            val[idx] = eps
                        end
                    elseif startswith(String(key), "sigma") && abs(xf - 1.0) < eps
                        val[idx] = 1.0 - 1.0e-4
                    elseif startswith(String(key), "lambda") && abs(xf) < eps
                        val[idx] = eps
                    elseif startswith(String(key), "alpha") && xf < 0.0
                        val[idx] = eps
                    end
                end
            end
        elseif val isa Real
            xf = float(val)
            if !isfinite(xf)
                if !(key in [:LndMAX, :eta_T, :omega_T, :omega_K])
                    PAR[key] = eps
                end
            elseif startswith(String(key), "sigma") && abs(xf - 1.0) < eps
                PAR[key] = 1.0 - 1.0e-4
            elseif startswith(String(key), "lambda") && abs(xf) < eps
                PAR[key] = eps
            end
        end
    end
    return PAR
end

function parameter_diagnostics(PAR; max_items::Int=25)
    rows = NamedTuple[]
    n_bad = 0
    n_sigma_one = 0
    for (key, val) in PAR
        if val isa Dict
            for (idx, x) in val
                if x isa Real
                    bad = !isfinite(float(x)) && !(key in [:omega_migr, :omega_F])
                    sig = startswith(String(key), "sigma") && abs(float(x) - 1.0) < 1.0e-8
                    n_bad += bad ? 1 : 0
                    n_sigma_one += sig ? 1 : 0
                    if (bad || sig) && length(rows) < max_items
                        push!(rows, (parameter=key, index=idx, value=x, issue=bad ? "nonfinite" : "sigma_equal_one"))
                    end
                end
            end
        elseif val isa Real
            bad = !isfinite(float(val)) && !(key in [:LndMAX, :eta_T, :omega_T, :omega_K])
            sig = startswith(String(key), "sigma") && abs(float(val) - 1.0) < 1.0e-8
            n_bad += bad ? 1 : 0
            n_sigma_one += sig ? 1 : 0
            if (bad || sig) && length(rows) < max_items
                push!(rows, (parameter=key, index=(), value=val, issue=bad ? "nonfinite" : "sigma_equal_one"))
            end
        end
    end
    return (bad_parameters=n_bad, sigma_equal_one=n_sigma_one, examples=rows)
end

function precompute_parameters(data::LinkageData)
    default_sets!(data)
    S = data.sets
    i=S[:i]; j=S[:j]; k=S[:k]; r=S[:r]; rp=S[:rp]; v=S[:v]; l=S[:l]; ul=S[:ul]; sl=S[:sl]; h=S[:h]; f=S[:f]; ins=S[:in]; t=S[:t]
    cr=S[:cr]; lv=S[:lv]; ip=S[:ip]; e=S[:e]; ft=S[:ft]; fd=S[:fd]; nnft=S[:nnft]; nnfd=S[:nnfd]
    PAR = Dict{Symbol,Any}()
    # scalar defaults
    for key in [:TY0,:FY0,:KY0,:GDEF0,:SAVE0,:INVEST0,:PNUM0,:WALRAS0]
        PAR[key] = 1.0
    end
    # one-dimensional tables
    for key in [:AT,:pi,:tau_p,:tau_y,:tau_f,:growth,:aeei,:popgrow,:savrate,:govshare,:invshare,:gdpdef]
        _fill!(PAR,key,i, key in [:pi,:tau_p,:tau_y,:tau_f] ? 0.0 : 1.0)
    end
    _fill!(PAR,:mu,h,0.25)
    # ELES subsistence quantities: 0 unless calibration supplies better information.
    _fill!(PAR,:theta,[(kk,hh) for kk in k for hh in h],0.0)
    _fill!(PAR,:lambda_w,[(rr,rrp,ii) for rr in r for rrp in rp for ii in i],1.0)
    _fill!(PAR,:wtr,[(rr,ii,inn) for rr in r for ii in i for inn in ins],0.0)
    _fill!(PAR,:trqrent,[(rr,rrp,ii) for rr in r for rrp in rp for ii in i],0.0)
    # common two/three-dimensional production parameters
    # Efficiency / technical-change indices (lambda_*) are 1 in the base year;
    # anything else rescales the CES duals and breaks benchmark replication.
    for key in [:alpha_nd,:alpha_va,:sigma_p,:alpha_l,:alpha_hktef,:sigma_v,:alpha_fert,:alpha_hkte,:sigma_f,:alpha_e,:alpha_hkt,:sigma_e,:alpha_h,:alpha_kt,:sigma_h,:alpha_k,:alpha_t,:alpha_ff,:lambda_k,:lambda_t,:lambda_f,:sigma_k,:alpha_ktel,:alpha_tfd,:alpha_feed,:sigma_feed,:alpha_hkte_liv]
        _fill!(PAR,key,[(ii,vv) for ii in i for vv in v], startswith(String(key),"lambda") ? 1.0 : 0.5)
    end
    for key in [:tau_Ap,:a_nd]
        _fill!(PAR,key,[(jj,ii) for jj in j for ii in i], key == :tau_Ap ? 0.0 : 0.1)
    end

    # Labor demand disaggregation parameters for production equations P-72--P-75.
    _fill!(PAR, :alpha_ul, [(ll,ii) for ll in ul for ii in i], length(ul) == 0 ? 0.0 : 1.0 / length(ul))
    _fill!(PAR, :alpha_sl, [(ll,ii) for ll in sl for ii in i], length(sl) == 0 ? 0.0 : 1.0 / length(sl))
    _fill!(PAR, :lambda_l, [(ll,ii) for ll in l for ii in i], 1.0)
    _fill!(PAR, :sigma_ul, i, 0.5)
    _fill!(PAR, :sigma_sl, i, 0.5)
    for key in [:alpha_ft,:lambda_ft,:sigma_ft,:alpha_ep,:lambda_ep,:sigma_ep,:alpha_fd,:lambda_fd,:sigma_fd]
        _fill!(PAR,key,[(jj,ii) for jj in i for ii in i], startswith(String(key),"lambda") ? 1.0 : 0.5)
    end
    for key in [:alpha_xd,:alpha_xm,:sigma_m,:alpha_w,:sigma_w,:alpha_eout,:alpha_dout,:sigma_x,:tau_m,:tau_e,:tmarg,:gamma,:rho,:kappa,:phi,:omega]
        _fill!(PAR,key,[(rr,rrp,ii) for rr in r for rrp in rp for ii in i], startswith(String(key),"tau") ? 0.0 : 0.5)
    end
    _fill!(PAR,:depr,[(ii,vv) for ii in i for vv in v],0.05)

    # Exogenous factor supplies used directly in JuMP equations.
    _fill!(PAR, :KSupply, [(ii,vv) for ii in i for vv in v], 100.0)
    _fill!(PAR, :LSupply, l, 1000.0)
    _fill!(PAR, :TSupply, i, 100.0)
    _fill!(PAR, :FSupply, i, 100.0)
    _fill!(PAR,:bop,r,0.0)

    # Income distribution parameters used in Y-5--Y-8.
    _fill!(PAR, :phi_T, h, 1.0 / length(h))
    _fill!(PAR, :phi_F, h, 1.0 / length(h))
    _fill!(PAR, :phi_K, h, 1.0 / length(h))
    _fill!(PAR, :phi_L, [(hh,ll) for hh in h for ll in l], 1.0 / length(h))
    _fill!(PAR, :TRG, h, 0.0)
    _fill!(PAR, :kappa_h, h, 0.0)
    PAR[:chi_kappa] = 1.0
    _fill!(PAR, :delta_f, [(ii,vv) for ii in i for vv in v], 0.05)
    _fill!(PAR, :WTR, [(rr,rrp,inn,"HH",hh) for rr in r for rrp in rp for inn in ins for hh in h], 0.0)

    # Factor market and dynamic parameters used by the complete paper-numbered equations.
    _fill!(PAR, :g_l, [(ll,gg) for ll in l for gg in S[:gz]], 0.0)
    _fill!(PAR, :LS0, [(ll,gg) for ll in l for gg in S[:gz]], 1000.0)
    # Labour-market closure (see the header of Factors.jl):
    #   :fixed_wage      — W = 1 and UE absorbs supply minus demand   (default)
    #   :full_employment — the wage TW clears sum_i LV = LS·(1-UE0)
    # :full_employment is implemented and square but does NOT converge yet
    # (PATH: SLOW_PROGRESS after ~440 major iterations of backtracking steps,
    # benchmark drift ~0.005-0.009 %); :fixed_wage stays the default because it
    # replicates the benchmark and solves in ~1 s.  See the Factors.jl header
    # and the hand-over notes for the remaining diagnosis.
    PAR[:labour_closure] = :fixed_wage
    # Numeraire under :full_employment (ignored under :fixed_wage, where W = 1):
    #   :pabs — PABS = 1 (default; benchmark drift 0.004 %)
    #   :cpi  — CPI[first(r)] = 1, PABS endogenous   (benchmark drift 0.60 %)
    # Both give the same start residual; :pabs is the more stable of the two.
    PAR[:numeraire] = :pabs
    _fill!(PAR, :UE0, l, 0.0)          # benchmark unemployment rate per skill
    _fill!(PAR, :chi_migr, l, 0.0)
    _fill!(PAR, :omega_migr, l, 0.5)
    _fill!(PAR, :phi_wage, [(ll,ii) for ll in l for ii in i], 1.0)
    _fill!(PAR, :tau_l, [(ll,ii) for ll in l for ii in i], 0.0)
    _fill!(PAR, :chi_wmin, [(ll,gg) for ll in l for gg in S[:gz]], 1.0)
    _fill!(PAR, :omega_ps, S[:gz], 0.0)
    _fill!(PAR, :omega_p, S[:gz], 0.0)
    _fill!(PAR, :omega_ue, S[:gz], 0.0)
    _fill!(PAR, :chi_T, [:land], 100.0)
    PAR[:eta_T] = 0.5; PAR[:omega_T] = 0.5; PAR[:LndMAX] = Inf; PAR[:gamma_ts] = 1.0; PAR[:PTLnd0] = 1.0
    _fill!(PAR, :gamma_T, i, 1.0/length(i))
    _fill!(PAR, :tau_t, i, 0.0)
    _fill!(PAR, :chi_F, i, 100.0)
    _fill!(PAR, :omega_F, i, 0.5)
    _fill!(PAR, :PF0, i, 1.0)
    _fill!(PAR, :gamma_K, i, 1.0/length(i))
    PAR[:omega_K] = 0.5; PAR[:eta_k] = 0.5; PAR[:delta] = 0.05; PAR[:nstep] = 1.0
    _fill!(PAR, :tau_k, [(ii,vv) for ii in i for vv in v], 0.0)
    _fill!(PAR, :K0, i, 100.0)
    _fill!(PAR, :g_pop, t, 0.0); _fill!(PAR, :g_L, t, 0.0); _fill!(PAR, :g_T, t, 0.0)
    _fill!(PAR, :g_F, [(ii,tt) for ii in i for tt in t], 0.0)
    _fill!(PAR, :g_y, t, 0.0); _fill!(PAR, :gamma_l, t, 0.0)
    _fill!(PAR, :pi_dyn, [(ii,tt) for ii in i for tt in t], 0.0)
    _fill!(PAR, :gamma_s, [(ii,tt) for ii in i for tt in t], 0.0)
    _fill!(PAR, :gamma_t, [(ii,tt) for ii in i for tt in t], 0.0)
    _fill!(PAR, :gamma_f, [(ii,tt) for ii in i for tt in t], 0.0)
    _fill!(PAR, :gamma_e, [(ee,ii,tt) for ee in e for ii in i for tt in t], 0.0)
    _fill!(PAR, :RGDPMP0, t, 100.0); _fill!(PAR, :Pop0, t, 100.0); _fill!(PAR, :ChiL0, t, 1.0); _fill!(PAR, :ChiT0, t, 1.0)
    _fill!(PAR, :ChiF0, [(ii,tt) for ii in i for tt in t], 1.0)
    _fill!(PAR, :alpha_p_share, [(ii,tt) for ii in i for tt in t], 0.4)
    _fill!(PAR, :eta_p, i, 0.0)

    _fill!(PAR, :lambda_l0, [(ll,ii,tt) for ll in l for ii in i for tt in t], 1.0)
    _fill!(PAR, :lambda_k0, [(ii,vv,tt) for ii in i for vv in v for tt in t], 1.0)
    _fill!(PAR, :lambda_t0, [(ii,tt) for ii in i for tt in t], 1.0)
    _fill!(PAR, :lambda_f0, [(ii,tt) for ii in i for tt in t], 1.0)
    _fill!(PAR, :lambda_ep0, [(ee,ii,tt) for ee in e for ii in i for tt in t], 1.0)
    _fill!(PAR, :aid_mu, k, 1.0/length(k)); _fill!(PAR, :aid_alpha, k, 0.2); _fill!(PAR, :aid_beta, k, 0.8)
    _fill!(PAR, :delta_c, [(kk,k2) for kk in k for k2 in k], 0.0)


    # Paper-exact final-demand and trade parameter defaults.
    _fill!(PAR, :PopH, h, 1.0)
    _fill!(PAR, :mu_c, [(kk,hh) for kk in k for hh in h], 1.0/length(k))
    _fill!(PAR, :PC0, [(kk,hh) for kk in k for hh in h], 1.0)
    _fill!(PAR, :GammaC, [(ii,kk,hh) for ii in i for kk in k for hh in h], 1.0/length(i))
    _fill!(PAR, :sigma_c, [(kk,hh) for kk in k for hh in h], 0.5)
    _fill!(PAR, :tau_Ac, [(ii,hh) for ii in i for hh in h], 0.0)
    _fill!(PAR, :a_f, [(ii,ff) for ii in i for ff in f], 1.0/length(i))
    _fill!(PAR, :tau_Af, [(ii,ff) for ii in i for ff in f], 0.0)
    _fill!(PAR, :alpha_dc, [(ii,hh) for ii in i for hh in h], 0.5)
    _fill!(PAR, :alpha_mc, [(ii,hh) for ii in i for hh in h], 0.5)
    _fill!(PAR, :sigma_mc, [(ii,hh) for ii in i for hh in h], 0.5)
    _fill!(PAR, :tau_Dc, [(ii,hh) for ii in i for hh in h], 0.0)
    _fill!(PAR, :tau_Mc, [(ii,hh) for ii in i for hh in h], 0.0)
    _fill!(PAR, :alpha_df, [(ii,ff) for ii in i for ff in f], 0.5)
    _fill!(PAR, :alpha_mf, [(ii,ff) for ii in i for ff in f], 0.5)
    _fill!(PAR, :sigma_mf, [(ii,ff) for ii in i for ff in f], 0.5)
    _fill!(PAR, :tau_Df, [(ii,ff) for ii in i for ff in f], 0.0)
    _fill!(PAR, :tau_Mf, [(ii,ff) for ii in i for ff in f], 0.0)
    _fill!(PAR, :beta_d, i, 0.5)
    _fill!(PAR, :beta_m, i, 0.5)
    _fill!(PAR, :sigma_top_m, i, 0.5)
    _fill!(PAR, :beta_1, [(rr,ii) for rr in r for ii in i], 1.0/length(r))
    _fill!(PAR, :beta_2, [(rr,ii) for rr in r for ii in i], 1.0/length(r))
    _fill!(PAR, :beta_w, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 1.0/length(r))
    _fill!(PAR, :sigma_w1, i, 0.5)
    _fill!(PAR, :sigma_w2, [(rr,ii) for rr in r for ii in i], 0.5)
    _fill!(PAR, :sigma_w3, [(rr,ii) for rr in r for ii in i], 0.5)
    _fill!(PAR, :beta_xd, i, 0.5)
    _fill!(PAR, :beta_es, i, 0.5)
    _fill!(PAR, :sigma_z, i, 0.5)
    # beta_z must sum to 1 over the FULL (rr,rrp) grid, not just over rrp: T_19
    # (Trade.jl) is a single flat CET dual price nest summing beta_z over both
    # rr and rrp together (sum(sum(...for rrp in rp) for rr in r)), and T_18's
    # bilateral CET supply split allocates one scalar ES[ii] across every
    # (rr,rrp) pair. 1/length(rp) (summing to length(r), not 1) left T_19's
    # price identity off by length(r)^(1/(1+sigma_z2)) at the benchmark.
    _fill!(PAR, :beta_z, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 1.0/(length(r)*length(rp)))
    _fill!(PAR, :sigma_z2, i, 0.5)
    _fill!(PAR, :tau_pr, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 0.0)
    _fill!(PAR, :tau_in, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 0.0)
    _fill!(PAR, :tau_out, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 0.0)
    _fill!(PAR, :tau_trq_share, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 0.0)
    _fill!(PAR, :zeta_t, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 0.0)
    _fill!(PAR, :lambda_w, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 1.0)
    _fill!(PAR, :alpha_TT, r, 1.0/length(r))
    PAR[:sigma_TT] = 0.5
    _fill!(PAR, :a_Mg, [(rr,ii) for rr in r for ii in i], 1.0/length(i))


    _fill!(PAR, :WTRgov_in, [(rr,inn) for rr in r for inn in ins], 0.0)
    _fill!(PAR, :WTRgov_out, [(rr,inn) for rr in r for inn in ins], 0.0)
    _fill!(PAR, :WTRinv_in, [(rr,inn) for rr in r for inn in ins], 0.0)
    _fill!(PAR, :WTRinv_out, [(rr,inn) for rr in r for inn in ins], 0.0)
    PAR[:chi_gov] = 0.2
    PAR[:chi_inv] = 0.2
    PAR[:FDInv0]  = 1.0
    _fill!(PAR, :Sfbar, r, 0.0)
    # Trade closure (see the header of Calibration.jl and README "Closures"):
    #   :bop      — small open economy.  World prices are exogenous in foreign
    #               currency (PWE0/PWM0) and converted with the real exchange
    #               rate ER; imports and exports are independent; the current
    #               account is pinned by Sfbar (C-BOP).  Default.
    #   :balanced — the legacy convention.  E-2 forces WTFd = lambda_w·WTFs, so
    #               world-price imports equal world-price exports good by good
    #               and Sfbar must be 0.
    PAR[:trade_closure] = :bop
    # Which variable clears the balance of payments (ignored under :balanced):
    #   :flex_er  — Sfbar exogenous, the real exchange rate ER adjusts (default,
    #               the standard LINKAGE closure).
    #   :fixed_er — ER = ER0 fixed, foreign saving Sf of the home region adjusts.
    PAR[:bop_closure] = :flex_er
    # Investment under :bop: :fixed (exogenous real investment, C-INV) or :savings
    # (C-9 kept, investment = savings; experimental, see README "Closures").
    PAR[:inv_closure] = :fixed
    PAR[:ER0] = 1.0
    # World prices in foreign currency, one per (origin, destination, good).
    _fill!(PAR, :PWE0, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 1.0)
    _fill!(PAR, :PWM0, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 1.0)
    _fill!(PAR, :WTF0, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 1.0)
    _fill!(PAR, :WPE0, [(rr,rrp,ii) for rr in r for rrp in rp for ii in i], 1.0)
    _fill!(PAR, :TR_region, r, 1.0)
    _fill!(PAR, :K_region, r, 1.0)


    _fill!(PAR, :aids_alpha, [(rrp,rr,ii) for rrp in rp for rr in r for ii in i], 0.0)
    _fill!(PAR, :aids_beta, [(rrp,rr,ii) for rrp in rp for rr in r for ii in i], 0.0)
    _fill!(PAR, :aids_alpha0, [(rr,ii) for rr in r for ii in i], 0.0)
    _fill!(PAR, :aids_gamma, [(rrp,rrpp,rr,ii) for rrp in rp for rrpp in rp for rr in r for ii in i], 0.0)
    _fill!(PAR, :alpha_dd, [(rr,ii) for rr in r for ii in i], 0.5)
    _fill!(PAR, :alpha_dm, [(rr,ii) for rr in r for ii in i], 0.5)
    _fill!(PAR, :sigma_b, [(rr,ii) for rr in r for ii in i], 0.5)

    # Overlay SAM-calibrated parameters. Keys in data.par replace defaults.
    for (kcal, vcal) in data.par
        PAR[kcal] = vcal
    end
    sanitize_parameters!(PAR)
    return PAR
end
