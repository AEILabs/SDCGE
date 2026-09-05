# Usage: add_production_equations!(model, data, PAR)
# Variables are declared in Variables.jl. This file only adds paper-numbered JuMP equations.

function add_production_equations!(model, data::LinkageData, PAR)
    S=data.sets; default_sets!(data)
    i=S[:i]; j=S[:j]; k=S[:k]; r=S[:r]; rp=S[:rp]; v=S[:v]; l=S[:l]; h=S[:h]; f=S[:f]; ins=S[:in]; t=S[:t]
    cr=S[:cr]; lv=S[:lv]; ip=S[:ip]; e=S[:e]; ft=S[:ft]; fd=S[:fd]; nf=S[:nf]; nnft=S[:nnft]; nnfd=S[:nnfd]; gz=S[:gz]; ul=S[:ul]; sl=S[:sl]

    XPv=model[:XPv]; XP=model[:XP]; ND=model[:ND]; VA=model[:VA]; PND=model[:PND]; PVA=model[:PVA]; UVCv=model[:UVCv]; UVC=model[:UVC]; AC=model[:AC]; Nfirm=model[:Nfirm]; PX=model[:PX]; PP=model[:PP]; PROFIT=model[:PROFIT]
    ULD=model[:ULD]; SLD=model[:SLD]; HKTEF=model[:HKTEF]; PHKTEF=model[:PHKTEF]; fert=model[:fert]; Pfert=model[:Pfert]; HKTE=model[:HKTE]; PHKTE=model[:PHKTE]; XEp=model[:XEp]; PEp=model[:PEp]; HKT=model[:HKT]; PHKT=model[:PHKT]; KT=model[:KT]; PKT=model[:PKT]; Kvd=model[:Kvd]; Td=model[:Td]; Fd=model[:Fd]
    XAp=model[:XAp]; PA=model[:PA]; PAp=model[:PAp]; Pfeed=model[:Pfeed]; feed=model[:feed]; KTEL=model[:KTEL]; PKTEL=model[:PKTEL]; TFD=model[:TFD]; PTFD=model[:PTFD]; W=model[:W]; UW=model[:UW]; SW=model[:SW]; LV=model[:LV]; LF_d=model[:LF_d]; KF_d=model[:KF_d]; PT=model[:PT]; PF=model[:PF]; R=model[:R]

    # Index sets of the ND (fixed-coefficient intermediate) bundle.  Energy is
    # bought through the CES energy nest (P-27/P-51/P-67), fertiliser through
    # P-25 (crops) and feed through P-49 (livestock), so those goods are NOT in
    # the ND bundle.  The quantity equations already used these sets; the price
    # indices P-30/P-54/P-70 summed over the wider sets, which made PND differ
    # from the true unit cost of the bundle and left a residual in P-3/P-1.
    cr_nd = [x for x in nnft if !(x in e)]
    lv_nd = [x for x in nnfd if !(x in e)]
    ip_nd = [x for x in i    if !(x in e)]

    # (P-1) Aggregate intermediate demand by vintage top nest
    @constraint(model, P_1[ii in i], (ND[ii]) - (PAR[:AT][ii] * sum(PAR[:alpha_nd][(ii,vv)] * XPv[ii,vv] * (UVCv[ii,vv] / PND[ii])^(PAR[:sigma_p][(ii,vv)]) for vv in v)) ⟂ ND[ii])
    # (P-2) Value added demand by vintage top nest
    @constraint(model, P_2[ii in i, vv in v], (VA[ii,vv]) - (PAR[:AT][ii] * PAR[:alpha_va][(ii,vv)] * XPv[ii,vv] * (UVCv[ii,vv] / PVA[ii,vv])^(PAR[:sigma_p][(ii,vv)])) ⟂ VA[ii,vv])
    # (P-3) Unit variable cost by vintage
    @constraint(model, P_3[ii in i, vv in v], (UVCv[ii,vv]) - ((PAR[:alpha_nd][(ii,vv)]*PND[ii]^(1-PAR[:sigma_p][(ii,vv)]) + PAR[:alpha_va][(ii,vv)]*PVA[ii,vv]^(1-PAR[:sigma_p][(ii,vv)]))^(1/(1-PAR[:sigma_p][(ii,vv)])) / PAR[:AT][ii]) ⟂ UVCv[ii,vv])
    # (P-4) Aggregate unit variable cost
    @constraint(model, P_4[ii in i], (UVC[ii] * sum(XPv[ii,vv] for vv in v)) - (sum(UVCv[ii,vv] * XPv[ii,vv] for vv in v)) ⟂ UVC[ii])
    # (P-5) Average cost including fixed costs
    @constraint(model, P_5[ii in i], (AC[ii]) - (UVC[ii] + Nfirm[ii]*(sum(W[ll,ii]*LF_d[ll,ii] for ll in l) + R[ii,"Old"]*KF_d[ii])/(XP[ii] + 1.0e-6)) ⟂ AC[ii])
    # (P-6) Net output price with markup
    @constraint(model, P_6[ii in i], (PX[ii]) - (UVC[ii]*(1 + PAR[:pi][ii])) ⟂ PX[ii])
    # (P-7) Gross producer price including output tax
    @constraint(model, P_7[ii in i], (PP[ii]) - (PX[ii]*(1 + PAR[:tau_p][ii])) ⟂ PP[ii])
    # (P-8) Gross profit
    @constraint(model, P_8[ii in i], (PROFIT[ii]) - (XP[ii]*(PX[ii] - AC[ii])) ⟂ PROFIT[ii])

    # Crop production P-9 to P-30
    @constraint(model, P_9[ii in cr], (ULD[ii]) - (sum(PAR[:alpha_l][(ii,vv)]*(PVA[ii,vv]/UW[ii])^(PAR[:sigma_v][(ii,vv)])*VA[ii,vv] for vv in v)) ⟂ ULD[ii])
    @constraint(model, P_10[ii in cr, vv in v], (HKTEF[ii,vv]) - (PAR[:alpha_hktef][(ii,vv)]*(PVA[ii,vv]/PHKTEF[ii,vv])^(PAR[:sigma_v][(ii,vv)])*VA[ii,vv]) ⟂ HKTEF[ii,vv])
    @constraint(model, P_11[ii in cr, vv in v], (PVA[ii,vv]) - ((PAR[:alpha_l][(ii,vv)]*UW[ii]^(1-PAR[:sigma_v][(ii,vv)]) + PAR[:alpha_hktef][(ii,vv)]*PHKTEF[ii,vv]^(1-PAR[:sigma_v][(ii,vv)]))^(1/(1-PAR[:sigma_v][(ii,vv)]))) ⟂ PVA[ii,vv])
    @constraint(model, P_12[ii in cr], (fert[ii]) - (sum(PAR[:alpha_fert][(ii,vv)]*(PHKTEF[ii,vv]/Pfert[ii])^(PAR[:sigma_f][(ii,vv)])*HKTEF[ii,vv] for vv in v)) ⟂ fert[ii])
    @constraint(model, P_13[ii in cr, vv in v], (HKTE[ii,vv]) - (PAR[:alpha_hkte][(ii,vv)]*(PHKTEF[ii,vv]/PHKTE[ii,vv])^(PAR[:sigma_f][(ii,vv)])*HKTEF[ii,vv]) ⟂ HKTE[ii,vv])
    @constraint(model, P_14[ii in cr, vv in v], (PHKTEF[ii,vv]) - ((PAR[:alpha_fert][(ii,vv)]*Pfert[ii]^(1-PAR[:sigma_f][(ii,vv)]) + PAR[:alpha_hkte][(ii,vv)]*PHKTE[ii,vv]^(1-PAR[:sigma_f][(ii,vv)]))^(1/(1-PAR[:sigma_f][(ii,vv)]))) ⟂ PHKTEF[ii,vv])
    @constraint(model, P_15[ii in cr, vv in v], (XEp[ii,vv]) - (PAR[:alpha_e][(ii,vv)]*(PHKTE[ii,vv]/PEp[ii,vv])^(PAR[:sigma_e][(ii,vv)])*HKTE[ii,vv]) ⟂ XEp[ii,vv])
    @constraint(model, P_16[ii in cr, vv in v], (HKT[ii,vv]) - (PAR[:alpha_hkt][(ii,vv)]*(PHKTE[ii,vv]/PHKT[ii,vv])^(PAR[:sigma_e][(ii,vv)])*HKTE[ii,vv]) ⟂ HKT[ii,vv])
    @constraint(model, P_17[ii in cr, vv in v], (PHKTE[ii,vv]) - ((PAR[:alpha_e][(ii,vv)]*PEp[ii,vv]^(1-PAR[:sigma_e][(ii,vv)]) + PAR[:alpha_hkt][(ii,vv)]*PHKT[ii,vv]^(1-PAR[:sigma_e][(ii,vv)]))^(1/(1-PAR[:sigma_e][(ii,vv)]))) ⟂ PHKTE[ii,vv])
    @constraint(model, P_18[ii in cr], (SLD[ii]) - (sum(PAR[:alpha_h][(ii,vv)]*(PHKT[ii,vv]/SW[ii])^(PAR[:sigma_h][(ii,vv)])*HKT[ii,vv] for vv in v)) ⟂ SLD[ii])
    @constraint(model, P_19[ii in cr, vv in v], (KT[ii,vv]) - (PAR[:alpha_kt][(ii,vv)]*(PHKT[ii,vv]/PKT[ii,vv])^(PAR[:sigma_h][(ii,vv)])*HKT[ii,vv]) ⟂ KT[ii,vv])
    @constraint(model, P_20[ii in cr, vv in v], (PHKT[ii,vv]) - ((PAR[:alpha_h][(ii,vv)]*SW[ii]^(1-PAR[:sigma_h][(ii,vv)]) + PAR[:alpha_kt][(ii,vv)]*PKT[ii,vv]^(1-PAR[:sigma_h][(ii,vv)]))^(1/(1-PAR[:sigma_h][(ii,vv)]))) ⟂ PHKT[ii,vv])
    @constraint(model, P_21[ii in cr, vv in v], (Kvd[ii,vv]) - (PAR[:lambda_k][(ii,vv)]^(PAR[:sigma_k][(ii,vv)]-1)*PAR[:alpha_k][(ii,vv)]*(PKT[ii,vv]/R[ii,vv])^(PAR[:sigma_k][(ii,vv)])*KT[ii,vv]) ⟂ Kvd[ii,vv])
    @constraint(model, P_22[ii in cr], (Td[ii]) - (sum(PAR[:lambda_t][(ii,vv)]^(PAR[:sigma_k][(ii,vv)]-1)*PAR[:alpha_t][(ii,vv)]*(PKT[ii,vv]/PT[ii])^(PAR[:sigma_k][(ii,vv)])*KT[ii,vv] for vv in v)) ⟂ Td[ii])
    @constraint(model, P_23[ii in cr], (Fd[ii]) - (sum(PAR[:lambda_f][(ii,vv)]^(PAR[:sigma_k][(ii,vv)]-1)*PAR[:alpha_ff][(ii,vv)]*(PKT[ii,vv]/PF[ii])^(PAR[:sigma_k][(ii,vv)])*KT[ii,vv] for vv in v)) ⟂ Fd[ii])
    @constraint(model, P_24[ii in cr, vv in v], (PKT[ii,vv]) - ((PAR[:alpha_k][(ii,vv)]*(R[ii,vv]/PAR[:lambda_k][(ii,vv)])^(1-PAR[:sigma_k][(ii,vv)]) + PAR[:alpha_t][(ii,vv)]*(PT[ii]/PAR[:lambda_t][(ii,vv)])^(1-PAR[:sigma_k][(ii,vv)]) + PAR[:alpha_ff][(ii,vv)]*(PF[ii]/PAR[:lambda_f][(ii,vv)])^(1-PAR[:sigma_k][(ii,vv)]))^(1/(1-PAR[:sigma_k][(ii,vv)]))) ⟂ PKT[ii,vv])
    @constraint(model, P_25[jj in ft, ii in cr], (XAp[jj,ii]) - (PAR[:lambda_ft][(jj,ii)]^(PAR[:sigma_ft][(jj,ii)]-1)*PAR[:alpha_ft][(jj,ii)]*(Pfert[ii]/((1+PAR[:tau_Ap][(jj,ii)])*PA[jj]))^(PAR[:sigma_ft][(jj,ii)])*fert[ii]) ⟂ XAp[jj,ii])
    # (an empty ft set — no such commodity in the economy — leaves the bundle price at 1)
    @constraint(model, P_26[ii in cr], (Pfert[ii]) - (isempty(ft) ? 1.0 : (sum(PAR[:alpha_ft][(jj,ii)]*((1+PAR[:tau_Ap][(jj,ii)])*PA[jj]/PAR[:lambda_ft][(jj,ii)])^(1-PAR[:sigma_ft][(jj,ii)]) for jj in ft))^(1/(1-PAR[:sigma_ft][(first(ft),ii)]))) ⟂ Pfert[ii])
    @constraint(model, P_27[jj in e, ii in cr], (XAp[jj,ii]) - (sum(PAR[:lambda_ep][(jj,ii)]^(PAR[:sigma_ep][(jj,ii)]-1)*PAR[:alpha_ep][(jj,ii)]*(PEp[ii,vv]/((1+PAR[:tau_Ap][(jj,ii)])*PA[jj]))^(PAR[:sigma_ep][(jj,ii)])*XEp[ii,vv] for vv in v)) ⟂ XAp[jj,ii])
    @constraint(model, P_28[ii in cr, vv in v], (PEp[ii,vv]) - ((sum(PAR[:alpha_ep][(jj,ii)]*((1+PAR[:tau_Ap][(jj,ii)])*PA[jj]/PAR[:lambda_ep][(jj,ii)])^(1-PAR[:sigma_ep][(jj,ii)]) for jj in e))^(1/(1-PAR[:sigma_ep][(first(e),ii)]))) ⟂ PEp[ii,vv])
    # P_29 restricted to nnft minus energy: energy inputs are handled by P_27 (CES).
    # Including energy here would give XAp[e,cr] two equations.
    @constraint(model, P_29[jj in cr_nd, ii in cr], (XAp[jj,ii]) - (PAR[:a_nd][(jj,ii)]*((PND[ii]/((1+PAR[:tau_Ap][(jj,ii)])*PA[jj]))^0.0)*ND[ii]) ⟂ XAp[jj,ii])
    @constraint(model, P_30[ii in cr], (PND[ii]) - (sum(PAR[:a_nd][(jj,ii)]*(1+PAR[:tau_Ap][(jj,ii)])*PA[jj] for jj in cr_nd)) ⟂ PND[ii])

    # Livestock production P-31 to P-54
    @constraint(model, P_31[ii in lv, vv in v], (KTEL[ii,vv]) - (PAR[:alpha_ktel][(ii,vv)]*VA[ii,vv]) ⟂ KTEL[ii,vv])
    @constraint(model, P_32[ii in lv, vv in v], (TFD[ii,vv]) - (PAR[:alpha_tfd][(ii,vv)]*VA[ii,vv]) ⟂ TFD[ii,vv])
    @constraint(model, P_33[ii in lv, vv in v], (PVA[ii,vv]) - (PAR[:alpha_ktel][(ii,vv)]*PKTEL[ii,vv] + PAR[:alpha_tfd][(ii,vv)]*PTFD[ii,vv]) ⟂ PVA[ii,vv])
    @constraint(model, P_34[ii in lv], (feed[ii]) - (sum(PAR[:alpha_feed][(ii,vv)]*(PTFD[ii,vv]/Pfeed[ii])^(PAR[:sigma_feed][(ii,vv)])*TFD[ii,vv] for vv in v)) ⟂ feed[ii])
    @constraint(model, P_35[ii in lv], (Td[ii]) - (sum(PAR[:lambda_t][(ii,vv)]^(PAR[:sigma_feed][(ii,vv)]-1)*PAR[:alpha_t][(ii,vv)]*(PTFD[ii,vv]/PT[ii])^(PAR[:sigma_feed][(ii,vv)])*TFD[ii,vv] for vv in v)) ⟂ Td[ii])
    @constraint(model, P_36[ii in lv, vv in v], (PTFD[ii,vv]) - ((PAR[:alpha_feed][(ii,vv)]*Pfeed[ii]^(1-PAR[:sigma_feed][(ii,vv)]) + PAR[:alpha_t][(ii,vv)]*(PT[ii]/PAR[:lambda_t][(ii,vv)])^(1-PAR[:sigma_feed][(ii,vv)]))^(1/(1-PAR[:sigma_feed][(ii,vv)]))) ⟂ PTFD[ii,vv])
    @constraint(model, P_37[ii in lv], (ULD[ii]) - (sum(PAR[:alpha_l][(ii,vv)]*(PKTEL[ii,vv]/UW[ii])^(PAR[:sigma_v][(ii,vv)])*KTEL[ii,vv] for vv in v)) ⟂ ULD[ii])
    @constraint(model, P_38[ii in lv, vv in v], (HKTE[ii,vv]) - (PAR[:alpha_hkte_liv][(ii,vv)]*(PKTEL[ii,vv]/PHKTE[ii,vv])^(PAR[:sigma_v][(ii,vv)])*KTEL[ii,vv]) ⟂ HKTE[ii,vv])
    @constraint(model, P_39[ii in lv, vv in v], (PKTEL[ii,vv]) - ((PAR[:alpha_l][(ii,vv)]*UW[ii]^(1-PAR[:sigma_v][(ii,vv)]) + PAR[:alpha_hkte_liv][(ii,vv)]*PHKTE[ii,vv]^(1-PAR[:sigma_v][(ii,vv)]))^(1/(1-PAR[:sigma_v][(ii,vv)]))) ⟂ PKTEL[ii,vv])
    @constraint(model, P_40[ii in lv, vv in v], (XEp[ii,vv]) - (PAR[:alpha_e][(ii,vv)]*(PHKTE[ii,vv]/PEp[ii,vv])^(PAR[:sigma_e][(ii,vv)])*HKTE[ii,vv]) ⟂ XEp[ii,vv])
    @constraint(model, P_41[ii in lv, vv in v], (HKT[ii,vv]) - (PAR[:alpha_hkt][(ii,vv)]*(PHKTE[ii,vv]/PHKT[ii,vv])^(PAR[:sigma_e][(ii,vv)])*HKTE[ii,vv]) ⟂ HKT[ii,vv])
    @constraint(model, P_42[ii in lv, vv in v], (PHKTE[ii,vv]) - ((PAR[:alpha_e][(ii,vv)]*PEp[ii,vv]^(1-PAR[:sigma_e][(ii,vv)]) + PAR[:alpha_hkt][(ii,vv)]*PHKT[ii,vv]^(1-PAR[:sigma_e][(ii,vv)]))^(1/(1-PAR[:sigma_e][(ii,vv)]))) ⟂ PHKTE[ii,vv])
    @constraint(model, P_43[ii in lv], (SLD[ii]) - (sum(PAR[:alpha_h][(ii,vv)]*(PHKT[ii,vv]/SW[ii])^(PAR[:sigma_h][(ii,vv)])*HKT[ii,vv] for vv in v)) ⟂ SLD[ii])
    @constraint(model, P_44[ii in lv, vv in v], (KT[ii,vv]) - (PAR[:alpha_kt][(ii,vv)]*(PHKT[ii,vv]/PKT[ii,vv])^(PAR[:sigma_h][(ii,vv)])*HKT[ii,vv]) ⟂ KT[ii,vv])
    @constraint(model, P_45[ii in lv, vv in v], (PHKT[ii,vv]) - ((PAR[:alpha_h][(ii,vv)]*SW[ii]^(1-PAR[:sigma_h][(ii,vv)]) + PAR[:alpha_kt][(ii,vv)]*PKT[ii,vv]^(1-PAR[:sigma_h][(ii,vv)]))^(1/(1-PAR[:sigma_h][(ii,vv)]))) ⟂ PHKT[ii,vv])
    @constraint(model, P_46[ii in lv, vv in v], (Kvd[ii,vv]) - (PAR[:lambda_k][(ii,vv)]^(PAR[:sigma_k][(ii,vv)]-1)*PAR[:alpha_k][(ii,vv)]*(PKT[ii,vv]/R[ii,vv])^(PAR[:sigma_k][(ii,vv)])*KT[ii,vv]) ⟂ Kvd[ii,vv])
    @constraint(model, P_47[ii in lv], (Fd[ii]) - (sum(PAR[:lambda_f][(ii,vv)]^(PAR[:sigma_k][(ii,vv)]-1)*PAR[:alpha_ff][(ii,vv)]*(PKT[ii,vv]/PF[ii])^(PAR[:sigma_k][(ii,vv)])*KT[ii,vv] for vv in v)) ⟂ Fd[ii])
    @constraint(model, P_48[ii in lv, vv in v], (PKT[ii,vv]) - ((PAR[:alpha_k][(ii,vv)]*(R[ii,vv]/PAR[:lambda_k][(ii,vv)])^(1-PAR[:sigma_k][(ii,vv)]) + PAR[:alpha_ff][(ii,vv)]*(PF[ii]/PAR[:lambda_f][(ii,vv)])^(1-PAR[:sigma_k][(ii,vv)]))^(1/(1-PAR[:sigma_k][(ii,vv)]))) ⟂ PKT[ii,vv])
    @constraint(model, P_49[jj in fd, ii in lv], (XAp[jj,ii]) - (PAR[:lambda_fd][(jj,ii)]^(PAR[:sigma_fd][(jj,ii)]-1)*PAR[:alpha_fd][(jj,ii)]*(Pfeed[ii]/((1+PAR[:tau_Ap][(jj,ii)])*PA[jj]))^(PAR[:sigma_fd][(jj,ii)])*feed[ii]) ⟂ XAp[jj,ii])
    # (an empty fd set — no such commodity in the economy — leaves the bundle price at 1)
    @constraint(model, P_50[ii in lv], (Pfeed[ii]) - (isempty(fd) ? 1.0 : (sum(PAR[:alpha_fd][(jj,ii)]*((1+PAR[:tau_Ap][(jj,ii)])*PA[jj]/PAR[:lambda_fd][(jj,ii)])^(1-PAR[:sigma_fd][(jj,ii)]) for jj in fd))^(1/(1-PAR[:sigma_fd][(first(fd),ii)]))) ⟂ Pfeed[ii])
    @constraint(model, P_51[jj in e, ii in lv], (XAp[jj,ii]) - (sum(PAR[:lambda_ep][(jj,ii)]^(PAR[:sigma_ep][(jj,ii)]-1)*PAR[:alpha_ep][(jj,ii)]*(PEp[ii,vv]/((1+PAR[:tau_Ap][(jj,ii)])*PA[jj]))^(PAR[:sigma_ep][(jj,ii)])*XEp[ii,vv] for vv in v)) ⟂ XAp[jj,ii])
    @constraint(model, P_52[ii in lv, vv in v], (PEp[ii,vv]) - ((sum(PAR[:alpha_ep][(jj,ii)]*((1+PAR[:tau_Ap][(jj,ii)])*PA[jj]/PAR[:lambda_ep][(jj,ii)])^(1-PAR[:sigma_ep][(jj,ii)]) for jj in e))^(1/(1-PAR[:sigma_ep][(first(e),ii)]))) ⟂ PEp[ii,vv])
    # P_53 restricted to nnfd minus energy: energy inputs handled by P_51 (CES).
    @constraint(model, P_53[jj in lv_nd, ii in lv], (XAp[jj,ii]) - (PAR[:a_nd][(jj,ii)]*ND[ii]) ⟂ XAp[jj,ii])
    @constraint(model, P_54[ii in lv], (PND[ii]) - (sum(PAR[:a_nd][(jj,ii)]*(1+PAR[:tau_Ap][(jj,ii)])*PA[jj] for jj in lv_nd)) ⟂ PND[ii])

    # Non-agricultural production P-55 to P-80
    @constraint(model, P_55[ii in ip], (ULD[ii]) - (sum(PAR[:alpha_l][(ii,vv)]*(PVA[ii,vv]/UW[ii])^(PAR[:sigma_v][(ii,vv)])*VA[ii,vv] for vv in v)) ⟂ ULD[ii])
    @constraint(model, P_56[ii in ip, vv in v], (HKTE[ii,vv]) - (PAR[:alpha_hkte][(ii,vv)]*(PVA[ii,vv]/PHKTE[ii,vv])^(PAR[:sigma_v][(ii,vv)])*VA[ii,vv]) ⟂ HKTE[ii,vv])
    @constraint(model, P_57[ii in ip, vv in v], (PVA[ii,vv]) - ((PAR[:alpha_l][(ii,vv)]*UW[ii]^(1-PAR[:sigma_v][(ii,vv)]) + PAR[:alpha_hkte][(ii,vv)]*PHKTE[ii,vv]^(1-PAR[:sigma_v][(ii,vv)]))^(1/(1-PAR[:sigma_v][(ii,vv)]))) ⟂ PVA[ii,vv])
    @constraint(model, P_58[ii in ip, vv in v], (XEp[ii,vv]) - (PAR[:alpha_e][(ii,vv)]*(PHKTE[ii,vv]/PEp[ii,vv])^(PAR[:sigma_e][(ii,vv)])*HKTE[ii,vv]) ⟂ XEp[ii,vv])
    @constraint(model, P_59[ii in ip, vv in v], (HKT[ii,vv]) - (PAR[:alpha_hkt][(ii,vv)]*(PHKTE[ii,vv]/PHKT[ii,vv])^(PAR[:sigma_e][(ii,vv)])*HKTE[ii,vv]) ⟂ HKT[ii,vv])
    @constraint(model, P_60[ii in ip, vv in v], (PHKTE[ii,vv]) - ((PAR[:alpha_e][(ii,vv)]*PEp[ii,vv]^(1-PAR[:sigma_e][(ii,vv)]) + PAR[:alpha_hkt][(ii,vv)]*PHKT[ii,vv]^(1-PAR[:sigma_e][(ii,vv)]))^(1/(1-PAR[:sigma_e][(ii,vv)]))) ⟂ PHKTE[ii,vv])
    @constraint(model, P_61[ii in ip], (SLD[ii]) - (sum(PAR[:alpha_h][(ii,vv)]*(PHKT[ii,vv]/SW[ii])^(PAR[:sigma_h][(ii,vv)])*HKT[ii,vv] for vv in v)) ⟂ SLD[ii])
    @constraint(model, P_62[ii in ip, vv in v], (KT[ii,vv]) - (PAR[:alpha_kt][(ii,vv)]*(PHKT[ii,vv]/PKT[ii,vv])^(PAR[:sigma_h][(ii,vv)])*HKT[ii,vv]) ⟂ KT[ii,vv])
    @constraint(model, P_63[ii in ip, vv in v], (PHKT[ii,vv]) - ((PAR[:alpha_h][(ii,vv)]*SW[ii]^(1-PAR[:sigma_h][(ii,vv)]) + PAR[:alpha_kt][(ii,vv)]*PKT[ii,vv]^(1-PAR[:sigma_h][(ii,vv)]))^(1/(1-PAR[:sigma_h][(ii,vv)]))) ⟂ PHKT[ii,vv])
    @constraint(model, P_64[ii in ip, vv in v], (Kvd[ii,vv]) - (PAR[:lambda_k][(ii,vv)]^(PAR[:sigma_k][(ii,vv)]-1)*PAR[:alpha_k][(ii,vv)]*(PKT[ii,vv]/R[ii,vv])^(PAR[:sigma_k][(ii,vv)])*KT[ii,vv]) ⟂ Kvd[ii,vv])
    @constraint(model, P_65[ii in ip], (Fd[ii]) - (sum(PAR[:lambda_f][(ii,vv)]^(PAR[:sigma_k][(ii,vv)]-1)*PAR[:alpha_ff][(ii,vv)]*(PKT[ii,vv]/PF[ii])^(PAR[:sigma_k][(ii,vv)])*KT[ii,vv] for vv in v)) ⟂ Fd[ii])
    @constraint(model, P_66[ii in ip, vv in v], (PKT[ii,vv]) - ((PAR[:alpha_k][(ii,vv)]*(R[ii,vv]/PAR[:lambda_k][(ii,vv)])^(1-PAR[:sigma_k][(ii,vv)]) + PAR[:alpha_ff][(ii,vv)]*(PF[ii]/PAR[:lambda_f][(ii,vv)])^(1-PAR[:sigma_k][(ii,vv)]))^(1/(1-PAR[:sigma_k][(ii,vv)]))) ⟂ PKT[ii,vv])
    @constraint(model, P_67[jj in e, ii in ip], (XAp[jj,ii]) - (sum(PAR[:lambda_ep][(jj,ii)]^(PAR[:sigma_ep][(jj,ii)]-1)*PAR[:alpha_ep][(jj,ii)]*(PEp[ii,vv]/((1+PAR[:tau_Ap][(jj,ii)])*PA[jj]))^(PAR[:sigma_ep][(jj,ii)])*XEp[ii,vv] for vv in v)) ⟂ XAp[jj,ii])
    @constraint(model, P_68[ii in ip, vv in v], (PEp[ii,vv]) - ((sum(PAR[:alpha_ep][(jj,ii)]*((1+PAR[:tau_Ap][(jj,ii)])*PA[jj]/PAR[:lambda_ep][(jj,ii)])^(1-PAR[:sigma_ep][(jj,ii)]) for jj in e))^(1/(1-PAR[:sigma_ep][(first(e),ii)]))) ⟂ PEp[ii,vv])
    # P_69 restricted to non-energy goods: energy inputs handled by P_67 (CES).
    @constraint(model, P_69[jj in ip_nd, ii in ip], (XAp[jj,ii]) - (PAR[:a_nd][(jj,ii)]*ND[ii]) ⟂ XAp[jj,ii])
    @constraint(model, P_70[ii in ip], (PND[ii]) - (sum(PAR[:a_nd][(jj,ii)]*(1+PAR[:tau_Ap][(jj,ii)])*PA[jj] for jj in ip_nd)) ⟂ PND[ii])
    # Labor-demand disaggregation P-72 to P-75.
    # These equations decompose the aggregate unskilled/skilled labor bundles into skill-specific labor demand LV.
    # λ_l captures labor-augmenting technical change; σ_ul and σ_sl govern within-bundle substitution.
    @constraint(model, P_72[ll in ul, ii in i], (LV[ll,ii]) - (PAR[:lambda_l][(ll,ii)]^(PAR[:sigma_ul][ii] - 1) *
                     PAR[:alpha_ul][(ll,ii)] *
                     (UW[ii] / (W[ll,ii] + 1.0e-9))^PAR[:sigma_ul][ii] * ULD[ii]) ⟂ LV[ll,ii])

    @constraint(model, P_73[ii in i], (UW[ii]) - ((sum(PAR[:alpha_ul][(ll,ii)] *
                       ((W[ll,ii] + 1.0e-9) / PAR[:lambda_l][(ll,ii)])^(1 - PAR[:sigma_ul][ii])
                       for ll in ul))^(1 / (1 - PAR[:sigma_ul][ii]))) ⟂ UW[ii])

    @constraint(model, P_74[ll in sl, ii in i], (LV[ll,ii]) - (PAR[:lambda_l][(ll,ii)]^(PAR[:sigma_sl][ii] - 1) *
                     PAR[:alpha_sl][(ll,ii)] *
                     (SW[ii] / (W[ll,ii] + 1.0e-9))^PAR[:sigma_sl][ii] * SLD[ii]) ⟂ LV[ll,ii])

    @constraint(model, P_75[ii in i], (SW[ii]) - ((sum(PAR[:alpha_sl][(ll,ii)] *
                       ((W[ll,ii] + 1.0e-9) / PAR[:lambda_l][(ll,ii)])^(1 - PAR[:sigma_sl][ii])
                       for ll in sl))^(1 / (1 - PAR[:sigma_sl][ii]))) ⟂ SW[ii])

    # Tax-inclusive intermediate price for diagonal elements (off-diagonal are fixed at PA[j] via PAp_wedge_offdiag in Other.jl).
    @constraint(model, P_aux_PAp_wedge[ii in i], (PAp[ii,ii]) - ((1 + PAR[:tau_Ap][(ii,ii)])*PA[ii]) ⟂ PAp[ii,ii])

    # Stub equations for increasing-returns firm variables.
    # In the static competitive benchmark these equal their calibrated values.
    # Nfirm = 1 (one representative firm per sector).
    @constraint(model, P_Nfirm[ii in i], (Nfirm[ii]) - (1.0) ⟂ Nfirm[ii])
    # Fixed-cost labor and capital per firm are zero in the competitive specification.
    @constraint(model, P_LFd[ll in l, ii in i], (LF_d[ll,ii]) - (0.0) ⟂ LF_d[ll,ii])
    @constraint(model, P_KFd[ii in i], (KF_d[ii]) - (0.0) ⟂ KF_d[ii])

    # ── Sector-specific variable stubs ───────────────────────────────────────
    # Variables declared for all i but only used in sector subsets.
    # Sectors outside the relevant subset get stub equations fixing them to 0 (quantities)
    # or 1 (prices) so the system remains square.

    # Crop nest variables (P_10..P_30): only active for cr; zero out for lv+ip.
    noncr = [x for x in i if !(x in cr)]
    @constraint(model, P_HKTEF_stub[ii in noncr, vv in v],  (HKTEF[ii,vv])  - (0.0) ⟂ HKTEF[ii,vv])
    @constraint(model, P_PHKTEF_stub[ii in noncr, vv in v], (PHKTEF[ii,vv]) - (1.0) ⟂ PHKTEF[ii,vv])
    @constraint(model, P_fert_stub[ii in noncr],             (fert[ii])      - (0.0) ⟂ fert[ii])
    @constraint(model, P_Pfert_stub[ii in noncr],            (Pfert[ii])     - (1.0) ⟂ Pfert[ii])

    # Livestock nest variables (P_31..P_54): only active for lv; zero out for cr+ip.
    nonlv = [x for x in i if !(x in lv)]
    @constraint(model, P_KTEL_stub[ii in nonlv, vv in v],  (KTEL[ii,vv])  - (0.0) ⟂ KTEL[ii,vv])
    @constraint(model, P_PKTEL_stub[ii in nonlv, vv in v], (PKTEL[ii,vv]) - (1.0) ⟂ PKTEL[ii,vv])
    @constraint(model, P_TFD_stub[ii in nonlv, vv in v],   (TFD[ii,vv])   - (0.0) ⟂ TFD[ii,vv])
    @constraint(model, P_PTFD_stub[ii in nonlv, vv in v],  (PTFD[ii,vv])  - (1.0) ⟂ PTFD[ii,vv])
    @constraint(model, P_feed_stub[ii in nonlv],            (feed[ii])     - (0.0) ⟂ feed[ii])
    @constraint(model, P_Pfeed_stub[ii in nonlv],           (Pfeed[ii])    - (1.0) ⟂ Pfeed[ii])

    return model
end
