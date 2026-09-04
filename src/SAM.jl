# Usage examples:
#   data = init_data(); default_sets!(data); setup_sam_accounts!(data)
#   build_default_large_sam!(data)       # creates an internally generated N-sector SAM
#   read_sam_csv!(data, "data/csv/sam.csv")
#   balance_sam_ras!(data)
#
#   data = init_data(); read_sets_csv!(data, "data/csv/sets.csv")   # N != 100 sectors
#
# SAM convention: rows receive payments from columns.
# The SAM has 2N + 16 accounts for N = length(data.sets[:i]) sectors: N activities,
# N commodities, 5 factors, 6 taxes, 4 institutions and 1 margin account.

"""Read the model sets from a two-column `set,item` CSV (see `data/csv/sets.csv`).

Call this *before* `default_sets!` (`prepare_data!(...; sets_path = ...)` does).  Every
group the file supplies is **assigned** into `data.sets`, so the `get!` calls in
`default_sets!` leave it alone; sets the file does not mention (`:r :v :l :h :f :t ...`)
keep their defaults.  Only `:i` is mandatory:

* `:j` and `:k` default to `:i`;
* every item of a sector set must be a member of `:i`;
* `ft` and `fd` must be disjoint from `e` — the `XAp[j,i]` columns are partitioned into
  fertiliser/feed, energy and "other" blocks, so an overlap gives one variable two
  equations and the model stops being square;
* `:ag = cr union lv`, `:ip = i without ag` and `:nf = i without ag` are derived when
  absent, and `:nnft = i without ft`, `:nnfd = i without fd` are always recomputed.
"""
function read_sets_csv!(data::LinkageData, path::AbstractString)
    raw = readdlm(path, ',', String, '\n')
    size(raw,1) >= 2 || error("Sets CSV must include a header row and at least one data row: $(path)")
    size(raw,2) >= 2 || error("Sets CSV must have two columns `set,item`: $(path)")
    header = lowercase.(strip.(String.(raw[1,1:2])))
    header == ["set", "item"] || error("Sets CSV header must be `set,item`, found `$(join(header, ","))`: $(path)")

    groups = Dict{Symbol,Vector{String}}()
    for r in 2:size(raw,1)
        g  = strip(String(raw[r,1]))
        it = strip(String(raw[r,2]))
        isempty(g) && isempty(it) && continue
        (isempty(g) || isempty(it)) &&
            error("Sets CSV row $(r) has an empty set name or item: $(path)")
        push!(get!(groups, Symbol(g), String[]), String(it))
    end

    haskey(groups, :i) || error("Sets CSV must define set `i` (the activity/product list): $(path)")
    allunique(groups[:i]) || error("Sets CSV: set `i` contains duplicate items: $(path)")

    S = data.sets
    for (g, items) in groups
        S[g] = items
    end
    haskey(groups, :j) || (S[:j] = copy(groups[:i]))
    haskey(groups, :k) || (S[:k] = copy(groups[:i]))

    # Every sector set must be drawn from :i.
    pool = Set(groups[:i])
    for g in (:j, :k, :cr, :lv, :ag, :ip, :e, :ft, :fd, :nf, :nnft, :nnfd)
        (haskey(groups, g) || g in (:j, :k)) || continue
        bad = [x for x in S[g] if !(x in pool)]
        isempty(bad) || error("Sets CSV: set `$(g)` contains items that are not in `i`: " *
                              join(bad, ", ") * " ($(path))")
    end

    # Squareness: fertiliser and feed must not also be energy goods.
    for (g, other) in ((:ft, :e), (:fd, :e))
        (haskey(S, g) && haskey(S, other)) || continue
        overlap = intersect(S[g], S[other])
        isempty(overlap) || error("Sets CSV: `$(g)` and `$(other)` overlap in " *
                                  join(overlap, ", ") *
                                  "; the model needs them disjoint to stay square ($(path))")
    end

    if haskey(S, :cr) && haskey(S, :lv)
        haskey(groups, :ag) || (S[:ag] = vcat(S[:cr], S[:lv]))
    end
    if haskey(S, :ag)
        agset = Set(S[:ag])
        haskey(groups, :ip) || (S[:ip] = [x for x in S[:i] if !(x in agset)])
        haskey(groups, :nf) || (S[:nf] = [x for x in S[:i] if !(x in agset)])
    end
    haskey(S, :ft) && (S[:nnft] = [x for x in S[:i] if !(x in Set(S[:ft]))])
    haskey(S, :fd) && (S[:nnfd] = [x for x in S[:i] if !(x in Set(S[:fd]))])
    return data
end

function setup_sam_accounts!(data::LinkageData)
    default_sets!(data)
    i = data.sets[:i]
    accounts = Dict{Symbol,Vector{String}}()
    accounts[:activities]  = ["ACT_" * x for x in i]
    accounts[:commodities] = ["COM_" * x for x in i]
    accounts[:factors]     = ["LAB_UNSK", "LAB_SK", "CAP", "LAND", "NRES"]
    accounts[:taxes]       = ["TAX_OUT", "TAX_INT", "TAX_IMP", "TAX_EXP", "TAX_FACT", "TAX_INC"]
    accounts[:institutions]= ["HH", "GOV", "INV", "ROW"]
    accounts[:margins]     = ["TRD_MRG"]
    accounts[:all] = vcat(accounts[:activities], accounts[:commodities], accounts[:factors],
                           accounts[:taxes], accounts[:institutions], accounts[:margins])
    data.sam_accounts = accounts
    data.sam_index = Dict(a => n for (n,a) in enumerate(accounts[:all]))
    data.metadata[:sam_account_count] = length(accounts[:all])
    return data
end

function _idx(data::LinkageData, account::String)
    return data.sam_index[account]
end

function set_sam!(data::LinkageData, accounts::Vector{String}, mat::AbstractMatrix{<:Real})
    n = length(accounts)
    n == size(mat,1) == size(mat,2) || error("SAM must be square and match account labels.")
    data.sam_accounts[:all] = copy(accounts)
    data.sam_index = Dict(a => i for (i,a) in enumerate(accounts))
    data.sam = Matrix{Float64}(mat)
    data.balanced_sam = copy(data.sam)
    return data
end

function build_default_large_sam!(data::LinkageData)
    setup_sam_accounts!(data)
    acc = data.sam_accounts[:all]
    n = length(acc)
    M = zeros(Float64, n, n)
    S = data.sets
    products = S[:i]

    for (pos, p) in enumerate(products)
        act = _idx(data, "ACT_" * p)
        com = _idx(data, "COM_" * p)
        output = 900.0 + 3.0*pos

        # Activity column: production costs paid by activity p.
        int_total = 0.42 * output
        for (qpos, q) in enumerate(products)
            # Sparse-ish deterministic IO: local block plus small diagonal; still complete enough for calibration.
            weight = (q == p ? 0.18 : (abs(qpos-pos) <= 2 ? 0.015 : 0.0))
            if weight > 0
                M[_idx(data, "COM_" * q), act] += int_total * weight
            end
        end
        # normalize intermediate column to int_total
        current_int = sum(M[_idx(data,"COM_"*q), act] for q in products)
        if current_int > 0
            for q in products
                M[_idx(data,"COM_"*q), act] *= int_total/current_int
            end
        end
        M[_idx(data,"LAB_UNSK"), act] += 0.20 * output
        M[_idx(data,"LAB_SK"),   act] += 0.08 * output
        M[_idx(data,"CAP"),      act] += 0.16 * output
        M[_idx(data,"LAND"),     act] += (p in S[:ag] ? 0.06 : 0.01) * output
        M[_idx(data,"NRES"),     act] += (p in S[:ip] ? 0.03 : 0.01) * output
        M[_idx(data,"TAX_OUT"),  act] += 0.03 * output
        M[_idx(data,"TAX_INT"),  act] += 0.02 * output

        # Commodity column: sales revenue paid to domestic activity, imports, tariffs, and margins.
        M[act, com] += output
        M[_idx(data,"ROW"), com] += 0.08 * output
        M[_idx(data,"TAX_IMP"), com] += 0.015 * output
        M[_idx(data,"TRD_MRG"), com] += 0.02 * output
    end

    # Commodity rows: final demands and exports to match column scale approximately.
    for p in products
        com = _idx(data, "COM_" * p)
        row_now = sum(M[com, :])
        col_now = sum(M[:, com])
        gap = max(col_now - row_now, 1.0)
        M[com, _idx(data,"HH")]  += 0.62 * gap
        M[com, _idx(data,"GOV")] += 0.12 * gap
        M[com, _idx(data,"INV")] += 0.18 * gap
        M[com, _idx(data,"ROW")] += 0.08 * gap
        M[_idx(data,"TAX_EXP"), _idx(data,"ROW")] += 0.002 * gap
    end

    # Factors pay households; taxes pay government; margin account pays transport/service commodities.
    for fac in data.sam_accounts[:factors]
        fi = _idx(data, fac)
        M[_idx(data,"HH"), fi] = sum(M[fi, :])
    end
    for tax in data.sam_accounts[:taxes]
        ti = _idx(data, tax)
        M[_idx(data,"GOV"), ti] = sum(M[ti, :])
    end
    mrg = _idx(data,"TRD_MRG")
    mrg_products = products[max(1, end-19):end]      # == products[81:100] at N = 100
    for p in mrg_products
        M[_idx(data,"COM_"*p), mrg] += sum(M[mrg,:]) / length(mrg_products)
    end

    # Close household/government/investment/ROW roughly; RAS will tighten.
    hh = _idx(data,"HH"); gov = _idx(data,"GOV"); inv = _idx(data,"INV"); row = _idx(data,"ROW")
    hh_gap = sum(M[hh,:]) - sum(M[:,hh])
    if hh_gap >= 0
        M[gov,hh] += 0.35*hh_gap; M[inv,hh] += 0.65*hh_gap
    else
        M[hh,gov] += -hh_gap
    end
    gov_gap = sum(M[gov,:]) - sum(M[:,gov])
    if gov_gap >= 0
        M[inv,gov] += gov_gap
    else
        M[gov,inv] += -gov_gap
    end
    inv_gap = sum(M[inv,:]) - sum(M[:,inv])
    if inv_gap >= 0
        M[row,inv] += inv_gap
    else
        M[inv,row] += -inv_gap
    end
    row_gap = sum(M[row,:]) - sum(M[:,row])
    if row_gap >= 0
        M[inv,row] += row_gap
    else
        M[row,inv] += -row_gap
    end

    data.sam = M
    data.balanced_sam = copy(M)
    validate_sam!(data; require_balanced=false)
    return data
end

"""Error unless the account labels read from a SAM file are exactly the ones
`setup_sam_accounts!` derived from `data.sets[:i]` (order is irrelevant — lookups go
through `sam_index` — but a missing or stray label silently breaks calibration)."""
function _assert_expected_sam_accounts(data::LinkageData, accounts::AbstractVector{<:AbstractString},
                                       path::AbstractString)
    expected = get(data.sam_accounts, :all, String[])
    isempty(expected) && return nothing            # no account list to check against yet
    got = Set(String.(accounts))
    Set(expected) == got && return nothing
    missing_acc = [a for a in expected if !(a in got)]
    extra_acc   = [a for a in accounts if !(a in Set(expected))]
    _show(v) = length(v) <= 10 ? join(v, ", ") : join(v[1:10], ", ") * ", ... ($(length(v)) total)"
    error("SAM file `$(path)` does not carry the expected $(length(expected)) accounts for " *
          "$(length(data.sets[:i])) sectors (2N+16).\n" *
          "  missing from the file: " * (isempty(missing_acc) ? "none" : _show(missing_acc)) * "\n" *
          "  not expected by the model: " * (isempty(extra_acc) ? "none" : _show(extra_acc)) * "\n" *
          "Check that the file matches `data.sets[:i]` (see read_sets_csv! / sets_path).")
end

function read_sam_csv!(data::LinkageData, path::AbstractString)
    raw = readdlm(path, ',', Any, '\n')
    size(raw,1) >= 2 || error("CSV SAM must include header row and at least one data row.")
    col_accounts = strip.(string.(raw[1,2:end]))
    row_accounts = strip.(string.(raw[2:end,1]))
    row_accounts == col_accounts || error("SAM row and column accounts differ.")
    _assert_expected_sam_accounts(data, row_accounts, path)
    mat = zeros(Float64, length(row_accounts), length(col_accounts))
    for r in 1:length(row_accounts), c in 1:length(col_accounts)
        x = raw[r+1,c+1]
        mat[r,c] = x isa Number ? Float64(x) : parse(Float64, strip(string(x)))
    end
    set_sam!(data, String.(row_accounts), mat)
    return data
end


function _excel_col(n::Int)
    s = ""
    while n > 0
        n, r = divrem(n-1, 26)
        s = string(Char(Int('A') + r)) * s
    end
    return s
end

function read_sam_excel!(data::LinkageData, path::AbstractString; sheet::AbstractString="SAM")
    # The account dimension (2N+16) comes from setup_sam_accounts!, i.e. from data.sets[:i];
    # the read range adapts to it, so an N != 100 workbook needs no change here.
    setup_sam_accounts!(data)
    n = length(data.sam_accounts[:all])
    lastcol = _excel_col(n + 1)
    range = "A1:" * lastcol * string(n + 1)
    raw = XLSX.readdata(path, sheet, range)
    col_accounts = strip.(string.(raw[1,2:end]))
    row_accounts = strip.(string.(raw[2:end,1]))
    row_accounts == col_accounts || error("Excel SAM row and column accounts differ.")
    _assert_expected_sam_accounts(data, row_accounts, path)
    mat = zeros(Float64, n, n)
    for r in 1:n, c in 1:n
        x = raw[r+1,c+1]
        mat[r,c] = x isa Number ? Float64(x) : parse(Float64, strip(string(x)))
    end
    set_sam!(data, String.(row_accounts), mat)
    return data
end

function validate_sam!(data::LinkageData; require_balanced::Bool=true, tol::Float64=1e-6)
    M = size(data.balanced_sam,1) > 0 ? data.balanced_sam : data.sam
    n = size(M,1)
    n == size(M,2) || error("SAM is not square: $(size(M)).")
    length(data.sam_accounts[:all]) == n || error("SAM account count does not match matrix dimension.")
    any(x -> x < -tol, M) && error("SAM contains negative values.")
    if require_balanced
        gap = maximum(abs.(sum(M,dims=1)' .- sum(M,dims=2)))
        gap <= tol || error("SAM is not balanced. Maximum row-column gap = $gap")
    end
    return true
end

function balance_sam_ras!(data::LinkageData; maxiter::Int=2_000, tol::Float64=1e-8)
    # Keep a copy of the raw, unbalanced SAM for diagnostics, then replace
    # the active SAM with the balanced SAM. This guarantees that calibration,
    # initialization, validation, and later model routines all use the same
    # balanced benchmark matrix.
    data.metadata[:raw_sam] = copy(data.sam)

    M = copy(data.sam)
    n = size(M,1)
    row0 = vec(sum(M,dims=2)); col0 = vec(sum(M,dims=1))
    target = 0.5 .* (row0 .+ col0)
    total = sum(target)
    target .*= sum(M) / total
    target .= max.(target, 1.0e-9)
    M .= max.(M, 1.0e-12)
    for iter in 1:maxiter
        M .*= target ./ vec(sum(M,dims=2))
        M .*= (target ./ vec(sum(M,dims=1)))'
        if maximum(abs.(vec(sum(M,dims=2)) .- vec(sum(M,dims=1)))) < tol
            data.metadata[:sam_balance_iterations] = iter
            break
        end
    end
    data.balanced_sam = M
    data.sam = copy(M)
    validate_sam!(data; require_balanced=true, tol=1e-5)
    return data
end


# ---------------------------------------------------------------------
# SAM balance diagnostics
# ---------------------------------------------------------------------

"""Return row sums, column sums, and row-minus-column gaps for the active SAM.

By default this uses `data.balanced_sam`, because calibration and initialization
must be based on the balanced benchmark SAM.
"""
function sam_balance_table(data::LinkageData; use_balanced::Bool=true)
    M = use_balanced && size(data.balanced_sam, 1) > 0 ? data.balanced_sam : data.sam
    n = size(M, 1)
    n == size(M, 2) || error("SAM is not square: $(size(M)).")
    length(data.sam_accounts[:all]) == n || error("SAM account count does not match matrix dimension.")

    accounts = data.sam_accounts[:all]
    row_sum = vec(sum(M, dims=2))
    col_sum = vec(sum(M, dims=1))
    gap = row_sum .- col_sum

    return DataFrame(
        account = accounts,
        row_sum = row_sum,
        column_sum = col_sum,
        row_minus_column = gap,
        abs_gap = abs.(gap),
    )
end

"""Return a compact dictionary with SAM balance diagnostics."""
function sam_balance_summary(data::LinkageData; use_balanced::Bool=true)
    tbl = sam_balance_table(data; use_balanced=use_balanced)
    total_abs_gap = sum(tbl.abs_gap)
    max_gap, max_pos = findmax(tbl.abs_gap)
    total_sam = sum(use_balanced && size(data.balanced_sam, 1) > 0 ? data.balanced_sam : data.sam)

    return Dict{Symbol,Any}(
        :balanced => max_gap <= 1.0e-6,
        :max_abs_gap => max_gap,
        :max_gap_account => tbl.account[max_pos],
        :total_abs_gap => total_abs_gap,
        :total_sam => total_sam,
        :relative_max_gap => total_sam <= 0 ? Inf : max_gap / total_sam,
    )
end

"""Throw an error if the active balanced SAM is not balanced within tolerance."""
function assert_balanced_sam!(data::LinkageData; tol::Float64=1.0e-6)
    validate_sam!(data; require_balanced=true, tol=tol)
    data.metadata[:sam_balance_summary] = sam_balance_summary(data)
    return true
end

"""Export SAM balance diagnostics to CSV files.

Files written:
- `sam_balance_table.csv`
- `sam_balance_summary.csv`
"""
function export_sam_balance_report!(data::LinkageData; outdir::AbstractString="results", tol::Float64=1.0e-6)
    mkpath(outdir)
    tbl = sam_balance_table(data)
    summary = sam_balance_summary(data)

    open(joinpath(outdir, "sam_balance_table.csv"), "w") do io
        println(io, "account,row_sum,column_sum,row_minus_column,abs_gap")
        for r in eachrow(tbl)
            println(io, join((r.account, r.row_sum, r.column_sum, r.row_minus_column, r.abs_gap), ","))
        end
    end

    open(joinpath(outdir, "sam_balance_summary.csv"), "w") do io
        println(io, "item,value")
        for key in sort(collect(keys(summary)); by=string)
            println(io, string(key), ",", summary[key])
        end
        println(io, "tolerance,", tol)
    end

    assert_balanced_sam!(data; tol=tol)
    return summary
end


"""Export the active balanced SAM to a CSV file with row and column labels."""
function export_balanced_sam!(data::LinkageData; path::AbstractString=joinpath("results", "balanced_sam.csv"))
    M = size(data.balanced_sam, 1) > 0 ? data.balanced_sam : data.sam
    validate_sam!(data; require_balanced=true, tol=1.0e-5)
    mkpath(dirname(path))
    accounts = data.sam_accounts[:all]
    open(path, "w") do io
        println(io, ",", join(accounts, ","))
        for (r, acc) in enumerate(accounts)
            println(io, acc, ",", join(M[r, :], ","))
        end
    end
    return path
end
