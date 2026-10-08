# Condense benchmark logs into error-vs-depth tables (error relative to max W0),
# with the gain of each closure over the hard cutoff at the same depth.
for file in ARGS
    case = ""
    rows = Dict{Tuple{String,Int},Float64}()
    cost = Dict{Tuple{String,Int},Float64}()
    function flush_case()
        isempty(rows) && return
        depths = sort(unique(last.(keys(rows))))
        println("\n", case)
        println("| depth | none | markov_full | nz_diag | nz_full | none at depth+1 |")
        println("|---:|---:|---:|---:|---:|---:|")
        for d in depths
            e(k, dd = d) = get(rows, (k, dd), NaN)
            f(x) = isnan(x) ? "—" : isfinite(x) ? string(round(x; sigdigits = 2)) : "blew up"
            g(k) = (x = e(k); isfinite(x) ? "$(f(x)) ($(round(e("none") / x; sigdigits = 2))×)" : f(x))
            println("| $d | $(f(e("none"))) | $(g("markov_full")) | $(g("nz_diag")) | $(g("nz_full")) | $(f(e("none", d + 1))) |")
        end
        middle(c) = sort(c)[cld(length(c), 2)]
        ratios(k) = [cost[(k, d)] / cost[("none", d)] for d in depths if haskey(cost, (k, d))]
        println("rhs cost relative to none (median over depths): ",
            join(["$k $(round(middle(ratios(k)); sigdigits = 2))×" for k in ("markov_full", "nz_diag", "nz_full")], ", "))
        empty!(rows); empty!(cost)
    end
    for line in eachline(file)
        if startswith(line, "==")
            flush_case(); case = strip(line, ['=', ' '])
        elseif (m = match(r"^(\w+)\s+(\d+)\s+\d+\s+\d+\s+([\d.]+)\s+\d+\s+[\d.]+\s+(\S+)$", line)) !== nothing
            rows[(m[1], parse(Int, m[2]))] = parse(Float64, m[4])
            cost[(m[1], parse(Int, m[2]))] = parse(Float64, m[3])
        end
    end
    flush_case()
end
