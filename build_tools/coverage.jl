using Coverage

cd(dirname(@__DIR__)) do
    # Count uncalled function bodies too, including files without a .cov report.
    coverage = withenv("DISABLE_AMEND_COVERAGE_FROM_SRC" => "no") do
        vcat(process_folder("src"), process_folder("ext"))
    end
    LCOV.writefile("lcov.info", coverage)
    covered, total = get_summary(coverage)
    if total == 0
        println(stderr, "No source coverage found. Run Pkg.test(coverage=true) first.")
        exit(1)
    end

    println(
        "Source line coverage: $covered/$total (",
        round(100 * covered / total; digits = 2),
        "%)",
    )
    for file in coverage
        missing = findall(count -> count === 0, file.coverage)
        if !isempty(missing)
            println(stderr, file.filename, ": uncovered lines ", join(missing, ", "))
        end
    end
    if covered != total
        println(stderr, "Source line coverage must be 100%. See lcov.info for details.")
        exit(1)
    end
end
