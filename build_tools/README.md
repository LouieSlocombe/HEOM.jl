# Development tools

This directory contains a Julia environment for JuliaFormatter and Coverage.
Keeping it separate from the root project avoids adding developer tools to the
package's runtime dependencies. Aqua and Test belong to the root project's
test target and are installed automatically by `Pkg.test()`.

Install Julia with [Juliaup](https://julialang.org/downloads/), then run all of
these commands from the repository root:

```bash
julia --project=build_tools -e 'using Pkg; Pkg.instantiate()'
julia --project=build_tools build_tools/format.jl --check
julia --project=. -e 'using Pkg; Pkg.test()'
```

Use `build_tools/format.jl --fix` to rewrite Julia files with the shared
`.JuliaFormatter.toml` settings. `--check` reports unformatted files and fails
without modifying them. The pre-commit hook uses the same check command.

For a fresh coverage report:

```bash
julia --project=build_tools -e 'using Coverage; foreach(clean_folder, ("src", "ext"))'
julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'
julia --project=build_tools build_tools/coverage.jl
```

Coverage is checked for the package source in `src/` and extensions in `ext/`.
The gate fails on missing coverage or uncovered executable lines and writes an
LCOV report to `lcov.info`.
The test suite and development scripts are not part of that coverage percentage.

To update the tools within their compatibility bounds:

```bash
julia --project=build_tools -e 'using Pkg; Pkg.update()'
```

Review formatting changes after tool updates, and rerun the checks. Dependabot
updates GitHub Actions; update Julia dependency bounds in the relevant
`Project.toml` files when adopting new major versions. This library template
does not commit generated manifests.

Git hooks are optional and need a separate installation of
[pre-commit](https://pre-commit.com/#installation). Once it is on your `PATH`:

```bash
pre-commit install
pre-commit run --all-files
```
