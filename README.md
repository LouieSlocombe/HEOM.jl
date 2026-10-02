# HEOM.jl

A fresh starting point for HEOM development, based on
[template_julia](https://github.com/LouieSlocombe/template_julia).
It uses `Project.toml` packaging, JuliaFormatter, Aqua, Julia's `Test` standard
library, full source-line coverage, pre-commit, and GitHub Actions.

The package currently contains the template's example numerical API and greeting
CLI. The previous HEOM implementation has been removed; HEOM functionality will
be developed from this scaffold. There are no external runtime dependencies.

## Requirements

- Julia 1.10 or newer in the 1.x series. Install it with
  [Juliaup](https://julialang.org/downloads/).
- Optional: [pre-commit](https://pre-commit.com/#installation) for Git hooks.

CI tests the minimum supported Julia version and the latest stable Julia on
Linux, macOS, and Windows.

## Quick start

From the repository root, instantiate and precompile the package:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Run the example:

```bash
julia --project=. bin/heom.jl Ada
julia --project=. bin/heom.jl --help
```

Omit the name to print `Hello, World!`. Quote names containing spaces, or use
`--` before a name beginning with a dash. The CLI returns status `2` for invalid
arguments.

Or start Julia with `julia --project=.` and use the library:

```julia
using HEOM

print_hello("Ada")
samples = line(-1.0, 1.0; num=5)
println(samples)
```

`line` returns a `Vector{Float64}` including both interval endpoints. A sample
count of zero returns an empty vector; a count of one returns only the starting
value. Negative counts raise an `ArgumentError`. Public functions have docstrings
available through Julia's help mode, for example `?line`.

To work on this package from another Julia project, use
`Pkg.develop(path="/path/to/HEOM.jl")` in that project's environment.

## Development

Install the separate formatting and coverage tools once:

```bash
julia --project=build_tools -e 'using Pkg; Pkg.instantiate()'
```

Run the local checks from the repository root:

```bash
julia --project=build_tools build_tools/format.jl --check
julia --project=. -e 'using Pkg; Pkg.test()'
```

`Pkg.test()` installs the test-only dependencies and runs behavior tests, CLI
subprocess tests, selected type-inference checks, and Aqua's package quality
checks. Aqua checks issues such as method ambiguities, undefined exports, stale
dependencies, and missing compatibility bounds.

Apply formatting changes with:

```bash
julia --project=build_tools build_tools/format.jl --fix
```

Run the coverage gate used by CI:

```bash
julia --project=build_tools -e 'using Coverage; clean_folder("src")'
julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'
julia --project=build_tools build_tools/coverage.jl
```

The gate requires 100% coverage of executable lines in `src/` and writes
`lcov.info`. Julia reports line coverage, so this is not the Python template's
branch-coverage metric. Clearing previous coverage first prevents stale runs
from hiding missing tests. Coverage reports do not require a hosted service or
secret token.

Install the optional Git hooks once, then run them on all files:

```bash
pre-commit install
pre-commit run --all-files
```

The hooks check file hygiene and Julia formatting. They require `julia` on your
`PATH` and the tooling environment instantiated as above. See
[build_tools/README.md](build_tools/README.md) for environment maintenance.

CI also tests loading this package from a fresh consumer environment. Julia
packages are distributed as source; there is no wheel-building step.

## Project layout

```text
.
├── .github/workflows/ci.yml   # formatting, tests, coverage, consumer smoke test
├── bin/heom.jl               # command-line entry point
├── build_tools/              # separate development tools and scripts
├── src/HEOM.jl               # package module and public exports
├── test/                     # behavior and package quality tests
├── .JuliaFormatter.toml      # formatting rules
└── Project.toml              # package metadata, compatibility, test dependencies
```

## Development starting point

The scaffold comes from `LouieSlocombe/template_julia` at commit
`98142bfc232b5efc2d19d7362a36c7d6a9b05b3e`, adapted to the `HEOM` package name
and existing package UUID.

Replace the example API and tests as HEOM functionality is implemented, keeping
the quality checks passing. Add runtime dependencies with `Pkg.add` and give them
explicit compatibility bounds. Update the supported Julia versions in `[compat]`
and the CI matrix together.

Manifests are ignored because this is a reusable library: its tests
resolve dependencies within the declared compatibility bounds. If you turn it
into an application that needs a fixed environment, commit the relevant
`Manifest.toml` and update `.gitignore` accordingly.

## License

Released under the [MIT License](LICENSE).
