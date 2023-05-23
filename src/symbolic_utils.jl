function make_symbolic_1d(x_vals, y_vals, x_symb)
    """
    This is a helper function to assist numerical -> symbolic
    """
    dx = x_vals[2] - x_vals[1]
    xl = minimum(x_vals)
    xh = maximum(x_vals)
    # Make the interpolation
    itp = interpolate(y_vals, BSpline(Cubic(Line(OnGrid()))))
    # Call it with the x vals and extrapolate
    sym_itp = extrapolate(scale(itp, xl:dx:xh), Line())
    # Register the function with the symbolic term
    return sym_itp(x_symb)
end

function make_symbolic_2d(x_vals, y_vals, z_vals, x_symb, y_symb)
    """
    This is a helper function to assist numerical -> symbolic for 2D
    """
    dx = x_vals[2] - x_vals[1]
    xl = minimum(x_vals)
    xh = maximum(x_vals)
    dy = y_vals[2] - y_vals[1]
    yl = minimum(y_vals)
    yh = maximum(y_vals)
    # Make the interpolation using the sampled data
    itp = interpolate(z_vals, BSpline(Cubic(Line(OnGrid()))))
    # Call it with the x and y vals
    sym_itp = extrapolate(scale(itp, xl:dx:xh, yl:dy:yh), Line())
    # Register the function with the symbolic terms
    return sym_itp(x_symb, y_symb)
end

function make_discretised_1d(f_target, sym_args, num_x, num_args)
    """
    This is a helper function to assist symbolic -> numerical
    """
    # Turn the symbolic expression into a callable function
    myf = eval(build_function(f_target, sym_args))
    # Loop over the x values and evaluate the function
    num_target = [Base.invokelatest(myf, [num_x[i], num_args...]) for i = 1:length(num_x)]
    # Return the result vector
    return num_target
end

function make_discretised_2d(f_target, sym_args, num_x, num_y, num_args)
    """
    This is a helper function to assist symbolic -> numerical for 2D
    """
    # Turn the symbolic expression into a callable function
    myf = eval(build_function(f_target, sym_args))
    # Loop over the x values and evaluate the function
    num_target = [Base.invokelatest(myf, [num_x[i], num_y[j], num_args...]) for i = 1:length(num_x), j = 1:length(num_y)]
    # Return the result vector
    return num_target
end

function latexify_nice(eq; f_compress=true, rm_params=true)
    expr = latexify(eq)
    # Replace the d with partial
    expr = replace(expr, "\\mathrm{d}" => "\\partial ")
    # Replace the h slash with h h_bar
    expr = replace(expr, "\\hslash" => "\\hbar")
    # Fix any problems with brackets
    expr = replace(expr, "\\left( q \\right)" => "(q)")
    expr = replace(expr, "\\left( p \\right)" => "(p)")
    expr = replace(expr, "\\left( t \\right)" => "(t)")

    # Replace L with the nice operator
    expr = replace(expr, "L" => "\\hat{\\mathcal{L}}_{\\mathrm{QM}}")

    # Fix fix_indexes
    expr = fix_indexes(expr)

    if rm_params
        expr = replace(expr, "(q, p, t)" => "")
        expr = replace(expr, "\\left( q, p, t \\right)" => "")
        expr = replace(expr, "(q, p)" => "")
        expr = replace(expr, "(q)" => "")
        expr = replace(expr, "(p)" => "")
        expr = replace(expr, "(t)" => "")
    end

    if f_compress
        # fix the partial potential
        expr = replace(expr, "\\frac{\\partial  V(q)}{\\partial q}" => "\\partial_q V(q)")
        expr = replace(expr, "\\frac{\\partial  V}{\\partial q}" => "\\partial_q V")

        # Replace the dt with partial t
        expr = replace(expr, "\\frac{\\partial }{\\partial t}" => "\\partial_t")
        # Replace the dq with partial q
        expr = replace(expr, "\\frac{\\partial }{\\partial q}" => "\\partial_q")
        # Replace the dp with partial p
        expr = replace(expr, "\\frac{\\partial }{\\partial p}" => "\\partial_p")

        # Replace the repeated partials
        expr = replace(expr, "\\partial_q \\partial_q \\partial_q" => "\\partial_q^3")
        expr = replace(expr, "\\partial_q \\partial_q" => "\\partial_q^2")
        expr = replace(expr, "\\partial_p \\partial_p \\partial_p" => "\\partial_p^3")
        expr = replace(expr, "\\partial_p \\partial_p" => "\\partial_p^2")
    end
    return expr
end

function fix_indexes(input::AbstractString; n=3)
    # https://en.wikibooks.org/wiki/Introducing_Julia/Strings_and_characters
    if n == 2
        pattern = r"}ˏ_"
        pattern1 = r"_{(\d+)_(\d+)"
        tmp1 = replace(input, pattern => "_")
        tmp2 = replace(tmp1, pattern1 => SubstitutionString("_{\\g<1>,\\g<2>}"))
    elseif n == 3
        pattern = r"ˏ_"
        pattern1 = r"_{(\d+)}_{(\d+)}_(\d+)"
        tmp1 = @. replace(input, pattern => "_")
        tmp2 = @. replace(tmp1, pattern1 => SubstitutionString("_{\\g<1>,\\g<2>,\\g<3>}"))
    end
    return tmp2
end

function eq_simple(eq)
    # Simplify the equation, expand the terms and remove large fractions
    return simplify(expand_derivatives(eq), simplify_fractions=false, expand=true)
end

function eq_inserter(eq, sub, expr)
    # Substitute the expression for the sub in the equation
    eq_out = substitute(eq, Dict(sub => expr))
    # Simplify the equation
    return eq_simple(eq_out)
end