function calculate_fwhm(x,y)
    y_hmax = maximum(y) * 0.5
    return x[findlast(y .>= y_hmax)] - x[findfirst(y .>= y_hmax)]
end