function plot_statevars(sim::AbstractDataFrame, vars::Vector{Symbol}; linecolor = :black, kwargs...)

    ncols = 4
    nrows = Int(ceil(length(vars)/4))

    plt = plot(
        layout = (nrows, ncols), 
        bottommargin = 5mm, leftmargin = 5mm, 
        size = (850,300*nrows); 
        kwargs...
        )

    for (i,var) in enumerate(vars)
        plot!(plt, subplot=i, sim.t, sim[:,var], xlabel = "t", ylabel = var, color = linecolor)
    end

    return plt
end
