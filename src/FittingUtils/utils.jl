using DataFrames 
using Latexify
using Distributions

import Plots:plot
plot(hyper::Hyperdist; kwargs...) = plot(hyper.dist; kwargs...)


function assign_value_by_label!(p, label, value)::Nothing

    labels = EcotoxSystems.ComponentArrays.labels(p)
    idx = findfirst(x -> x == label, labels)
    
    p[idx] = value

    return nothing

end


function assign_values_from_file!(
    p::ComponentVector, 
    posterior_summary::AbstractDataFrame; 
    exceptions = OrderedDict()
    )

    for (label,value) in zip(posterior_summary.param, posterior_summary.best_fit)
        if !(label in keys(exceptions))
            assign_value_by_label!(p, label, value)
        else
            exceptions[label](p, label, value)
        end
    end

end

function assign_values_from_file!(
    p::ComponentVector, 
    file::AbstractString; 
    exceptions = OrderedDict()
    )::Nothing

    posterior_summary = CSV.read(file, DataFrame)
    assign_values_from_file!(p, posterior_summary; exceptions)

    return nothing
end

"""
    plot_metam_phase(sim::DataFrame)

Plot state varariables with indication of life stages
"""
function plot_metam_phase(sim::DataFrame)

    plt = @df sim plot(
        plot(:t, :S, color = :black, lw = 1.5, ylabel = "S"), 
        plot(:t, :E_mt, color = :black, lw = 1.5, ylabel = "E_mt"), 
        xlabel = "Time (d)", xlim = (0,100), leg = false
    )

    for (i,var) in enumerate([:S, :E_mt])
        @df sim plot!(subplot = i, :t, :larva .* maximum(sim[:,var]), fill = true, fillalpha = .2, lw = 0, color = :purple, linetype = :stepmid)
        @df sim plot!(subplot = i, :t, :metamorph .* maximum(sim[:,var]), fill = true, fillalpha = .2, lw = 0, color = :steelblue, linetype = :stepmid)
        @df sim plot!(subplot = i, :t, :juvenile .* maximum(sim[:,var]), fill = true, fillalpha = .2, lw = 0, color = :teal, linetype = :stepmid)
    end

    return plt
end

macro h(x)
    quote
        display("text/markdown", @doc $x)
    end    
end


"""
 clean(df::AbstractDataFrame)

Removes all rows with any non-finite and missing values from dataframe. Mostly to clean up data before plotting.
"""
function clean(df::AbstractDataFrame)
    
    valid_idxs = [sum(.!isfinite.(Vector(row)))==0 for row in eachrow(df)]

    return dropmissing(df[valid_idxs,:])

end

"""
Match the order of elements in `a` to the order of elements in `b`, assuming that `a` is a subset of `b`.
"""
function match_order(a::Vector, b::Vector)
    
    index_map = Dict(elem => idx for (idx, elem) in enumerate(a))

    return [a[index_map[elem]] for elem in b]
end


"""
    fround(x; sigdigits=2)
Formatted rounding to significant digits (omitting decimal point when appropriate). 
Returns rounded number as string.

"""
function fround(x; sigdigits=2)
    xround = string(round(x, sigdigits = sigdigits))
    if xround[end-1:end]==".0"
        xround = string(xround[1:end-2])
    end
    return xround
end


function _df_to_tex(df::AbstractDataFrame, fname::AbstractString; colnames::Union{Nothing,Vector{AbstractString}} = nothing)::Nothing

    tex_table = @chain df begin
        !isnothing(colnames) ? rename(_, colnames) : _
        latexify(env = :table, booktabs = true, latex = false, fmt = FancyNumberFormatter(3))   
    end 
    
    open(fname, "w") do f
        write(f, tex_table)
    end

    @info "Writing latex table to $fname"

    return nothing
end

"""
    fdist(dist::Truncated{Normal{Float64}})

Format truncated Normal distribution as String.
"""
function fdist(dist::Truncated{Normal{Float64}})
    return "TN($(fround(dist.untruncated.μ)), $(fround(dist.untruncated.σ)), $(dist.lower), $(dist.upper))"
end

"""
    fdist(dist::Dirac)

Format Dirac distribution as String.
"""
function fdist(dist::Dirac)
    return "Dirac($(fround(dist.value)))"
end

"""
Convert parameter object to table (`DataFrame`).
"""
function _as_table(p::EcotoxSystems.ComponentVector; printtable = true)

    df = DataFrame(
        param = EcotoxSystems.ComponentArrays.labels(p), 
        value = vcat(p...)
    )

    if printtable
        show(df, allrows = true)
    end

    return df
end