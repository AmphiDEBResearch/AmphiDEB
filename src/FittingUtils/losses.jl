# losses.jl
# all loss functions in this file are designed to take a single pair of time-series,
# or pair of scalar daat values as arguments.
# for example, if we have observations from multiple time-resolved treatments and multiple response variables, 
# `a` and `b` are the observed and predicted time-series for a single response variable in a single treatment
# the data weights are (currently) assumed to be assigned on the level of the response variable, rather than for individual observations

# loss functions apply a penalty if the length of the prediction does not match the length of the data
# we cannot simply use length(a) because the predictions have already been matched with data at this point, 
# dropping entries for which not both exist. 
# the penalty applies if we encounter instabilities during ODE solving or a certain stage of the simulation cannot be reached
# by default, the penalty is `Inf`, which will lead to rejection of the parameter sample in the ABC algorithm
# any other definitions of missing_values_penalty should be considered very carefully

missing_values_penalty(nominal_length, actual_length) = actual_length < nominal_length ? Inf : 1 #(((nominal_length)+1)/(actual_length+1))^2

"""
    loss_mse(a::Vector{Float64}, b::Vector{Float64}, weight = 1, nominal_length::Int = length(b))::Float64

Mean squared error with missing values penalty.
"""
function loss_mse(a::Vector{Float64}, b::Vector{Float64}, weight = 1, nominal_length::Int = length(b))::Float64
    return missing_values_penalty(nominal_length, length(b)) * sum(weight * (a .- b).^2)/length(a)
end

function loss_mse_logtransform(a::Vector{Float64}, b::Vector{Float64}, weight = 1, nominal_length::Int = length(b))::Float64

    # negative values are replcaced with NaN
    # NaN will lead to rejection of the particle, and we should only be getting negtaive values for parameter vectors which are unreasonable to begin with, 
    # so we don't expect this to have undesired side-effects

    b[b .< 0] .= NaN 

    return missing_values_penalty(nominal_length, length(b)) * sum(weight * (log.(a .+ 1) .- log.(b .+ 1)).^2)/length(a)
end

"""
    loss_symmbound(a::Vector{Float64}, b::Vector{Float64}, weight = 1, nominal_length::Int = length(b))::Float64

Symmetric bounded loss with missing values penalty.
"""
function loss_symmbound(a::Vector{Float64}, b::Vector{Float64}, weight = 1, nominal_length::Int = length(b))::Float64

    squared_err = (a .- b).^2
    denom = mean(a)^2 + mean(b)^2
    n = length(a)

    return missing_values_penalty(nominal_length, length(b)) * (weight/n) * sum((squared_err ./ denom))
end


#function loss_dtw(a::Vector{Float64}, b::Vector{Float64}, nominal_length::int = length(b))::Float64
#
#
#end


# log mean relative error

function loss_logratio(a::Vector{Float64}, b::Vector{Float64}, weight = 1, nominal_length::Int = length(b)) 
    
    missing_values_penalty(nominal_length, length(b)) * weight/length(a) * sum(abs.(log.((a .+ 1) ./ (b .+ 1))))

end


"""
    nrmsd(a, b)

Normalized root mean squared deviation. 
Normalization is done by the mean of `a` to avoid division by 0 if `length(a)=1`.
"""
function nrmsd(a, b)
    return sqrt(sum((a .- b) .^2)/length(a))/mean(a)
end



"""
    euclidean_distance_fixed_scale(
        a::Vector{Float64}, 
        b::Vector{Float64}, 
        weight::Real = 1, 
        nominal_length::Int = length(b)
    )::Float64


Computes euclidean distance, 
assuming that `a` and `b` are already scaled. 
"""
function euclidean_distance_fixed_scale(
    a::Vector{Float64}, 
    b::Vector{Float64}, 
    weight::Real = 1, 
    nominal_length::Int = length(b)
    )::Float64

    penalty = missing_values_penalty(nominal_length, length(b))
    ed = sum((weight .* sqrt.((a .- b) .^2)))

    return penalty * ed
end


"""
    euclidean_distance_adaptive_scale(
        a::Vector{Float64}, 
        b::Vector{Float64}, 
        scale::Float64,
        weight::Real = 1, 
        nominal_length::Int = length(b)
    )::Float64

Computes euclidean distance with adaptive scaling factor `scale`.
"""
function euclidean_distance_adaptive_scale(
    a::Vector{Float64}, 
    b::Vector{Float64}, 
    scale::Float64,
    weight::Real = 1, 
    nominal_length::Int = length(b)
    )::Float64

    return missing_values_penalty(nominal_length, length(b)) * sum((weight .* ((a .- b)).^2)/scale) 

end