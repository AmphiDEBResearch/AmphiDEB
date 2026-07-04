
"""
Smooth absolute function. 
Parameter λ controls the "smoothness".
"""
function sabs(x; λ = 1e-6)
    return sqrt(x^2 + λ)
end

"""
Smooth maximum function.
Parameter λ controls the "smoothness" via sabs.
"""
function smax(a, b; λ = 1e-6)
    return (a+b+sabs(a-b; λ = λ))/2
end
