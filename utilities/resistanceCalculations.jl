"""
Factored axial capacity of section
"""
function axialCapacity(fc′::Real, Atotal::Real, steelAreas::Vector{<:Real}, fpe::Vector{<:Real}; reductionFactor = 1e3, ϕ_c = 0.65, ϕ_s = 0.8)
    
    #individual forces
    fconcrete = 0.85 * fc′ * Atotal / reductionFactor
    fsteel = sum(fpe .* steelAreas) / reductionFactor

    #unfactored
    Pr = fconcrete - fsteel

    #Factored
    return ϕ_c * ϕ_s * Pr
end
