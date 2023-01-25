"""
Return the stress-at-failure value of post-tensioned steel
"""
function getfps(fₚₑ::Real, fc::Real, ρ::Real, fₚᵤ::Real)

    #possible stress values
    temp1 = fₚₑ + 70 + fc / (100 * ρ)
    temp2 = fₚₑ + 420
    temp3 = 1300

    #take minimum
    fₚₛ = min(temp1, temp2, temp3)

    #check value
    if fₚₑ >= 0.5*fₚᵤ && fₚₛ == temp1
        error("Invalid Assumption for Eq.1")
    end

    #return value
    return fₚₛ
end

"""
Get the elastic modulus of concrete
"""
function getEc(fc′::Real, ρ::Real, highStrength::Bool) 
    if highStrength
        Ec = (3300 * sqrt(fc′ + 1000) + 6900) * (ρ / 2000)^1.5
    else
        Ec = 4700sqrt(fc′)
    end

    #check strain at Maximum
    strainMax = 0.85 * fc′ / Ec
    if strainMax > 1e-3
        println("WARNING: Section compression failure at pure axial")
    end

    return Ec
end

"""
Get ϕ₁ factor
"""
function phi1(fc′::Real) 
    ϕ1 = 0.85 - 0.05 * (fc′ - 28) / 7

    return clamp(0.65, ϕ1, 0.85)
end

"""
Get ϕ₂ factor
"""
function phi2(strains::Vector{Float64})
    factors = 0.65 .+ 0.25 .* (strains .- 2e-3) ./ 1e-3
    return clamp.(factors, 0.65, 0.9)
end

"""
Get compression area/depth\\
targetArea: required area for concrete compression region to match\\
sectionMap: vector of cumulative section area from top\\
centroidMap: vector of compression centroids at each slice
"""
function compressioncentroid(targetArea::Real, areaMap::Vector{<:Real}, centroidMap::Vector{<:Real}; tolerance = 3e-3) 
    @assert length(areaMap) == length(centroidMap) "maps must be of same length"

    #compare area map to required area
    comparison = (areaMap .- targetArea) ./ targetArea

    if all(comparison .< 0)
        error("Insufficient concrete area.")
    end
    
    #get index
    i = findfirst(0 .< comparison .< tolerance)

    #verify a result is found
    if typeof(i) != Nothing
        return centroidMap[i]
    else
        ibest = findfirst(comparison .> 0)
        err = comparison[ibest]
        println("Solution within tolerance not found; returning best solution with error $err")
        return centroidMap[ibest]
    end

end

"""
Checks the stresses of top/bottom fibres for transfer state and service state to check if they exceed limits
"""
function servicecheck(ecc::Vector{<:Real}, steelForces::Vector{<:Real}, Atotal::Real, Scompression::Real, Stension::Real, Mdead::Real, Mtotal::Real, fc′::Real) 

    #limits
    compression_limit = -0.5fc′
    tension_limit = 0.25sqrt(fc′)
    
    #global properties
    ext_moment = sum(steelForces .* ecc)
    steelForce = sum(steelForces)

    #stress states
    #transfer states
    TS_top = -steelForce / Atotal + (ext_moment - Mdead) / Scompression
    TS_bot = -steelForce / Atotal - (ext_moment - Mdead) / Stension

    #service state
    SS_top = -steelForce / Atotal + (ext_moment - Mtotal) / Scompression

    #collect
    states = [TS_top, TS_bot, SS_top]

    #check
    check = compression_limit .< states .< tension_limit

    if any(.!check)
        println("Serviceability state not met in one or more criteria")
    end

    return check
end