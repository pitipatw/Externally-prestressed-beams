"""
Factored axial capacity of section
"""
function axialCapacity(fc′::Real, Atotal::Real, steelAreas::Vector{<:Real}, fpe::Vector{<:Real}; 
        reductionFactor = 1e3, 
        ϕ_c = 0.65, 
        ϕ_s = 0.8)
    
    #individual forces
    fconcrete = 0.85 * fc′ * Atotal / reductionFactor
    fsteel = sum(fpe .* steelAreas) / reductionFactor

    #unfactored
    Pr = fconcrete - fsteel

    #Factored
    return ϕ_c * ϕ_s * Pr
end

# a bit unsure about clamp in the case where min > max.
"""
    Adopted from fib2010 Section 7.7.3.2.2: Beams without shear reinforcement
    Notations: 
    γc: partial safety factor for the concrete without fibres;
    k : factor that takes into account the size effect
    d : effective depth of the cross-sectin [mm]
    ρl: longitudinal reinforcement ratio, Area of the longitudinal steel / gross area of concrete
    PTsteel_area is Asl in the code : the longitudinal steel, (non-prestressed + prestressed)
    fFtuk : characteristic value of the ultimate residual tensile strength for Fiber-reinforced Concrete
        when considering wᵤ = 1.5mm, accoring to eq 5.6-6
    fctk   : characteristic value of the tensile strength for the concrete without fibres -- default = 2 MPa (conservative)
    wᵤ : parameter for fiber testing, default = 1.5mm
    CMOD₃ : parameter fore fiber testing, default = 1.5
    fck    : fc' (characteristic value of cylindrical compressive strength)
    σcp :average stress acting on the concrete section for an axial force, due to loadings or prestressing actions
        σcp = NEd/Ac <= 0.2fcd (Ned > 0 for compression)
        fcd : design compressive strength = fc′/1.2
    
    Ashear : area that resist shear, could be a percentage of the total area (%Aconcrete).
"""
function shearCapacity(fc′::Float64, fR1::Float64 , fR3:: Float64 , Ac::Float64,d::Float64, ρl::Float64, ratio::Float64, Ned::Float64;
        fctk = 2., 
        wᵤ = 1.5, 
        CMOD₃ = 1.5, 
        γc = 1.2)

    # calculation starts
    fck = fc′ # fc′ parameter
    fFts = 0.45*fR1  # a fiber parameter
    fFtuk = fFts - wᵤ/CMOD₃*(fFts - 0.5*fR3 + 0.2*fR1) #a fiber parameter
    fFtuk < 0 ? fFtuk = 0 : nothing  # check fFtuk not less than 0
    
    fcd = fc′ / γc #fc' with safety factor
    σcp = min(Ned / Ac, 0.2fcd)

    k = clamp(1+ sqrt(200. / d), 0., 2.) # k <=2

    # Shear per unit area [N/mm²]
    unitV = (0.18 / γc * k * (100 * ρl * (1 + 7.5fFtuk / fctk) * fck)^(1/3) + 0.15σcp) #shear stress (before multiply by the shear area
    # unitV will not be smaller than unitV_min = 0.035*k^1.5*fck^0.5 + 0.15*σcp (below)
    unitV = clamp(unitV, 0.035k^1.5 * fck^0.5 + 0.15σcp , unitV)
    # Shear [N], refer in code as V_Rd,F
    Ashear = ratio * Ac
    V = Ashear * unitV

    Vu = 0.75 * V / 1000 #[kN] %factored shear
    return Vu
end
