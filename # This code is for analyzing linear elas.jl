# This code is for analyzing linear elastic cracked section and nonlinear creacked section

# Setting up packages
using CSV
using DataFrames
using UnPack
using Makie, GLMakie

# Setting up the data

"""
Need Mcr equation
"""
# This state is when the concrete's stress is = 0.4fc′
# Since sections are usually under-reinforced, the behavior will govern by the steel yielding. 
# Therefore, the nonlinear behavior of the concrete is neglected.
# Inputs\
# ..........Notes..........\
# Use Ld = Ls (this test only) \
# Eccentricities measured from the neutral axis\
# M is the moment in the constant region\
# Mg = moment due to the selfweight\
# M(x) is the moment equation due to the load\
# Units N, mm, MPa
# Material Properties
struct Material
    fc′::Float64 # Concrete strength [MPa] ****Should update on the test day using cylinder test***
    Ec::Float64 # MPa  ACI fc-> Concrete modulus relationship [MPa]
    Eps::Float64 #Post tensioning steel modulus [MPa]
    fpy::Float64 #MPa  
    #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
    # is ~ 150 MPa. Currently 140 MPa :)
end

@show fc′= 30. #(*) Concrete strength [MPa] ****Should update on the test day using cylinder test***
@show Ec = 4700.0*sqrt(fc′) # MPa  ACI fc-> Concrete modulus relationship [MPa]
@show Eps = 70000.0 #Post tensioning steel modulus [MPa]
@show fpy = 0.002*Eps #MPa  
#Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
# is ~ 150 MPa. Currently 140 MPa :)

struct Section
    em::Float64 # Eccentricity at the middle of the member [mm]
    es::Float64 # Eccentricity at the support of the member   [mm]
    em0::Float64 # Initial eccentricity at the midspan        [mm]
    Ls::Float64 # Distance from support to the first load point [mm]
    Ld::Float64 # Distance from support to the first deviator [mm]
    L::Float64 # Total length of the member [mm]
    # two 1/4" bars with 1200 lb capacity
    Aps::Float64 # Total area of the steel in the section [mm^2]
    Atr::Float64 # Transformed area of the cross section [mm^2]
    Itr::Float64 # Moment of inertia of the transformed cross section [mm^4]
    Zb::Float64 # Section modulus of the concrete section from the centroid to extreme tension fiber [mm^3]

end

# PixelFrame section/element properties
# Eccentricity is measured from the centroid of the concrete crossection to the centroid of the steels
@show em = 225.76 # Eccentricity at the middle of the member [mm]
# Since the ropes at the supports are at the centroid of the concrete section
@show es = 4.24    # Eccentricity at the support of the member[mm]
@show em0 = em   # Initial eccentricity at the midspan      [mm]
@show Ls = 502.7 # Distance from support to the first load point [mm]
@show Ld = Ls    # Distance from support to the first deviator [mm]
@show L = 2000.0 # Total length of the member [mm]

# Steel properties
    # two 1/4" bars with 1200 lb capacity
@show Aps = 2.0*(0.25*25.4)^2*pi/4.0 # Total area of the post tensioned steel [mm2]
# If there are multiple materials, transformed section geometry is needed for Zb (and everything related to section area)

@show Atr = 18537.69 + 347.96  # Transformed area of the cross section (Concrete + Steel) [mm2]

# v -> possible error.
#previous value : 8.9795e7
@show Itr  = 1.0788e8  # Moment of inertia of the transformed cross section [mm4]
# # Section modulus of the concrete section from the centroid to extreme tension fiber [mm3]
@show c = 137.51 # Distance from the centroid of the entire section to the centroid of the steel section (extreme tension) [mm]
@show Zb = Itr/c

# Its was = 6.4198e+07 #moment of inertia [mm4]
# Zb  was = 452894.24

struct Loads
    w::Float64 # Selfweight [N/mm]
    mg::Float64 # Moment due to selfweight [Nmm]
    fr::Float64 # Concrete cracking strenght [MPa]
    r::Float64 # Radius of gyration [mm]
    #ps_force::Float64 # Post tensioning force [N]
    fpe::Float64 # Effective post tensioning stress [MPa]
end




"""
(19)
"""
function getFps()
    first_term = Eps*(ϵpe + Ωc*ϵc)
    second_term = Ωc*fc*Eps/Ec*(dps/c-1)
    if first_term + second_term > fpy
        return fpy
    else
        return first_term + second_term
    end
end

"""
(20)
"""
function getLc()
    Lc = L - 2*Ls*Mcr/M 
    if Lc < L-2*Ls
        return L-2*Ls
    else
        return Lc
    end
end

"""
(21)
"""
function getΩ2(Ω::Float64, L::Float64,Lc::Float64, Ls::Float64, Ld::Float64,Icr::Float64,Itr::Float64
     if Ld < Ls 
        if (L - 2*Ls) < Lc < (L - 2*Ld)
            Ωc = Ω*Icr/Itr + (1-Icr/Itr)*
                (1 - L/(4*Ls)  + Lc/(2*Ls) - Lc^2/(4*L*Ls) - Ls/L)
        elseif Lc >= L-2*Ld
            Ωc = Ω*Icr/Itr + (1-Icr/Itr)*
                (1 - L/(4*Ls)  + Ld/(2*Ls) - Ld^2/(4*L*Ls) - Ls/L)
        else
            Ωc = Ω
        end
    elseif Ld >= Ls
        Ωc = Ω
    end
    return Ωc
"""