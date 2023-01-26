"""
    Stress at Ultimate in unbonded post tensioning tendons for simply supported beams
    A State-of-the-Art Review
The included models are 
    (1) ACI 318-19 (in the paper was originally ACI 318-83)
    (2) Pannell (1969) # No non-prestressing steel
    (3) Tam and Pannell (1976) # Partially prestressed beams span/depth = 20 to 45
    (4) Du and Tao (1985) # 26 beams under third point loading with span/depth = 19.1
    (5) Harajli (1990)
    (6) Harajli and Kanj (1991) # was proposed to replace ACI318-83
    (7) Harajli and Hijazi (1991)
    (8) Naaman and Alkhairi (1991) *** this one is famouse, a lot of papers use this eq.
    (9) Chakrabarti (1995)
    (10)Li-Hyung Lee et al (1999)
    (11)Au and Du (2004)
"""
# input variables

# todo 
# sort these function, 
# rename as the numbers above 
# specify inputs
# unit conversion from ... to psi world. and psi world back!

function getfps1(fₚₑ::Float64, fc::Float64, ρ::Float64, fₚᵤ::Float64, fpy::Float64, r::Float64)
    # ref 20.3.2.4 in ACI 318M-19
    #possible stress values
    if r <= 35
        temp1 = fₚₑ + 70 + fc / (100 * ρ)
        temp2 = fₚₑ + 420
    else
        temp1 = fₚₑ + 70 + fc / (300 * ρ)
        temp2 = fₚₑ + 210
    end
    #take minimum
    fps = min(temp1, temp2, fpy)

    #check value
    if fₚₑ >= 0.5*fₚᵤ && fps == temp1
        println("Invalid Assumption for Eq.1")
    end

    #return value
    return fps
end

function getfps()
    #MPa
    #Pannell
    Lₚ = 10*c
    Ψ = 10
    λ  = Ψ*ρₚϵcu*Eps*dₚₛ/(L*fc′)
    qₑ = Aps*fpe/(b*dₚ*fc′)
    qᵤ = (qₑ + λ)/(1+2*λ)
    fps = qᵤ/ρₚ*fc′
    return fps
end

function getfps2(fₚₑ::Float64, fc::Float64, ρ::Float64, fₚᵤ::Float64, fpy::Float64, r::Float64, strain::Float64)
    # Tam and Pannell 1976
    # span to depth ratio of 20 to 45
    # Modified Pannell 1969
    α = 0.85*β₁ # β₁ depends on strain
    qₑ = Aps*fpe/(b*dₚ*fc′)
    qₛ = As*fy/(b*dₚ*fc′)
    λ = ψ*r*ϵᵤ*Es*d/L/fcu
    Lp = 10.5*c
    
    fps = fc′*((qₑ+λ)/(1+λ/α)-qₛλ/(α+λ))/ρₚ
    return fps

end

function getfps3()
    # Du and Tao 1985
    q₀  = As/Ac
    fps = fₚₑ+785-1920*q₀
    fps = clamp(fps , 0 , fpy)
return fps
end

function getfps4
    # this is in psi 
    # ACI-ASCE 2002
    fps = fpe + Ωᵤ*Eps*ϵcu*(dₚ/c-1)*L1/L2 # psi

    return fps
end

function getfps5
    # this is in MPa

    cpe = (Aps*fpe+As*fy)/(0.85*β₁*fc′*b)
    fpd = fpe + 0.0279*Eps*(dp-cpe)/Ie

    return
end

function getfps 
    # Chakrabarti (1995) 
    # psi
    A = clamp(fc′/(100.*ρₛ)*dₚ/dₛ*60000./fy*(1+ρₛ/0.025) , 0 , 20000.)

    if S/d <= 33 
        r = 1.0
    else
        r = 0.8
    end

    B = clamp(r*fc′/(100*ρₚ*fpe) , 0. , 0.25

    if ρₛ == 0 && S/d > 33.
        fps = fpe + ((fpe + 10000 + A)/(1-B) - fpe) * 0.65
    else
        fps = fpe + (fpe + 10000. + A)/(1-B) #psi
    end

    if S/d <= 33 
        fps = clamp(fps,0, fpe + 60000)
    else
        fps = clamp(fps,0, fpe + 40000)
    end

return fps
end

function getfps 
    # Naaman and Alkhairi (1991) 
    # L1 = length of loaded span or sum of lengths of loaded spans, influenced by the same tendon
    # L2 = leng of tendon between end enchorages.
    if p ==1  # one point loading
        Ωᵤ = 1.5*dₚₛ/L
    elseif p == 3 # third point loading
        Ωᵤ = 3.*dₚₛ/L
    end

    fps = clamp(fpe + Ωᵤ*Eps*ϵcu*(dₚₛ/c-1)*L1/L2, 0 , 0.94*fpy) # psi

    return fps
end


function getfps 
    #Harajli an Kanj (1991) 
    γ₀ = n₀/n/(0.12+2.5/(S/dₚ))
    A = clamp((Aps*fpe+As*fy)/(b*dₚ*fc′) , 0 , 0.23)
    fps = fpe + γ₀*fₚᵤ*(1.0-3.0*A) #psi
    fps = clamp(fps, fpe + 0.3 * γ₀*fpu , fps)

return fps
end

function getfps
    # Harajli and Hijazi (1991)
    # Partially prestressed members
    # Non linear analysis
if point == 0 #uniform loading
    γₛ = n₀/n*(0.25+ 1.2/(S/dₚ))
    β₀ = 1.75
elseif point == 1 #single point loading
    γₛ = n₀/n*(0.1 +   2/(S/dₚ))
    β₀ = 1.8
elseif point == 3 
    γₛ = n₀/n*(0.4 + 1.1/(S/dₚ))
    β₀ = 1.75
end

    fps = clamp(fpe + γₛ*fpu*(1-β₀*c/dₚ) , 0. , fpy) #psi

function getfps
    # Au, F.T.K and DU, J.S (2004)
    cₚₑ = clamp( (Aps*fpe + As*fy)/(0.85*β₁*fc′*b) , 0. , fpy) #very weird cap value, might have to recheck.
    fps = clamp( fpe + 0.0279*Eps*(dₚ-cₚₑ)/lₑ , 0., fpy) #MPa
    return fps 
end

function getfps 
    # Li-Hyung Lee et al. (1999)
    # what is 'f' ? 
    fps = 10000. + 0.8*fse + (As′-As)*fy/(15*Aps) + 80*sqrt(ds/dp*fc′/ρₚ*(1/f+ dp/L)) #psi
return fps
end




function getfps 
    # Harajli (1990) 
    fps = minimum(fpe + (10000. + fc′/(100*ρₚ))*(0.4+8/(S/dₚ)) , fpe + 60000. , fpy) #psi
    return fps
end
 
    