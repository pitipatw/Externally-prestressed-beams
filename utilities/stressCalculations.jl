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

function getfps1(fpe::Float64, fc′::Float64, ρ::Float64, fpu::Float64, fpy::Float64, r::Float64)
    # MPa
    # ref 20.3.2.4 in ACI 318M-19
    #possible stress values
    # r: span/depth ratio
    if r <= 35
        temp1 = fpe + 70 + fc′ / (100 * ρ)
        temp2 = fpe + 420
    else
        temp1 = fpe + 70 + fc′ / (300 * ρ)
        temp2 = fpe + 210
    end
    #take minimum
    fps = min(temp1, temp2, fpy)

    #check assumption
    if fpe >= 0.5*fpu && fps == temp1
        println("Invalid Assumption for Eq.1")
    end

    return fps
end

function getfps2(fpe::Float64, fc′::Float64, Aps::Float64, ρₚ::Float64, ϵcu::Float64, Eps::Float64, b::Float64, dₚₛ::Float64, L::Float64, )
    # MPa
    # Pannell 1969
    # Lₚ = 10*c this is to check the assumption, 
    # plastic zone (Lp) = 10*compression depth (c).
    # b could be 2x the width of the PixelFrame
    Ψ = 10
    λ  = Ψ*ρₚ*ϵcu*Eps*dₚₛ/(L*fc′)
    qₑ = Aps*fpe/(b*dₚₛ*fc′)
    qᵤ = (qₑ + λ)/(1+2*λ)
    fps = qᵤ/ρₚ*fc′
    return fps
end

function getfps3(fpe::Float64, fc′::Float64, Aps::Float64, ρₚ::Float64, ϵcu::Float64, β₁::Float64, Eps::Float64, b::Float64, dₚₛ::Float64, L::Float64, )
    # MPa
    # Tam and Pannell 1976 
    # Modified Pannel 1969
    # span to depth ratio of 20 to 45
    # Modified Pannell 1969
    Ψ = 10
    α = 0.85*β₁ # β₁ depends on strain
    qₑ = Aps*fpe/(b*dₚₛ*fc′)
    # qₛ = As*fy/(b*dₚ*fc′)
    qₛ = 0 # for no non-prestressed steel
    λ  = Ψ*ρₚ*ϵcu*Eps*dₚₛ/(L*fc′)
    Lp = 10.5*c #this is for checking the assumption
    
    fps = fc′*((qₑ+λ)/(1+λ/α)-qₛ*λ/(α+λ))/ρₚ
    return fps

end

function getfps4(fpe::Float64, ρₚ::Float64, fpy::Float64)
    # MPa    
    # Du and Tao 1985
        q₀  = clamp(ρₚ, 0., 0.3)
        fps = fpe + 785. -1920. *q₀
        fps = clamp(fps , 0. , fpy)
        if fpe < 0.55*fpy || fpe > 0.65*fpy 
            println("Du and Tao Assumption failed")
        end
    return fps
end

function getfps5(fpe::Float64, fc′::Float64, ρₚ::Float64, S::Float64, fpy::Float64, dₚ::Float64)
    # psi
    # Harajli (1990) 
    fps = min(fpe + (10000. + fc′/(100. *ρₚ))*(0.4+8. /(S/dₚ)) , fpe + 60000. , fpy) #psi
    return fps
end

function getfps6(fpe::Float64, fc′::Float64, ρ::Float64, fpu::Float64, fpy::Float64, r::Float64)
    # psi
    #Harajli an Kanj (1991) 
    γ₀ = n₀/n/(0.12+2.5/(S/dₚ))
    A = clamp((Aps*fpe+As*fy)/(b*dₚ*fc′) , 0 , 0.23)
    fps = fpe + γ₀*fpu*(1.0-3.0*A) #psi
    fps = clamp(fps, fpe + 0.3 * γ₀*fpu , fps)

return fps
end

function getfps7(fpe::Float64, fpu::Float64, fpy::Float64, S::Float64, dₚ::Float64, c::Float64)
    # psi
    # Harajli and Hijazi (1991)
    # Partially prestressed members
    # Non linear analysis

    """
    Have to check n0 and n in their paper.
    """
    if point == 0 #uniform loading
        γₛ = n₀/n*(0.25+ 1.2/(S/dₚ)) #S/dₚ => span to depth ratio
        β₀ = 1.75
    elseif point == 1 #single point loading
        γₛ = n₀/n*(0.1 +   2/(S/dₚ))
        β₀ = 1.8
    elseif point == 3 
        γₛ = n₀/n*(0.4 + 1.1/(S/dₚ))
        β₀ = 1.75
    end

    fps = clamp(fpe + γₛ*fpu*(1-β₀*c/dₚ) , 0. , fpy) #psi
    return fps
end

function getfps8(fpe::Float64, dₚₛ::Float64, L::Float64, Eps::Float64, fpy::Float64, ϵcu::Float64, c::Float64, L1::Float64, L2::Float64)
    # Naaman and Alkhairi (1991) 
    # L1 = length of loaded span or sum of lengths of loaded spans, influenced by the same tendon
    # L2 = leng of tendon between end enchorages.
    if p == 1.  # one point loading
        Ωᵤ = 1.5*dₚₛ/L
    elseif p == 3. # third point loading
        Ωᵤ = 3. *dₚₛ/L
    end
    fps = clamp(fpe + Ωᵤ*Eps*ϵcu*(dₚₛ/c-1)*L1/L2, 0 , 0.94*fpy) # psi

    return fps
end

function getfps9(fpe::Float64, fc′::Float64, ρₚ::Float64, ρₛ::Float64, dₚ::Float64, dₛ::Float64, fy::Float64)
    # psi
    # Chakrabarti (1995) 
    # non-prestressing steel-> ds is undefined. -> use dp/ds = 1 ? 
    # sub script "s" is for steel, "p" is for prestressed
    A = clamp(fc′/(100. *ρₛ)*dₚ/dₛ*60000. /fy*(1+ρₛ/0.025) , 0 , 20000.)

    if S/d <= 33 
        r = 1.0
    else
        r = 0.8
    end

    B = clamp(r*fc′/(100. *ρₚ*fpe) , 0. , 0.25)

    if ρₛ == 0. && S/d > 33.
        fps = fpe + ((fpe + 10000. + A)/(1. -B) - fpe) * 0.65
    else
        fps = fpe + (fpe + 10000. + A)/(1. -B) #psi
    end

    if S/d <= 33. 
        fps = clamp(fps,0., fpe + 60000.)
    else
        fps = clamp(fps,0., fpe + 40000.)
    end

return fps
end








function getfps11fpe::Float64, fc′::Float64, ρ::Float64, fpu::Float64, fpy::Float64, r::Float64)
    # Au, F.T.K and DU, J.S (2004)
    cₚₑ = clamp( (Aps*fpe + As*fy)/(0.85*β₁*fc′*b) , 0. , fpy) #very weird cap value, might have to recheck.
    fps = clamp( fpe + 0.0279*Eps*(dₚ-cₚₑ)/lₑ , 0., fpy) #MPa
    return fps 
end

function getfps10fpe::Float64, fc′::Float64, ρ::Float64, fpu::Float64, fpy::Float64, r::Float64)
    # Li-Hyung Lee et al. (1999)
    # what is 'f' ? 
    fps = 10000. + 0.8*fse + (As′-As)*fy/(15*Aps) + 80*sqrt(ds/dp*fc′/ρₚ*(1/f+ dp/L)) #psi
return fps
end



"""
outdated
function getfpsfpe::Float64, fc′::Float64, ρ::Float64, fpu::Float64, fpy::Float64, r::Float64)
    # this is in psi 
    #need to look at Ωᵤ definition, which has many version
    # ACI-ASCE 2002
    fps = fpe + Ωᵤ*Eps*ϵcu*(dₚ/c-1)*L1/L2 # psi

    return fps
end
"""