"""
    Stress at Ultimate in unbonded post tensioning tendons for simply supported beams
    A State-of-the-Art Review
The included models are 
    (1) ACI 318-19 (in the paper was originally ACI 318-83)
        [ready]
    (2) Pannell (1969) # No non-prestressing steel
    (3) Tam and Pannell (1976) # Partially prestressed beams span/depth = 20 to 45
    (4) Du and Tao (1985) # 26 beams under third point loading with span/depth = 19.1
    (5) Harajli (1990)
    [n,n₀ def](6) Harajli and Kanj (1991) # was proposed to replace ACI318-83
    [n,n₀ def](7) Harajli and Hijazi (1991)
    (8) Naaman and Alkhairi (1991) *** this one is famouse, a lot of papers use this eq.
    (9) Chakrabarti (1995)
    [debugging, f](10)Li-Hyung Lee et al (1999)
    [check Original paper](11)Au and Du (2004)
"""
# input variables
# r : S/dₚₛ or dₚₛ
# todo 
# sort these function, 
# rename as the numbers above 
# specify inputs
# unit conversion from ... to psi world. and psi world back!
begin

"""
ACI 318M-18 ref 20.3.2.4

r : span/depth ratio
"""
function getfps1(fpe::Float64, fc′::Float64, ρ::Float64, fpu::Float64, fpy::Float64, r::Float64)
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

"""
Pannell 1969

Lₚ = 10*c this is to check the assumption\\
plastic zone (Lp) = 10*compression depth (c)\\
b could be 2x the width of the PixelFrame
"""
function getfps2(fpe::Float64, fc′::Float64, Aps::Float64, ρₚ::Float64, ϵcu::Float64, Eps::Float64, b::Float64, dₚₛ::Float64, L::Float64, )
    Ψ = 10
    λ  = Ψ * ρₚ * ϵcu * Eps * dₚₛ / (L * fc′)
    qₑ = Aps * fpe / (b * dₚₛ * fc′)
    qᵤ = (qₑ + λ) / (1 + 2λ)
    fps = qᵤ / ρₚ * fc′
    return fps
end

"""
Tam and Pannell (1976), Modified Pannell (1969)\\
-span to depth ratio of 20-45\\
-Modified Pannell 1969
"""
function getfps3(fpe::Float64, fc′::Float64, Aps::Float64, ρₚ::Float64, ϵcu::Float64, β₁::Float64, Eps::Float64, b::Float64, dₚₛ::Float64, L::Float64, )

    Ψ = 10
    α = 0.85 * β₁ # β₁ depends on strain
    qₑ = Aps * fpe / (b * dₚₛ * fc′)
    # qₛ = As*fy/(b*dₚₛ*fc′)
    qₛ = 0 # for no non-prestressed steel
    λ  = Ψ * ρₚ * ϵcu * Eps * dₚₛ / (L * fc′)
    #Lp = 10.5*c #this is for checking the assumption
    
    fps = fc′ * ((qₑ + λ) / (1 + λ / α) - qₛ * λ / (α + λ)) / ρₚ
    return fps

end

"""
Du and Tao 1985
"""
function getfps4(fpe::Float64, ρₚ::Float64, fpy::Float64)
    # MPa    
    # Du and Tao 1985
        q₀  = clamp(ρₚ, 0., 0.3)
        fps = fpe + 785. -1920. *q₀
        fps = clamp(fps , 0. , fpy)
        if fpe < 0.55fpy || fpe > 0.65fpy 
            println("Du and Tao Assumption failed")
        end
    return fps
end

"""
Harajli (1990)
psi
"""
function getfps5(fpe::Float64, fc′::Float64, ρₚ::Float64, r::Float64, fpy::Float64)
    fps = min(fpe + (10000. + fc′ / 100ρₚ) * (0.4 + 8r) , fpe + 60000. , fpy) #psi
    return fps
end

"""
Harajli and Kanj (1991)
"""
function getfps6(fpe::Float64, fc′::Float64, As::Float64, fpu::Float64, b::Float64, r::Float64, dₚₛ::Float64)
    # psi
    γ₀ = n₀ / n / (0.12 + 2.5r)
    A = clamp((Aps * fpe + As * fy) / (b * dₚₛ * fc′), 0. , 0.23)
    fps = fpe + γ₀ * fpu * (1.0 - 3A) #psi
    fps = clamp(fps, fpe + 0.3γ₀ * fpu, fps)

return fps
end

"""
Harajli and Hijazi (1991)\\
Partially prestressed members\\
Non-linear analysis
"""
function getfps7(fpe::Float64, fpu::Float64, fpy::Float64, S::Float64, dₚₛ::Float64, c::Float64)
    ###Have to check n0 and n in their paper.

    if point == 0 #uniform loading
        γₛ = n₀ / n * (0.25 + 1.2 / (S / dₚₛ)) #S/dₚₛ => span to depth ratio
        β₀ = 1.75
    elseif point == 1 #single point loading
        γₛ = n₀ / n * (0.1 + 2 / (S / dₚₛ))
        β₀ = 1.8
    elseif point == 3 
        γₛ = n₀ / n * (0.4 + 1.1 / (S / dₚₛ))
        β₀ = 1.75
    end

    fps = clamp(fpe + γₛ * fpu * (1 - β₀ * c / dₚₛ), 0., fpy) #psi
    return fps
end


"""
Naaman and Alkhairi (1991)

L₁ : Length of loaded span or sum of lengths of loaded spans ,influenced by same tendon
L₂ : length of tendon between end anchorages
"""
function getfps8(fpe::Float64, dₚₛ::Float64, L::Float64, Eps::Float64, fpy::Float64, ϵcu::Float64, c::Float64, L1::Float64, L2::Float64,p::Int8)

    if p == 1.  # one point loading
        Ωᵤ = 1.5dₚₛ / L
    elseif p == 3. # third point loading
        Ωᵤ = 3dₚₛ / L
    end
    fps = clamp(fpe + Ωᵤ * Eps * ϵcu * (dₚₛ / c - 1) * L1 / L2, 0., 0.94fpy) # psi

    return fps
end

"""
Chakrabarti (1995)

Non-prestressing steel -> ds is undefined -> use dp/ds = 1 ?
subscript "s" is for steel, "p" is for prestressed
"""
function getfps9(fpe::Float64, fc′::Float64, ρₚ::Float64, ρₛ::Float64, dₚₛ::Float64, dₛ::Float64, fy::Float64, r::Float64)
    # psi
    A = clamp(fc′ / (100. * ρₛ) * dₚₛ / dₛ * 60000. / fy * (1 + ρₛ / 0.025),0. , 20000.)

    #define R
    r <= 33 ? R = 1.0 : R = 0.8

    B = clamp(R * fc′ / (100. * ρₚ * fpe), 0., 0.25)

    if ρₛ == 0. && r > 33.
        fps = fpe + ((fpe + 10000. + A) / (1. - B) - fpe) * 0.65
    else
        fps = fpe + (fpe + 10000. + A) / (1. - B) #psi
    end

    if r <= 33. 
        fps = clamp(fps, 0., fpe + 60000.)
    else
        fps = clamp(fps, 0., fpe + 40000.)
    end

return fps
end

"""
Li-Hyung Lee et al. (1990)

what is 'f'?
"""
function getfps10(fpe::Float64, fc′::Float64, As′::Float64, As::Float64, fy::Float64,Aps::Float64,dₛ::Float64, dₚₛ::Float64, ρₚ::Float64,f::Float64,L::Float64)
    # psi
    fps = 10000. + 0.8fpe + (As′ - As) * fy / (15. * Aps) + 80. * sqrt(dₛ / dₚₛ * fc′ / ρₚ * (1 / f + dₚₛ / L)) #psi
return fps
end

"""
Au, F.T.K and Du, J.S (2004)
Le : span length between the anchorage divided by number of plastic hinges
"""
function getfps11(fpe::Float64, fc′::Float64, Aps::Float64, As::Float64,fpy::Float64,fy::Float64,Eps::Float64,β₁::Float64,b::Float64, dₚₛ::Float64,Lₑ::Float64)
    # MPa
    cₚₑ = clamp((Aps * fpe + As * fy) / (0.85β₁ * fc′ * b), 0. , fpy) #very weird cap value, might have to recheck.
    fps = clamp(fpe + 0.0279Eps * (dₚₛ - cₚₑ) / Lₑ, 0., fpy) #MPa
    return fps 
end



"""
outdated
function getfpsfpe::Float64, fc′::Float64, ρ::Float64, fpu::Float64, fpy::Float64, r::Float64)
    # this is in psi 
    #need to look at Ωᵤ definition, which has many version
    # ACI-ASCE 2002
    fps = fpe + Ωᵤ*Eps*ϵcu*(dₚₛ/c-1)*L1/L2 # psi

    return fps
end
"""
end

function callAllmodel(fc′::Float64, fpe::Float64, fpy::Float64, fpu::Float64,fy::Float64,r::Float64,Ac::Float64, Aps::Float64,As::Float64,As′::Float64,Eps::Float64,dₚₛ::Float64,dₛ::Float64,c::Float64,β₁::Float64,b::Float64,L::Float64,L1::Float64,L2::Float64)
    ϵcu = 0.003
    ρ = (Aps+As+As′)/Ac
    ρₚ = Aps/Ac
    ρₛ = As/Ac
    #working on these
    n = 1.
    n₀ = 1.
    p::Int8 = 1
    f = 1.
    Lₑ = 200.
    f6 = 0.
    f7 = 0.
    f1 = getfps1(fpe::Float64, fc′::Float64, ρ::Float64, fpu::Float64, fpy::Float64, r::Float64)
            # MPa
            # ref 20.3.2.4 in ACI 318M-19
            #possible stress values
            # r: span/depth ratio
    f2 = getfps2(fpe::Float64, fc′::Float64, Aps::Float64, ρₚ::Float64, ϵcu::Float64, Eps::Float64, b::Float64, dₚₛ::Float64, L::Float64, )
            # MPa
            # Pannell 1969
            # Lₚ = 10*c this is to check the assumption, 
            # plastic zone (Lp) = 10*compression depth (c).
            # b could be 2x the width of the PixelFrame
    f3 = getfps3(fpe::Float64, fc′::Float64, Aps::Float64, ρₚ::Float64, ϵcu::Float64, β₁::Float64, Eps::Float64, b::Float64, dₚₛ::Float64, L::Float64, )
            # MPa
            # Tam and Pannell 1976 
            # Modified Pannel 1969
            # span to depth ratio of 20 to 45
            # Modified Pannell 1969
    f4 = getfps4(fpe::Float64, ρₚ::Float64, fpy::Float64)
            # MPa    
            # Du and Tao 1985
    
    # turn variables into US units

    f5 = getfps5(fpe::Float64, fc′::Float64, ρₚ::Float64, r::Float64, fpy::Float64)
            # psi
            # Harajli (1990) 
    #f6 = getfps6(fpe::Float64, fc′::Float64, As::Float64, fpu::Float64, b::Float64, r::Float64, dₚₛ::Float64)
            # psi
            #Harajli an Kanj (1991) 
    #f7 = getfps7(fpe::Float64, fpu::Float64, fpy::Float64, S::Float64, dₚₛ::Float64, c::Float64)
            # psi
            # Harajli and Hijazi (1991)
            # Partially prestressed members
            # Non linear analysis
    f8 = getfps8(fpe::Float64, dₚₛ::Float64, L::Float64, Eps::Float64, fpy::Float64, ϵcu::Float64, c::Float64, L1::Float64, L2::Float64,p::Int8)
            # psi
            # Naaman and Alkhairi (1991) 
            # L1 = length of loaded span or sum of lengths of loaded spans, influenced by the same tendon
            # L2 = leng of tendon between end enchorages.
    f9 = getfps9(fpe::Float64, fc′::Float64, ρₚ::Float64, ρₛ::Float64, dₚₛ::Float64, dₛ::Float64, fy::Float64, r::Float64)
            # psi
            # Chakrabarti (1995) 
            # non-prestressing steel-> ds is undefined. -> use dp/ds = 1 ? 
            # sub script "s" is for steel, "p" is for prestressed
    f10 = getfps10(fpe::Float64, fc′::Float64, As′::Float64, As::Float64, fy::Float64,Aps::Float64,dₛ::Float64, dₚₛ::Float64, ρₚ::Float64,f::Float64,L::Float64)
            # psi
            # Li-Hyung Lee et al. (1999)
            # what is 'f' ? 
    f11 = getfps11(fpe::Float64, fc′::Float64, Aps::Float64, As::Float64,fpy::Float64,fy::Float64,Eps::Float64,β₁::Float64,b::Float64, dₚₛ::Float64,Lₑ::Float64)
            # MPa
            # Au, F.T.K and Du, J.S (2004)
            # Le: span length between the anchorage divided by number of plastic hinges

    (f5,f8,f9,f10) = 0.00689476.*(f5,f8,f9,f10)
   return [f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11] 
end
