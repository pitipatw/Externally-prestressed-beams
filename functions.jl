"""
(3) 
Only calculate once
"""
function getMcr(Mat::Material, Sec::Section, f::Loads, Ω::Float64 ) 
    @unpack fc′, Ec, Eps, fpy = Mat
    @unpack Aps, Atr, Itr, Zb, em = Sec
    @unpack w, mg, fr, r, fpe = f

    
    @show mcre = Aps*fpe*(em + Zb/Atr) + (fr * Zb) # Cracking moment due to initial effective prestress (mcre)
    @assert mcre > 0 # mcre should be positive
    @show dmcr = (Aps*em*(em + Zb/Atr)*(mcre - mg)) / ((1/Ω*Itr*Ec/Eps) + Aps*(r^2-em*Zb/Atr)) # Moment due to stress increase in external tendons.
    @assert dmcr > 0
    mcr = mcre + dmcr
    return mcr
end
"""
(12)
"""
function getDelta(Mat::Material, Sec::Section, f::Loads, Itr::Float64, M::Float64, e::Float64, fps::Float64)

    @unpack fc′, Ec, fpy = Mat
    @unpack Aps, Atr, Zb, em, es, Ls, Ld, L = Sec
    @unpack w, mg, fr, r, fpe = f
    
    #displacement 
    # due to the PS force
    #at mid span
    δ_mid⁻ = fps*Aps/(Ec*Itr) * (em * L^2 / 8 - (em-es)*Ls^2/6)
    # at deviator 
    δ_dev⁻ = fps*Aps/(Ec*Itr)*(es*Ls^2/6 + em*(L*Ls/2-2/3*Ls^2))

    # due to the applied force
    #at the mid span
    δ_mid⁺ = M*L^2/(6*Ec*Itr)*(3/4-(Ls/L)^2)
    # at a deviator
    δ_dev⁺ = M*L^2/(6*Ec*Itr)*( 3*(Ls/L)*(1-Ls/L)-(Ls/L)^2)

    δ_mid = δ_mid⁺ - δ_mid⁻

    δ_mid_cal = M*L^2/(6*Ec*Itr)*(3/4-(Ls/L)^2) - fps*Aps/(Ec*Itr) * (e * L^2 / 8 - (e-es)*Ls^2/6)
    # @show Itr
    # @show δ_mid, δ_mid_cal
    @assert abs(δ_mid - δ_mid_cal) < 1e-6
    Δ = δ_mid - (δ_dev⁺ - δ_dev⁻)
    K1 = Ls/L-1 
    K2 = 0.0
    Δcalc = M*L^2/(6*Ec*Itr)*(3 * (Ls/L) * K1 + 3/4 + K2) - fps*Aps*e/(Ec*Itr)*(L^2/8 - L*Ls/2 + Ls^2/2)

    # @show Δ - Δcalc
    @assert abs(Δ - Δcalc) < 1e-9
    e = (em + M*L^2/(6*Ec*Itr)*(3/4-(Ls/L)^2)) / (1 - fps*Aps/(Ec*Itr) * (L^2/8 - L*Ld/2 +Ld^2/2))
    return Δ, δ_mid , e 
end

"""
(19)
"""
function getFps1(Mat::Material , Sec::Section, f::Loads , Ωc::Float64, c::Float64, dps::Float64)

    @unpack fc′, Ec, Eps, fpy = Mat
    @unpack em, es, em0, dps0, Aps, Atr, Itr, Zb = Sec
    @unpack w, mg, fr, r, fpe, ϵpe, ϵce = f

    first_term = Eps*(ϵpe + Ωc*ϵce)
    second_term = Ωc*fc′*Eps/Ec*(dps/c-1)
    fps = first_term + second_term
    if fps > fpy
        return fpy
    else
        return fps
    end
end

"""
(23)
"""
function getFps2(Mat::Material , Sec::Section, f::Loads , Ωc::Float64, c::Float64, dps::Float64, fc::Float64)

    @unpack fc′, Ec, Eps, fpy = Mat
    @unpack em, es, em0, Aps, Atr, Itr, Zb = Sec
    @unpack w, mg, fr, r, fpe, ϵpe, ϵce = f
    #fc is suppose to be the stress in top concrete fiber, 
    # use constitution equation to get the stress in the top concrete fiber
    first_term = Eps*(ϵpe + Ωc*ϵce)
    second_term = Ωc*fc*Eps/Ec*(dps/c-1)
    fps = first_term + second_term
    if fps > fpy
        return fpy
    else
        return fps
    end
end

"""
(20)
"""
function getLc(Sec::Section , Mcr::Float64, M::Float64)
    @unpack Ls, L = Sec
    Lc = L - 2*Ls*Mcr/M 
    if Lc < L-2*Ls
        return L-2*Ls
    else
        return Lc
    end
end

function getOmega(Sec::Section)
    @unpack em, es, Ls, Ld, L = Sec
    Ω = 1.0 - (es / em)*(Ls/L) + (es - em)/ em * (Ls^2/(3*L*Ld) +Ld/L)
    return Ω
end

"""
(21)
"""
function getΩc( Ω::Float64, Icr::Float64, Lc::Float64, Sec::Section)
    @unpack em, es, L, Ld, Ls, Itr = Sec
    # will have to add more variable to each structure.
    # println("In getOmega_c")
    # println("Ld = $Ld")
    # println("Ls = $Ls")
    # println("L = $L")

    if Ld < Ls 
        if (L - 2*Ls) < Lc < (L - 2*Ld)
            Ωc = Ω*Icr/Itr + (1-Icr/Itr)*
                (1 - L/(4*Ls)  + Lc/(2*Ls) - Lc^2/(4*L*Ls) - Ls/L)
        elseif Lc >= L-2*Ld
            Ωc = Ω*Icr/Itr + (1-Icr/Itr)*
                (1 - Ls/L - Ld^2/(L*Ls) + 
                (1 - es/em) * ( L*Lc/(4*Ld*Ls) - Lc^2/(4*Ld*Ls) + 
                Lc^3/(12*L*Ld*Ls) - L^2/(12*Ld*Ls) + 2*Ld^2/(3*L*Ls) )+
                es/em*(Lc/(2*Ls) - L/(4*Ls) - Lc^2/(4*L*Ls) + Ld^2/(L*Ls)))
        else
            println("Warning: Lc is out of range")
        end
    elseif Ld >= Ls
        Ωc = Ω*Icr/Itr + (1-Icr/Itr)*
            (1 - 2*Ls/L + 
            (1 - es/em)*( L*Lc/(4*Ld*Ls) - Lc^2/(4*Ld*Ls) +
            Lc^3/(12*L*Ld*Ls) - L^2/(12*Ld*Ls) + Ld/L - Ls^2/(3*L*Ld)) +
            es/em*(Lc^2/(4*L*Ls) - L/(4*Ls) + 2*Ld/L - Ls/Ls ))
        @assert Ωc > 0
    else
        println("Warning: Ld is out of range")
    end
    # println("Icr", Icr)
    # println("Itr", Itr)
    @assert Ωc > 0 
    return Ωc
end

"""
(25)
"""
function getDeltamid()
    first_term = M*L^2/(6*Ec*Ie)*(3/4 - (Ls/L)^2)
    second_term = fps*Aps/(Ec*Ie)*( e*L^2/8 - (e-es)*Ld^2/6)
    return first_term + second_term
end 

"""
(26)
Mdec is decompression moment
is the moment from externally post tension 
Mdec  = fpe*Aps*em
"""
function getIe(Mcr::Float64, Mdec::Float64, M::Float64, Icr::Float64, Itr::Float64)
    first_term = ((Mcr- Mdec)/(M-Mdec))^3*Itr
    second_term = 1 - ((Mcr- Mdec)/(M-Mdec))^3*Icr
    # @show first_term + second_term
    return clamp( abs(first_term + second_term) , 0, Itr)
end

# """
# (27)
# This function is similar to linear uncrack
# """

"""
(28)
dps0 : initial effective ost tension dendon depth
"""
function getDps(dps0::Float64, Δ::Float64)
    K1 = Ls/L-1 
    K2 = 0.0

    dps = dps0 + M*L^2/(6*Ec*Ie)*(3*Ld/L*(-K1) -3/4 - K2) +
        fps*Aps/(Ec*Ie)* e *(L^2/8 - L*Ld/2 + Ld^2/2)
    @assert dps == dps0 + Δ
    return dps
end