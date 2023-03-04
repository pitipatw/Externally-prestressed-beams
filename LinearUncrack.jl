"""
Adopted from 
    "Flexural behavior of externally prestressed beams Part 1: Analytical models"
    Chee Khoon Ng, Kiang Hwee Tan.
"""
#This function plot the flexural behavior of externally prestressed beams
# Pseudo-section analysis

#Inputs
fc′ = 35 #MPa 
fr = 0.7*sqrt(fc′) #MPa
Ec = 
Eps = 

em = 
es = 
Ls = 
Ld = Ls
Aps = 
fpe = 
Zb = # elastic modulus
Atr = #Transformed area of crosssection.
Itr = 
Mg = 
#for Ld == Ls (our test) 
#eccentricity measured from the neutral axis
# es is eccentricity at the support
# em is eccentricity at the midspan at start 
# L is effective span of the beam
# M is the moment in the constant region
# Mg = moment due to the selfweight
# M(x) is the moment equation due to the load
# Itr = moment of inertia of the transformed section
# r = radius of gyration of the transformed section

#create function omega that calculate omega bases on the line 38 
#
#
#
# which type is faster, Real or Float64?
#

function getOmega(Ls::Float64, L::Float64, es::Float64, em::Float64)
    omega = 1-(es/em)+(es-em)/em*(Ls/L*4/3)
    return omega
end

#calculate the moment at the critical section
function getMcr(em::Float64, Mg::Float64, omega::Float64, Itr::Float64, Ec::Float64, Eps::Float64, Aps::Float64, r::Float64, Zb::Float64, Atr::Float64)
    mcre = Aps*fpe*(em + Zb/Atr) + fr*Zb
    dmcr = (Aps*em*(em + Zb / Atr)*(Mcre - Mg)) / ((1/omega*Itr*Ec/Eps) + Aps*(r^-em*Zb/Atr))
    mcr = mcre + dmcr
    return mcr
end

#pass these 2 lines first
omega = getOmega(Ls, L, es, em)
mcr = getMcr(em, Mg, omega, Itr, Ec, Eps, Aps, r, Zb, Atr)

##
# Under Mcr fps will follow this equation 
# M here is the Moment at the critical section (variable)
# so this is a function between fps and M,e (e is the middle of the beam)

function getDelta()
    #displacement 
    # due to the PS force
    #at mid span
    delta_mid_neg = fps*Aps/(Ec*Itr)*(eL^2/8 - (e-es)*Ls^2/6)

    # at deviator 
    delta_dev_neg = fps*Aps/(Ec*Itr)*(es*Ls^2/6 + e*(L*Ls/2-2/3*Ls^2))

    # due to the applied force
    delta_mid_pos = M*L^2/(6*Ec*Itr)*(3/4-(Ls/L)^2)
    delta_dev_pos = M*L^2/(6*Ec*Itr)*( 3*(Ls/L)*(1-Ls/L)-(Ls/L)^2)

    delta_mid = delta_mid_pos - delta_mid_neg
    Δ = delta_mid - delta_dev_pos + delta_dev_neg
    Δcalc = M*L^2/(6*Ec*Itr)*(3/4) - fps*Aps*e/(Ec*Itr)*(L^2/8 - L*Ls/2 + Ls^2/2)
    e = em - ( delta_mid - delta_dev_pos + delta_dev_neg )
    @assert Δ == Δcalc
    return Δ
end



### START HERE ###

#Evaluate the deflection at the midspan
#Given a value of M(t) as an input (increments) -> get this from load graph. 
P = 0.:0.1:100.  # Dummy load as from 0 to 100 with 0.1 increment
Ls = 50. # distance from support to the loading point.

M = P*Ls/2.
displacements = zeros(length(M))
fps = fpe # initial guess of the stress in the tendons.
omega = getOmega(Ls, L, es, em)
fps_old = fpe #will have to update at the end as the current fps and use that in the next loop.
#loop M
for i in eachindex(M)
    Mi = M[i] #value of the current applied moment
    if Mi <= Mcr

    #M0 = 0.1 
    #Assume fps = fpe as an initial guess
    
        check_term = 1.
        tol = 1e-6
        while check_term > tol
            #calculate the total deflection
            Δ = getDelta(Mi, fps0_old, Itr, Ec, Eps, Aps, r)

            #find the actual eccentricity of the tendon
            e = em - Δ 

            #calculate the stress in the tendons, limited by fpy
            fps_new = clamp (fpe +(omega*Mi*e)/(Itr*Ec/Eps+ Aps*(r^2+e^2)*omega) , fpe , fpy)

            check_term = abs(fps_new - fps)/fps

            fps_old = fps_new # this one used for the next loop.
        end

        delta_mid_neg = fps*Aps/(Ec*Itr)*(eL^2/8 - (e-es)*Ls^2/6)
        delta_mid_pos = M*L^2/(6*Ec*Itr)*(3/4-(Ls/L)^2)
        delta_mid = delta_mid_pos - delta_mid_neg
        displacements[i] = delta_mid
    elseif Mi > Mcr && Mi <= My
        println("Exceeds the cracking moment")
        println("Using Linear Crack scheme")
    elseif Mi > My
        println("Exceeds the yielding moment")
        println("Beam reaching Ultimate Moment capacity (Mu) and will fail")
    else
        println("Something is wrong")
    end
    
end

#repeat the process

#After it converges,
#Calculate the deflection

