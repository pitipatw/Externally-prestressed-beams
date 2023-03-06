"""
Adopted from 
    "Flexural behavior of externally prestressed beams Part 1: Analytical models"
    Chee Khoon Ng, Kiang Hwee Tan.
"""

using ProgressBars
#This function plot the flexural behavior of externally prestressed beams
# Pseudo-section analysis

#Inputs
#Units N, mm, MPa
#   Material Properties
    fc′= 30. # Concrete strength [MPa] 
    Ec = 4700.0*sqrt(fc′) # MPa  ACI fc-> Concrete modulus relationship [MPa]
    Eps = 70000.0 #Post tensioning steel modulus [MPa]
    fpy = 0.002*Eps #MPa  
    #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
    # is ~ 150 MPa.
# PixelFrame section/element properties

    em = 228.9 #mm
    es = 0.
    em0 = em # Initial eccentricity at the midspan
    Ls = 502.7 #mm
    Ld = Ls
    L = 2000.0 #mm
    # two 1/4" bars with 1200 lb capacity
    Aps = 2.0*(0.25*25.4)^2*pi/4.0 #mm2
    fpe = 0.1*fpy # MPa will input the one on the test day.
    Zb = 452894.24 #mm3 elastic modulus
    Atr = 18537.69 # mm2 Transformed area of crosssection.
    Itr = 6.4198e+07 #mm4 moment of inertia 

    #forces
    w = Atr/10^9*2400.0*9.81 # N/mm
    Mg = w*2000.0^2/8.0 # Nmm
    fr = 0.7*sqrt(fc′) #MPa
    r  = sqrt(Itr/Atr) #mm

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
    omega = 1.0 - (es / em)*(Ls/L) + (es - em)/ em * (Ls/L*4/3)
    return omega
end

#calculate the moment at the critical section
function getMcr(em::Float64, Mg::Float64, omega::Float64, Itr::Float64, Ec::Float64, Eps::Float64, Aps::Float64, r::Float64, Zb::Float64, Atr::Float64)
    @show mcre = Aps*fpe*(em + Zb/Atr) + (fr * Zb)
    @assert mcre > 0
    @show dmcr = (Aps*em*(em + Zb / Atr)*(mcre - Mg)) / ((1/omega*Itr*Ec/Eps) + Aps*(r^-em*Zb/Atr))
    @assert dmcr > 0
    mcr = mcre + dmcr
    return mcr
end

#pass these 2 lines first
omega = getOmega(Ls, L, es, em)
mcr = getMcr(em, Mg, omega, Itr, Ec, Eps, Aps, r, Zb, Atr) # Nmm
##
# Under Mcr fps will follow this equation 
# M here is the Moment at the critical section (variable)
# so this is a function between fps and M,e (e is the middle of the beam)

function getDelta(fps::Float64, M::Float64, L::Float64, Ls::Float64, es::Float64, em::Float64, Itr::Float64, Ec::Float64, Aps::Float64)
    #displacement 
    # due to the PS force
    #at mid span
    delta_mid_neg = fps*Aps/(Ec*Itr) * (em * L^2 / 8 - (em-es)*Ls^2/6)
    # at deviator 
    delta_dev_neg = fps*Aps/(Ec*Itr)*(es*Ls^2/6 + em*(L*Ls/2-2/3*Ls^2))

    # due to the applied force
    #at the mid span
    delta_mid_pos = M*L^2/(6*Ec*Itr)*(3/4-(Ls/L)^2)
    # at a deviator
    delta_dev_pos = M*L^2/(6*Ec*Itr)*( 3*(Ls/L)*(1-Ls/L)-(Ls/L)^2)

    delta_mid = delta_mid_pos - delta_mid_neg
    delta_mid_calc = M*L^2/(6*Ec*Itr)*(3/4-(Ls/L)^2) - fps*Aps/(Ec*Itr) * (em * L^2 / 8 - (em-es)*Ls^2/6)
    @assert abs(delta_mid - delta_mid_calc) < 1e-3
    Δ = delta_mid - (delta_dev_pos - delta_dev_neg)
    K1 = Ls/L-1 
    K2 = 0.0
    Δcalc = M*L^2/(6*Ec*Itr)*(3 * (Ls/L) * K1 + 3/4 + K2) - fps*Aps*em/(Ec*Itr)*(L^2/8 - L*Ls/2 + Ls^2/2)

    # @show Δ - Δcalc
    @assert abs(Δ - Δcalc) < 1e-6
    return Δ
end



### START HERE ###

#Evaluate the deflection at the midspan
#Given a value of M(t) as an input (increments) -> get this from load graph. 

# Dummy load as from 0 to 6600 N (1.5kips) with 0.1 increment
P = 0.:10:10600.  # N
M = P*Ls/2.
displacements = zeros(length(M))
displacements_mid_pos = zeros(length(M))
iteration_exceeded = zeros(length(M)) 
fps = fpe # initial guess of the stress in the tendons.
omega = getOmega(Ls, L, es, em)
fps_old = fpe #will have to update at the end as the current fps and use that in the next loop.
max_it = 1000
fps_sub_hist = zeros(max_it)
fps_history = zeros(length(M))
#loop M
for i in ProgressBar(eachindex(M))
    Mi = M[i] #value of the current applied moment
    if Mi <= mcr
        #Better format here
        #print mcr with long format
        println("Moment is less than the critical moment")
        println("M = $Mi, Mcr = $mcr")

    
        #Assume fps = fpe as an initial guess
        check_term = 1.
        tol = 1e-6
   
        it = 0
        while check_term > tol
            it += 1 
            if it > max_it
                println("Maximum number of iterations reached for M = $Mi")
                iteration_exceeded[i] = 1
                break
            end
            #calculate the total deflection
            Δ = getDelta(fps_old,Mi, L, Ls, es, em, Itr, Ec, Aps)

            #find the actual eccentricity of the tendon
            e = em - Δ 

            #calculate the stress in the tendons, limited by fpy
            fps_new = fpe +(omega*Mi*e)/(Itr*Ec/Eps + Aps*(r^2+e^2)*omega)
            if fps_new > fpy 
                println("Exceeds the yielding stress")
                println("Beam reaching Ultimate Stress capacity (fpy) and will fail")
                break
            end
            check_term = abs(fps_new - fps_old)/fps_old
            fps_sub_hist[it] = fps_new
            fps_old = fps_new # this one used for the next loop.
        end

        displacements[i] = getDelta(fps_old,Mi, L, Ls, es, em, Itr, Ec, Aps)
        displacements_mid_pos[i] = Mi*L^2/(6*Ec*Itr)*(3/4-(Ls/L)^2)
        fps_history[i] = fps_old
    elseif Mi > mcr 
        println(Mi) 
        Pi = P[i]
        Pi_kips = Pi/4.44822/1000.
        println("Exceeds the cracking moment at load = $Pi N ($Pi_kips kips)")
        println("Using Linear Crack scheme")
        break
    elseif Mi > My
        println("Exceeds the yielding moment")
        println("Beam reaching Ultimate Moment capacity (Mu) and will fail")
    else
        println("Something is wrong")
    end
end
fps_sub_hist
#plot the deflection
using Plots

plot(P, displacements_mid_pos, label = "Deflection")

plot(P, fps_history, label = "Stress")