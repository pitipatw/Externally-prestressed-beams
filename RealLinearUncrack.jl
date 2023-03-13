"""
Adopted from 
    "Flexural behavior of externally prestressed beams Part 1: Analytical models"
    Chee Khoon Ng, Kiang Hwee Tan.
"""
# To update
# fc' from actual test
# fpe from the test


using ProgressBars
using CSV
using DataFrames
include("utilities\\stressCalculations.jl")
#This function plot the flexural behavior of externally prestressed beams
# Pseudo-section analysis


#Inputs
#..........Notes..........
# Use Ld = Ls (this test only) 
# Eccentricities measured from the neutral axis
# M is the moment in the constant region
# Mg = moment due to the selfweight
# M(x) is the moment equation due to the load
#Units N, mm, MPa

#   Material Properties
    fc′= 30. # Concrete strength [MPa] ****Should update on the test day using cylinder test***
    Ec = 4700.0*sqrt(fc′) # MPa  ACI fc-> Concrete modulus relationship [MPa]
    Eps = 70000.0 #Post tensioning steel modulus [MPa]
    fpy = 0.002*Eps #MPa  
    #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
    # is ~ 150 MPa. Currently 140 MPa :)

# PixelFrame section/element properties
    em = 228.9 # Eccentricity at the middle of the member [mm]
    es = 0. # Eccentricity at the support of the member   [mm]
    em0 = em # Initial eccentricity at the midspan        [mm]
    Ls = 502.7 # Distance from support to the first load point [mm]
    Ld = Ls    # Distance from support to the first deviator [mm]
    L = 2000.0 # Total length of the member [mm]
    # two 1/4" bars with 1200 lb capacity
    Aps = 2.0*(0.25*25.4)^2*pi/4.0 # Total area of the post tensioned steel [mm2]
    Zb = 452894.24 # Elastic modulus of the concrete section from the centroid to extreme tension fiber [mm3]
    # If there are multiple materials, transformed section geometry is needed for Zb (and everything related to section area)

    Atr = 18537.69 # Transformed area of the cross section [mm2]
    Itr = 6.4198e+07 #moment of inertia [mm4]

    #forces
    w = Atr/10^9*2400.0*9.81 # Selfweight [N/mm]
    mg = w*L^2/8.0 # Moment due to selfweight [Nmm]
    fr = 0.7*sqrt(fc′) # Concrete cracking strenght [MPa]
    r  = sqrt(Itr/Atr) # Radius of gyration [mm]
    ps_force = 890 # Post tensioning force [N]
    fpe = 0.0#ps_force/Aps # Effective post tensioning stress [MPa] ***will input the one on the test day***


"""
Bond reduction coefficient for the linear elastic uncracked regime, Naaman's
"""
function getOmega(Ls::Float64,Ld::Float64, L::Float64, es::Float64, em::Float64)
    omega = 1.0 - (es / em)*(Ls/L) + (es - em)/ em * (Ls^2/(3*L*Ld) +Ld/L)
    return omega
end

"""
#calculate the cracking moment of the beam
"""
function getMcr(Aps::Float64,fpe:: Float64, em::Float64,Zb::Float64, Atr::Float64, mg::Float64, omega::Float64, Itr::Float64, Ec::Float64, Eps::Float64,  r::Float64)
    @show mcre = Aps*fpe*(em + Zb/Atr) + (fr * Zb) # Cracking moment due to initial effective prestress (mcre)
    @assert mcre > 0 # mcre should be positive
    @show dmcr = (Aps*em*(em + Zb/Atr)*(mcre - mg)) / ((1/omega*Itr*Ec/Eps) + Aps*(r^2-em*Zb/Atr)) # Moment due to stress increase in external tendons.
    @assert dmcr > 0
    mcr = mcre + dmcr
    return mcr
end

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
    @assert abs(delta_mid - delta_mid_calc) < 1e-9
    Δ = delta_mid - (delta_dev_pos - delta_dev_neg)
    K1 = Ls/L-1 
    K2 = 0.0
    Δcalc = M*L^2/(6*Ec*Itr)*(3 * (Ls/L) * K1 + 3/4 + K2) - fps*Aps*em/(Ec*Itr)*(L^2/8 - L*Ls/2 + Ls^2/2)

    # @show Δ - Δcalc
    @assert abs(Δ - Δcalc) < 1e-9
    return Δ
end



### START HERE ###

#Evaluate the deflection at the midspan
#Given a value of M(t) as an input (increments) -> get this from load graph. 

# Dummy load as from 0 to 6600 N (1.5kips) with 0.1 increment
# P = 0.:10:10600.  # N
# M = P*Ls/2.

#Actual data
P_lb = 0:10:4000  #[lb]
P_N  = 4.448*P_lb # [N]

P = P_N

M = P*Ls/2.
#Run these 2 lines first, thye dont change during the calculation.
omega = getOmega(Ls,Ld, L, es, em) # [unitless]
mcr = getMcr(Aps, fpe, em, Zb, Atr, mg, omega, Itr, Ec, Eps, r) # [Nmm]

displacements = zeros(length(M))
displacements_mid_pos = zeros(length(M))
iteration_exceeded = zeros(length(M)) 
fps = fpe # initial guess of the stress in the tendons.
fps_old = fpe #will have to update at the end as the current fps and use that in the next loop.
max_it = 10000
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


plot(displacements,P, label = "Deflection", ylabel = "Deflection (mm)", xlabel = "Load (N)")


dis_in = displacements
plot(dis_in,P_lb, ylabel = "Load [lb]", xlabel = "Displacement [in]")
df = CSV.File(joinpath(@__DIR__,"pixelframe_beam1.csv"))
df = DataFrame(df)
test_P = df[!,2]
test_d = df[!,3]
plot(test_d,test_P, label = "Test Data")
plot!(dis_in, P_lb, label = "Model")
plot(P, fps_history, label = "Stress")