# This code is for analyzing linear elastic "cracked" section 
# and potentially nonlinear creacked section as well.
"""
Adopted from 
    "Flexural behavior of externally prestressed beams Part 1: Analytical models"
    Chee Khoon Ng, Kiang Hwee Tan.
"""
# Setting up packages
using CSV
using DataFrames
using UnPack
using Makie, GLMakie

# Setting up the data

include("input_data.jl")
include("Interpolate.jl")
"""
Need Mcr equation
"""


# This state is when the concrete's stress is = 0.4fc′
# Since sections are usually under-reinforced, the behavior will govern by the steel yielding. 
# Therefore, the nonlinear behavior of the concrete is neglected.

# ..........Notes..........
# Use Ld = Ls (this test only) 
# Eccentricities measured from the neutral axis
# M is the moment in the constant region
# Mg = moment due to the selfweight
# M(x) is the moment equation due to the load
# Units N, mm, MPa




#iteration procedure starts here. 

P_lb = 0:10:4000  #[lb]
P_N  = 4.448*P_lb # [N]
P = P_N
#given M
M = P*Ls/2.0

#Assume
Icr_old = Itr
fps_old = fpe
dps_old = dps0
# workflow follows fig 7 in the paper.
conv1 = 1
conv2 = 1  
counter = 0 
while conv1 > 1e-6
    #assume value of Itr and fps
    Icr = Icr_old
    fps = fps_old
    dps = dps_old

    while conv2 > 1e-6
        counter += 1 
        if counter > 1000
            println("Warning: iteration did not converge")
            break
        end
        Ωc = getΩc(Mat, Sec, Lc, Icr)
        
        #Neutral axis depth (c) calculation 
        # this might be wrong, need to check
        ps_force = Aps*fps
        Ac_req = ps_force / 0.85/fc′
        c = get_C(Ac_req)
        #calculate Icr
        Icr_calc = getIcrack(c)
        
        conv2 = abs(Icr_calc - Icr)/Icr_calc

        




Icr_calc = 20 #calculate Icr 

#check convergence of Icr and Icr_calc
conv1 = abs(Icr_calc -Icr) / Icr

#eq 26 27 28
Ie = 
e = 
dps = 

fps_calc = getFps2()

conv2 = abs(fps_calc - fps) / fps
fps = fps_calc
end
end

δmid = getDeltamid()

println("δmid = ", δmid)



