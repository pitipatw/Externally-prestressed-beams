# This code is for analyzing linear elastic "cracked" section 
# and potentially nonlinear creacked section as well.
"""
Adopted from 
    "Flexural behavior of externally prestressed beams Part 1: Analytical models"
    Chee Khoon Ng, Kiang Hwee Tan.
"""
# Setting up packages
using CSV, DataFrames
using UnPack
using Makie, GLMakie

# Setting up the data
begin
include("input_data.jl")
include("functions.jl")
include("Interpolate.jl")
end


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

#setup test values
st = 50.0
P_lb = 0:st:4000  #[lb]
P_N  = 4.448*P_lb # [N]
P = P_N
#given M
M = P*Ls/2.0

#set up history containers
begin
fps_history = zeros(length(P))
dps_history = zeros(length(P))
Icr_history = zeros(length(P))
Ie_history = zeros(length(P))
c_history = zeros(length(P))
dis_history = zeros(length(P))
end
#Assume
begin
Icr = Itr
fps = fpe
dps = dps0
Ω =  getOmega(Sec)
#we could do Mcr = 0 , becuase we crack at the begining anyway. 
Mcr = getMcr(Mat, Sec, f, Ω)
# Mcr = 10.0
Ie = Icr
end

#These lines just to make the variables global
Ωc = 0
c  = 0
Ac_req  = 0 
Lc = 0
fc = 0.0
δ_mid = 0

begin
title_name = [ "dps", "fps", "DisMid", "c", "Inertia(s)"]
fig_monitor = Figure(resolution = (1200, 2000))
axs_monitor = [Axis(fig_monitor[i, 1], title = title_name[i]) for i in 1:5]
# workflow follows fig 7 in the paper.

conv1 = 1
counter1 = 0
counter2 = 0
for i in eachindex(M)
    Mi = M[i] 

    println("M = ", Mi)
    Lc = getLc(Sec,Mcr,Mi)
    # println(Lc)
    # break
    counter1 = 0
    conv1 = 1
    while conv1 > 1e-6
        counter1 += 1 
        if counter1 > 1000
            println("Warning: 1st iteration did not converge")
            break
        end
        # println("HI")
        #assume value of Itr and fps

        conv2 = 1
        counter2 = 0
        while conv2 > 1e-6
            # println("counter")
            counter2 += 1 
            if counter2 > 1000
                println("Warning: 2nd iteration did not converge")
                break
            end

            Ωc = getΩc(Ω, Icr, Lc, Sec)
            ps_force_i = Ωc*Aps*fps
            Ac_req = ps_force_i /0.85/fc′
            c = get_C(Ac_req)
            #calculate Icr
            Icr_calc = get_Icrack(c)

            conv2 = abs(Icr_calc - Icr)/Icr_calc

            Icr = Icr_calc
            
        end
       
        # println("Icr = ", Icr)
        # println("Ac_req ", Ac_req)
        # println("c: ", c)
        # @show Mcr , Mdec, Mi , Icr, Itr
        Ie = getIe(Mcr, Mdec, Mi, Icr, Itr)
        # println("Ie/Icr" , Ie/Icr)
        Δ, δ_mid , e = getDelta(Mat, Sec, f, Ie, Mi, em,fps)
        dps = dps0 - Δ
        fc = fps/Eps/Ωc*c/(dps-c)
        println(fc)
        @assert fc <= 0.003
        @show fps_calc = getFps2(Mat, Sec, f , Ωc, c, dps, fc)
        conv1 = abs(fps_calc - fps) / fps
        fps = fps_calc
        #plot convergence of fps, icr and dps using Makie

    end

    
    # δmid = getDeltamid()
    #record the history
    fps_history[i] = fps
    dps_history[i] = dps
    Icr_history[i] = Icr
    Ie_history[i] = Ie
    c_history[i] = c
    dis_history[i] = δ_mid
end

scatter!(axs_monitor[1], P, dps_history, color = :red)
scatter!(axs_monitor[2], P, fps_history, color = :red)
scatter!(axs_monitor[3], P, dis_history, color = :red)
scatter!(axs_monitor[4], P, c_history, color = :red)
scatter!(axs_monitor[5], P, Ie_history, color = :red ,label = "Ie")
scatter!(axs_monitor[5], P, Icr_history, color = :blue, label= "Icr")
axislegend()
display(fig_monitor)

end




#compare the result with the test data.
begin
df = CSV.File(joinpath(@__DIR__,"pixelframe_beam1.csv"))
df = DataFrame(df)
test_P = df[!,2]
test_d = df[!,3]

dis_in = dis_history/25.4

figure2 = Figure(resolution = (800, 600))
ax1 = Axis(figure2[1, 1], ylabel = "Load [lb]", xlabel = "Displacement [in]")
ax2 = Axis(figure2[2, 1], ylabel = "fps[MPa]", xlabel = "Displacement [in]")
plot!(ax1,dis_in[1:end-30],P_lb[1:end-30], label = "calc", color = :blue)
plot!(ax1,test_d,test_P, label = "test", color = :red)
plot!(ax2, dis_in, fps_history, label = "fps", color = :blue)
display(figure2)
axislegend()
#plot 

end