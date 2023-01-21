# using Plots
# x = range(0,10000,length=100)
# y = get_fps.(x,1000.,0.02,1000)
# plot(x,y)

# inputs from Grasshopper
# Probably going to use JSON, (Dictionary mapping?)
# Concrete Properties
fc′= 56 #MPa
ρ_concrete = 23 #kg/m³
high_strength = false
ϵ_cu = 0.003 # Maximum compressive strain for Concrete
ϵ_ce = 0.005 # Maximum tensile strain for Concrete :Obsolete

# Steel Properties
# use tuple for mutiple rods
PTsteelᵣ= (9,0,0) #dia mm 
PTsteel_area = pi.*(PTsteelᵣ.^2)./4

# Postion of the vertical axis relative section's centroid
# aka e : eccentricity [ratio]
PTpos   = (80,-5,-5) # mm 
# "Should" be symmetrical such that e2 = e3 = e1*cos60
# Indexing
# ................
# ..x(2)...x(3)...
# ................
# ......x(1)......
fpe     = (-800,-0,-0) #MPa
fpu     = (1860,1860,1860) #MPa
steel_modulus = 200000 #MPa

# load detail
DL = 0.0070 # [N/mm²] = 7kPa
LL = 0.0048 # [N/mm²] = 4.8kPa
# buidling detail
BayDepth = 2000 # [mm]
SpanLength = 3000 # [mm]

if high_strength
    temp1 = (3300*(fc′+1000)^0.5 +6900)*(ρ_concrete/2300)^1.5
else
    temp1 = 4700*(fc′^0.5)
end
concrete_modulus = temp1

# section properties
Total_area = 200000 #mm2
centroid = 0 #origin in y axis
I = 2000000 #mm4
ct = 20   ; cb = 40  #mm
st = I/ct ; sb = I/cb #mm



###### NO INPUTs AFTER THIS LINE ######
## Calculation process
# Start of Pure Axial Capacity Calculation
Ø₁ = 0.85 - 0.05*(fc′-28)/7
if Ø₁ > 0.85 ; Ø₁ = 0.85
elseif Ø₁ <0.65 ; Ø₁ = 0.65 ;end
Force_concrete = 0.85*fc′*Total_area/1000 #[kN]
Force_steel = sum(fpe.*PTsteel_area)/1000 #[kN]
Axial_capacity = Force_concrete + Force_steel #[kN]
Factored_AxialCap = 0.65*0.8*Axial_capacity #[kN]

concrete_strain = 0.85*fc′/concrete_modulus #[-]
if concrete_strain > 0.003
    println("WARNING Section Compression Failure at Pure Axial")
end

println("Factored Axial Capacity: $Factored_AxialCap kN")

# End of Pure Axial Capacity Calculation

# Start of Moment Capacity Calculation
# input section divisions started from the top.
section_division = range(1,10, length = 901)
section_centroid = range(1,50, length = 901)
fps₁ = 1000  ; fps₂ = 300 ; fps₃ = 300 
fps = (fps₁,fps₂,fps₃) 
Ext_steel_force = sum(fps.*PTsteel_area)

target = Ext_steel_force
n = size(section_division, 1)
println(n)
tol = 1
i = 0
target_area = Aps*fps/0.85/fc′
check = false
for i = 1:n
    i +=1 
    comp_area = section_division[i]
    tol = abs(comp_area - target_area)/target_area
    if tol<0.003
        check = true
        break
    end
end

if !check
    println("Warning Insufficient Xsection Size")
end


centroid = section_centroid[i] 
arm = PTpos .- centroid # [mm] tuple of lenght of moment arm
strain = fps/E
temp1 = 0.65 .+ 0.25 .*(strain.-0.002)./0.003
n = size(temp1,1)
for i in 1:n
    temp2 = temp1[i]
    if temp2 >0.9
        Ø₂[i]= 0.9
    elseif temp2<0.65
        Ø₂[i]= 0.65
    else
        Ø₂[i]= temp2
    end
end

Factored_Moment = sum(Ø₂.*arm.*fps.*PTsteel_area)/1e6

# 2 and 3 might be discarded, i.e. fps₁ and fps₂ = 0 

# End of Moment capacity Calculation

# Constraints check

#Transfer state (TS)
comp_limit  = -0.5*fc′
ten_limit = 0.25*fc′^0.5
TS_top = -Force_steel/Aconcrete + (P1*e1  + P2e2+P3e3 -Mdead)/st
TS_bot = -Force_steel/Aconcrete + (-P1*e1 + P2e2+P3e3 +Mdead)/sb
#Service state (SS)
SS_top = -Force_steel/Aconcrete + (P1*e1  + P2e2+P3e3 -Mtotal)/st

constraint = (TS_top, TS_bot, SS_top)
c_check = Array{Float64}(undef,3,1)
if i in constraints 
    if comp_limit <= i <= ten_limit
        c_check[i] = true
    else
        c_check[i] = false
    end
end

        
