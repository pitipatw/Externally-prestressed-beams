# variable assignments
begin
    #Inputs
    #..........Notes..........
    # Use Ld = Ls (this test only) 
    # Eccentricities measured from the neutral axis
    # M is the moment in the constant region
    # Mg = moment due to the selfweight
    # M(x) is the moment equation due to the load
    #Units N, mm, MPa

    #   Material Properties
    fc′ = 36.0 # Concrete strength [MPa] ****Should update on the test day using cylinder test***
    # Ec = 4700.0*sqrt(fc′) # MPa  ACI fc-> Concrete modulus relationship [MPa]
    Ec = 58000.0 # MPa  from the cylinder test
    Eps = 70000.0 #Post tensioning steel modulus [MPa]
    fpy = 0.002 * Eps #MPa  
    #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
    # is ~ 150 MPa. Currently 140 MPa :)

    # PixelFrame section/element properties
    centroid_to_top = 91.5 #[mm]
    em = 228.9 # Eccentricity at the middle of the member [mm]
    es = 0.0 # Eccentricity at the support of the member   [mm]
    em0 = em # Initial eccentricity at the midspan        [mm]
    dps0 = centroid_to_top + em0 # Initial distance from the top to the point of application of the load [mm]
    Ls = 502.7 # Distance from support to the first load point [mm]
    Ld = Ls    # Distance from support to the first deviator [mm]
    L = 2000.0 # Total length of the member [mm]
    # two 1/4" bars with 1200 lb capacity
    Aps = 2.0 * (0.25 * 25.4)^2 * pi / 4.0 # Total area of the post tensioned steel [mm2]
 
    Atr = 18537.69 # Transformed area of the cross section [mm2]
    Itr = 6.4198e+07 #moment of inertia [mm4]
    # Itr = 1.082e+8 

    Zb = Itr/centroid_to_top # Elastic modulus of the concrete section from the centroid to extreme tension fiber [mm3]
    # If there are multiple materials, transformed section geometry is needed for Zb (and everything related to section area)


    #forces
    w = Atr / 10^9 * 2400.0 * 9.81 # Selfweight [N/mm]
    mg = w * L^2 / 8.0 # Moment due to selfweight [Nmm]
    fr = 0.7 * sqrt(fc′) # Concrete cracking strenght [MPa]
    r = sqrt(Itr / Atr) # Radius of gyration [mm]
    ps_force = 890.0/sin(24.0*pi/180.0)
         # Post tensioning force [N]
    concrete_force = ps_force*cos(24.0*pi/180.0)
    fpe = ps_force/Aps # Effective post tensioning stress [MPa] ***will input the one on the test day***
    ϵpe = fpe / Eps # Effective post tensioning strain [mm/mm]
    #find moment due to the applied force.
    ϵce = ps_force*em/Zb/Ec - concrete_force/Atr/Ec # effetive strain in the concrete [mm/mm]
end

begin
    # Constructing types
    # Material Properties
    struct Material
        fc′::Float64 # Concrete strength [MPa] ****Should update on the test day using cylinder test***
        Ec::Float64 # MPa  ACI fc-> Concrete modulus relationship [MPa]
        Eps::Float64 #Post tensioning steel modulus [MPa]
        fpy::Float64 #MPa  
        #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
        # is ~ 150 MPa. Currently 140 MPa :)
    end

    struct Section
        em::Float64 # Eccentricity at the middle of the member [mm]
        es::Float64 # Eccentricity at the support of the member   [mm]
        em0::Float64 # Initial eccentricity at the midspan        [mm]
        dps0::Float64 # Initial distance from the top to the point of application of the load [mm]
        Ls::Float64 # Distance from support to the first load point [mm]
        Ld::Float64 # Distance from support to the first deviator [mm]
        L::Float64 # Total length of the member [mm]
        # two 1/4" bars with 1200 lb capacity
        Aps::Float64 # Total area of the steel in the section [mm^2]
        Atr::Float64 # Transformed area of the cross section [mm^2]
        Itr::Float64 # Moment of inertia of the transformed cross section [mm^4]
        Zb::Float64 # Section modulus of the concrete section from the centroid to extreme tension fiber [mm^3]


    end

    struct Loads
        w::Float64 # Selfweight [N/mm]
        mg::Float64 # Moment due to selfweight [Nmm]
        fr::Float64 # Concrete cracking strenght [MPa]
        r::Float64 # Radius of gyration [mm]
        #ps_force::Float64 # Post tensioning force [N]
        fpe::Float64 # Effective post tensioning stress [MPa]
        ϵpe::Float64 # Effective post tensioning strain [-]
        ϵce::Float64 # Effective concrete strain [-]

    end

end

# Create structs
Mat = Material(fc′, Ec, Eps, fpy)
Sec = Section(em, es, em0, dps0, Ls, Ld, L, Aps, Atr, Itr, Zb)
f   = Loads(w, mg, fr, r, fpe, ϵpe, ϵce)

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
function getFps2(Mat::Material , Sec::Section, f::Loads , Ωc::Float64, c::Float64, dps::Float64)

    @unpack fc′, Ec, Eps, fpy = Mat
    @unpack em, es, em0, Aps, Atr, Itr, Zb = Sec
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
(20)
"""
function getLc()
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
function getΩc( Mat::Material, Sec::Section, Lc::Float64)
    @unpack fc′, Ec, Eps, fpy = Mat
    @unpack em, es, em0, Aps, Atr, Itr, Zb = Sec
    @unpack w, mg, fr, r, fpe = Loads
    # will have to add more variable to each structure.

    Ω =  getOmega(Sec)
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
            (1 - 2*Ls*L + (1 - es/em)*( L*Lc/(4*Ld*Ls) - Lc^2/(4*Ld*Ls) +
            Lc^3/(12*L*Ld*Ls) - L^2/(12*Ld*Ls) + Ld/L - 2*Ls^2/(3*L*Ld)) +
            es/em*(Lc^2/(4*L*Ls) - L/(4*Ls) + 2*Ld/L -Ls/Ls ))
    else
        println("Warning: Ld is out of range")
    end
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

    return clamp( first_term + second_term, 0, Itr)
end

"""
(27)
This function is similar to linear uncrack
"""

"""
(28)
dps0 : initial effective ost tension dendon depth
"""
function getDps(dps0::Foloat64, Δ::FLoat64)
    K1 = Ls/L-1 
    K2 = 0.0

    dps = dps0 + M*L^2/(6*Ec*Ie)*(3*Ld/L*(-K1) -3/4 - K2) +
        fps*Aps/(Ec*Ie)* e *(L^2/8 - L*Ld/2 + Ld^2/2)
    @assert dps == dps0 + Δ
    return dps
end