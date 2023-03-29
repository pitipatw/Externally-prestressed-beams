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
    Mdec = ps_force*em
    concrete_force = ps_force*cos(24.0*pi/180.0)
    fpe = ps_force/Aps # Effective post tensioning stress [MPa] ***will input the one on the test day***
    ϵpe = fpe / Eps # Effective post tensioning strain [mm/mm]
    #find moment due to the applied force.
    ϵce = ps_force*em/Zb/Ec - concrete_force/Atr/Ec # effetive strain in the concrete [mm/mm]
end

begin
    # Constructing types
    # Material Properties
    mutable struct Material
        fc′::Float64 # Concrete strength [MPa] ****Should update on the test day using cylinder test***
        Ec::Float64 # MPa  ACI fc-> Concrete modulus relationship [MPa]
        Eps::Float64 #Post tensioning steel modulus [MPa]
        fpy::Float64 #MPa  
        #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
        # is ~ 150 MPa. Currently 140 MPa :)
    end

    mutable struct Section
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

    mutable struct Loads
        w::Float64 # Selfweight [N/mm]
        mg::Float64 # Moment due to selfweight [Nmm]
        fr::Float64 # Concrete cracking strenght [MPa]
        r::Float64 # Radius of gyration [mm]
        #ps_force::Float64 # Post tensioning force [N]
        fpe::Float64 # Effective post tensioning stress [MPa]
        ϵpe::Float64 # Effective post tensioning strain [-]
        ϵce::Float64 # Effective concrete strain [-]
        Mdec::Float64 # decompression moment [Nmm]

    end

end

# Create structs
Mat = Material(fc′, Ec, Eps, fpy)
Sec = Section(em, es, em0, dps0, Ls, Ld, L, Aps, Atr, Itr, Zb)
f   = Loads(w, mg, fr, r, fpe, ϵpe, ϵce, Mdec)