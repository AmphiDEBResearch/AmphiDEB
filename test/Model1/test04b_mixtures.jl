using Pkg; Pkg.activate("test")
using Test
using Distributions
using OrdinaryDiffEq

using Plots, StatsPlots, Plots.Measures
default(leg = false, lw = 1.5)

include("testutils.jl")

using DataFrames, DataFramesMeta
using StatsBase

using EcotoxSystems
import EcotoxSystems: constrmvec

using Revise

using AmphiDEB
@time import AmphiDEB: defaultparams
import AmphiDEB: ComponentVector

@testset "Re-construction of parameter vector with increased number of stressors" begin

    # TODO: should not be necessary to do the re-construction manually. should be possible with macros?
    glb = ComponentVector(
        t_max = 56., # 8 - week simulation
        N0 = 1., # start with single value [] - only possible setting for ODE_simulator
        dX_in = [20., 20.], # food input rate [mg d^-1] - assuming *ad libitum* feeding per default
        k_V = [0., 0.], # dilution rate in the aquatic medium [V] - does not matter for ad libitum conditions
        V_patch = [1., 1.], # simulated volume [d^-1] - does not matter for ad libitum conditions
        T = 293.15, # ambient temperature [K]
        C_W = [0., 0.], # exposure concentrations - currently only constant; given as a matrix. each row is a treatment, each column is a compound
        pathogen_inoculation_dose = 0., # amount of pathogen spores added to aquatic medium [# spores]
        pathogen_inoculation_time = 30., # time-point of pathogen inoculation [t]
        medium_renewals = [0.] # time-points on which media renewals occur; results in removal of spores
    )

    spc = ComponentVector(

        Z = Dirac(1.), # zoom factor
        propagate_zoom = ComponentVector( # parameters affected by zoom factor, and their scaling components
            dI_max_emb = 1/3, # dI ∝ Z^(1/3)
            dI_max_lrv = 1/3, 
            dI_max_juv = 1/3, 
            X_emb_int = 1., # X_emb ∝ Z 
            H_j1 = 1., # H ∝ Z
            H_p = 1., 
            K_X_lrv = 1/3, # K_X = dI/F_max => K_X ∝ dI ∝ Z^(1/3) (F_max is maximum area-specific searching rate)
            K_X_juv = 1/3
        ),
        
        #=
        Physiological baseline (DEB) parameters
        =#

        X_emb_int = 1, # Initial vitellus (≈ dry mass of an egg)
        K_X_lrv = 1.,  # larval half-saturation constant for food uptake
        K_X_juv = 1., # juvenile and adult half-saturation constant for food uptake
        dI_max_emb = 1, # embryonic maximum specific ingestion rate
        dI_max_lrv = 1, # larval maximum specific ingestion rate
        dI_max_juv = 1, # juvenile and adult maximum specific ingestion rate
        kappa_emb = 0.8, # embryonic to metamorph allocation fraction; default value is suggested value from DEBkiss book
        kappa_juv = 0.8, # juvenile and adult allocation fraction to soma; default value is suggested value from DEBkiss book
        gamma = 0.5, # larval allocation fraction to metamorphic reserves (split occurs downstream of κ-split)
        eta_IA = 0.54, # assimilation efficiency; default value from Pfab et al. (2020)
        eta_AS_emb = 0.4, # embryonic to metamorph growth efficiency; default value from Pfab et al. (2020)
        eta_AS_juv = 0.4, # juvenile and adult growth efficiency; initially assumed eval to eta_AS_emb
        eta_AR = 0.95, # reproduction efficiency [-]; default value is suggested value from DEBkiss book and also the default in add-my-pet
        eta_SA = 0.8, # shrinking efficiency
        k_M_emb = 0.11, # embryonic to metamorph somatic maintenance rate constant; default value from Pfab et al. (2020)
        k_M_juv = 0.11, # juvenile and adult somatic maintenance rate constant; initially assumed equal to k_M_emb
        delta_k_M_mt = 1., # metamorph somatic maintenance rate constant, relative to that of embryos nd larvae
        k_J_emb = 0.027, # embryonic to metamorph maturity maintenance rate constant; initially calculated based on (1-kappa)/kappa-ratio
        k_J_juv = 0.027, # juvenile and adult maturity maintenance rate constant; initially assumed equal to k_J_emb
        H_j1 = 1, # maturity at the start of metamorphosis (decline of feeding)
        H_p = 55., # maturity at puberty
        delta_E = 1., # energy density of E_mt, relative to remaining dry mass

        T_A = 8000., # Arrhenius temperature (K) 
        T_ref = 293.15, # reference 
        b_T = 40., # effect strength of temperature on resource allocation

        #=
        TKTD parameters    
        =#

        # TK feedbacks (1.0/0.0 for on/off) 
        # only including growth feedback for now
        fb_G = 0., 

        h_b = 0., # background mortality
        
        KD = [0. 0. 0. 0. 0. 0. 0.; 0. 0. 0. 0. 0. 0. 0.], # k_D - value per PMoA (G,M,A,R,H,kap) and stressor (1 row = 1 stressor)
        B = [2. 2. 2. 2. 2. 2. 2.; 2. 2. 2. 2. 2. 2. 2.], # slope parameters
        E = [1e10 1e10 1e10 1e10 1e10 1e10 1e10; 1e10 1e10 1e10 1e10 1e10 1e10 1e10], # sensitivity parameters (thresholds)
        KD_h = [0.; 0.], # k_D - value for GUTS-Sd module (1 row = 1 stressor)
        E_h = [1e10; 1e10], # sensitivity parameter (threshold) for GUTS-SD module
        B_h = [1.; 1.], # slope parameter for GUTS-SD module 
        C_h = [1.; 1.], # proportionality constant to convert relative response to hazard rate 

        # these are curently only used in an individual-based context, but could find application in the pure-ODE implementation 
        # for example by triggering emptying of the reproduction buffer through callbacks

        S_rel_crit = 0.66, # initial guess on how much body mass can be lost - X. laevis can lose up to 45% of mass, but accompanied with considerable starvation mortality (Merkle & Hanke (1987) Comp. Biochem. Physiol.)
        h_S = 0.6, # hazard rate below critical mass - ca. 50% daily survival probability ()
        a_max = truncated(Normal(15 * 365, 1.5 * 365), 0, Inf), # maximum age [d]
        tau_R = 365., # reproduction period [d]
        
        #=
        Pathogen dynamics and effect parameters
        =#

        Chi = LogNormal(log(1)+1^2, 1), # killing rate modifier, log-normal distribution with mode 1 and sigma 1, (mu = ln(mode)+sigma^2).
            E_P = [Inf, Inf, Inf, Inf], # sensitivity parameter (threshold) for GUTS-SD module
            B_P = [2., 2., 2., 2.], # slope parameter for GUTS-SD module 
    )

    p = ComponentVector(
        glb = glb, 
        spc = spc, 
        pth = AmphiDEB.pth
    )


    sim = AmphiDEB.ODE_simulator(p)

    # check that we have output columns for the second stressor
    @test "C_W_2" in names(sim)
    @test "D_j_2_1" in names(sim)

end

