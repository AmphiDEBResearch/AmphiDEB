abstract type AmphibianEnenergyBudget <: EcotoxSystems.AbstractEnergyBudget end

Base.@kwdef mutable struct AmphibianEnergyBudget1  <: AmphibianEnenergyBudget
    parameters::ComponentVector = parameters_M1
end


generate_individual_params(p) = EcotoxSystems.generate_individual_params(p; pth = p.pth)

"""
Simulate M1.
"""
function simulate(
    model::AmphibianEnergyBudget1;
    additional_callbacks = [],
    kwargs...
    )

<<<<<<< HEAD:src/AmphiDEB_M1/AmphiDEB_M1.jl
    return EcotoxSystems.ODE_simulator(
        model.parameters;
        model = M1_complete_ODE_with_loglogistic_TD!,
        statevars_init = initialize_statevars_M1,
        tstops = [model.parameters.glb.pathogen_inoculation_time, model.parameters.glb.medium_renewals...],
        gen_ind_params = M1_generate_individual_params,
        callbacks = M1_callbacks(additional_callbacks...),
=======
    EcotoxSystems.ODE_simulator(
        p;
        model = model,
        statevars_init = statevars_init,
        tstops = [p.glb.pathogen_inoculation_time, p.glb.medium_renewals...],
        gen_ind_params = gen_ind_params,
        callback = callbacks,
>>>>>>> 92b859b (restructured embryo and larva derivatives):src/simulators.jl
        kwargs...
    )


end
 
#"""
#    ODE_simulator(
#        p::ComponentVector; 
#        model = AmphiDEB_ODE_M1!, 
#        callbacks = AmphODE_callbacks(), 
#        statevars_init = initialize_statevars,
#        kwargs...
#        )
#A wrapper around the `EcotoxSystems.jl` ODE_simulator 
#for application with the AmphiDEB model.
#"""
#function ODE_simulator(
#    p::ComponentVector; 
#    model = M1_complete_ODE_with_loglogistic_TD!, 
#    callbacks = AmphODE_callbacks(), 
#    statevars_init = initialize_statevars,
#    gen_ind_params = p -> EcotoxSystems.generate_individual_params(p; pth = p.pth),
#    kwargs...
#    )
#
#    EcotoxSystems.ODE_simulator(
#        p;
#        model = model,
#        statevars_init = statevars_init,
#        tstops = [p.glb.pathogen_inoculation_time, p.glb.medium_renewals...],
#        gen_ind_params = gen_ind_params,
#        callbacks = callbacks,
#        kwargs...
#    )
#end

"""
    IBM_simulator(
        p::ComponentVector;

        init_global_statevars = initialize_global_statevars,
        global_ode! = AmphiDEB_global!,
        global_rules! = AmphiDEB_global_rules!,
        
        init_individual_statevars = initialize_individual_statevars,
        individual_ode! = M1_individual_ODE_with_loglogistic_TD!, 
        individual_rules! = default_individual_rules!,
        gen_ind_params = p -> EcotoxSystems.generate_individual_params(p; pth = p.pth),

        kwargs...
        )

Compose and execute an individual-based simulation from the given components. <br>
This function is a simple wrapper around `EcotoxSystems.IBM_simulator`. <br>
Its main purpose is to provide a default configuration for the `AmphiDEB` model. 

## Positional arguments

- `p`: Parameters for all components: global (`p.glb`), Amphibian species-specific (`p.spc`) and pathogen-specific (`p.pth`)

## Keyword arguments

**Specification of global part of the system**

- `init_global_statevars`: Function that initializes global state variables as `ComponentVector`, e.g. food abundance, temperature. 
- `global_ode!`: Defintion of global ODE-based portion of the system (e.g. growth or experimental addition of food resource).
- `global_rules!`: Definition of global rule-based portion of the system (e.g. computing population aggregates).

**Specification of species/individual-specific part of the system**

- `init_individual_statevars`: Function that initializes individual-level state variables as `ComponentVector`, e.g. structural mass, scaled adamge, etc.
- `individual_ode!`: Definition of individual-level ODE-based portion of the system (e.g. derivatives of the DEB-TKTD module).
- `individual_rules!`: Definition of individual-level rule-based portion of the system (e.g. death, discrete reproduction events).
- `gen_ind_params`: A function that converts species-specific parameters to individual-specific parameters. For parameters which are provided as distributions, this takes a random sample from the distribution each time an individual is generated. 
This argument probably does not have to be changed, unless additional components (other than `glb`, `spc` and `pth`) are added to the system.
"""
function IBM_simulator(
    # parameters
    p::ComponentVector;

    # default global model
    init_global_statevars = initialize_global_statevars,
    global_ode! = AmphiDEB_global!,
    global_rules! = AmphiDEB_global_rules!,
    
    # default individual model
    init_individual_statevars = initialize_individual_statevars,
    individual_ode! = M1_individual_ODE_with_loglogistic_TD!, 
    individual_rules! = default_individual_rules!,
    gen_ind_params = p -> EcotoxSystems.generate_individual_params(p; pth = p.pth),

    kwargs...
    )

    EcotoxSystems.IBM_simulator(
        p;

        # global model
        init_global_statevars = init_global_statevars,
        global_ode! = global_ode!, 
        global_rules! = global_rules!,
        
        # individual model
        individual_ode! = individual_ode!,
        init_individual_statevars = init_individual_statevars,
        individual_rules! = individual_rules!,
        gen_ind_params = gen_ind_params,

        kwargs...
    )

end