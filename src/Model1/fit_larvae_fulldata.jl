# ModelFitting_01_Discoglossus_larvae.jl
# Estimating AmphiDEB parameters from UCLM D. galganoi

using Base.Threads
using DataFrames, DataFramesMeta
using DataStructures
using CSV
using Chain
using LaTeXStrings
using Suppressor
using Logging
using StatsBase
using Distributions

using StatsPlots, Plots.Measures

default(leg = false)
theme(:default)

# custom packages

using EcotoxSystems, EcotoxModelFitting, AmphiDEB
import EcotoxModelFitting: Hyperdist
import AmphiDEB: ComponentVector
using AmphiDEB.Model1

# constants

const EGG_DRYMASS_MG = 1 # dry mass of an egg, if not fitted
const WETMASS_AT_REPRO_MEASURE_MG = 16530 # wet mass of adult female at time of clutch size measurement (not including egg mass)
const REPRODUCTION_PERIOD = 365 # reproduction period in days

# source files

using Revise
include("losses.jl")
include("traits.jl") # functions to infer traits from ODE solutions
include("utils.jl") # utility functions

#=
## Loading & plotting data
=#

function load_data(
    input_path; # datadir("exp_raw", "Discoglossus_UCLM")
    paths::OrderedDict = OrderedDict(
        :aquatic => [joinpath(input_path, "01_aquatic.csv"), 1], # number indicates row where data header is located (omitting metadata)
        :metamorphs => [joinpath(input_path, "02_metamorphs.csv"), 1],
    ))

    data = OrderedDict()

    for (key,info) in pairs(paths)
        path, header = info
        data[key] = CSV.read(path, DataFrame, header = header)
    end

    dropmissing!(data[:aquatic])

    return data
end

function plot_data(;kwargs...)

    plt_aqua = @df f.data[:aquatic] plot(
        lineplot(
            :t_exp, :wetmass_mg, 
            ylabel = "Wet mass (mg)", 
            leg = :bottomright,
            marker = true, fillalpha = .2, color = :black, label = "Observed mean (P₅-P₉₅)", 
            leftmargin = 5mm, lw = 0,
            title = "Larvae",
            ylim = (100,400),
            xlim = (-1,20)
        ), 
        lineplot(
            :t_exp, :num_tadpoles, 
            xlabel = "Time since start of experiment (d)", ylabel = "Individuals <G42 \n per aquarium", 
            marker = true, fillalpha = .2, color = :black, leg = false,
            leftmargin = 5mm, lw = 0,
            xlim = (-1,20)
        ), 
        layout = (2,1), linecolor = :black, fillcolor = :gray, fillalpha = .5, lw = 1
    )

    @df f.data[:aquatic] scatter!(plt_aqua, :t_exp, :wetmass_mg, color = :gray, label = "")
    
    plt_meta = plot(
        plot(ylabel = "Time since start \n of experiment (d)", title = "Metamorphs"), 
        plot(xlabel = "Gosner stage", ylabel = "Wet mass (mg)"), 
        layout = (2,1), linecolor = :black, fillcolor = :gray, fillalpha = .5, lw = 1
    )

    # violin plots for metam timing
    for (lab, y) in zip(["G42", "G46"], [:t_exp_G42, :t_exp_G46])
        violin!(
            plt_meta, subplot = 1, 
            [lab], f.data[:metamorphs][:,y], 
            xrotation = 45, 
            side = :left, color = :gray, fillalpha = .25, label = y == :t_exp_G42 ? "Observed (KDE)" : ""
            )

    end

    # violin plots for metam growth
    for (lab, y) in zip(["G42", "G46"], [:wetmass_G42_mg :wetmass_G46_mg])
        violin!(
            plt_meta, subplot = 2, 
            [lab], f.data[:metamorphs][:,y], 
            xrotation = 45,
            side = :left, color = :gray, fillalpha = .25, leg = false
            )

    end
    
    @df f.data[:metamorphs] scatter!(plt_meta, markersize = 3, repeat(["G42"], length(:t_exp_G42)), :t_exp_G42, color = :black, subplot = 1, label = "Observed")
    @df f.data[:metamorphs] scatter!(plt_meta, markersize = 3, repeat(["G46"], length(:t_exp_G42)), :t_exp_G46, color = :black, subplot = 1, label = "")
    @df f.data[:metamorphs] scatter!(plt_meta, markersize = 3, repeat(["G42"], length(:t_exp_G42)), :wetmass_G42_mg, color = :black, subplot = 2, label = "")
    @df f.data[:metamorphs] scatter!(plt_meta, markersize = 3, repeat(["G46"], length(:t_exp_G42)), :wetmass_G46_mg, color = :black, subplot = 2, label = "")

    plt = plot(
        plt_aqua, plt_meta, #plt_ad_growth, plt_ad_repro, 
        bottommargin = 10mm, leftmargin = 5mm, topmargin = 5mm, 
        layout = (1,2), size = (800,700); 
        kwargs...
        )

    return plt
end


#=
Plotting simulations
=#

function plot_sims!(plt, predictions::AbstractVector; label = "Prior predictive check")::Nothing

    predictions_valid = filter(x -> !isnothing(x), predictions)
    df_aqua = sort(vcat([df[:aquatic] for df in predictions_valid]...), :t_exp) 

    @df df_aqua lineplot!(
        :t_exp, :wetmass_mg, estimator = mean, color = :steelblue, label = label, leg = :topleft, 
        lw = 2, fillalpha = .2, subplot = 1; 

        )
    
    @df df_aqua lineplot!(
        :t_exp, :num_tadpoles, estimator = mean, color = :steelblue, 
        lw = 2, fillalpha = .2, subplot = 2
    )

    df_metam = vcat([df[:metamorphs] for df in predictions_valid]...) |> clean
    
    for (lab, y) in zip(["G42", "G46"], [:t_exp_G42, :t_exp_G46])
        violin!(
            plt, subplot = 3, 
            [lab], df_metam[:,y], 
            xrotation = 45, 
            side = :right, color = :steelblue, leg = false, 
            fillalpha = .2
            )
    
    end
    
    for (lab, y) in zip(["G42", "G46"], [:wetmass_G42_mg, :wetmass_G46_mg])
        violin!(
            plt, subplot = 4, 
            [lab], df_metam[:,y], 
            xrotation = 45, 
            side = :right, color = :steelblue, leg = false,
            fillalpha = .2
            )
    
    end

    #df_ad = vcat([df[:adults] for df in predictions]...) |> clean
    #
    #@df df_ad density!(plt, subplot = 5, :post_metam_weightchange, fill = true, fillalpha = .25, color = :steelblue; label = label) 
    #@df df_ad density!(plt, subplot = 6, :max_wetmass_mg, fill = true, fillalpha = .25, color = :steelblue; label = "")
    #@df df_ad density!(plt, subplot = 7, :age_at_maturity_assum ./ 365, fill = true, fillalpha = .25, color = :steelblue, label = "")
    #@df df_ad density!(plt, subplot = 8, :mean_clutchsize ./ (1e-3 * WETMASS_AT_REPRO_MEASURE_MG), fill = true, fillalpha = .25, color = :steelblue; label = "")

    return nothing
end

#=
## Setting up the default parameter set
=#

function define_defaultparams()::ComponentVector

    p = EcotoxSystems.ComponentVector(
        glb = params.glb, 
        pth = params.pth,
        spc = ComponentVector(
            params.spc; 
            H_j1_prime = 1., 
            H_p_prime = 1.,
            watercontent_larvae = 0.93, 
            watercontent_juveniles = 0.85,
            time_since_birth = 20., # time since birth at the start of the experiment
            emb_dev_time = 2. 
    ))

    # setting global parameters

    p.glb.t_max = 365. 
    p.glb.pathogen_inoculation_time = Inf
    p.glb.dX_in = [1e10, 1e10] # ad libitum feeding conditions

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    # propagation of zoom factor to H_j1 is turned off => we want variability in the transition to metamorphs
    p.spc.propagate_zoom.H_j1 = 0.

    p.spc.X_emb_int = 1. # ≈ initial dry mass of an egg (mg)

    return p
end

"""
    simulate_preexperiment(p::ComponentVector)::ComponentVector

Simulate individual life history before the start of the experiment. <br>
This is used to simulate different conditions before and during the experiment. 
"""
function simulate_preexperiment(p::ComponentVector)::ComponentVector

    sim_emb, u0lrv, p_ind = sim_embryo(p) # simulate embryo 
    p_ind.glb.t_max = p.ind.time_since_birth # set max. simulation time to time since birth
    _, u_t0, p_ind = sim_larva(p_ind, u0lrv, p_ind) # simulate larvae from birth to `time_since_birth`

    return u_t0, p_ind
end


#=
## Simulating the experiment
=#

"""
    simulator(
        p::EcotoxSystems.ComponentVector; 
        return_raw::Bool = false, 
        param_links::NamedTuple = (ind = link_ind_params!,),
        kwargs... 
        )

Simulate life history of *D. galganoi* for calibration with UCLM dataset.
"""
function simulator(
    p::EcotoxSystems.ComponentVector; # parameters and forcings
    return_raw::Bool = false, # return raw simulation output? if false, converts output to format of the data
    kwargs... # additional arguments for ODE_simulator
    )

    let Z = deepcopy(p.spc.Z)

        p.spc.time_since_birth = ceil(p.spc.time_since_birth)

        # add parameter links
        
        p.spc.dI_max_emb = p.spc.dI_max_lrv
        p.spc.k_M_juv = p.spc.k_M_emb
        
        # estimate population mean of the embryonic development time

        p.spc.Z = Dirac(1.) # turn off individual variabiilty to get estimate based on popmean 
        p.spc.emb_dev_time = estimate_emb_dev_time(p) # assign embryonic development time
        p.spc.Z = Z # re-assign original zoom factor

        # solve the model for a number of replicates
        # in the context of simulation, "replicates" are the number of tadpoles per aquarium (rather than the number of aquaria)
        
        u_t0, p_ind = simulate_preexperiment(p)

        try
            sim = @replicates sim_all(p) 10

            sim[!,:drymass_mg] = sim.S .+ sim.E_mt
            sim[!,:wetmass_mg] = [calc_wetmass(r, p.spc.watercontent_larvae, p.spc.watercontent_juveniles) for r in eachrow(sim)]
            sim[!,:t_exp] = sim.t #.- age_at_birth(sim) .- p.spc.time_since_birth # conversion is not needed if simulate_preexperiment is used in call to AmphiDEB.ODE_simulator

            # optionally, return the simulation output
            if return_raw
                return sim 
            end

            # convert raw simulation output to dataset
            sim_data = OrderedDict(
                :aquatic => extract_aquatic_data(sim),
                :metamorphs => extract_metamorph_data(sim, 0.), # second argument set to 0 because we have already previously corrected for time_since_birth and emb_dev_time
            )

            #sim_data[:metamorphs][!,:H_j1] .= p.spc.H_j1
            #sim_data[:metamorphs][!,:H_eq_lrv] .= H_eq_lrv

            return sim_data
        catch e
            # Log a task-specific message
            with_logger(logger) do
                @info("
                Encountered error in simulator: 
                $e 
                Error ocurred for the following parameter sample: 
                $p
                ")
            end

            # write buffered message to logging file
            flush(io)
            
            return nothing
        end
    end

end

"""
    extract_aquatic_data(sim::AbstractDataFrame)::DataFrame

Convert raw simulation output to data key "aquatic". <br>
Selects larval life stages, computes average wet mass and number of tadpoles (∝ 1 - time-resolved probability to reach metamorphosis). 
"""
function extract_aquatic_data(sim::AbstractDataFrame)::DataFrame
 
    aquatic = @chain sim begin
        @select(:t_exp, :larva, :wetmass_mg, :replicate)
        @subset(:larva .== 1.) # only larvae 
        groupby(:t_exp) # for every time point
        combine(
            [:wetmass_mg,:replicate] => 
            ((m, r) -> 
                (wetmass_mg = mean(m), # calculate average dry mass
                num_tadpoles = length(unique(r))) # count how many tadpoles we have
                ) => AsTable
            )
    end

    return aquatic
end

"""
    extract_metamorph_data(sim::AbstractDataFrame, time_since_birth)::DataFrame

Convert raw simulation output to metamorphs data key "metamorphs". <br>
Computes metamorphosis traits recorded in fixed-stage designs.
"""
function extract_metamorph_data(sim::AbstractDataFrame, time_since_birth)::DataFrame

    metamorphs = combine(groupby(sim, :replicate)) do df
        DataFrame(
            t_exp_G42 = metamorphosis_timing(df)[1] .- age_at_birth(df) .- time_since_birth,
            t_exp_G46 = metamorphosis_timing(df)[2] .- age_at_birth(df) .- time_since_birth,
            dt_mt = metamorphosis_duration(df),
            wetmass_G42_mg = wetmass_at_G42(df),
            wetmass_G46_mg = wetmass_at_G46(df)
        )
    end |> df -> DataFrame(
        t_exp_G42 = mean(df.t_exp_G42),
        t_exp_G46 = mean(df.t_exp_G46),
        dt_mt = mean(df.dt_mt),
        wetmass_G42_mg = mean(df.wetmass_G42_mg),
        wetmass_G46_mg = mean(df.wetmass_G46_mg)
    )
    
    return metamorphs
end

function extract_adult_data(sim::AbstractDataFrame)::DataFrame
    return @chain sim begin
        groupby(:replicate)
        combine(_) do df
            DataFrame(
                age_at_maturity_assum = EcotoxSystems.robustmin(@subset(df, :adult .> 0.5).t),
                clutchsize = clutchsize_at_specific_wetmass(df), 
                max_wetmass_mg = maximum(df.wetmass_mg),
                post_metam_weightchange = post_metamorphic_weightchange(df),
            )
        end
        DataFrame(
            mean_clutchsize = EcotoxSystems.robustmean(_.clutchsize),
            age_at_maturity_assum = EcotoxSystems.robustmean(_.age_at_maturity_assum),
            max_wetmass_mg = EcotoxSystems.robustmean(_.max_wetmass_mg),
            post_metam_weightchange = EcotoxSystems.robustmean(_.post_metam_weightchange)
        )
    end
end

function link_ind_params!(ind::ComponentVector)::Nothing
    ind.dI_max_emb = ind.dI_max_lrv # ingestion rate for embryos assumed to be same as for larave
    return nothing
end


"""
    setup_modelfit(;
        σ_factor = 1., 
        loss_function = loss_symmbound
    )

Instantiate a `ModelFit` instance, including definition of priors, data weights 
and response variable configuration. 

Note that the output of this function may be modified in the notebook, in order to generate multiple alternative fits.

kwargs

- `σ_factor`: multiplies the prior SD of `dI_max_lrv` and `k_M_emb` with a common factor
- `loss_function`: error model applied to all response variables. `loss_symmboud`, `loss_mse` and `loss_mse_logtransform` are currently built-in choices.
"""
function setup_modelfit(;
    σ_factor = .5, 
    loss_function = loss_mse_logtransform
    )
    
    if !isdir(datadir("sims", savetag))
        mkdir(datadir("sims", savetag))
    end

    # set up logger
    # FIXME: this does not save the log file to the intended directory, because the setup_modelfit does not know the composite "savetag" passed on to run_PMC1
    global io = open(datadir("sims", savetag, "log.txt"), "w+")
    global logger = SimpleLogger(io)

    data = load_data()
    defparams = define_defaultparams()

    ### ---- settting up ModelFit object ---- ####
    
    f = ModelFit( 
        prior = Prior(
            "spc.Z" => Hyperdist(
                σ -> truncated(Normal(1, σ), 0, Inf), # distribution of Z as a function of sigma
                truncated(Normal(0.1, 0.1), 0, Inf) # prior distribution of sigma
                ),  
            "spc.dI_max_lrv" => truncated(Normal(1.6, 0.3*σ_factor), 0, Inf), #truncated(Normal(1.6, 0.3*σ_factor), 0, Inf), 
            "spc.k_M_emb" => truncated(Normal(0.2, 0.1*σ_factor), 0, Inf),     
            "spc.eta_AS_emb" => truncated(Normal(0.75, 0.75), 0, 1),
            "spc.H_j1" => truncated(Normal(9, 4.5), 0, Inf),
            # H_j1_prime is relevant for an alternative parameterization of the model; this is currently turned off
            #"spc.H_j1_prime" => truncated(Normal(0.2, 0.02), 0, 1),
            "spc.gamma" => truncated(Normal(0.5, 0.5), 0, 1), 
            "spc.k_J_emb" => truncated(Normal(0.05, 0.05), 0, Inf),
            # prior for kappa_emb has a lower limit >0  because extremely low values result in instability of the ODE
            "spc.kappa_emb" => truncated(Normal(0.8, 0.4), 0.2, 1),
            "spc.watercontent_larvae" => truncated(Normal(0.94, 0.004), 0, 1),
            "spc.watercontent_juveniles" => truncated(Normal(0.84, 0.042), 0, 0.93), 
            "spc.time_since_birth" => truncated(Normal(20,1),15,25)
            ),
        defaultparams = defparams, 
        simulator = simulator, 
        data = data, 
        response_vars = [
            [:wetmass_mg, :num_tadpoles], 
            [:wetmass_G42_mg, :wetmass_G46_mg, :dt_mt]
        ],
        data_weights = [
            [6., 6.],
            [1., 1., 1.]
        ],
        time_resolved = [true, false],
        time_var = :t_exp,
        plot_data = plot_data,
        loss_functions = loss_function
    )

    return f
end

# misc

paramlabels = OrderedDict(
    "spc.eta_AS_emb" => L"\eta_{AS}^{emb}",
    "spc.eta_AS_juv" => L"\eta_{AS}^{juv}",
    "spc.eta_AR" => L"\eta_{AR}",
    "spc.dI_max_lrv" => L"\{ \dot{I} \}_{max}^{lrv}",
    "spc.dI_max_juv" => L"\{ \dot{I} \}_{max}^{juv}",
    "spc.k_M_emb" => L"k_M^{emb}",
    "spc.k_J_emb" => L"k_J^{emb}",
    "spc.Z" => L"\sigma_{Z,M}",
    "spc.H_j1" => L"H^{j}",
    "spc.H_j1_prime" => L"H^{j'}",
    "spc.H_p_prime" => L"H^{p'}",
    "spc.gamma" => L"\gamma",
    "spc.kappa_emb" => L"\kappa",
    "spc.watercontent_larvae" => L"Larval\ water\ content",
    "spc.watercontent_juveniles" => L"Juvenile\ water\ content",
    "spc.time_since_birth" => L"Time\ since\ birth",
    "spc.delta_E" => L"\delta_E", 
    "spc.delta_k_M_mt" => L"\delta_{k_M^{mt}}"
)

"""
    fit_model!(
        f; 
        n = n, 
        n_init = n_init, 
        q_dist = q_dist, 
        t_max = t_max, 
        savetag = nothing, 
        paramlabels = paramlabels, 
        continue_from = nothing,
        n_posterior_check = 100,
        kwargs...
        )

Model fit for *D. galganoi* larave/metamorphs to UCLM data.
"""
function fit_model!(
    f; 
    n = n, 
    n_init = n_init, 
    q_dist = q_dist, 
    t_max = t_max, 
    savetag = nothing, 
    paramlabels = paramlabels, 
    continue_from = nothing,
    kwargs...
    )

    @info "Using savetag $(savetag)"

    let pmchist, posterior_check
            
        # run the calibration
        pmchist = run_PMC!(
            f; 
            n = n, 
            n_init = n_init, 
            q_dist = q_dist, 
            t_max = t_max, 
            evals_per_sample = 3,
            savedir = datadir("sims"),
            savetag = savetag, 
            paramlabels = paramlabels, 
            continue_from = continue_from,
            logweights = true
        )

        closeall() # reset plots

        posterior_check = run_diagnostics(f, pmchist, savetag)

        return pmchist, posterior_check
    end
end

function update_data_weights!(f::ModelFit, w)

    f.data_weights = w |> x -> x ./ sum(vcat(x...))
    f.loss = generate_loss_function(f)

end

function VPC_bestfit(f::ModelFit, savetag::AbstractString)

    @info "#### ---- Best fit ---- ####"

    let plt = f.plot_data()

        p_opt = f.accepted[:,argmin(vec(f.losses))]
        sim_opt = [f.simulator(p_opt) for _ in 1:100]

        plot_sims!(plt, sim_opt, label = "Best fit")
        display(plt)
        
        if !isnothing(savetag)
            savefig(plot(plt, dpi = 400), datadir("sims", "$(savetag)", "VPC_posterior_bestfit.png"))
        end

    end
end

function VPC_posterior(f::ModelFit; n_posterior_check = 100)
    
    @info "#### ---- Posterior retrodictions ---- ####"

    let posterior_check
        @suppress posterior_check = posterior_predictions(f, n_posterior_check) 

        let plt = f.plot_data()
            

            plot_sims!(plt, posterior_check.predictions, label = "Retrodictions")

            if !isnothing(savetag)
                savefig(plot(plt, dpi = 400), datadir("sims", "$(savetag)", "VPC_posterior_samples.png"))
            end

            display(plt)
        end

        return posterior_check
    end
end

function run_diagnostics(
    f::ModelFit, 
    pmchist, 
    savetag::Union{AbstractString,Nothing}
    )


    posterior_check = VPC_posterior(f)
    VPC_bestfit(f, savetag)

    @info "#### ---- Marginal posteriors ---- ####"

    let plt, num_params = length(f.prior.dists), num_cols = 4, num_rows = Int(ceil(num_params / num_cols))

        plt = plot(
            plot.(f.prior.dists, color = :black)..., layout = (num_rows, num_cols), 
            leg = hcat(vcat(true, repeat([false], num_params-1))...), 
            label = "Prior", 
            size = (1200,200*num_rows), bottommargin = 5mm
            )
    
            
        for (i,param) in enumerate(f.prior.labels)
    
            histogram!(
                plt, subplot = i, 
                f.accepted[i,:], weights = Weights(vec(f.weights)), 
                xlabel = param in keys(paramlabels) ? paramlabels[param] : param, 
                normalize = :pdf, label = "Posterior", color = :gray, lw = 0.5, fillalpha = .5
                )
        end
    
        display(plt)
        
        if !isnothing(savetag)
            savefig(plot(plt, dpi = 400), datadir("sims", "$(savetag)", "marginal_posteriors.png"))
        end
    end

    @info "#### ---- Convergence behaviour ---- ####"

    let plt
        plt = plot(eachindex(pmchist.dists) .- 1, map(median, pmchist.dists), marker = true, lw = 1.2, xlabel = "PMC step", ylabel = "Loss", label = "Median", xticks = eachindex(pmchist.dists) .- 1)
        plot!(plt, eachindex(pmchist.dists) .- 1, map(minimum, pmchist.dists), marker = true, lw = 1.2, label = "Minimum")
        display(plt)
        savefig(plot(plt, dpi = 400), datadir("sims", savetag, "loss.png"))
    end

    begin
        posterior_medians = [mapslices(x -> median(x, Weights(pmchist.weights[i])), pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
        x -> hcat(x...)

        q25 = [mapslices(x -> quantile(x, Weights(pmchist.weights[i]), 0.25), pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
        x -> hcat(x...)
        
        q75 = [mapslices(x -> quantile(x, Weights(pmchist.weights[i]), 0.75), pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
        x -> hcat(x...)
        
        num_params = length(f.prior.dists) 
        num_cols = 4 
        num_rows = Int(ceil(num_params / num_cols))

        plt = plot(
            eachindex(pmchist.particles) .-1,
            posterior_medians', layout = size(posterior_medians)[1], 
            size = (1200,200*num_rows), marker = true,
            leg = hcat(vcat(:topleft, repeat([false], length(f.prior.dists)-1))...),
            label = "Median",
            ylabel = hcat(f.prior.labels...), titlefontsize = 12,
            bottommargin = 5mm, leftmargin = 5mm, 
            xlabel = "PMC step"
            )

        plot!(
            eachindex(pmchist.particles) .-1, 
            q25', fillrange = q75', 
            lw = 0, fillcolor = :gray, fillalpha = .25,
            markersize = 0, label = "IQR"
            )


        for (i,dist) in enumerate(f.prior.dists)

            # set ylim based in prior limits

            q1 = quantile(dist, 0.01)
            q2 = quantile(dist, 0.99)

            # indicate IQR of priors 

            l = repeat([quantile(dist, 0.25)], length(pmchist.particles))
            u = repeat([quantile(dist, 0.75)], length(pmchist.particles))

            plot!(plt, subplot = i, ylim = (q1,q2))

            hline!(plt, subplot = i, [quantile(dist, 0.25)], color = :gray, linestyle = :dash, label = "")
            hline!(plt, subplot = i, [quantile(dist, 0.75)], color = :gray, linestyle = :dash, label = "")

        end
        display(plt)
    end

    @info "#### ---- Posterior summary ---- ####"

    generate_posterior_summary(
        f; 
        tex = false,
        paramlabels = paramlabels,
        savetag = nothing
    ) |> display

    return posterior_check
end