using Colors
using CSV 
using DataFrames
using Dates
using Distributions
using HypothesisTests
using MCMCChains
using BridgeStan
using Plots
using Plots.PlotMeasures: px
using Random
using SplittableRandoms: SplittableRandom
using Statistics
using StanSample
using StatsPlots

using Pigeons
try
    using autoHMC
catch
    using Pkg
    Pkg.instantiate()
    Pkg.precompile()
end

###############################################################################
# ESS and friends
###############################################################################

include("batch_means.jl")

const PigeonsSample = Union{Pigeons.SampleArray{<:AbstractVector}, Vector{<:AbstractVector}}

get_component_samples(samples::PigeonsSample, idx_component) =
    collect(hcat((s[idx_component] for s in samples)...)')
get_component_samples(samples::PigeonsSample, idx_component::Int) =
    [s[idx_component] for s in samples]
get_component_samples(samples::DataFrame, idx_component) = Array(samples[:,idx_component])

# ESS
function margin_ess(samples, model_exact_means=nothing, margin_idx=1, args...)
    margin = get_component_samples(samples,margin_idx)
    isnothing(model_exact_means) && return batch_means_ess(margin)
    batch_means_ess(margin, special_margin_mean_std(model_exact_means, args...)...) 
end
special_margin_idx(model::AbstractString, dim::Int) = 
    startswith(model, "two_component_normal") ? dim : one(dim)

special_margin_mean_std(model::String, args...) = 
    if startswith(model,"funnel")
        (0.,3.)
    elseif startswith(model,"banana")
        (0., sqrt(10))
    elseif startswith(model,"eight_schools")    # approx using variational PT (see 8schools_pilot.jl)
        (3.574118538746056, 3.1726880307401455)
    elseif startswith(model, "normal") 
        (0., 1.)
    elseif startswith(model, "two_component_normal")
        (0., last(two_component_normal_stdevs(args...))) # the margin with the largest stdev 
    else
        error("unknown model $model")
    end
two_component_normal_stdevs(e::Real,args...) = (10. ^(-e), 10. ^(e))
n_vars(samples::PigeonsSample) = length(first(samples))
n_vars(samples::DataFrame) = size(samples,2)
min_ess_batch(samples) = minimum(1:n_vars(samples)) do i
        batch_means_ess(get_component_samples(samples, i))
    end
to_chains(samples::PigeonsSample) = Chains(samples)
to_chains(samples::DataFrame) = Chains(Array(samples))
min_ess_chains(samples) = minimum(ess(to_chains(samples),kind=:basic).nt.ess) # :basic is actual ess of the vars. default is :bulk which computes ess on rank-normalized vars

# TODO: should we discard the warmup samples?
function MSJD(sample::PigeonsSample) 
    n = length(sample) 
    msjd = 0.0
    for i in 2:n 
       msjd += sum((sample[i] .- sample[i-1]) .^ 2) / (n-1)
    end 
    return msjd
end 

function MSJD(sample::DataFrame) 
    sample_vec = df_to_vec(sample) 
    return MSJD(sample_vec) 
end 

function Statistics.mean(sample::PigeonsSample, margin)
    x = [sample[i][margin] for i in eachindex(sample)]
    return mean(x) 
end

function Statistics.mean(sample::DataFrame, margin)
    sample_vec = df_to_vec(sample)
    return mean(sample_vec, margin)
end

function Statistics.var(sample::PigeonsSample, margin)
    x = [sample[i][margin] for i in eachindex(sample)]
    return var(x) 
end

function Statistics.var(sample::DataFrame, margin)
    sample_vec = df_to_vec(sample)
    return var(sample_vec, margin)
end

#=
Kolmogorov-Smirnov test.  
=#
function KS_statistic(sample::PigeonsSample, d::Distribution, margin = 1)
    x = [sample[i][margin] for i in 1:length(sample)]
    t = HypothesisTests.ApproximateOneSampleKSTest(x, d)
    return sqrt(t.n)*t.δ, pvalue(t) 
end 

function KS_statistic(sample::DataFrame, d::Distribution, margin = 1) 
    sample_vec = df_to_vec(sample) 
    return KS_statistic(sample_vec, d, margin) 
end

function df_to_vec(df::DataFrame) 
    n = size(df,1)
    df_vec = Vector{Vector{Float64}}(undef, n) 
    for i in 1:n 
        df_vec[i] = Vector(df[i, :]) 
    end 
    return df_vec
end 

###############################################################################
# sampling
###############################################################################

function pt_sample_from_model(
    target, seed, explorer, explorer_type, args...;
    n_chains = 1, kwargs...
    )
    inp = Inputs(
        target      = target, 
        seed        = seed,
        n_rounds    = 1,
        n_chains    = n_chains, 
        explorer    = explorer, 
        show_report = true
    )
    pt_sample_from_model(inp, args...; kwargs...)
end

function pt_sample_from_model(inp::Inputs, n_rounds; keep_traces=true)
    # build pt
    recorders = [record_default(); Pigeons.explorer_acceptance_pr] 
    keep_traces && (recorders = vcat(recorders, Pigeons.traces))
    inp.record = recorders
    pt = PT(inp)

    # iterate rounds
    n_leapfrog = 0
    for _ in 1:n_rounds
        pt = pigeons(pt)
        n_leapfrog += Pigeons.explorer_n_steps(pt)[1]
        pt = Pigeons.increment_n_rounds!(pt, 1)
    end
    time   = sum(pt.shared.reports.summary.last_round_max_time)
    if keep_traces 
        sample = get_sample(pt) # compute summaries of last round (0.5 of all samples)
        @assert length(sample) == 2^n_rounds
    else 
        sample = nothing 
    end
    return pt, time, sample, n_leapfrog
end

function nuts_sample_from_model(model, seed, n_rounds; kwargs...)
    stan_model = model_string(model; kwargs...)
    sm = SampleModel(model, stan_model) 
    sm.num_threads      = 1
    sm.num_julia_chains = 1
    sm.num_chains       = 1
    sm.num_samples      = 2^n_rounds 
    sm.num_warmups      = 2^n_rounds - 2
    sm.save_warmup      = true
    sm.seed             = seed
    sm.show_logging     = true

    data = stan_data(model; kwargs...)
    time = @elapsed begin 
        rc = stan_sample(sm; data)
    end
    sample = DataFrame(read_samples(sm))[(sm.num_warmups+1):end,:] # discard warmup
    @assert nrow(sample) == 2^n_rounds == sm.num_samples "nrow(sample) = $(nrow(sample)), sm.num_samples=$(sm.num_samples)"
    info = DataFrame(CSV.File(joinpath(sm.tmpdir, model * "_chain_1.csv"), comment = "#"))
    @assert size(info,1) == sm.num_samples + sm.num_warmups
    n_leapfrog = sum(info.n_leapfrog__) # count leapfrogs during warmup
    return time, sample, n_leapfrog
end

get_step_size(explorer) = explorer.step_size
get_step_size(explorer::Mix) = first(explorer.explorers).step_size

###############################################################################
# loading data
###############################################################################

function base_dir()
    base_folder = dirname(dirname(Base.active_project()))
    endswith(base_folder, "autoHMC-mev") || error("please activate the autoHMC-mev julia-environment")
    return base_folder
end

function get_summary_df(experiment::String)
    base_folder = base_dir()
    csv_path    = joinpath(base_folder, "deliverables", experiment, "aggregated", "summary.csv")
    return DataFrame(CSV.File(csv_path))
end

###############################################################################
# sampling utilities
###############################################################################

function model_string(model; dataset=nothing, kwargs...)
    if model == "normal" # dont have the standard normal on Pigeons examples
        return "data {
          int<lower=1> dim;
        }
        parameters {
          vector[dim] x;
        }
        model {
          x ~ std_normal();
        }"
    end
    if startswith(model, "two_component_normal")
        return read(joinpath(
            base_dir(), "stan", "two_component_normal.stan"), String)
    end
    if startswith(model, "horseshoe")
        is_logit = any(Base.Fix1(startswith,dataset), ("prostate", "ionosphere", "sonar"))
        return read(joinpath(
            base_dir(), "stan", "horseshoe_" * (is_logit ? "logit" : "linear") * ".stan"
        ), String)
    end
    if model == "mRNA"
        return read(joinpath(base_dir(), "stan", "mRNA.stan"), String)
    end
    pigeons_stan_dir = joinpath(dirname(dirname(pathof(Pigeons))),"examples","stan")
    if startswith(model, "eight_schools_") 
        return read(joinpath(pigeons_stan_dir,"$model.stan"), String)
    end
    model_class = first(split(model,"_"))
    if model_class in ("banana","funnel") 
        return read(joinpath(pigeons_stan_dir,"$model_class.stan"), String)
    end
    error("model_string: model $model unknown")
end

function stan_data(model::String; dataset=nothing, dim=nothing, scale=nothing) 
    if model in ("funnel", "banana")
        Dict("dim" => dim-1, "scale" => scale)
    elseif model in ("funnel_scale", "banana_scale") 
        Dict("dim" => 1, "scale" => inv(dim))
    elseif model == "normal"
        Dict("dim" => dim) 
    elseif model == "two_component_normal_scale"
        s_lo, s_hi = two_component_normal_stdevs(dim)
        Dict("n" => 1, "s_lo" => s_lo, "s_hi" => s_hi)
    elseif model == "horseshoe"
        x,y = isnothing(dim) ? make_HSP_data(dataset) : make_HSP_data(dataset,dim) # interpret dim as n_obs
        Dict("n" => length(y), "d" => size(x,2), "x" => x, "y" => y)
    elseif startswith(model,"eight_schools")
        Dict("J" => 8, "y" => [28, 8, -3, 7, -1, 1, 18, 12],
        "sigma" => [15, 10, 16, 11, 9, 11, 10, 18])
    elseif model == "mRNA"
        dta = DataFrame(CSV.File(joinpath(base_dir(), "data", "transfection.csv")))
        Dict("N" => nrow(dta), "ts" => dta[:,1], "ys" => dta[:,3])
    else
        error("stan_data: unknown model $model") 
    end 
end 

# Two component normal for testing preconditioner
function make_2_comp_norm_target(n, exponent)
    s_lo, s_hi = two_component_normal_stdevs(exponent)
    json_str = Pigeons.json(; n=n, s_lo=s_lo, s_hi=s_hi)
    StanLogPotential(joinpath(
        base_dir(), "stan", "two_component_normal.stan"
    ), json_str)
end

# build the horseshoe prior target with varying number of observations
load_HSP_df(dataset::String) = 
    DataFrame(CSV.File(
        joinpath(base_dir(), "data", dataset * ".csv") ))
make_HSP_data(dataset::String, n_obs::Int=typemax(Int)) = 
    make_HSP_data(dataset,load_HSP_df(dataset),n_obs)
function make_HSP_data(dataset::String, df::DataFrame, n_obs::Int)
    iszero(n_obs) && return (zeros( ( n_obs,size(df,2)-1 ) ), Int64[])
    n = min(n_obs, size(df, 1))
    if startswith(dataset,"prostate")
        x = Matrix(df[1:n,2:end])
        y = df[1:n,1]
    elseif startswith(dataset,"ionosphere")
        x = Matrix(hcat(df[1:n,1], df[1:n,3:(end-1)])) # col 2 is constant
        y = Int.(df[1:n,end] .== "g")
    elseif startswith(dataset,"sonar")
        x = Matrix(df[1:n,1:(end-1)])
        y = Int.(df[1:n,end] .== "Mine")
    end
    x,y
end
function make_HSP_target(dataset::String, n_obs::Int=typemax(Int))
    xmat,y = make_HSP_data(dataset, n_obs)
    d = size(xmat,2)
    json_str = if iszero(n_obs)
        Pigeons.json(; n=n_obs,d=d,x="[[]]",y="[]")
    else
        x = [copy(r) for r in eachrow(xmat)]
        Pigeons.json(; n=length(y), d=d, x=x, y=y)
    end
    is_logit = any(Base.Fix1(startswith,dataset), ("prostate", "ionosphere", "sonar"))
    StanLogPotential(joinpath(
        base_dir(), "stan", "horseshoe_" * (is_logit ? "logit" : "linear") * ".stan"
    ), json_str)
end

###############################################################################
# plotting utilities
###############################################################################

# Boxplots
const JULIA_AUTO = theme_palette(:auto).colors.colors
const COLORBLIND_IBM = [colorant"#785EF0", colorant"#DC267F", colorant"#FE6100", colorant"#FFB000", colorant"#648FFF"] # IBM palette from https://davidmathlogic.com/colorblind/
const COLORBLIND_WONG = [
    colorant"#000000", colorant"#E69F00", colorant"#56B4E9", colorant"#009E73", 
    colorant"#F0E442", colorant"#0072B2", colorant"#D55E00", colorant"#CC79A7"] # Wong palette from https://davidmathlogic.com/colorblind/
get_palette(n_levels::Int) = 
    ((n_levels <= 2 || n_levels > 8) ? JULIA_AUTO : n_levels <= 5 ? COLORBLIND_IBM : COLORBLIND_WONG)[1:n_levels]
dataset_nickname(d::AbstractString) = d=="sonar" ? "Sonar" : (d=="prostate_small" ? "Prostate" : "Ion.")
function make_boxplots(df::DataFrame; model = first(df.model), fn_end = ".pdf")
    n_samplers = length(unique(df.sampler))
    only_two_samplers = n_samplers == 2
    colors = get_palette(n_samplers)

    # preprocessing
    sort!(df)
    is_funnel = occursin("funnel",model)
    is_highdim = occursin("highdim",model)
    is_banana = !is_highdim && occursin("banana",model) # check if its banana_scale
    is_log2_x = is_banana || is_highdim
    is_log2_x && (df.dim .= log2.(df.dim)) # boxplot not working with xaxis=:log 
    is_hsp = occursin("horseshoe",model) && hasproperty(df, :dataset)
    is_2_comp_norm = occursin("two_component_normal",model)
    xvar = :dim
    if is_hsp
        xvar = :data_nobs
        df[!,xvar] = map(
            t -> (dataset_nickname(t[1]) * "\n" * (t[2] > 100 ? "full" : string(t[2]))),
            zip(df.dataset,df.dim)
        )
        sort!(df, xvar)
    end
    df[!,:min_ess] = min.(df.ess, df.ess_exact)
    df[!,:nleap_to_min_ess] = df.n_leapfrog ./ df.min_ess
    df_means = combine(
        groupby(df, [xvar, :sampler]),
        :margin1_mean => mean,
        :margin1_var => mean,
        :n_leapfrog => mean,
        :min_ess => mean,
        :nleap_to_min_ess => mean, 
        renamecols=false)
    sort!(df_means)

    # common properties
    size   = (650,300)
    xlab   = is_hsp ? "Dataset" : (is_highdim ? "Dimension" : "Inverse scale") * (is_log2_x ? " (log₂)" : "")
    mar    = 15px
    path   = joinpath(base_dir(), "deliverables", model)
    n_samplers = length(unique(df.sampler))

    # plots for known margins
    if !is_hsp
        # margin1 mean
        margin_idx = is_2_comp_norm ? 2 : 1
        p=@df df groupedboxplot(
            :dim, 
            :margin1_mean, 
            group=:sampler,
            bar_position = :dodge,
            size=size,
            xlab=xlab,
            palette=colors,
            ylab="Margin $margin_idx mean",
            left_margin = mar, bottom_margin = mar,
        )
        if !only_two_samplers
            @df df_means plot!(p,
                :dim,
                :margin1_mean,
                group=:sampler,
                palette=colors,
                linewidth=2,
                label=""
            )
        end
        savefig(p,joinpath(path, "boxplots-margin-mean" * fn_end))

        # margin1 var
        p=@df df groupedboxplot(
            :dim, 
            :margin1_var, 
            group=:sampler,
            bar_position = :dodge, 
            size=size,
            xlab=xlab,
            legend = is_2_comp_norm ? :topleft : :best,
            yaxis= is_2_comp_norm ? :log : :identity,
            palette=colors,
            ylab="Margin $margin_idx var",
            left_margin = mar, bottom_margin = mar,
        )
        if !only_two_samplers
            @df df_means plot!(p,
                :dim,
                :margin1_var,
                group=:sampler,
                yaxis= is_2_comp_norm ? :log : :identity,
                palette=colors,
                linewidth=2,
                label=""
            )
        end
        if is_2_comp_norm
            dim_vals = sort(unique(df.dim))
            plot!(dim_vals, 10. .^ (2*dim_vals), label = "true",
            linestyle=:dash, color=colorant"#648FFF")
        end
        savefig(p,joinpath(path, "boxplots-margin-var" * fn_end))
    end

    # n_leapfrog
    p=@df df groupedboxplot(
        (is_hsp ? :data_nobs : :dim), 
        :n_leapfrog, 
        group=:sampler,
        bar_position = :dodge, 
        legend = is_hsp || is_2_comp_norm ? :outerright : :best,
        yaxis= :log,
        palette=colors,
        size=size,
        xlab=xlab,
        ylab="Total number of leapfrog steps",
        left_margin = mar, bottom_margin = mar,
    )
    if !only_two_samplers && is_hsp
        plot!(p,
            repeat(first(first(xticks(p))),inner=n_samplers),
            df_means.n_leapfrog,
            group=df_means.sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    elseif !only_two_samplers
        @df df_means plot!(p,
            :dim,
            :n_leapfrog,
            group=:sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    end
    savefig(p,joinpath(path, "boxplots-n_leapfrog" * fn_end))

    # min ess
    p=@df df groupedboxplot(
        (is_hsp ? :data_nobs : :dim), 
        :min_ess,
        group=:sampler,
        bar_position = :dodge, 
        legend = is_2_comp_norm ? :outerright : :best,
        yaxis= :log,
        palette=colors,
        size=size,
        xlab=xlab,
        ylab="minESS",
        left_margin = mar, bottom_margin = mar
    )
    if !only_two_samplers && is_hsp
        plot!(p,
            repeat(first(first(xticks(p))),inner=n_samplers),
            df_means.min_ess,
            group=df_means.sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    elseif !only_two_samplers
        @df df_means plot!(p,
            :dim,
            :min_ess,
            group=:sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    end
    savefig(p,joinpath(path, "boxplots-miness" * fn_end))

    # nleap to miness
    p=@df df groupedboxplot(
        (is_hsp ? :data_nobs : :dim), 
        :nleap_to_min_ess, 
        group=:sampler,
        bar_position = :dodge, 
        legend = (only_two_samplers || is_funnel) ? :bottomright : (is_2_comp_norm ? :topleft : (is_hsp ? :outerright : :best)),
        yaxis= :log,
        palette=colors,
        size=size,
        xlab=xlab,
        ylab="Leapfrogs per minESS",
        left_margin = mar, bottom_margin = mar,
    )
    if !only_two_samplers && is_hsp
        plot!(p,
            repeat(first(first(xticks(p))),inner=n_samplers),
            df_means.nleap_to_min_ess,
            group=df_means.sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    elseif !only_two_samplers
        @df df_means plot!(p,
            :dim,
            :nleap_to_min_ess,
            group=:sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    end
    savefig(p,joinpath(path, "boxplots-nleap_to_min_ess" * fn_end))
end


