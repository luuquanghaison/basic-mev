include { crossProduct; collectCSVs; setupPigeons; head; pow; deliverables; checkGitUpdated; commit } from './utils.nf'
params.dryRun = false

def variables = [
    dim: (1..20).collect{0.2*it},
    seed: (1..20),
    model: ["funnel_scale"],
    nleaps: [1],// currently not used // (10..30).step(10), // collect{1<<it}, // == 2^it but using bit shift
    sampler: ["AM","AH_fix", "AH_unif", "AH_exp","NUTS"]
]

model_string = [
    funnel_scale: "Pigeons.stan_funnel(1, 1/dim)", 
]

sampler_string = [ 
    AM: "AutoMALA(base_n_refresh=1)",
    AH_fix: "SimpleAHMC()", // base_n_refresh=1 by default on pkg autoHMC. also jitter = Dirac(1.0)
    AH_unif: "SimpleAHMC(jitter_n_leaps=Uniform(0.,2.))",
    AH_exp: "SimpleAHMC(jitter_n_leaps=Exponential(1.))",
    NUTS: "Pigeons.MALA()", // ignored, just use it to compile
]

n_rounds = params.dryRun ? 4 : 20
PT_n_chains = 10
def julia_env_dir = file("julia-environment")
def julia_depot_dir = file(".depot")
def deliv = deliverables(workflow)

workflow {
    args = crossProduct(variables, params.dryRun)
        .filter { it.sampler.startsWith("AH") || it.nleaps == variables.nleaps.first() } // nleaps only relevant to AHMC
    	//.collect()
    	//.view()    
    julia_env = setupPigeons(julia_depot_dir, julia_env_dir)
    agg_path = runSimulation(julia_depot_dir, julia_env, args) | collectCSVs
}

process runSimulation {
    memory { 5.GB * task.attempt }
    time { 2.hour * task.attempt }
    errorStrategy 'retry'
    maxRetries '1'
    input:
        env JULIA_DEPOT_PATH
        path julia_env
        val arg
    output:
        tuple val(arg), path('csvs')
  script:
    template 'scale_main.jl'
}

