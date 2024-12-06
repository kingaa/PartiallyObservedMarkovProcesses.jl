var documenterSearchIndex = {"docs":
[{"location":"#POMP.jl","page":"Home","title":"POMP.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is a Julia implementation of the pomp package for R.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Implementation of POMP models\nSimulation\nParticle filter\nWorkhorses (low-level interface to basic model components)\nHelper functions\nExamples","category":"page"},{"location":"#Function-Documentation","page":"Home","title":"Function Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"mkpath(\"assets/figures\")","category":"page"},{"location":"#Implementation-of-POMP-models","page":"Home","title":"Implementation of POMP models","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pomp","category":"page"},{"location":"#POMP.pomp","page":"Home","title":"POMP.pomp","text":"pomp is the constructor for the PompObject class.\n\npomp(\n    data;\n    t0, times, timevar,\n    params,\n    accumvars,\n    rinit, rprocess,\n    rmeasure, logdmeasure\n    )\n\nArguments\n\ndata: observations. The default constructor takes a vector of NamedTuples as data. One can also supply a DataFrame.\nt0: zero time, t₀.\ntimes: observation times. If data is supplied as a DataFrame, times should be a Symbol which is the time variable in the DataFrame.\ntimevar: optional symbol.  Name of the time variable.\nparams: parameters. A NamedTuple or vector of NamedTuples.\naccumvars: a NamedTuple of state variables to be reset (usually to zero) immediately before each simulation stage.\nrinit: simulator of the latent-state distribution at t₀. This component should be a function that takes parameters and, optionally, t0, the initial time.\nrprocess: simulator of the latent-state process. This component should be a plugin.\nrmeasure: simulator of the measurement process. This component should be a function that takes states, parameters, and, optionally, t, the current time.\nlogdmeasure: log pdf of the measurement process. This component should be a function that takes data, states, parameters, and, optionally, t, the current time.\n\n\n\n\n\nGiven an AbstractPompObject, object, pomp(object) returns the underlying concrete PompObject. Calling pomp(object, args...) returns a copy of object, modified according to args....\n\n\n\n\n\npomp(object::AbstractPompObject; params=missing, accumvars=missing, rinit=missing, rprocess=missing, rmeasure=missing, logdmeasure=missing)\n\nThis form returns a modified version of object. Individual basic components can be modified or removed. The default is to leave them unchanged.\n\n\n\n\n\n","category":"function"},{"location":"#Simulation","page":"Home","title":"Simulation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"simulate","category":"page"},{"location":"#POMP.simulate","page":"Home","title":"POMP.simulate","text":"simulate(object; nsim = 1, params, rinit, rprocess, rmeasure, args...)\n\nsimulate simulates the POMP.  At least the rinit, rprocess, and rmeasure basic components, are needed.\n\n\n\n\n\n","category":"function"},{"location":"#Particle-filter","page":"Home","title":"Particle filter","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pfilter","category":"page"},{"location":"#POMP.pfilter","page":"Home","title":"POMP.pfilter","text":"pfilter(object; Np = 1, params, rinit, rprocess, logmeasure, args...)\n\npfilter runs a basic particle filter. At least the rinit, rprocess, and logdmeasure basic components are needed. args... can be used to modify or unset additional fields.\n\n\n\n\n\npfilter(object; Np = object.Np, args...)\n\nRunning pfilter on a PfilterdPompObject re-runs the particle filter. One can adjust the parameters, number of particles (Np), or pomp model components.\n\n\n\n\n\n","category":"function"},{"location":"#Workhorses","page":"Home","title":"Workhorses","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"rinit\nrinit!","category":"page"},{"location":"#POMP.rinit","page":"Home","title":"POMP.rinit","text":"rinit(object; t0=timezero(object), params=coef(object), nsim=1)\n\nrinit is the workhorse for the simulator of the initial-state distribution.\n\nArguments\n\nobject: the PompObject\nparams: a NamedTuple of parameters or vector of NamedTuples\nt0: the time at which rinit is to be simulated. This should be a single scalar.\nnsim: the number of simulations desired.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.rinit!","page":"Home","title":"POMP.rinit!","text":"rinit!(object, x0; t0=timezero(object), params = coef(object))\n\nrinit! is the in-place version of the rinit workhorse.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"rprocess\nrprocess!","category":"page"},{"location":"#POMP.rprocess","page":"Home","title":"POMP.rprocess","text":"rprocess(object; x0, t0 = timezero(object), times=times(object), params = coef(object))\n\nrprocess is the workhorse for the simulator of the process\n\nIf there is no user-supplied rprocess component, the dynamics are trivial.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.rprocess!","page":"Home","title":"POMP.rprocess!","text":"rprocess!(object, x; x0 = init_state(object), t0 = timezero(object), times=times(object), params = coef(object))\n\nrprocess! is the in-place version of the rprocess workhorse.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"rmeasure","category":"page"},{"location":"#POMP.rmeasure","page":"Home","title":"POMP.rmeasure","text":"rmeasure(object; x, times=times(object), params=coef(object))\n\nrmeasure is the workhorse for the simulator of the measurement distribution.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"logdmeasure\nlogdmeasure!","category":"page"},{"location":"#POMP.logdmeasure","page":"Home","title":"POMP.logdmeasure","text":"logdmeasure(object; times=times(object), y=obs(object), x=states(object), params=coef(object))\n\nlogdmeasure is the workhorse for the evaluator of the log measurement density.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.logdmeasure!","page":"Home","title":"POMP.logdmeasure!","text":"logdmeasure!(object, ell; times=times(object), y=obs(object), x=states(object), params=coef(object))\n\nlogdmeasure! is the in-place version of the logdmeasure workhorse.\n\n\n\n\n\n","category":"function"},{"location":"#Helper-functions","page":"Home","title":"Helper functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"coef\nobs\nstates\ninit_state\ntimes\ntimezero","category":"page"},{"location":"#POMP.coef","page":"Home","title":"POMP.coef","text":"coef(object)\n\ncoef extracts the parameter vector of a PompObject.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.obs","page":"Home","title":"POMP.obs","text":"obs(object)\n\nobs extracts the vector of observables from a PompObject.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.states","page":"Home","title":"POMP.states","text":"states(object)\n\nstates extracts the latent state trajectory of a PompObject.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.init_state","page":"Home","title":"POMP.init_state","text":"init_state(object)\n\ninit_state extracts the latent state at time t0.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.times","page":"Home","title":"POMP.times","text":"times(object)\n\ntimes extracts the time vector from a PompObject.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.timezero","page":"Home","title":"POMP.timezero","text":"timezero(object)\n\ntimezero extracts the zero-time (t0) from a PompObject.\n\n\n\n\n\n","category":"function"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using POMP, RCall\nR\"\"\"\noptions(tidyverse.quiet=TRUE)\nlibrary(tidyverse)\n\"\"\"","category":"page"},{"location":"#The-Gompertz-model","page":"Home","title":"The Gompertz model","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"gompertz","category":"page"},{"location":"#POMP.gompertz","page":"Home","title":"POMP.gompertz","text":"gompertz()\n\ngompertz is a PompObject containing Parus major data and a simple Gompertz population model. The population model has a single scalar state variable, X_t, which obeys\n\nX_t = X_t-1^SK^1-Svarepsilon_t\n\nwhere S = e^-rdeltat and varepsilon_t sim mathrmLogNormal(0sigma_p). The time-step is one unit: deltat=1. The data are assumed to be drawn from a log-normal distribution. In particular,\n\nmathrmpop_t sim mathrmLogNormal(logX_tsigma_m)\n\nParameters\n\nr: the growth rate\nK: the equilibrium population density\nX₀: the initial population density\nσₚ: process noise s.d.\nσₘ: measurement noise s.d.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"View the Parus data:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using POMP, RCall\nP = gompertz()\nd = melt(P)\nR\"\"\"\nsvg(\"assets/figures/gompertz1.svg\",width=7,height=4) #hide\n$d |>\n  pivot_longer(-year) |>\n  ggplot(aes(x=year,y=value))+\n  geom_line()+\n  facet_wrap(~name,scales=\"free_y\",ncol=1)+\n  labs(y=\"\")+\n  theme_bw() -> pl\nprint(pl)\ndev.off() #hide\n\"\"\"\nnothing #hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: gompertz_data)","category":"page"},{"location":"","page":"Home","title":"Home","text":"View a few representative simulations:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using POMP, RCall\nP = gompertz()\nQ = simulate(P;params=(r=4.5,K=210.0,σₚ=0.7,σₘ=0.1,X₀=150.0),nsim=5)\nd = melt(Q,:rep,:parset)\nR\"\"\"\nsvg(\"assets/figures/gompertz2.svg\",width=7,height=5) #hide\n$d |>\n  pivot_longer(-c(year,rep,parset)) |>\n  ggplot(aes(x=year,y=value,group=rep,color=factor(rep)))+\n  geom_line()+\n  geom_point()+\n  facet_wrap(~name,scales=\"free_y\",ncol=1)+\n  labs(y=\"\",color=\"replicate\")+\n  theme_bw() -> pl\nprint(pl)\ndev.off() #hide\n\"\"\"\nnothing #hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: gompertz_sims)","category":"page"},{"location":"#A-simple-SIR-model","page":"Home","title":"A simple SIR model","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"sir","category":"page"},{"location":"#POMP.sir","page":"Home","title":"POMP.sir","text":"sir(\n    β = 0.5, γ = 0.25, N = 10000,\n    ρ = 0.3, k = 10,\n    S₀ = 0.9, I₀ = 0.01, R₀ = 0.1,\n    δt = 0.1, t₀ = 0.0,\n    times = range(start=1.0,stop=90,step=1.0)\n   )\n\nsir returns a PompObject containing simulated SIR data.\n\nParameters\n\nβ: transmission rate\nγ: recovery rate\nN: population size\nρ: reporting rate\nk: overdispersion coefficient (negative binomial size parameter)\nS₀, I₀, R₀: relative proportions of susceptible, infected, recovered (respectively) in the population at t=t₀.\nδt: Euler stepsize\nt₀: zero-time\ntimes: vector of observation times\n\n\n\n\n\n","category":"function"},{"location":"#The-Rosenzweig-MacArthur-model","page":"Home","title":"The Rosenzweig-MacArthur model","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"rmca","category":"page"},{"location":"#POMP.rmca","page":"Home","title":"POMP.rmca","text":"rmca(\n     r = 1, K = 1e4, A = 1e3,\n     b = 1e-3, c = 1, m = 0.8,\n     V = 100, σ = 0.01,\n     N₀ = 3000, P₀ = 4, t₀ = 0.0,\n     δt = 0.01,\n     times=range(start=0,stop=500,step=0.2)\n   )\n\nParameters\n\nr: intrinsic growth rate of prey\nK: carrying capacity for prey\nA: half-saturation prey density\nc: predator foraging rate\nb: predator yield (predators born per prey item killed)\nm: predator death rate\nV: system size\nσ: measurement noise magnitude\nN₀, P₀: initial densities\nt₀: zero-time\nδt: Euler stepsize\ntimes: vector of observation times\n\nObservables\n\nn: prey density\np: predator density\n\nState variables\n\nX = log(N)\nY = log(P)\n\nDetails\n\nrmca returns a PompObject containing simulated data from a Rosenzweig-MacArthur model implemented as an Itô diffusion. Specifically, if N and P are prey and predator densities, respectively, then dN = dG - dC - dS and dP = b dS - dM, where\n\nbeginaligned\ndG = r N dt + sqrtfrac1V r N dW_1 \ndC = fracr N^2K dt + sqrtfrac1V fracr N^2K dW_2 \ndS = fracc N P1+NA dt + sqrtfrac1V fracc N P1+NA dW_3 \ndM = m P dt + sqrtfrac1V m P dW_4 \nendaligned\n\nHere, the dW_i are increments of independent standard Wiener processes. Thus, the process noise scales demographically. Specifically, the system size, V, converts the densities N, P into numbers. It controls the relative magnitude (coefficient of variation) of the demographic process noise. Moreover, V determines a lower threshold on the population sizes, such that if ever N V  1 or P V  1, the population is taken to be extinct. Otherwise, it plays no role in the dynamics. The measurement error is assumed to scale environmentally:\n\nbeginaligned\nn sim mathrmLogNormal(logNsigma) \np sim mathrmLogNormal(logPsigma) \nendaligned\n\nNote that, in the limit Vtoinfty, the Itô diffusion becomes the ordinary differential equation\n\nbeginaligned\nfracdNdt = r N left(1-fracNKright) - fracc N P1+NA \nfracdPdt = fracb c N P1+NA - m P \nendaligned\n\nwhich is the classical Rosenzweig-MacArthur model.\n\nIn this system, the predator is inviable unless R = fracbcAm  1. Even if the predator is viable, the environment is too impoverished to support predators unless R1+fracAK. If the environment is rich enough, and if moreover Rfrac1+fracAK1-fracAK, then the nontrivial equilibrium of the system is unstable. For the default parameters, we have R = 125 and fracAK = 01, so the latter condition holds.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"A sample simulation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using POMP, RCall #hide\nP = rmca(σ=0.1,times=range(0,400.0,step=1.0))\nd = melt(P)\nR\"\"\"\nsvg(\"assets/figures/rmca1.svg\",width=7,height=6) #hide\n$d |>\n  pivot_longer(-time) |>\n  ggplot(aes(x=time,y=value))+\n  geom_path()+\n  facet_wrap(~name,scales=\"free_y\",ncol=1)+\n  labs(y=\"\")+\n  theme_bw() -> pl\nprint(pl)\ndev.off() #hide\n\"\"\"\nnothing #hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: rmca_dynamics)","category":"page"}]
}
