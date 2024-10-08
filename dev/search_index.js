var documenterSearchIndex = {"docs":
[{"location":"#POMP.jl","page":"Home","title":"POMP.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Summary.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Implement POMP models.\nSimulate models with process noise and measurement error.\nParticle filters.","category":"page"},{"location":"#Function-Documentation","page":"Home","title":"Function Documentation","text":"","category":"section"},{"location":"#Basic-constructor","page":"Home","title":"Basic constructor","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pomp\npomp!","category":"page"},{"location":"#POMP.pomp","page":"Home","title":"POMP.pomp","text":"pomp is the constructor for the PompObject class.\n\npomp(\n    data;\n    t0, times,\n    accumvars,\n    rinit, rprocess,\n    rmeasure, logdmeasure\n    )\n\nArguments\n\ndata: observations. The default constructor takes a vector of NamedTuples as data. One can also supply a DataFrame.\nt0: zero time, t₀.\ntimes: observation times. If data is supplied as a DataFrame, times should be a Symbol which is the time variable in the DataFrame.\naccumvars: a NamedTuple of state variables to be reset (usually to zero) immediately following each observation.\nrinit: simulator of the latent-state distribution at t₀.\nrprocess: simulator of the latent-state process.\nrmeasure: simulator of the measurement process.\nlogdmeasure: log pdf of the measurement process.\n\n\n\n\n\nAlternatively, one can construct a PompObject from a DataFrame. In this case, times should be the Symbol of the time-variable.\n\n\n\n\n\nGiven an AbstractPompObject, object, pomp(object) returns the underlying concrete PompObject. Calling pomp(object, args...) returns a copy of object, modified according to args....\n\n\n\n\n\n","category":"function"},{"location":"#POMP.pomp!","page":"Home","title":"POMP.pomp!","text":"pomp!(object, args...) modifies the PompObject object in place. Individual basic components can be modified, set, or unset.\n\n\n\n\n\n","category":"function"},{"location":"#Simulation","page":"Home","title":"Simulation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"simulate\nsimulate!","category":"page"},{"location":"#POMP.simulate","page":"Home","title":"POMP.simulate","text":"simulate(object; params, nsim = 1, args...)\n\nsimulate simulates the POMP. args... can be used to modify or unset fields.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.simulate!","page":"Home","title":"POMP.simulate!","text":"simulate!(object; args...)\n\nsimulate! simulates in place. args... can be used to modify or unset fields.\n\n\n\n\n\n","category":"function"},{"location":"#Particle-filter","page":"Home","title":"Particle filter","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pfilter","category":"page"},{"location":"#POMP.pfilter","page":"Home","title":"POMP.pfilter","text":"pfilter(object; params, Np = 1, args...)\n\npfilter runs a basic particle filter. args... can be used to modify or unset fields.\n\n\n\n\n\n","category":"function"},{"location":"#Workhorses","page":"Home","title":"Workhorses","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"rinit","category":"page"},{"location":"#POMP.rinit","page":"Home","title":"POMP.rinit","text":"rinit(object; t0=timezero(object), params=coef(object), nsim=1)\n\nrinit is the workhorse for the simulator of the initial-state distribution.\n\nThe user can supply an rinit component as a function that takes parameters and, optionally, t0, the initial time.\n\nArguments\n\nobject: the PompObject\nparams: a NamedTuple of parameters\nt0: the time at which rinit is to be simulated. This should be a single scalar.\nnsim: the number of simulations desired.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"rmeasure","category":"page"},{"location":"#POMP.rmeasure","page":"Home","title":"POMP.rmeasure","text":"rmeasure(object; x, times=times(object), params=coef(object))\n\nrmeasure is the workhorse for the simulator of the measurement distribution.\n\nThe user can supply an rmeasure component as a function that takes states, parameters, and, optionally, t, the current time.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"rprocess!","category":"page"},{"location":"#POMP.rprocess!","page":"Home","title":"POMP.rprocess!","text":"rprocess(object; x0, t0 = timezero(object), times=times(object), params=coef(object))\n\nrprocess is the workhorse for the simulator of the process\n\nThe user can supply an rprocess component as a function that takes states, parameters, and current time (t) and returns the updated time and state.\n\nIf there is no user-supplied rprocess component, the dynamics are trivial.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"logdmeasure!","category":"page"},{"location":"#POMP.logdmeasure!","page":"Home","title":"POMP.logdmeasure!","text":"logdmeasure(object; times=times(object), y = obs(object), x = states(object), params)\n\nlogdmeasure is the workhorse for the evaluator of the log measurement density.\n\nThe user can supply a logdmeasure component as a function that takes data, states, parameters, and, optionally, t, the current time.\n\n\n\n\n\n","category":"function"},{"location":"#Helper-functions","page":"Home","title":"Helper functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"coef\nstates\nobs\ntimes\ntimezero","category":"page"},{"location":"#POMP.coef","page":"Home","title":"POMP.coef","text":"coef(object::SimPompObject)\n\ncoef extracts the parameters stored in a SimPompObject.\n\n\n\n\n\ncoef(object::PfilterdPompObject)\n\ncoef extracts the parameters stored in a PfilterdPompObject.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.states","page":"Home","title":"POMP.states","text":"states(object::AbstractSimPompObject)\n\nstates extracts the state trajectory of a SimPompObject.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.obs","page":"Home","title":"POMP.obs","text":"obs(object)\n\nobs extracts the time vector of a PompObject.\n\n\n\n\n\nobs(object::AbstractSimPompObject)\n\nobs extracts the observations of a SimPompObject.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.times","page":"Home","title":"POMP.times","text":"times(object)\n\ntimes extracts the time vector of a PompObject.\n\n\n\n\n\n","category":"function"},{"location":"#POMP.timezero","page":"Home","title":"POMP.timezero","text":"timezero(object)\n\ntimezero extracts the zero-time (t0) of a PompObject.\n\n\n\n\n\n","category":"function"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"sir\ngompertz","category":"page"},{"location":"#POMP.sir","page":"Home","title":"POMP.sir","text":"sir(\n    β = 0.5, γ = 0.25, N = 10000,\n    ρ = 0.3, k = 10,\n    S₀ = 0.9, I₀ = 0.01, R₀ = 0.1,\n    δt = 0.1, t₀ = 0.0,\n    times = [t for t ∈ range(start=1.0,stop=90,step=1.0)]\n   )\n\nsir returns a SimPompObject containing simulated SIR data.\n\nParameters\n\nβ: transmission rate\nγ: recovery rate\nN: population size\nρ: reporting rate\nk: overdispersion coefficient (negative binomial size parameter)\nS₀, I₀, R₀: relative proportions of susceptible, infected, recovered (respectively) in the population at t=t₀.\nδt: Euler stepsize\nt₀: zero-time\ntimes: vector of observation times\n\n\n\n\n\n","category":"function"},{"location":"#POMP.gompertz","page":"Home","title":"POMP.gompertz","text":"gompertz()\n\ngompertz is a PompObject containing Parus major data and a simple Gompertz population model. The population model has a single scalar state variable, X_t, which obeys\n\nX_t = X_t-1^SK^1-Svarepsilon_t\n\nwhere S = e^-rdeltat and varepsilon_t sim mathrmLogNormal(0sigma_p). The time-step is one unit: deltat=1. The data are assumed to be drawn from a log-normal distribution. In particular,\n\nY_t sim mathrmLogNormal(logX_tsigma_m)\n\nParameters\n\nr: the growth rate\nK: the equilibrium population density\nX₀: the initial population density\nσₚ: process noise s.d.\nσₘ: measurement noise s.d.\n\n\n\n\n\n","category":"function"}]
}
