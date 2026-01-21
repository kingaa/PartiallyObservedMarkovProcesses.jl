struct PfilterdPompObject{
    T,P,A,X0,X,Y,F,U,
    Z <: NamedTuple,
    W <: AbstractFloat
    } <: AbstractPompObject{T,P,A,X0,X,Y,F,U}
    pompobj::PompObject{T,P,A,X0,X,Y,F,U}
    Np::Integer
    x0::Array{Z,1}
    filt::Array{Z,2}
    pred::Array{Z,2}
    traj::Vector{Z}
    weights::Array{W,2}
    eff_sample_size::Array{W,1}
    cond_logLik::Array{W,1}
    logLik::W
end


pomp(object::PfilterdPompObject) = object.pompobj


"""
    pfilter(object; Np = 1, params, rinit, rprocess, logmeasure, args...)

`pfilter` runs a basic particle filter.
At least the `rinit`, `rprocess`, and `logdmeasure` basic components are needed.
`args...` can be used to modify or unset additional fields.
"""
pfilter(
    object::ValidPompData;
    Np::Integer = 1,
    params::P = coef(object),
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{PompPlugin,Nothing,Missing} = missing,
    logdmeasure::Union{Function,Nothing,Missing} = missing,
    args...,
) where {P<:NamedTuple} = begin
    object = pomp(
        object;
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        logdmeasure=logdmeasure,
        args...,
    )
    t0 = timezero(object)
    t = times(object)
    y = obs(object)
    x0 = POMP.rinit(object;t0=t0,nsim=Np)
    @assert(size(x0)==(1,Np))
    xf = similar(x0,length(t),Np)
    xp = similar(x0,length(t),Np)
    xt = similar(x0,length(t))
    w = Array{LogLik}(undef,length(t),Np)
    cond_logLik = similar(w,length(t))
    eff_sample_size = similar(w,length(t))
    perm = Array{Int64}(undef,length(t),Np)
    pfilter_internal!(
        object,
        x0,
        reshape(xf,length(t),1,Np),
        reshape(xp,length(t),1,Np),
        reshape(w,length(t),1,Np,1),
        t0,t,
        reshape(y,length(t),1,1),
        eff_sample_size,
        cond_logLik,
        perm
    )
    trace_ancestry!(xt,xf,perm)
    PfilterdPompObject(
        object,Np,
        vec(x0),xf,xp,xt,w,
        eff_sample_size,
        cond_logLik,
        sum(cond_logLik)
    )
end

"""
    pfilter(object; Np = object.Np, args...)

Running `pfilter` on a `PfilterdPompObject` re-runs the particle filter.
One can adjust the parameters, number of particles (`Np`), or pomp model components.
"""
pfilter(
    object::PfilterdPompObject;
    Np::Integer = object.Np,
    args...,
) =
    pfilter(pomp(object;args...);Np=Np)

pfilter(_...) = error("Incorrect call to `pfilter`.")

pfilter_internal!(
    object::AbstractPompObject,
    x0::AbstractArray{X,2},
    xf::AbstractArray{X,3},
    xp::AbstractArray{X,3},
    w::AbstractArray{LogLik,4},
    t0::T,
    t::AbstractVector{T},
    y::AbstractArray{Y,3},
    eff_sample_size::AbstractVector{LogLik},
    cond_logLik::AbstractVector{LogLik},
    perm::AbstractArray{Int64,2}
) where {T,X,Y} = begin
    for k ∈ eachindex(t)
        rprocess!(
            object,
            @view(xp[[k],:,:]);
            x0=x0,
            t0=t0,
            times=@view(t[[k]])
        )
        logdmeasure!(
            object,
            @view(w[[k],:,:,:]);
            times=@view(t[[k]]),
            y=@view(y[[k],:,:]),
            x=@view(xp[[k],:,:])
        )
        pfilt_step_comps!(
            @view(cond_logLik[k]),
            @view(eff_sample_size[k]),
            @view(w[k,1,:,1]),
            @view(perm[k,:]),
            @view(xp[k,1,:]),
            @view(xf[k,1,:]),
        )
        t0 = t[k]
        x0 = view(xf,k,:,:)
    end
end

pfilt_step_comps!(
    logLik::AbstractArray{W,0},
    ess::AbstractArray{W,0},
    w::AbstractVector{W},
    p::AbstractVector{I},
    xp::AbstractVector{X},
    xf::AbstractVector{X},
    n::Int64 = length(w),
) where {W<:AbstractFloat,I<:Integer,X<:NamedTuple} = begin
    wmax::W = -Inf
    s::W = 0
    ss::W = 0
    for k ∈ eachindex(w)
        @inbounds wmax = (w[k] > wmax) ? w[k] : wmax
    end
    if isfinite(wmax)
        for k ∈ eachindex(w)
            @inbounds v::W = exp(w[k]-wmax)
            s += v
            ss += v*v
            @inbounds w[k] = s
        end
        ess[] = s*s/ss
        logLik[] = wmax+log(s/n)
        du::W = s/n
        u::W = -du*rand(LogLik)
        i::I = 1
        for j ∈ eachindex(p)
            u += du
            @inbounds while (u > w[i] && i < n)
                i += 1
            end
            @inbounds p[j] = i
        end
        @inbounds xf[:] = xp[p]
    else
        s = 0
        ss = 0
        wmax = 0
        ess[] = 0
        logLik[] = -Inf
        p[:] = collect(eachindex(p))
        @inbounds xf[:] = xp[:]
    end
end

trace_ancestry!(
    traj::AbstractVector{X},
    filt::AbstractArray{X,2},
    perm::AbstractArray{I,2},
) where {X,I<:Integer} = begin
    @assert size(traj,1)==size(perm,1)
    @assert size(filt)==size(perm)
    j = rand(axes(perm,2))
    for i ∈ Iterators.reverse(axes(perm,1))
        @inbounds traj[i] = filt[i,j]
        @inbounds j = perm[i,j]
    end
end
