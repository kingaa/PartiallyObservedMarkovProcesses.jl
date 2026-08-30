"""
    paramsymbs(object)

Attempt to determine the names of parameters.  Return the result as a
`Vector{Symbol}`.  If `object` is an `AbstractPompObject` without
`init_state` or `states`, this function may mistakenly include state
variables.
"""
paramsymbs(m::Method) = Base.kwarg_decl(m)

paramsymbs(f::Function) = begin
    m = methods(f)
    @assert length(m)==1 "In specifying basic model components, do not use functions with more than one method!"
    paramsymbs(first(m))
end

paramsymbs(f::Nothing) = Symbol[]

paramsymbs(f::NamedTuple) = Symbol[keys(f)...]

paramsymbs(object::AbstractPompObject) = paramsymbs(pomp(object))

paramsymbs(
    object::PompObject{T,X,A,Y,F,U},
) where {T,X,A,Y,F,U} = begin
    components = [:rinit, :rprocess, :rmeasure, :logdmeasure, :rprior, :logdprior, :params]
    symbs = map(components) do c
        paramsymbs(getfield(object,c))
    end
    excls = [
        Symbol[object.timevar,Symbol("_...")],
        fieldnames(A),fieldnames(X),fieldnames(U),fieldnames(Y)
    ]
    setdiff(union(symbs...),union(excls...))
end

paramsymbs(f::VectorfieldPlugin) = begin
    setdiff(paramsymbs(f.vf),[:t,f.statenames...])
end

argnames(m::Method) = Base.rest(Base.method_argnames(m),2)

argnames(f::Function) = begin
    m = methods(f)
    @assert length(m)==1 "In specifying basic model components, do not use functions with more than one method!"
    argnames(first(m))
end

