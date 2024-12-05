using Distributions

rin = function(;x₀,_...)
    d = Poisson(x₀)
    (x=rand(d),)
end

rlin = function (;t,a,x,_...)
    d = Poisson(a*x)
    (x=rand(d),)
end

rmeas = function (;x,k,_...)
    d = NegativeBinomial(k,k/(k+x))
    (y=rand(d),)
end

logdmeas = function (;x,y,k,_...)
    d = NegativeBinomial(k,k/(k+x))
    logpdf(d,y)
end

theta = (a=1.0,k=7.0,x₀=5.0);

nothing
