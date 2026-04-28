using Distributions

rin = function(;x0,_...)
    d = Poisson(x0)
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

logdpri = function (;x0,_...)
    d = LogNormal(1,1)
    logpdf(d,x0)
end

rpri = function (;x0,args...)
    d = LogNormal(1,1)
    (;x0=rand(d),args...)
end

theta = (a=1.0,k=7.0,x0=5.0);

nothing
