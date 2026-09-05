import ThreadPools: bmap
import Base.Threads: nthreads

"""
    flexmap!(f, itrs...; kwargs...)

Map `f` over the collections `itrs...` using parallel threads if there
is more than one thread or serially if necessary.  Nothing is
returned: the result must be achieved through side-effects of `f`.
"""
flexmap!(
    f, itrs...;
    kwargs...,
) = begin
    if nthreads() > 1
        bmap(f, itrs...; kwargs...)
    else
        map(f, itrs...)
    end
    nothing
end
