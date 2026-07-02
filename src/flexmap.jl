import Distributed: nprocs, nworkers, pmap
import ThreadPools: bmap
import Base.Threads: nthreads

"""
    flexmap(f, itrs...; prefer_threads = false, kwargs...)

Map `f` over the collections `itrs...`
using parallel workers if there is more than one worker,
using threads if there is more than one thread,
or serially if neither of these conditions holds.
If there are both multiple workers and multiple threads,
the threads are preferred if and only if `prefer_threads == true`.
"""
flexmap!(
    f, itrs...;
    prefer_threads = false,
    batch_size = 1,
    kwargs...,
) = begin
    if nthreads() > 1 && (prefer_threads || nprocs() == 1)
        bmap(f, itrs...; kwargs...)
    elseif nprocs() > 1
        pmap(f, itrs...; batch_size=batch_size, kwargs...)
    else
        map(f, itrs...)
    end
    nothing
end

num_workers(;prefer_threads = false,) = begin
    if nthreads() > 1 && (prefer_threads || nprocs() == 1)
        nthreads()::Int
    elseif nprocs() > 1
        nworkers()::Int
    else
        one(Int)
    end
end
