import Serialization: serialize, deserialize
import MacroTools: striplines
import Random

"""
    @bake dish recipe

A facility for caching results of computations in a file and
retrieving them when needed. `@bake` checks to see whether the code
used to produce the results matches that of the call, and recomputes
as needed.

`@bake` first parses, then computes a hash of the expression
`recipe`. It then checks to see if the file at path `dish` exists. If
it does, and if the hash stored in this file matches the computed
hash, the results saved in `dish` are loaded and returned. If either
`dish` does not exist, or if the digests do not match, it evaluates
`recipe` and stores the result (and the digest) in `dish`.

`@bake` uses a serialized representation to store the results.
"""
macro bake(dish, recipe)
    digest = hash(striplines(Meta.parse("$recipe")))
    quote
        local file = $(esc(dish))
        local reload = false
        local res, result
        if isfile(file)
            res = deserialize(file)
            reload = (res.digest == $digest)
            if !reload
                @warn "in `bake`: recomputing '$file'."
            end
        end
        if reload
            result = res.result # COV_EXCL_LINE (false positive)
        else
            result = $(esc(recipe))
            serialize(file,(;digest=$digest,result=result))
        end
        result
    end
end

"""
    @freeze [rng,] seed, code

Fixes the state of the pseudorandom number generator stream `rng` to
`seed` for the duration of the evaluation of `code`, then restores the
state of `rng` to its original value.  By default `rng = Random.default_rng()`.
"""
macro freeze(rng, seed, code)
    quote
        local rng_state = copy($(esc(rng)))
        Random.seed!($(esc(rng)), $seed)
        local result = $(esc(code))
        copy!($(esc(rng)), rng_state)
        result
    end
end

macro freeze(seed, code)
    quote
        @freeze($(esc(Random.default_rng())), $seed, $(esc(code)))
    end
end
