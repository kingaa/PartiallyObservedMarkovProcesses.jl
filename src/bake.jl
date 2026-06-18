import Serialization: serialize, deserialize
import MacroTools: striplines
import Random

"""
    @bake file code

A facility for caching results of computations in a file and
retrieving them when needed. `@bake` checks to see whether the code
used to produce the results matches that of the call, and recomputes
as needed.

`@bake` first parses, then computes a hash of the expression
`code`. It then checks to see if the file at path `file` exists. If it
does, and if the hash stored in this file matches the computed hash,
the results saved in `file` are loaded and returned. If either `file`
does not exist, or if the digests do not match, it evaluates `code`
and stores the result (and the digest) in `file`.

`@bake` uses a serialized representation to store the results.
"""
macro bake(file, code)
    digest = hash(striplines(Meta.parse("$code")))
    expr = quote
        local reuse = false
        local result
        if isfile($file)
            local res = deserialize($file)
            reuse = (res.digest == $digest)
            if !reuse
                @warn "in `bake`: recomputing...."
            end
        end
        if reuse
            result = res.result # COV_EXCL_LINE (false positive)
        else
            result = $code
            serialize($file,(;result=result,digest=$digest))
        end
        result
    end
    esc(expr)
end

"""
    @freeze [rng,] seed, code

Fixes the state of the pseudorandom number generator stream `rng` to
`seed` for the duration of the evaluation of `code`, then restores the
state of `rng` to its original value.  By default `rng = Random.default_rng()`.
"""
macro freeze(rng, seed, code)
    expr = quote
        local rng_state = copy($rng)
        Random.seed!($rng,$seed)
        local result = $code
        copy!($rng,rng_state)
        result
    end
    esc(expr)
end

macro freeze(seed, code)
    expr = quote
        local rng = Random.default_rng()
        @freeze(rng,$seed,$code)
    end
    esc(expr)
end
