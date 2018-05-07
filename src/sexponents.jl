export SExponents

const __EXPONENTS = Dict{UInt64, Matrix}()

struct SExponents{E}
    function SExponents{E}() where { E }
        @assert haskey(__EXPONENTS, E) "No exponent matrix with key $E exists."
        new()
    end
end

function SExponents(exponents::Matrix{<:Integer})
    E = hash(exponents)
    if !haskey(__EXPONENTS, E) || __EXPONENTS[E] != E
        _all_addexponent(E, exponents)
    end
    return SExponents{E}()
end

function _all_addexponent(id, M)
    n = nprocs()
    @sync begin
        for p=1:n
            @async begin
                remotecall_fetch(_addexponent, p, id, M)
            end
        end
    end
end
_addexponent(id, M) = push!(__EXPONENTS, id => M)

"""
    exponents(::SExponents)

Converts exponents stored in a `SExponents` to a matrix.
"""
function exponents(::Type{SExponents{E}}, nvars) where {E}
    __EXPONENTS[E]
end
exponents(::S, nvars) where {S<:SExponents} = exponents(S, nvars)

function Base.show(io::IO, ::Type{SExponents{E}}) where {E}
    print(io, "SExponents{$(E)}")
end
Base.show(io::IO, S::SExponents) = print(io, typeof(S), "()")
