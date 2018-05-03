export Polynomial, coefficients, exponents, nvariables, coefficienttype


"""
    Polynomial(f::MP.AbstractPolynomial, [variables])

Construct a Polynomial from `f`.
"""
struct Polynomial{T, NVars, E<:SExponents}
    coefficients::Vector{T}
    variables::SVector{NVars, Symbol}
    deboor_coefficients::Vector{T}
    deboor_invperm::Vector{Int}

    function Polynomial{T, NVars, SExponents{E}}(coefficients::Vector{T}, variables::SVector{NVars, Symbol}) where {T, NVars, E}
        @assert length(coefficients) == div(length(E), NVars) "Coefficients size does not match exponents size"
        deboor_coefficients = coefficients[deboor_perm(exponents(SExponents{E}, NVars))]
        new(coefficients, variables, deboor_coefficients)
    end
end

function deboor_perm(E)
    exponents = [E[:,j] for j=1:size(E, 2)]
    D = maximum(sum, exponents)
    @show D
    by_degree = [Vector{Vector{Int}}() for i=0:D]
    coefficients = Dict{Vector{Int}, Int}()
    for (i, e) in enumerate(exponents)
        d = sum(e)
        push!(by_degree[d + 1], e)
        coefficients[e] = i
    end

    perm = Int[]
    for d=D:-1:0
        d_exponents = sort!(by_degree[d + 1], lt=lexless, rev=true)
        for e in d_exponents
            push!(perm, coefficients[e])
        end
    end
    perm
end


function Polynomial(coefficients::Vector{T}, nvars, exponents::E, variables) where {T, E<:SExponents}
    return Polynomial{T, nvars, E}(coefficients, variables)
end

function Polynomial(coefficients::Vector{T}, exponents::Matrix{<:Integer}, variables=SVector((Symbol("x", i) for i=1:size(exponents, 1))...)) where {T}
    nvars = size(exponents, 1)
    @assert length(coefficients) == size(exponents, 2) "Coefficients size does not match exponents size"
    E, p = revlexicographic_cols(exponents)
    return Polynomial(coefficients[p], nvars, SExponents(E), variables)
end

# Implementation from Base.sort adapted to also reorder an associated vector
"""
    revlexicographic_cols(A, v)

Sorts the columns of `A` in reverse lexicographic order and returns the permutation vector
to obtain this ordering.
"""
function revlexicographic_cols(A::AbstractMatrix; kws...)
    inds = Compat.axes(A,2)
    T = Base.Sort.slicetypeof(A, :, inds)
    cols = map(i -> (@view A[end:-1:1, i]), inds)
    if VERSION <= v"0.6.2"
        p = sortperm(cols; kws..., order=Base.Sort.Lexicographic)
    else
        p = sortperm(cols; kws...)
    end
    A[:,p], p
end

function Polynomial(p::MP.AbstractPolynomialLike{T}, vars = MP.variables(p)) where T
    terms = MP.terms(p)
    nterms = length(terms)
    nvars = length(vars)

    exponents = Matrix{Int}(undef, nvars, nterms)
    coefficients = Vector{T}(undef, nterms)
    for (j, term) in enumerate(terms)
        coefficients[j] = MP.coefficient(term)
        for (i, var) in enumerate(vars)
            exponents[i, j] = MP.degree(term, var)
        end
    end
    Polynomial(coefficients, exponents, SVector((Symbol(var) for var in vars)...))
end

"""
    coefficients(f)

Return the coefficients of `f`.
"""
coefficients(f::Polynomial) = f.coefficients

"""
    exponents(f)

Return the exponents of `f` as an matrix where each column represents
the exponents of a monomial.
"""
function exponents(::Polynomial{T, NVars, E}) where {T, NVars, E<:SExponents}
    exponents(E, NVars)
end

"""
    nvariables(f::Polynomial)

Return the number of variables `f`.
"""
nvariables(::Polynomial{T, NVars}) where {T, NVars} = NVars

"""
    coefficienttype(f::Polynomial)

Return the type of the coefficients of `f`.
"""
coefficienttype(::Polynomial{T, NVars}) where {T, NVars} = T


function Base.:(==)(f::Polynomial{T, NVars, E}, g::Polynomial{T, NVars, E}) where {T, NVars, E<:SExponents}
    coefficients(f) == coefficients(g)
end
