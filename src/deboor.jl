@generated function evaluate_deboor9(f::Polynomial{T, NVars, E}, x::AbstractVector) where {T, NVars, E}
    evaluate_deboor_impl(f)
end

function evaluate_deboor_impl(f::Type{Polynomial{T, NVars, E}}) where {T, NVars, E<:SExponents}
    quote
        @boundscheck length(x) ≥ NVars
        c = f.deboor_coefficients
        @inbounds out = $(generate_evaluate_deboor(exponents(E, NVars), T))
        out
    end
end

# function coefficient_permutation(E)
#     exponents = [E[:,j] for j=1:size(E, 2)]
#     D = maximum(sum, E)
#     by_degree = [Vector{Vector{Int}}() for i=0:D]
#     coefficients = Dict{Vector{Int}, Int}()
#     for (i, e) in enumerate(exponents)
#         d = sum(e)
#         push!(by_degree[d + 1], e)
#         coefficients[e] = i
#     end
#
#     perm = Int[]
#     for d=D:-1:0
#         d_exponents = sort!(by_degree[d + 1], lt=lexless, rev=true)
#         for e in d_exponents
#             push!(perm, coefficients[e])
#         end
#     end
#     perm
# end


function generate_evaluate_deboor(E, ::Type{T}) where T
    exprs = Expr[]
    perm = invperm(deboor_perm(E))

    exponents = [E[:,j] for j=1:size(E, 2)]
    D = maximum(sum, exponents)

    coefficients = Dict{Vector{Int}, Union{Expr, Symbol}}()
    by_degree = [Vector{Vector{Int}}() for i=0:D]
    for (i, e) in enumerate(exponents)
        d = sum(e)

        push!(by_degree[d + 1], e)
        coefficients[e] = :(c[$(perm[i])])
    end

    cmp = Compat.Fix2(>, 0)

    for d=D:-1:1
        d_exponents = sort!(by_degree[d + 1], lt=lexless, rev=true)

        for e in d_exponents
            i = findfirst(cmp, e)

            @gensym c

            new_e = copy(e)
            new_e[i] -= 1
            coeff = pop!(coefficients, e)

            if haskey(coefficients, new_e)
                ncoeff = pop!(coefficients, new_e)
                push!(exprs, :($c = muladd(x[$i], $coeff, $ncoeff)))
            else
                push!(by_degree[d], new_e)
                push!(exprs, :($c = muladd(x[$i], $coeff, zero($T))))
            end

            coefficients[new_e] = c
        end
    end
    Expr(:block, exprs...)
end


#
# D = 4
# exponents = [[0, 0, 2, 2], [0, 2, 0, 2], [0, 2, 2, 0]]
#
#
#
#
# exprs
# by_degree
# coefficients
#
#
#
#
#
#
#
#
#
# foreach_exponent()
#
#
# function update_exponent!(exponent, d)
#     i = N = length(exponent)
#     while i > 0
#         if exponent[i] == d
#             for j=1:N-1
#                 if j == i - 1
#                     exponent[j] = 1
#                 else
#                     exponent[j] = 0
#                 end
#             end
#             exponent[N] = d - 1
#             return exponent
#         elseif 0 < exponent[i] # < iter.d
#             i₀ = i - 1
#             exponent[N] = exponent[i] - 1
#             if i₀ ≥ 1
#                 exponent[i₀] += 1
#             end
#             for j=max(1, i):N - 1
#                 exponent[j] = 0
#             end
#
#             return exponent
#         end
#
#         i -= 1
#     end
#
#     exponent
# end
#
# x = [0, 1, 0, 1, 0, 2, 0, 0, 0, 2, 0]
# sum(x)
# @btime $x == $x
#
# @btime update_exponent!($x, $6)
#
# struct MonomialIterator{N}
#     d::Int
# end
#
#
# function Base.start(iter::MonomialIterator{N}) where N
#     (ntuple(i -> i == N ? iter.d: 0, Val{N}), false)
# end
# Base.eltype(::MonomialIterator{N}) where N = NTuple{N, Int}
# Base.done(iter::MonomialIterator, state) = state.last[1] == iter.d
# function Base.next(iter::MonomialIterator{N}, state) where N
#     last, first
#     if !state.first
#         return last, MonomialIteratorState(last, true)
#     end
#
#     i = N
#     while i > 0
#         if last[i] == iter.d
#             let i=i
#                 next::NTuple{N, Int} = ntuple(N) do j
#                     if j == i - 1
#                         1
#                     elseif j == N
#                         iter.d - 1
#                     else
#                         0
#                     end
#                 end
#                 return next, MonomialIteratorState(next, true)
#             end
#         elseif 0 < last[i] # < iter.d
#             let i=i
#                 next::NTuple{N, Int} = ntuple(N) do j
#                     if j < i - 1
#                         last[j]
#                     elseif j == i - 1
#                         last[j] + 1
#                     elseif j == N
#                         last[i] - 1
#                     else
#                         0
#                     end
#                 end::NTuple{N, Int}
#                 return next, MonomialIteratorState(next, true)
#             end
#         end
#
#         i -= 1
#     end
#     last, state
# end
#
#
# struct MonomialIterator{N}
#     d::Int
# end
# struct MonomialIteratorState2
#     last::Vector{Int}
#     first::Bool
# end
#
# function Base.start(iter::MonomialIterator{N}) where N
#     MonomialIteratorState2([i == N ? iter.d : 0 for i=1:N], false)
# end
# Base.eltype(::MonomialIterator{N}) where N = Vector{Int}
# Base.done(iter::MonomialIterator, state) = state.last[1] == iter.d
# function Base.next(iter::MonomialIterator{N}, state) where N
#     last, first = state.last, state.first
#     if !first
#         return last, MonomialIteratorState2(last, true)
#     end
#
#     i = N
#     while i > 0
#         if last[i] == iter.d
#             let i=i
#                 for j=1:N
#                     if j == i - 1
#                         last[j] = 1
#                     elseif j == N
#                         last[j] = iter.d - 1
#                     else
#                         last[j] = 0
#                     end
#                 end
#                 return last, state
#             end
#         elseif 0 < last[i] # < iter.d
#             let i=i
#                 i₀ = i - 1
#                 last[N] = last[i] - 1
#                 if i₀ ≥ 1
#                     last[i₀] += 1
#                 end
#                 for j=max(1, i):N - 1
#                     last[j] = 0
#                 end
#
#                 return last, state
#             end
#         end
#
#         i -= 1
#     end
#     last, state
# end
#
# Base.length(iter::MonomialIterator{N}) where N = binomial(iter.d + N - 1, N - 1)
# iter = MonomialIterator{4}(2)
# for x in MonomialIterator{4}(2)
#     println(x)
# end
#
# function iterall(iter)
#     for i in iter
#     end
#     true
# end
#
# iter = MonomialIterator{4}(5)
# state =  start(iter)
# _, state = next(iter, state)
#
# iter = MonomialIterator{10}(5)
# @btime iterall($iter)
# @trace next(iter, state)
#
# @code_warntype next(iter, state)
# @btime next($iter, $state)
#
# v = @SVector zeros(5)
#
# setindex(v, 2, 3)
#
# function next(v)
#     d = sum(v)
#
#     for k=length(v):-1:2
#         if v[k] > 0
#     end
# end
#
#
# #
# # function Base.start(iter::MonomialIterator{N}) where N
# #     last = ntuple(i -> i == N ? iter.d : 0, Val{N})
# #     k = N
# #     last, k, false
# # end
# # Base.eltype(::MonomialIterator{N}) where N = NTuple{N, Int}
# # Base.done(iter::MonomialIterator, state) = state[1][1] == iter.d
# # function Base.next(iter::MonomialIterator{N}, state) where N
# #     # @show state
# #     last, k, first = state
# #     !first && return last, (last, k, true)
# #
# #     for i=N:-1:1
# #         if last[i] == iter.d
# #             next = ntuple(Val{N}) do j
# #                 if j == i - 1
# #                     1
# #                 elseif j == N
# #                     iter.d - 1
# #                 else
# #                     0
# #                 end
# #             end
# #             return next, (next, k, true)
# #         elseif 0 < last[i] # < iter.d
# #             next = ntuple(Val{N}) do j
# #                 if j < i - 1
# #                     last[j]
# #                 elseif j == i - 1
# #                     last[j] + 1
# #                 elseif j == N
# #                     last[i] - 1
# #                 else
# #                     0
# #                 end
# #             end
# #             return next, (next, k, true)
# #         end
# #     end
# #     @assert false
# #     last, state
# # end
