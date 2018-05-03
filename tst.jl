using StaticArrays, BenchmarkTools
using DynamicPolynomials: @polyvar
using StaticPolynomials
using TestSystems


g = Polynomial(equations(katsura10())[3])

v = @SVector rand(11)

@btime StaticPolynomials.evaluate_deboor8($g, $v)
@btime StaticPolynomials.evaluate($g, $v)

@polyvar x y z

f = Polynomial(x^2+y^2+z^2)

w = @SVector rand(3)

StaticPolynomials.evaluate(f, w)
StaticPolynomials.exponents(f)
StaticPolynomials.evaluate_deboor_impl(typeof(f))

StaticPolynomials.evaluate_deboor(f, w)

function f10(x)
   f  = 48*x[1]^3 + 72*x[1]^2*x[2] + 72*x[1]^2*x[3] + 72*x[1]^2*x[4] + 72*x[1]^2*x[5] + 72*x[1]^2*x[7]
   f += 72*x[1]^2*x[8] + 72*x[1]*x[2]^2 + 144*x[1]*x[2]*x[4] + 144*x[1]*x[2]*x[7] + 72*x[1]*x[3]^2
   f += 144*x[1]*x[3]*x[5] + 144*x[1]*x[3]*x[8] + 72*x[1]*x[4]^2 + 144*x[1]*x[4]*x[7] + 72*x[1]*x[5]^2
   f += 144*x[1]*x[5]*x[8] + 72*x[1]*x[7]^2 + 72*x[1]*x[8]^2 + 48*x[2]^3 + 72*x[2]^2*x[3]
   f += 72*x[2]^2*x[4] + 72*x[2]^2*x[6] + 72*x[2]^2*x[7] + 72*x[2]^2*x[9] + 72*x[2]*x[3]^2
   f += 144*x[2]*x[3]*x[6] + 144*x[2]*x[3]*x[9] + 72*x[2]*x[4]^2 + 144*x[2]*x[4]*x[7] + 72*x[2]*x[6]^2
   f += 144*x[2]*x[6]*x[9] + 72*x[2]*x[7]^2 + 72*x[2]*x[9]^2 + 48*x[3]^3 + 72*x[3]^2*x[5]
   f += 72*x[3]^2*x[6] + 72*x[3]^2*x[8] + 72*x[3]^2*x[9] + 72*x[3]*x[5]^2 + 144*x[3]*x[5]*x[8]
   f += 72*x[3]*x[6]^2 + 144*x[3]*x[6]*x[9] + 72*x[3]*x[8]^2 + 72*x[3]*x[9]^2 + 48*x[4]^3
   f += 72*x[4]^2*x[5] + 72*x[4]^2*x[6] + 72*x[4]^2*x[7] + 72*x[4]^2*x[10] + 72*x[4]*x[5]^2
   f += 144*x[4]*x[5]*x[6] + 144*x[4]*x[5]*x[10] + 72*x[4]*x[6]^2 + 144*x[4]*x[6]*x[10] + 72*x[4]*x[7]^2
   f += 72*x[4]*x[10]^2 + 48*x[5]^3 + 72*x[5]^2*x[6] + 72*x[5]^2*x[8] + 72*x[5]^2*x[10]
   f += 72*x[5]*x[6]^2 + 144*x[5]*x[6]*x[10] + 72*x[5]*x[8]^2 + 72*x[5]*x[10]^2 + 48*x[6]^3
   f += 72*x[6]^2*x[9] + 72*x[6]^2*x[10] + 72*x[6]*x[9]^2 + 72*x[6]*x[10]^2 + 48*x[7]^3
   f += 72*x[7]^2*x[8] + 72*x[7]^2*x[9] + 72*x[7]^2*x[10] + 72*x[7]*x[8]^2 + 144*x[7]*x[8]*x[9]
   f += 144*x[7]*x[8]*x[10] + 72*x[7]*x[9]^2 + 144*x[7]*x[9]*x[10] + 72*x[7]*x[10]^2 + 48*x[8]^3
   f += 72*x[8]^2*x[9] + 72*x[8]^2*x[10] + 72*x[8]*x[9]^2 + 144*x[8]*x[9]*x[10] + 72*x[8]*x[10]^2
   f += 48.0*x[9]^3 + 72*x[9]^2*x[10] + 72*x[9]*x[10]^2
   return f
end

@polyvar x[1:10]
p10 = Polynomial(f10(x))

w = @SVector rand(Complex128, 10)
p10(w)
StaticPolynomials.evaluate_deboor(p10, w)

@btime evaluate($p10, $w)
@btime StaticPolynomials.evaluate_deboor6($p10, $w)

@code_native StaticPolynomials.evaluate_deboor(p10, w)
StaticPolynomials.evaluate_deboor_impl(typeof(p10))

@btime gradient($p10, $w)








function katsura(n)
  @polyvar x[0:n] # This creates variables x0, x1, ...

  return [
    (sum(x[abs(l)+1]*x[abs(m-l)+1] for l=-n:n if abs(m-l)<=n) -
    x[m+1] for m=0:n-1)...,
    x[1] + 2sum(x[i+1] for i=1:n) - 1
  ]
end

F = system(katsura(15))

f1 = F.f1
xx = @SVector rand(16)

@btime StaticPolynomials.gradient($f1, $xx)


@btime StaticPolynomials.gradient($f1, $xx)

StaticPolynomials._val_gradient_impl(typeof(f1))


@time jacobian(F, xx)
f = P.f15
@btime StaticPolynomials.evaluate($F, $xx)
@btime StaticPolynomials.jacobian($F, $xx) #   300.358 ns (0 allocations: 0 bytes)
