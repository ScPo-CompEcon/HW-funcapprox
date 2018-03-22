module funcapp

using Plots
using FastGaussQuadrature
using ApproXD
using Distributions
using ApproxFun

# use chebyshev to interpolate this:
function q1(n)
f(x::LinSpace{Float64})= x + 2 * x^2 - exp(-x)
f(x::Vector{Float64}) = x + 2.*x.^2 - exp(-x)
lb = -3
ub = 3
n = 100
deg = n-1

nodes = gausschebyshev(n)[1]
x_nodes = 0.5 * (lb + ub) + 0.5 * (ub - lb) * nodes
chubbychef = Float64[cos((n - i + 0.5) * (j - 1)* pi / n) for i = 1:n, j = 1:n]
c = \(chubbychef,f(x_nodes))

n_new = 100
nodes_new = linspace(lb,ub,n_new)
x_new = [lb + (i - 1)/(n_new - 1)*(ub - lb) for i in 1:n_new]
z = 2 * (x_new - lb)/(ub - lb) - 1
chubbychef_new = [cos(acos(z[i])*j) for i = 1:n_new, j = 0:deg]
f_hat = chubbychef_new*c
true_y = f.(nodes_new)
err = true_y - f_hat
return Dict(:err => 1.0)
end

function q2(n)
function f(x_afun)
return(x_afun + 2 * x_afun^2 - exp(-x_afun))
end
f_afun = Fun(f, Chebyshev(-3..3))
n_afun = 100
x_afun = linspace(-3,3,n_afun)
true_y = f.(x_afun)
cheby_y = f_afun.(x_afun)
error = true_y - cheby_y
plot1 = plot(x_afun, [true_y, cheby_y], label = ["f(x)" "f_afun(x)"])
plot2 = plot(x_afun, true_y - cheby_y, label = "f(x) - f_afun(x)")
end

# plot the first 9 basis Chebyshev Polynomial Basisi Fnctions
function q3()
n = 100
m = 9
nodes = gausschebyshev(n)[1]
chubbychef = [cos((n - i + 0.5) * (j - 1)* π / n) for i = 1:n, j = 1:n]
p = Any[]
push!(p, plot(chubbychef[:, 1], label = "Basis 1", ylim=(-1,1)))
push!(p, plot(chubbychef[:, 2], label = "Basis 2", ylim=(-1,1)))
push!(p, plot(chubbychef[:, 3], label = "Basis 3", ylim=(-1,1)))
push!(p, plot(chubbychef[:, 4], label = "Basis 4", ylim=(-1,1)))
push!(p, plot(chubbychef[:, 5], label = "Basis 5", ylim=(-1,1)))
push!(p, plot(chubbychef[:, 6], label = "Basis 6", ylim=(-1,1)))
push!(p, plot(chubbychef[:, 7], label = "Basis 7", ylim=(-1,1)))
push!(p, plot(chubbychef[:, 8], label = "Basis 8", ylim=(-1,1)))
push!(p, plot(chubbychef[:, 9], label = "Basis 9", ylim=(-1,1)))
plot(p...)
end

ChebyT(x,deg) = cos(acos(x)*deg)
unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]
type ChebyType
 f::Function # function to approximate
 nodes::Union{Vector,LinSpace} # evaluation points
 basis::Matrix # basis evaluated at nodes
 coefs::Vector # estimated coefficients
 deg::Int 	# degree of chebypolynomial
 lb::Float64 # bounds
 ub::Float64
# constructor
 function ChebyType(_nodes::Union{Vector,LinSpace},_deg,_lb,_ub,_f::Function)
  n = length(_nodes)
  y = _f.(_nodes)
  _basis = Float64[ChebyT(unitmap(_nodes[i],_lb,_ub),j) for i=1:n,j=0:_deg]
  _coefs = _basis \ y  # type `?\` to find out more about the backslash operator. depending the args given, it performs a different operation
  # create a ChebyType with those values
  new(_f,_nodes,_basis,_coefs,_deg,_lb,_ub)
 end
end
# function to predict points using info stored in ChebyType
function predict(Ch::ChebyType,x_new)
 true_new = Ch.f.(x_new)
 basis_new = Float64[ChebyT(unitmap(x_new[i],Ch.lb,Ch.ub),j) for i=1:length(x_new),j=0:Ch.deg]
 basis_nodes = Float64[ChebyT(unitmap(Ch.nodes[i],Ch.lb,Ch.ub),j) for i=1:length(Ch.nodes),j=0:Ch.deg]
 preds = basis_new * Ch.coefs
 preds_nodes = basis_nodes * Ch.coefs
 return Dict("x"=> x_new,"truth"=>true_new, "preds"=>preds, "preds_nodes" => preds_nodes)
end

function q4a(deg=(5,9,15),lb=-5,ub=5)
# Couldn't solve
#r(x) = 1./(1 + 25 .* x.^2)
#Eqspaced = Dict()
#Chebynodes = Dict()
print("could not solve")
end

function q4b()
r(x) = 1./(1 + 25 .* x.^2)
b = BSpline(13, 3, -5.0, 5.0)
Basis = full(getBasis(collect(linspace(-5,5,65)),b))
y = r.(linspace(-5,5,65))
c = Basis \ y
x_new = linspace(-5,5,100)
true_y = r.(x_new)
Basis1 = full(getBasis(collect(x_new),b))
Bs_y1 = Basis1*c
my_knots = [-2.5, -2, -1.5, -1, -0.5, -0.25, 0., 0.25, 0.5, 1, 1.5, 2, 2.5]
b2 = BSpline(my_knots, 3)
Basis2 = full(getBasis(collect(x_new),b2))
y = r.(x_new)
c2 = Basis2 \ y
Basis3 = full(getBasis(collect(x_new),b2))
Bs_y2 = Basis3*c2
plot1 = plot(x_new, [true_y, Bs_y1, Bs_y2],
labels = ["Runge" "Uniform knots" "Clustered at zero"])
plot2 = scatter(knots, zeros(15), labels = "knot positions")
plot2 = plot!(x_new, [true_y - Bs_y1, true_y-Bs_y2], labels = ["Uniform knots" "knots"])
end

function q5()
f(x) = abs(x).^.5
knots = 13
deg = 5
lb = -1
ub = 1
y = linspace(lb,ub,65)
bs_unif = BSpline(knots,3,lb,ub)
B_unif = full(getBasis(collect(y),bs_unif))
c_unif = B_unif\f(y)
approx_g_unif = vec(B_unif*c_unif)
knots_mult1 = vec([-1 -0.8 -0.6 -0.4 -0.2 -0.01 0 0.01 0.2 0.4 0.6 0.8 1])
bs_mult = BSpline(knots_mult1,3)
B_mult = full(getBasis(collect(y),bs_mult))
c_mult = B_mult\f(collect(y))
approx_g_mult = vec(B_mult*c_mult)
p = plot(y,f(collect(y)),layout=3)
plot!(p[2],approx_g_unif)
plot!(p[3],approx_g_mult)
end

  # function to run all questions
function runall()
println("running all questions of HW-funcapprox:")
q1(15)
q2(15)
q3()
q4a()
q4b()
q5()
end
end
