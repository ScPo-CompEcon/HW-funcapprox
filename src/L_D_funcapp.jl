module funcapp

using Plots
using FastGaussQuadrature
using ApproXD
using Distributions
using ApproxFun

# use chebyshev to interpolate this:
function q1(n)
  lb = -3
  ub = 3
  deg = n-1

  nodes = gausschebyshev(n)[1]
  x_nodes = 0.5 * (lb + ub) + 0.5 * (ub - lb) * z
  values = [f(x_nodes[i]) for i in 1:length(x_nodes)]

  chubbychef = [cos((n - i + 0.5) * (j - 1)* π / n) for i = 1:n, j = 1:n]
  c = \(chubbychef, values)
  function fhat(x)
    z = 2 * (x - lb)/(ub - lb) -1
    chubbychef = [cos(acos(nodes)*j) for j = 0:n-1]
    return(transpose(c)*chubbychef)
  end

  n_new = 100
  new_x = [a + (i-1)/(n_new -1)*(ub-lb) for i in 1:n_new] #new equally spaced points
  predict = [fhat(new_x[i]) for i in 1:length(new_x)] #approximation
  new_val = [f(new_x[i]) for i in 1:length(new_x)] #true value
  error = new_values-prediction #error between prediction and true value
  return Dict(:err => 1.0)
end

function q2(n)
  f_afun = Fun(f_afun, Interval(-3,3))
  n_afun = 100
  x_afun = linspace(-3, 3, n_afun)
  true_y = map(f, x_afun)
  cheby_y = map(f_afun, x_afun)
  plot1 = plot(x_afun, [true_y, cheby_y], label = ["f(x)" "f_afun(x)"])
  plot2 = plot(x_afun, true_y - cheby_y, yformatter = :scientific, label = "f(x) - f_afun(x)")
  return plot1, plot2
end

# plot the first 9 basis Chebyshev Polynomial Basisi Fnctions
function q3()
    lb = -3
    ub = 3
    z = 2 * (x - lb)/(ub - lb) -1
    chubbychef = cos(acos(z)j)
	p = Any[]
	push!(p, plot(chubbychef[:, 1], label = "1", ylim=(-1, 1)))
    push!(p, plot(chubbychef[:, 2], label = "2", ylim=(-1, 1)))
	push!(p, plot(chubbychef[:, 3], label = "3", ylim=(-1, 1)))
	push!(p, plot(chubbychef[:, 4], label = "4", ylim=(-1, 1)))
	push!(p, plot(chubbychef[:, 5], label = "5", ylim=(-1, 1)))
	push!(p, plot(chubbychef[:, 6], label = "6", ylim=(-1, 1)))
	push!(p, plot(chubbychef[:, 7], label = "7", ylim=(-1, 1)))
	push!(p, plot(chubbychef[:, 8], label = "8", ylim=(-1, 1)))
	push!(p, plot(chubbychef[:, 9], label = "9", ylim=(-1, 1)))
end

function q4a()
    unodes = Array(Vector,3)
          i = 1
          for n in [5, 10, 15]
              unodes[i] = linspace(-1,1,n)
              i = i+1
          end
     chnodes = Array(Vector,3)
          i = 1
          for n in [5, 10, 15]
              chnodes[i] = gausschebyshev(n)[1]
              i = i+1
          end
    c = get_coeffs(-5.0, 5.0, unodes[1],r)
    y = Array{Vector}(3,2)
    j=1
        for nodes in [unodes, chnodes]
    i = 1
        for n in [5, 10, 15]
            y[i,j] = predict(100, n, get_coeffs(-5.0,5.0,nodes[i],r), -5.0, 5.0,r)[1]
            i = i+1
            end
        j=j+1
        end
        x_new = linspace(-5,5,100)
          true_y = r.(x_new)
          plot1 = plot(x_new, [true_y, y[1,1], y[2,1], y[3,1]],
          labels = ["Runge" "5" "9" "15"], title = "Uniform")
          plot2 = plot(x_new, [true_y, y[1,2], y[2,2], y[3,2]],
          labels = ["Runge" "5" "9" "15"], title = "Chebyshev")
          return plot1, plot2
    end

function q4b()
  b = BSpline(13, 3, -5.0, 5.0)
  Basis = full(getBasis(collect(linspace(-5,5,65)),b))
  y = r.(linspace(-5,5,65))
  c = Basis \ y
  x_new = linspace(-5,5,65)
  true_y = r.(x_new)
  Basis1 = full(getBasis(collect(x_new),b))
  Bs_y1 = Basis1*c
  knots = vcat( collect(linspace(-5,-1.75,3)), collect(linspace(-1,1,7)), collect(linspace(1.75,5,3)) )
  b2 = BSpline(knots, 3)
  Basis2 = full(getBasis(collect(x_new),b2))
  y = r.(x_new)
  c2 = Basis2 \ y
  Basis3 = full(getBasis(collect(x_new),b2))
  Bs_y2 = Basis3*c2
  plot1 = plot(x_new, [true_y, Bs_y1, Bs_y2],
    linecolor = [:auto "purple" "red"], line = [:solid :dot :dot],
    labels = ["Runge" "Uniform knots" "Clustered at zero"])
  plot2 = scatter(knots, zeros(15), labels = "knot positions")
  plot2 = plot!(x_new, [true_y - Bs_y1, true_y-Bs_y2], labels = ["Uniform knots" "knots"])
  return plot1, plot2
  	end

function q5()
    knots=13
    deg=5
    lb=-1
    ub=1
    y=linspace(a,b,65)

    bs_unif = BSpline(knots,3,lb,ub)
    B_unif = full(getBasis(collect(y),bs_unif))
    c_unif=B_unif\g(collect(y))
    approx_g_unif=vec(B_unif*c_unif)

    knots_mult1=vec([-1 -0.8 -0.6 -0.4 -0.2 -0.01 0 0.01 0.2 0.4 0.6 0.8 1])
    bs_mult= BSpline(knots_mult1,3)
    B_mult = full(getBasis(collect(y),bs_mult))
    c_mult=B_mult\g(collect(y))
    approx_g_mult=vec(B_mult*c_mult)
    p=plot(y,g(collect(y)),layout=3)
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
