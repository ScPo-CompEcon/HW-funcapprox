

module funcapp


using FastGaussQuadrature
using ApproxFun
using PyPlot
using DataFrames
#Pkg.clone("https://github.com/floswald/ApproXD.jl")
using ApproXD

	# use chebyshev to interpolate this:
	function q1(n)
	   # Define the function
	   f(x) = x + 2x^2 - exp(-x)
	   # Specify degree of approximation
	   deg = n-1
	   # Define bounds of domain
	   ub = 3.0
	   lb = -3.0
	   # Get interpolation nodes
	   nodes = gausschebyshev(n)[1]
	   tnodes = (nodes+1)*(ub-lb)/2 + lb
	   # Get function value at x, y=f(x)
	   ytnodes = map(f, tnodes)
	   # Evaluate the chebyshev basis matrix at Z
	   basmat = [cos.(acos.(nodes[i])*j) for i=1:n,j=0:(deg)]
	   #Compute approx coeffs c by matrix inversion
	   coef = basmat \ ytnodes

	   # Evaluation at new points

	   n_new = 100
	   x_new = linspace(lb,ub,n_new)
	   ## Evaluate the basis matrix at x_new
	   ## We need to map x_new to -1,1 to find the base matrix
	   x_new_t = 2(x_new-lb)/(ub-lb) -1
	   basmat_new = [cos.(acos.(x_new_t[i]) * j) for i=1:n_new,j=0:(deg)]
	   y_new = basmat_new * coef
	   ytrue = map(f, x_new)
	   error = y_new - ytrue
	   figure("Question 1", figsize=(10, 8))
	   subplot(121)
	   plot(x_new,ytrue,color="blue",label="True Values")
	   scatter(x_new,y_new,label="Approximation")
	   title("Plain Vanilla Chebyshev approximation")
	   subplot(122)
	   plot(x_new, error, color="red", label="Approximation error")
        return(DataFrame(error = abs.(error)))
end

	function q2(n)
	     ff = x -> x + 2*x.^2 - exp.(-x)
	     F = Fun(ff,Chebyshev(Interval(-3.0,3.0)))
	     x = collect(linspace(-3.0, 3.0, n))
	     y = []
	     for i in 1:n
	      push!(y, F(x[i]))
	     end
	     error = ff(x) - y

	     figure("Question 2", figsize = (10,8))
	     subplot(121)
	     scatter(x, y, label = "Function Approximation",  s =  50, c = "red", marker = "+", alpha=0.5)
	     plot(x, ff(x), label = "True Function")
	     legend()
	     subplot(122)
	     plot(x, error)
	     title("Approximation Error")

	end

	ChebyT(x,deg) = cos.(acos.(x)*deg)

	# plot the first 9 basis Chebyshev Polynomial Basisi Fnctions
	function q3()
		x = collect(linspace(0,1.0,500))

		data = [x]
		for i in 0:8
		     data = hcat(data, [ChebyT(x, i)])
		end

		figure("Question 3", figsize = (8, 10))
		for i in 331:339
		     subplot(i)
		     plot(data[1], data[i-329])
		     title("Polynomial of order $(i-331)", fontsize = 10)
		     xticks(size = 6)
		     yticks(size = 6)
		end

	end


	unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]

	type ChebyType
		f::Function # fuction to approximate
		nodes::Union{Vector,LinSpace} # evaluation points
		basis::Matrix # basis evaluated at nodes
		coefs::Vector # estimated coefficients

		deg::Int 	# degree of chebypolynomial
		lb::Float64 # bounds
		ub::Float64

		# constructor
		function ChebyType(_nodes::Union{Vector,LinSpace},_deg,_lb,_ub,_f::Function)
			n = length(_nodes)
			y = _f(_nodes)
			_basis = Float64[ChebyT(unitmap(_nodes[i],_lb,_ub),j) for i=1:n,j=0:_deg]
			_coefs = _basis \ y  # type `?\` to find out more about the backslash operator. depending the args given, it performs a different operation
			# create a ChebyType with those values
			new(_f,_nodes,_basis,_coefs,_deg,_lb,_ub)
		end
	end

	# function to predict points using info stored in ChebyType
	function predict(Ch::ChebyType,x_new)

		true_new = Ch.f(x_new)
		basis_new = Float64[ChebyT(unitmap(x_new[i],Ch.lb,Ch.ub),j) for i=1:length(x_new),j=0:Ch.deg]
		basis_nodes = Float64[ChebyT(unitmap(Ch.nodes[i],Ch.lb,Ch.ub),j) for i=1:length(Ch.nodes),j=0:Ch.deg]
		preds = basis_new * Ch.coefs
		preds_nodes = basis_nodes * Ch.coefs

		return Dict("x"=> x_new,"truth"=>true_new, "preds"=>preds, "preds_nodes" => preds_nodes)
	end

	function q4a(deg=(5,9,15),lb=-1.0,ub=1.0)

        h(x) = 1 ./ (1+ 25 .* x .^2)
        x_new = linspace(lb, ub, 100)
        y_new = h(x_new)
        #Chebyshev nodes

        nodes = gausschebyshev(deg[1]+1)[1]
        d5 = predict(ChebyType(nodes, deg[1], lb, ub, h),
                            x_new)
        nodes = gausschebyshev(deg[2]+1)[1]
        d9 = predict(ChebyType(nodes, deg[2], lb, ub, h),
                            x_new)
        nodes = gausschebyshev(deg[3]+1)[1]
        d16 = predict(ChebyType(nodes, deg[3], lb, ub, h),
                            x_new)

		# Uniform nodes

	    grid = collect(linspace(lb, ub, (deg[1]+1)))
	    u5 = predict(ChebyType(grid, deg[1], lb, ub, h),
	   					 x_new)
	    grid = collect(linspace(lb, ub, (deg[2]+1)) )
	    u9 = predict(ChebyType(grid, deg[2], lb, ub, h),
	   					 x_new)
	    grid = collect(linspace(lb, ub, (deg[3]+1)))
	    u16 = predict(ChebyType(grid, deg[3], lb, ub, h),
	   					 x_new)

		#Plots

		figure("Question 4-a", figsize=(10, 8))
		subplot(121)
        plot(x_new, y_new, label="True function")
        title("Function Approximation\nwith Chebyshev nodes")
        plot(d5[:"x"], d5[:"preds"], label = "k=5")
        plot(d9[:"x"], d9[:"preds"], label = "k=9")
        plot(d16[:"x"], d16[:"preds"], label = "k=16")
        subplot(122)
        plot(x_new, y_new, label="True function")
        title("Function Approximation\nwith Uniform nodes")
        plot(u5[:"x"], u5[:"preds"], label = "k=5")
        plot(u9[:"x"], u9[:"preds"], label = "k=9")
        plot(u16[:"x"], u16[:"preds"], label = "k=16")

	end

	function q4b()
	 #Spaced knots
	 nknots = 13
	 deg = 5
	 lb = -5.0
	 ub = 5.0
	 h(x) = 1 ./ (1 + 25 .* x .^2)
	 b1 = BSpline(nknots,deg,lb,ub)
	 nevals = 5 * nknots
	 evals = collect(linspace(lb, ub, nevals))
	 mat1 = full(getBasis(evals,b1))
	 coef1 = mat1 \ h(evals)
	 y1 = mat1 * coef1
	 error1 = h(evals) - y1

	 #Knots around zero

	 t = gausschebyshev(12)[1]
	 tt = (t+1)*(ub-lb)/2 + lb
	 t1 = tt[1:6] + 1
	 t2 = tt[7:12] - 1
	 T = sort(vcat(t1, 0.0, t2))
	 b2 = BSpline(T, deg)
	 mat2 = full(getBasis(evals, b2))
	 coef2 = mat2 \ h(evals)
	 y2 = mat2 * coef2
	 error2 = h(evals) - y2

	 # Plots

	 x_new = collect(linspace(lb, ub, 200))
	 ytrue = h(x_new)
	 figure("Question 4-b", figsize=(10, 8))
	 subplot(121)
	 plot(x_new, ytrue)
	 title("Runge's function")
	 subplot(122)
	 plot(evals, error1, label =  "spaced knots")
	 plot(evals, error2, label = "concentrated knots")
	 scatter(T, zeros(13), label = "knots location", s =  50, c = "green", marker = "x", alpha=0.5 )
	 legend()
	 title("Error in Runge's function")
	end

	function q5()

	 v(x) = sqrt.(abs.(x))
	 lb = -1.0
	 ub = 1.0
	 nknots = 13
      deg=3
	 nevals = 5 * nknots
	 evals = collect(linspace(lb, ub, nevals))

	 # Uniform Knots
	 b1 = BSpline(nknots,deg,lb,ub)
	 mat1 = full(getBasis(evals,b1))
	 coef1 = mat1 \ v(evals)
	 y1 = mat1 * coef1
	 error1 = v(evals) - y1

	 # Knots with multiplicity at 0

	 k = vcat(collect(linspace(-1, -0.1, 5)), 0, 0, 0, collect(linspace(0.1, 1, 5)) )
	 b2 = BSpline(k, deg)
	 mat2 = full(getBasis(evals, b2))
	 coef2 = mat2 \ v(evals)
	 y2 = mat2 * coef2
	 error2 = v(evals) - y2

	 x_new = collect(linspace(-1.0, 1.0, 200))
	 y1_new = full(getBasis(x_new, b2)) * coef2
	 y2_new = full(getBasis(x_new, b1)) * coef1
	 ytrue = v(x_new)

	 figure("Question 5", figsize = (10,10))
	 subplot(221)
	 plot(x_new, ytrue)
	 title("True function")
	 subplot(222)
	 plot(x_new, y1_new, label = "Uniform knots")
	 plot(x_new, y2_new, label = "Knots with multiplicity at 0")
	 legend()
	 title("Function Approximations")
	 subplot(223)
	 plot(evals, error1, label = "Uniform knots")
	 plot(evals, error2, label = "Knots with multiplicity at 0")
	 legend()
	 title("Approximation errors")

	end


		# function to run all questions
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
