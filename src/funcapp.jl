module funcapp

	using FastGaussQuadrature
	using Plots
	using LaTeXStrings
	using ApproxFun
	using ApproXD
	using Base.Test

	"""
		Approximate `f(k)	= k + 2k^2 - exp(-k)` using Chebysev approximation.

		#### Fields
		- `n::Integer`: Number of Chebyshev nodes

		#### Returns
		None if the test failed^[1]. Otherwise, plots two labels
		 	- the true and the approximated function on a `linspace(-3,3,100)`,
			- the error of the approximation on the same grid.

		^[1]: test whether the maximal deviation between the true and function and the
		approxiamtion is smaller than `1e-9`
	"""
	function q1(n)
		f(k)	= k + 2k^2 - exp(-k)
		a,b 	= -3.,3.
		deg		= n-1
		# get the chebyshev nodes & adjust to the appropriate interval
		z			= gausschebyshev(n)[1]
		x			= 0.5 .* (a + b) .+ 0.5 .* (b-a) .* z
		y			= f.(x)

		# checking the position of the nodes
		Plots.plot(linspace(-3,3,100), f, line = 3, xlab = "x", label = L"f(x)", legendfont = font(12),
							grid = false)
		Plots.scatter!(x, f, label = "Nodes", legendfont = font(12), markersize=10)

		# build the Chebyshev basis matrix
		V			= Array{Float64}(n,deg+1)
		for d in 0:deg
			V[:, d+1] = cos(acos.(z) * d)
		end
		c			= V \ y

		# construct the basis matrix off the chebyshev grid
		x2		= linspace(-3.,3.,100)
		z2		= 2(x2 - a)/(b - a) - 1
		V2		= Array{Float64}(length(x2),deg+1)
		for d in 0:deg
			V2[:, d+1] = cos(acos.(z2) * d)
		end

		# approximate the function
		y_app	= V2 * c
		y2		= f.(x2)
		err		= y2 .- y_app

		@test maxabs(err) < 1e-9

		plot1 = Plots.plot(x2, [y2, y_app], line = 1, label = [L"f(x)" L"\hat{f}(x)"],
						xlab = L"x", title = "Function Approximation", legendfont = font(12))
		plot2	= Plots.plot(x2, err, line = 3, label = L"$f(x) - \hat{f}(x)$", xlab = L"x",
						title = "Approximation error", legendfont = font(12))

		return Plots.plot(plot1, plot2, layout = 2)

	end

	"""
		Approximate `f(k)	= k + 2k^2 - exp(-k)` using Chebysev approximation, using
		the ApproxFun package.

		#### Fields
		- `n::Integer`: Number of Chebyshev nodes

		#### Returns
		None if the test failed^[1]. Otherwise, plots two labels
		 	- the true and the approximated function on a `linspace(-3,3,100)`,
			- the error of the approximation on the same grid.

		^[1]: test whether the maximal deviation between the true and function and the
		approxiamtion is smaller than `1e-9`
	"""
	function q2(n)
		# Note that the ApproxFun package use by default the Chebyshev space.
		f(k)	= k + 2k^2 - exp(-k)
		S 		= Chebyshev(-3..3)
		x			= points(S, n)
		v			= f.(x)
		ff		= Fun(S, ApproxFun.transform(S,v))

		x2 		= linspace(-3.,3.,100)
		y_app = ff.(x2)
		y		 	= f.(x2)
		err 	= y .- y_app

		@test maxabs(err) < 1e-9

		plot1 = Plots.plot(x2, [y, y_app], line = 1, label = [L"f(x)" L"\hat{f}(x)"],
						xlab = L"x", title = "Function Approximation", legendfont = font(12))
		plot2	= Plots.plot(x2, err, line = 3, label = L"$f(x) - \hat{f}(x)$", xlab = L"x",
						title = "Approximation error", legendfont = font(12))

		return Plots.plot(plot1, plot2, layout = 2)
	end

	"""
	Function plotting the first 9 basis Chebyshev Polynomial basis functions

	#### Fields

	None.

	#### Returns

	9 panels plot of the first 9 basis Chebyshev Polynomial basis functions.
	"""
	function q3()
		n			= 100
		x			= linspace(-1,1,n)
		V			= Array{Float64}(n,9)
		titles = Array{String}(9)
		for d in 0:8
			V[:, d+1] = cos(acos.(x) * d)
			titles[d+1] = "Degree $d"
		end
		Plots.plot(x, V, title = titles', line = 3, layout = 9)
	end

	##NOTE: thanks Florian for ChebyType and predict(), saved me one hour ;)
	ChebyT(x,deg) = cos(acos(x)*deg)
	unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]

	type ChebyType
		f::Function # function to approximate
		nodes::Union{Vector,LinSpace} # evaluation points
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

		true_new 		= Ch.f.(x_new)
		basis_new 	= Float64[ChebyT(unitmap(x_new[i],Ch.lb,Ch.ub),j) for i=1:length(x_new),j=0:Ch.deg]
		basis_nodes = Float64[ChebyT(unitmap(Ch.nodes[i],Ch.lb,Ch.ub),j) for i=1:length(Ch.nodes),j=0:Ch.deg]
		preds 			= basis_new * Ch.coefs
		preds_nodes = basis_nodes * Ch.coefs

		return Dict("x"=> x_new,"truth"=>true_new, "preds"=>preds, "preds_nodes" => preds_nodes)
	end

	"""
	Function approxiamting the Runge's function via Chebyshev polynomials, with different degrees and
	different node placement.

	#### Fields

	- `deg::Tuple`: Different degree to be used in the splines.
	- `a::Float64`: Lower bound of the function
	- `b::Float64`: Upper bound of the function

	#### Returns

	2 panels plot: the first show the approximation with uniformly spaced grid,
	the second with Chebyshev nodes -- for all the degree passed in parameters.
	"""
	function q4a(deg=(5,9,15),a=-5.0,b=5.0)
		f(x) 	= 1 / (1 + 25 * x * x)
		x_n 	= linspace(a, b, 100)

		# use the predict() function to approximate f(x) at different degree d
		# can use either chebyshev nodes or a uniformly spaced grid
		function appro(d::Integer, cheby::Bool)
			n			= d + 1

			if cheby
				nod 	= gausschebyshev(n)[1]
				nod		= 0.5 .* (a + b) .+ 0.5 .* (b-a) .* nod
				che 	= ChebyType(nod, d, a, b, f)
				pre 	= predict(che, x_n)
			else
				nod		= linspace(a,b,n)
				che 	= ChebyType(nod, d, a, b, f)
				pre 	= predict(che, x_n)
			end

			return pre
		end

		app_g 		= Array{Vector}(4)
		app_c 		= Array{Vector}(4)
		labels 		= Array{String}(4)
		# True function
		app_g[1] 	= app_c[1] = f.(x_n)
		labels[1] = L"f(x)"
		# Aproximations
		for d in 1:length(deg)
			app_c[d+1] = appro(deg[d], true)["preds"]
			app_g[d+1] = appro(deg[d], false)["preds"]
			labels[d+1] = "k = $(deg[d])"
		end

		plot1 = Plots.plot(x_n, app_g, label = labels', line = 2, title = "Uniformly spaced grid", xlab = L"x")
		plot2 = Plots.plot(x_n, app_c, label = labels', line = 2, title = "Chebyshev interpolation points", xlab = L"x")

		return Plots.plot(plot1, plot2, layout = 2)
	end

	"""
	Function approxiamting the Runge's function via splines, with different node placement, using ApproXD.jl

	#### Fields

	None.

	#### Returns

	2 panels plot: the first show the Runge's function, the second the approximation error for the
	two set of nodes.
	"""
	function q4b()
		f(x) 	= 1 / (1 + 25 * x * x)
		a, b 	= -5., 5
		x_n		= linspace(a, b, 200)

		# build our own knots, concentrated toward 0
		k				= linspace(-pi/2 + .2,pi/2 - .2,13)
		own_k		= tan.(k)
		own_k		= (b - a) .* (own_k .- minimum(own_k)) ./ (maximum(own_k) .- minimum(own_k)) .+ a
		own_k		= sort(own_k)

		# Approximate the function with BSpline
		function appro(my_knot::Bool)
			n 		= 65
			x 		= linspace(a, b, n)				# grid
			y			= f.(x)

			println(my_knot)

			# NOTE: nothing was specified regarding the degree of the spline, so picked 5
			if my_knot
				# pick knots concentrated around 0
				bs	= ApproXD.BSpline(own_k, 5)
			else
				# ApproXD pickes the knots automatically
				bs 	= ApproXD.BSpline(13, 5, a, b)
			end
			d 		= full(getBasis(collect(x), bs))				# Basis function, evaluated at x
			c			= d \ y																	# Approximate coefficients

			# use the finer grid to evaluate the approximated function
			d_1 	= full(getBasis(collect(x_n), bs))			# Basis function, evaluated at x_n
			y_r	 	= d_1 * c

			return y_r
		end

		y_n			= f.(x_n)
		y_app 	= Array{Vector}(2)
		y_app[1]= appro(false)
		y_app[2]= appro(true)
		err 		= Array{Vector}(2)
		err[1] 	= y_n .- y_app[1]
		err[2] 	= y_n .- y_app[2]

		plot1 = Plots.plot(x_n, y_n, line = 2, xlab = L"x", title = "Runge's function")
		plot2 = Plots.plot(x_n, err, line = 2, xlab = L"x", title = "Approximation error",
						label = ["Version 1" "Version 2"])
		Plots.scatter!(own_k, zeros(13), markersize = 6, label = "Own knots")

		return plot(plot1, plot2, layout = 2)

	end

	function q5()


	end


		# function to run all questions
	function runall()
		println("running all questions of HW-funcapprox:")
		q1(15)
		q2(15)
		q3()
		q4a()
		q4b()
		# q5()
	end


end
