

module funcapp


using FastGaussQuadrature
using ApproxFun
using Plots
using DataFrames
Pkg.clone("https://github.com/floswald/ApproXD.jl")
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
	    basmat = [cos(acos(nodes[i])*j) for i=1:n,j=0:(deg)]
	    #Compute approx coeffs c by matrix inversion
	    coef = Phi \ ytnodes

	    # Evaluation at new points

	    n_new = 100
	    x_new = linspace(lb,ub,n_new)
	    ## Evaluate the basis matrix at x_new
	    ## We need to map x_new to -1,1 to find the base matrix
	    x_new_t = 2(x_new-lb)/(ub-lb) -1
	    basmat_new = [cos(acos(x_new_t[i]) * j) for i=1:n_new,j=0:(deg)]
	    y_new = basmat_new * coef
	    ytrue = map(f, x_new)
	    error = y_new - ytrue
	    D = Dict()
	    D[:approxerror] = error
	    plot1 = plot(x_new,ytrue,color="blue",label="True Values")
	    scatter!(x_new,y_new,label="Approximation")
	    plot2 = plot(x_new, error, color="red", label="Approximation error")
	    return plot(plot1, plot2)

	end

	function q2(n)
		x = Fun(Interval(-3.0,3.0))
		f(x) = x + 2x^2 - exp(-x)
		plot(f)
	end


	# plot the first 9 basis Chebyshev Polynomial Basisi Fnctions
	function q3()

		cheb1(x) = 1
		cheb2(x) = x
		cheb3(x) = 2*x*cheb2(x)-cheb1(x)
		cheb4(x) = 2*x*cheb3(x)-cheb2(x)
		cheb5(x) = 2*x*cheb4(x)-cheb3(x)
		cheb6(x) = 2*x*cheb5(x)-cheb4(x)
		cheb7(x) = 2*x*cheb6(x)-cheb5(x)
		cheb8(x) = 2*x*cheb7(x)-cheb6(x)
		cheb9(x) = 2*x*cheb8(x)-cheb7(x)
		p1 = plot(cheb1)
		p2 = plot(cheb2)
		p3 = plot(cheb3)
		p4 = plot(cheb4)
		p5 = plot(cheb5)
		p6  = plot(cheb6)
		p7 = plot(cheb7)
		p8 = plot(cheb8)
		p9 = plot(cheb9)
		plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, layout=9, label="")
	end

	ChebyT(x,deg) = cos(acos(x)*deg)
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


	end

	function q4b()


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
		q5()
	end


end
