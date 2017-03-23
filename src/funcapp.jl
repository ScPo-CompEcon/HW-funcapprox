module funcapp

using FastGaussQuadrature
using Plots
import ApproXD: getBasis, BSpline
using Distributions
using ApproxFun
using CompEcon
using Base.Test

	ChebyT(x,deg) = cos(acos(x)*deg)
	unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]

	function q1(n)
	ChebyT(x,deg) = cos(acos(x)*deg)
	unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]
	f(x) = x .+ 2x.^2 - exp(-x)
	 deg = n-1
	 lb,ub = (-3.0,3.0)
	 z = gausschebyshev(n)
	 x_nodes = Float64[(z[1][i]+1)*(ub-lb) / 2 + lb for i=1:n]  # map z[-1,1] to x[a,b]
	 y = f(x_nodes)
	 Phi = Float64[cos((n-i+0.5)*(j-1)*pi/n) for i=1:n,j=1:(deg+1)]
	 c = Phi \ y
	 n_new = 100
	 x_new = linspace(lb,ub,n_new)
	 # evaulate Cheby basis at new poitns
	 Phi_x_new = Float64[ChebyT(unitmap(x_new[i],lb,ub),j) for i=1:n_new,j=0:deg]
	 y_new = Phi_x_new * c
	 y_true = f(x_new)
	 err = y_new - y_true
	 # plot
	 p = Any[]
	 true_approx = plot(1:n_new, y_true, label="True f(x) = x+2x^2-exp(-x)",
	 title="True vs Apprx_Cheb")
	 plot!(y_new, label="Approximated with Cheby")
	 push!(p,true_approx )
	 push!(p, plot(1:n_new, err,
	 label="Deviation from true",
	 title="Chebyshev error"))
	 @test_approx_eq_eps maximum(err) 0 1e-9
	 p1 = plot(p...)
	 display(p1)
	end

	function q2(n)
		f(x) = x .+ 2x.^2 - exp(-x)
		lb,ub = (-3.0,3.0)
		S = Chebyshev(lb..ub)
		x1 = points(S,n)
		#x1 = linspace(lb,ub,n)
		v = [f(k) for k in x1]   # values at the non-default grid
		V = Array(Float64,n,n) # Create a Vandermonde matrix by evaluating the basis at the grid
		for k = 1:n
    	V[:,k] = Fun(S,[zeros(k-1);1]).(x1)
		end
		ff = Fun(S,V\v)

		n_new = 100
		x_new = linspace(lb,ub,n_new);   # a non-default grid
		f_app = [ff(xx) for xx in x_new]
		y_true = f(x_new)
		err = f_app - y_true
		q = Any[]
		true_approxfun = plot(1:n_new, y_true, label="True f(x)", title="True vs ApproxFun")
		plot!(f_app, label="ApproxFun")
		push!(q, true_approxfun)
		push!(q, plot(1:n_new, err,
		label="Deviation from true",
		title="Chebyshev error"))
		q1=plot(q...)
		display(q1)
	end

	# plot the first 9 basis Chebyshev Polynomial Basisi Fnctions
	function q3(n)
		m=9
		Phi = Float64[cos((n-i+0.5)*(j-1)*pi/n) for i=1:n,j=1:m]
		t = Any[]
		push!(t, plot(Phi[:, 1], label = "Basis 1", ylim=(-1.1, 1.1)))
		push!(t, plot(Phi[:, 2], label = "Basis 2", ylim=(-1.1, 1.1)))
		push!(t, plot(Phi[:, 3], label = "Basis 3", ylim=(-1.1, 1.1)))
		push!(t, plot(Phi[:, 4], label = "Basis 4", ylim=(-1.1, 1.1)))
		push!(t, plot(Phi[:, 5], label = "Basis 5", ylim=(-1.1, 1.1)))
		push!(t, plot(Phi[:, 6], label = "Basis 6", ylim=(-1.1, 1.1)))
		push!(t, plot(Phi[:, 7], label = "Basis 7", ylim=(-1.1, 1.1)))
		push!(t, plot(Phi[:, 8], label = "Basis 8", ylim=(-1.1, 1.1)))
		push!(t, plot(Phi[:, 9], label = "Basis 9", ylim=(-1.1, 1.1)))
		t1 = plot(t...)
		display(t1)
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

	function q4a(deg,lb,ub)
		runge(x) = 1.0./(1+25.*x.^2)
		@assert length(deg)==3
		xnew = linspace(lb,ub,1000)
		colors   = ["orange","green","purple"]
		w1 = plot(xnew, runge(xnew),color="red", lw=1.5, label="runge",
		title="Uniform nodes")
		for j in 1:length(deg)
 			nodes = linspace(lb,ub,deg[j]+1)
 			C = ChebyType(nodes,deg[j],lb,ub,runge)
 			p = predict(C,xnew)
			plot!(xnew, p["preds"],label="deg=$(deg[j])",color=colors[j],lw=1.5)
 		end
		w2 = plot(xnew, runge(xnew),color="red", lw=1.5, label="runge",
		title="Cheby Nodes")
		for j in 1:length(deg)
			n = deg[j] + 1
 			z = gausschebyshev(n)
 			nodes = Float64[(z[1][i]+1)*(ub-lb) / 2 + lb for i=1:n]
 			C = ChebyType(nodes,deg[j],lb,ub,runge)
 			p = predict(C,xnew)
			plot!(xnew, p["preds"],label="deg=$(deg[j])",color=colors[j],lw=1.5)
 		end
		r = Any[]
		push!(r, w1)
		push!(r, w2)
		r1 = plot(r...)
		display(r1)
	end


	function q4b()
		runge(x) = 1.0./(1+25.*x.^2)
		deg = 3 #cubic spline
		ub,lb = (5,-5)
		#Equidistant knots
		nknots = 13
		b1 = BSpline(nknots,deg,lb,ub)
		nevals = 5 * b1.numKnots #65

		#Knots concentrated around zero
		 G(k,s) = GeneralizedPareto(k,s,0)
		 pf(k,s) = quantile(GeneralizedPareto(k,s,0),linspace(0.05,cdf(G(0.5,1),5),6))
		 my_knots = vcat(-reverse(pf(0.5,1)),0.0,pf(0.5,1))
		 b2 = BSpline(my_knots,deg)

		#Coefficients, get our approximations
		points = collect(linspace(lb,ub,nevals))
  	c1 = getBasis(points,b1) \ runge(points)
  	c2 = getBasis(points,b2) \ runge(points)
		#Errors
		xnew = collect(linspace(lb,ub,1000))
 		runge_true = runge(xnew)
 		err1 = getBasis(xnew,b1) * c1 - runge_true;
 		err2 = getBasis(xnew,b2) * c2 - runge_true;

		#Plots
		s1 = plot(xnew, runge(xnew),color="red", lw=1.5, label="runge",
		title="Runge's function")
		s2 = plot(xnew, err1, color="blue", lw=1.5, label="approx1")
		title!("Error in Runge's function")
		plot!(xnew, err2, color="yellow", lw=1.5, label="approx2")
		#Show knots
		plot!(unique(b1.knots),zeros(nknots),color="pink", marker=(5,:ellipse),
		label = "Equidistant knots")
	 	plot!(my_knots,zeros(my_knots),color="green", marker=(5,:xcross),
		label = "Concentrated knots")
		g = Any[]
		push!(g, s1)
		push!(g, s2)
		g1 = plot(g...)
		display(g1)
	end

	function q5()
		f(x) = abs(x).^0.5
	 	lb,ub = (-1.0,1.0)
		deg = 3
	 	nknots = 13
		#Trying with ApproXD
		#b1 = ApproXD.BSpline(nknots,3,lb,ub)
		b1 = SplineParams(linspace(lb,ub,nknots),0,deg)
		nevals = 65
		my_knots = vcat(linspace(-1,-0.1,5),0,0,0, linspace(0.1,1,5))
	 	b2 = SplineParams(my_knots,0,deg)
		#Coeffs
		points = collect(linspace(lb,ub,nevals))
		#c1 = ApproXD.getBasis(points, b1) \ f(points)
		c1 = CompEcon.evalbase(b1,points) \ f(points)
		c2 = CompEcon.evalbase(b2,points) \ f(points)
		#Errors
		xnew = collect(linspace(lb,ub,1000))
 		true_f = f(xnew)
		#p1 = ApproXD.getBasis(xnew, b1) * c1
		p1 = CompEcon.evalbase(b1,xnew) * c1
		p2 = CompEcon.evalbase(b2,xnew) * c2
		err1 = p1 - true_f
 		err2 = p2 - true_f
		#Plots
		h1 = plot(xnew, true_f,color="black", lw=1.5, label="abs(x)^0.5",
		title="True function")
		h2 = plot(xnew, p1, color="green", lw=1.5,
		label="Uniform knot vector",
		legend= (:top),
		title="Two Spline approximations")
		plot!(xnew, p2, color="purple", lw=1,
		label="Knot vector with knot multiplicity at 0")
		h3 = plot(xnew, err1, color="green", lw=1.5,
		label="Uniform knot vector")
		title!("Errors")
		plot!(xnew, err2, color="purple", lw=1.5,
		label="Knot vector with knot multiplicity at 0")
		k = Any[]
		push!(k, h1)
		push!(k, h2)
		push!(k, h3)
		k1 = plot(k...)
		display(k1)
	end

	function runall()
		println("running all questions of HW-funcapprox:")
		println("Question 1")
		q1(15)
		println("Question 2")
		q2(15)
		println("Question 3")
		q3(15)
		println("Question 4a:")
		q4a((5,9,15),-1,1)
		println("Regarding q4a plots: Runge's phenomenon is a problem of oscillation at the edges of an interval that occurs when using polynomial interpolation with polynomials of high degree over a set of equispaced interpolation points.")
		println("Going to higher degrees does not always improve accuracy.")
		println("Question 4b")
		q4b()
		println("Couldn't handle knot multiplicity with ApproXD, and other issues too")
		println("Used CompEcon instead for question 5")
		q5()
		println("Setting more than one knot equal to zero made the plots disappear")
		println("With CompEcon, 3 zero knots do the job.")
	end

end
