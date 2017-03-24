module funcapp

using FastGaussQuadrature
using PyPlot
using Base.Test
Pkg.rm("ApproXD")
Pkg.clone("https://github.com/floswald/ApproXD.jl")
import ApproXD: getBasis, BSpline
using Distributions
using Plots
Pkg.add("ApproxFun")
using ApproxFun
using CompEcon



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


 	 Phi_x_new = Float64[ChebyT(unitmap(x_new[i],lb,ub),j) for i=1:n_new,j=0:deg]
 	 y_new = Phi_x_new * c
	 y_true = f(x_new)
	 err = y_new - y_true
	 @test_approx_eq_eps maximum(err) 0 1e-9
	 p = Any[]
	 true_approx = plot(1:n_new, y_true, title="True vs Apprx_Cheb")
	 plot!(y_new, label="Approximatd_y_Cheb")
	 push!(p,true_approx )
	 push!(p, plot(1:n_new, err, title="Chebyshev error"))
	 plot(p...)
	 savefig("Nesrine-q1.png")
   end


	function q2(n)
		f(x) = x .+ 2x.^2 - exp(-x)
		lb,ub = (-3.0,3.0)
		S = Chebyshev(lb..ub)
		x1 = points(S, n)
		n_new = 100
		v = [f(k) for k in x1];
		V = Array(Float64,n,n);
			for k = 1:n
		   V[:,k] = Fun(S,[zeros(k-1);1]).(x1)
			end
			ff = Fun(S,V\v);

			x_new = linspace(lb,ub,n_new);
			f_app = [ff(xx) for xx in x_new]
			y_true = f(x_new)
			err = f_app - y_true

			q = Any[]
			true_approxfun = plot(1:n_new, y_true, label="True function", title="True vs ApproxFun")
			plot!(f_app, label="ApproxFun")
			push!(q, true_approxfun)
 		push!(q, plot(1:n_new, err, title="Chebyshev error"))
 		plot(q...)
    savefig("Nesrine-q2.png")
	end


	function q3(n)
		m=9
		Phi = Float64[cos((n-i+0.5)*(j-1)*pi/n) for i=1:n,j=1:m]
		p = Any[]
		push!(p, plot(Phi[:, 1], label = "Basis 1", ylim=(-1.1, 1.1)))
		push!(p, plot(Phi[:, 2], label = "Basis 2", ylim=(-1.1, 1.1)))
		push!(p, plot(Phi[:, 3], label = "Basis 3", ylim=(-1.1, 1.1)))
		push!(p, plot(Phi[:, 4], label = "Basis 4", ylim=(-1.1, 1.1)))
		push!(p, plot(Phi[:, 5], label = "Basis 5", ylim=(-1.1, 1.1)))
		push!(p, plot(Phi[:, 6], label = "Basis 6", ylim=(-1.1, 1.1)))
		push!(p, plot(Phi[:, 7], label = "Basis 7", ylim=(-1.1, 1.1)))
		push!(p, plot(Phi[:, 8], label = "Basis 8", ylim=(-1.1, 1.1)))
		push!(p, plot(Phi[:, 9], label = "Basis 9", ylim=(-1.1, 1.1)))
		plot(p...)
		gui()
		savefig("Nesrine-q3.png")
	end


	ChebyT(x,deg) = cos(acos(x)*deg)
	unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]

	type ChebyType
		f::Function
		nodes::Union{Vector,LinSpace}
		basis::Matrix
		coefs::Vector

		deg::Int
		lb::Float64
		ub::Float64

		# constructor
		function ChebyType(_nodes::Union{Vector,LinSpace},_deg,_lb,_ub,_f::Function)
			n = length(_nodes)
			y = _f(_nodes)
			_basis = Float64[ChebyT(unitmap(_nodes[i],_lb,_ub),j) for i=1:n,j=0:_deg]
			_coefs = _basis \ y

			new(_f,_nodes,_basis,_coefs,_deg,_lb,_ub)
		end
	end

	function predict(Ch::ChebyType,x_new)

		true_new = Ch.f(x_new)
		basis_new = Float64[ChebyT(unitmap(x_new[i],Ch.lb,Ch.ub),j) for i=1:length(x_new),j=0:Ch.deg]
		basis_nodes = Float64[ChebyT(unitmap(Ch.nodes[i],Ch.lb,Ch.ub),j) for i=1:length(Ch.nodes),j=0:Ch.deg]
		preds = basis_new * Ch.coefs
		preds_nodes = basis_nodes * Ch.coefs

		return Dict("x"=> x_new,"truth"=>true_new, "preds"=>preds, "preds_nodes" => preds_nodes)
	end

	function q4a(deg=(5,9,15),lb=-5.0,ub=5.0)
	runge(x) = 1.0./(1+25.*x.^2)
    xnew = linspace(lb,ub,500)

    colors = ["orange","purple","green"]
     p1=plot(xnew,runge(xnew),color="black",lw=2,label="runge",title=("Uniform Nodes"))
       for j in 1:length(deg)
        nodes = linspace(lb,ub,deg[j]+1)
        C = ChebyType(nodes,deg[j],lb,ub,runge)
        p = predict(C,xnew)
        plot!(xnew,p["preds"],label="deg=$(deg[j])",color=colors[j],lw=2)
       end
     p2=plot(xnew,runge(xnew),color="black",lw=2,label="runge",title=("Cheby Nodes"))
       for j in 1:length(deg)
        n = deg[j] + 1
        z = gausschebyshev(n)
        nodes = Float64[(z[1][i]+1)*(ub-lb) / 2 + lb for i=1:n]
        C = ChebyType(nodes,deg[j],lb,ub,runge)
        p = predict(C,xnew)
        plot!(xnew,p["preds"],label="deg=$(deg[j])",color=colors[j],lw=2)
       end
    r=Any[]
		push!(r,p1)
		push!(r,p2)
		plot(r...)
    gui()
    savefig("Nesrine-q4a.png")
	end

	function q4b()
		ub,lb = (5,-5)
		runge(x) = 1.0./(1+25.*x.^2)

		nknots = 13
		deg = 3
		bs1 = BSpline(nknots,deg,lb,ub)
		nevals = 5 * bs1.numKnots

		G(k,s) = GeneralizedPareto(k,s,0)
		pf(k,s) = quantile(GeneralizedPareto(k,s,0),linspace(0.05,cdf(G(0.5,1),5),6))
		myknots = vcat(-reverse(pf(0.5,1)),0.0,pf(0.5,1))
    bs2 = BSpline(myknots,deg)

		eval_points = collect(linspace(lb,ub,nevals))
		c1 = getBasis(eval_points,bs1) \ runge(eval_points)
		c2 = getBasis(eval_points,bs2) \ runge(eval_points)
		test_points = collect(linspace(lb,ub,1000));
		truth = runge(test_points);
		e1 = getBasis(test_points,bs1) * c1 - truth;
		e2 = getBasis(test_points,bs2) * c2 - truth;
		s1 = plot(test_points, runge(test_points),color="red", lw=2, label="runge",
		title="Runge's function")
		s2 = plot(test_points, e1, color="blue", lw=2, label="approx1")
		title!("Error in Runge's function")
		plot!(test_points, e2, color="yellow", lw=2, label="approx2")
		#Show knots
		plot!(unique(bs1.knots),zeros(nknots),color="pink", marker=(5,:ellipse),
		label = "Equidistant knots")
	 	plot!(my_knots,zeros(my_knots),color="green", marker=(5,:xcross),
		label = "Concentrated knots")
		g = Any[]
		push!(g, s1)
		push!(g, s2)
		plot(g...)
    savefig("Nesrine-q4b.png")
	end



	function q5()
		f(x) = abs(x).^0.5
		lb,ub = (-1.0,1.0)
		nknots = 13
		deg = 3
		b1 = SplineParams(linspace(lb,ub,nknots),0,deg)
		nevals = 75

		myknots = vcat(linspace(-1,-0.1,5),0,0,0,    linspace(0.1,1,5))
		b2 = SplineParams(myknots,0,deg)


		eval_points = collect(linspace(lb,ub,nevals))
		c1 = CompEcon.evalbase(b1,eval_points) \ f(eval_points)
		c2 = CompEcon.evalbase(b2,eval_points) \ f(eval_points)


		test_points = collect(linspace(lb,ub,1000));
 truth = f(test_points);
 p1 = CompEcon.evalbase(b1,test_points) * c1;
 p2 = CompEcon.evalbase(b2,test_points) * c2;
 e1 = p1 - truth;
 e2 = p2 - truth;


l1=plot(test_points, f(test_points), color="black", lw=2, label="Fun")
l2=plot(test_points, p1, color="red", lw=2, label="uniform vector", title="Two Spline Approx")
plot!(test_points, p2, color="green", lw=2, label="knot multiplicity at 0")
l3=plot(test_points,e1,color="blue",label="uniform vector", title="Errors")
plot!(test_points,e2,color="yellow", label="knot multiplicity at 0")

p=Any[]
push!(g,l1)
push!(g,l2)
push!(g,l3)
plot(g...)
savefig("Nesrine-q5.png")
end

		# function to run all questions
	function runall()
		println("running all questions of HW-funcapprox:")
		display(q1(15))
		display(q2(15))
		dispaly(q3(15))
		display(q4a(deg=(5,9,15),lb=-5.0,ub=5.0)))
		display(q4b())
		display(q5())
	end


end
