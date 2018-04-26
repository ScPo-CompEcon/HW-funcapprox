
module funcapp

	using Plots
	using FastGaussQuadrature
	using ApproxFun
	using Base.Test
	using ApproXD

	function f(x)
		return(x + 2 * x^2 - exp(-x))
	end

	# use chebyshev to interpolate this:
	function q1(n)
		#deg = n-1
		a = -3
		b = +3
		n = 15
		nodes = gausschebyshev(n)[1]
		z = 0.5 * (a+b) + 0.5 * (b-a) * nodes
		values = f.(z)
		Phi = [cos((n-i+0.5) * (j-1) * π/n) for i = 1:n, j = 1:n]
		c = Phi\values
		function prediction(x)
			node = 2 * (x - a)/(b - a) -1
			Phi = [cos(acos(node)*j) for j = 0:n-1]
			return(transpose(c)*Phi)
		end
		new_n = 50
		new_nodes = [a + (i-1)/(new_n -1)*(b-a) for i in 1:new_n]
		y = f.(new_nodes)
		yhat = prediction.(new_nodes)
		global e = y - yhat
		global q1_plot = plot(new_nodes, e, title="Q1 - Deviation in approximation from true f")
	end

	function q2(n)
		nodes = Chebyshev(-3..3)
		grid = points(nodes,n)
		values = f.(grid)
		predictions = Fun(nodes,ApproxFun.transform(nodes,values))
		new_n =100
		x = linspace(-3,3,new_n)
		yhat = predictions.(x)
		y = f.(x)
		e = y - yhat
		global q2_plot = plot(x,e,title="Q2 - Deviation in approximation from true f using ApproxFun")
	end

	ChebyT(x,deg) = cos(acos(x)*deg)
	unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]
	# plot the first 9 basis Chebyshev Polynomial Basisi Fnctions
	function q3()
		plot_total = Dict()
		x = [-1 + (i-1)/(100 -1)*2 for i in 1:100]
		for deg in 0:8
			y = ChebyT.(x, deg)
			plot_total[deg] = plot(x, y, ylim = [-1.0, 1.0],labels = ["degree = $deg"])
		end
		global q3_plot = plot(title = "Q3 - Chebyshev  basis functions", plot_total[0], plot_total[1], plot_total[2], plot_total[3], plot_total[4], plot_total[5], plot_total[6], plot_total[7],  plot_total[8], layout = (3,3))
	end


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

	function runge(x)
		return(1/(1+25*x^2))
	end

	function q4a(deg=(5,9,15),lb=-5.0,ub=5.0)
		deg=(5,9,15)
		lb=-5.0
		ub=5.0
		x = [lb + (i-1)/(50 -1)*(ub-lb) for i in 1:50]
		unif_y = Array[]
		for j in 1:3
			n = 1 + deg[j]
			unif_nodes = [lb + (i-1)/(n -1)*(ub-lb) for i in 1:n]
			unif_cheby = ChebyType(unif_nodes, n, lb, ub, runge)
			push!(unif_y, predict.preds(unif_cheby, x))
		end
		for el in deg
			n = el + 1
			cheby_nodes = gausschebyshev(n)
			cheby_cheby = ChebyType(cheby_nodes, n, lb, ub, runge)
			cheby_y[el] = predict(cheby_nodes,x)
		end

	end
	#
	# function q4b()
	#
	#
	# end
	#
	# function q5()
	#
	#
	# end


		# function to run all questions
	function runall()
		println("running all questions of HW-funcapprox:")
		q1(15)
		display(q1_plot)
		# @test maximum(e) < 1e-9
		q2(15)
		display(q2_plot)
		q3()
		display(q3_plot)
		# q4a()
		# q4b()
		# q5()
		end



end
