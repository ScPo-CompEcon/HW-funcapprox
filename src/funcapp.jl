module funcapp

	using Plots
	using FastGaussQuadrature
	using ApproxFun
	using ApproXD

	f(x) = x .+ 2x.^2 - exp(-x)

	# use chebyshev to interpolate this:
	function q1(n)

		deg = n - 1
		a = -3
		b = +3
		z = gausschebyshev(n)[1]
		x = 0.5 * (a + b) + 0.5 * (b - a) * z
		y = f(x[i])
		cheby = [cos((n - i + 0.5) * (j - 1)* π / n) for i = 1:n, j = 1:n]
		c = cheby \ y
		fcheby(x, deg) = [cos(acos(x) * deg)]
		ztransf(x,a,b) = 2 * (x - a)/(b - a) - 1
		n2 = 100
		x2 = linspace(a,b,n2)
		cheby2 = [fcheby(ztransf(x2[i],a,b),j) for i in 1:n2, j = 0:n2 - 1]
		y2 = cheby2 * c
		y_control = f(x2[i])
		error = y2 - y_control
		plot1 = plot(x2, [y2, y_control], labels = ["approximation", "truth"], ylim = [-5,5], title = "Plain Vanilla")
		plot2 = plot(x2, error, labels = ["error"])
		plot_total = plot(plot1, plot2)
		return Dict(:err => 1.0)

	end

	function q2(n)

		a = -3
		b = 3
		n_new = 100
		estimations = Fun(f,Chebyshev(Interval(a,b)))
		x = linspace(-3,3,n_new)
		y_control = estimations.(x)
		y = f.(x)
		error = y - y_control
		plot(x,error,title="Deviation using ApproxFun")

	end


	# plot the first 9 basis Chebyshev Polynomial Basisi Fnctions
	function q3()

		plot_arr = Dict()
		x = linspace(0,1.0,100)
		fcheby(x,deg) = [cos(acos(x) * deg)]
		for deg in 0:8
			y = fcheby.(x, deg)
			plot_arr[deg+1] = plot(x, y, ylim = [-1.0, 1.0],labels = ["degree = $deg"])
		end
		Plot_total = plot(title = "Chebyshev  basis functions", plot_arr[1], plot_arr[2], plot_arr[3], plot_arr[4], plot_arr[5], plot_arr[6], plot_arr[7], plot_arr[8], plot_arr[9], layout = (3,3))

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
		display(plot_total)
		q2(15)
		q3()
		display(Plot:total)
		#q4a()
		#q4b()
		#q5()
	end
end



end
