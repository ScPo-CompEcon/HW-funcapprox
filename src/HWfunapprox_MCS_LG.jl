module funcapp
	using Plots
	using FastGaussQuadrature
	using ApproxFun
	using ApproXD

	function f1(x)
		return(x + 2 * x^2 - exp(-x))
	end

# use chebyshev to interpolate this:
function q1(n)
	deg = n - 1 #degree of approximation
	a = -3 #lower bound
	b = +3
	z = gausschebyshev(n)[1]
	x = 0.5 * (a + b) + 0.5 * (b - a) * z
	values = [f1(x[i]) for i in 1:length(x)]
	Cheby = [cos((n - i + 0.5) * (j - 1)* π / n) for i = 1:n, j = 1:n]
	c = \(Cheby, values)
	function fhat(x)
		z = 2 * (x - a)/(b - a) -1
		Cheby = [cos(acos(z)*j) for j = 0:n-1]
		return(transpose(c)*Cheby)
	end
	n_new = 100
	new_x = [a + (i-1)/(n_new -1)*(b-a) for i in 1:n_new] #new equally spaced points
	prediction = [fhat(new_x[i]) for i in 1:length(new_x)] #approximation
	new_values = [f1(new_x[i]) for i in 1:length(new_x)] #true value
	error = new_values-prediction #error between prediction and true value
	global plot1 = plot(new_x, [prediction, new_values], labels = ["approximation", "truth"], ylim = [-5,5], title = "Plain Vanilla")
	global plot2 = plot(new_x, error, labels = ["deviation"])
	global plot_total1 = plot(plot1, plot2)
	return Dict(:error => maximum(error))
	display(plot_total1)
end

function q2(n)
	x = Fun(identity,-3..3)
	fhat2 = x + 2 * x^2 - exp(-x)
	n_new = 100
	new_x = [-3 + (i-1)/(n_new -1)*6 for i in 1:n_new] #new equally spaced points
	prediction = [fhat2(new_x[i]) for i in 1:length(new_x)] #approximation
	new_values = [f1(new_x[i]) for i in 1:length(new_x)] #true value
	error = new_values-prediction #error between prediction and true value
	#return Dict(:err => 1.0)
	plot1 = plot(new_x, [prediction, new_values], labels = ["approximation", "truth"])
	plot2 = plot(new_x, error, labels = ["deviation"])
	global plot_total2 = plot(plot1, plot2)
end

ChebyT(x,deg) = cos(acos(x)*deg)
# plot the first 9 basis Chebyshev Polynomial Basis Fnctions
function q3()
#We want to plot the following function for x = -1,1
#for deg in 0,1,2,...,8
	plot_total = Dict()
	for deg in 0:8
		global x = [-1 + (i-1)/(100 -1)*2 for i in 1:100]
		global y = ChebyT.(x, deg)
		global plot_total[deg] = plot(x, y, ylim = [-1.01, 1.01], labels = ["degree = $deg"])
	end
	global plot_total3 = plot(plot_total[0], plot_total[1], plot_total[2], 	plot_total[3], plot_total[4], plot_total[5], plot_total[6], plot_total[7],  plot_total[8], layout = (3,3))
 #new equally spaced points
	global plot_eqsp = plot(x, y, ylim = [-1,1])
end

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

function Runge(x)
	return(1/(1+25*x^2))
end

function q4a(deg=(5,9,15),lb=-5,ub=5)
	Eqspaced = Dict()
	Chebynodes = Dict()
	x_new = [lb + (j-1)/(100 -1)*(ub-lb) for j in 1:100] #we will compare the functions for 100 points
	for i in deg
		n = i + 1
		nodes = [lb + (j-1)/(n -1)*(ub-lb) for j in 1:n]
		Cheb = ChebyType(nodes, n, lb, ub, Runge)
		Eqspaced[i] = predict(Cheb, x_new)
	end
	for i in deg
		n = i +1
		z = gausschebyshev(n)[1]
		nodes = 0.5 * (lb + ub) + 0.5 * (ub - lb) * z
		Cheb = ChebyType(nodes, n, lb, ub, Runge)
		Chebynodes[i] = predict(Cheb, x_new)
	end
	#maybe change the plots to show in the x-axis the value of the xs and not their ordinal
	#in the list
	plot1 = plot(x = Chebynodes[5]["x"], [Chebynodes[5]["truth"], Chebynodes[5]["preds"], Chebynodes[9]["preds"], Chebynodes[15]["preds"]], title = "Cheby Nodes", labels = ["true", "deg = 5", "deg = 9", "deg = 15"])
	plot2 = plot(x = Eqspaced[5]["x"], [Eqspaced[5]["truth"], Eqspaced[5]["preds"], Eqspaced[9]["preds"], Eqspaced[15]["preds"]], title = "Equidistanced Nodes", labels = ["true", "deg = 5", "deg = 9", "deg = 15"])
	global plot_total4a = plot(plot1, plot2)
end

function q4b()
	b1 = BSpline(13,3,-5,5)
	my_knots = [-5.0, -2.5, -1, -0.5, -0.2, -0.1,
	0.0, 0.1, 0.2, 0.5, 1, 2.5, 5.0]
	b2 = BSpline(my_knots,3)
	nevals=65
	nodes = collect(linspace(-5,5,nevals))
  	coef1 = getBasis(nodes,b1) \ Runge.(nodes)
  	coef2 = getBasis(nodes, b2) \ Runge.(nodes)
	new_x = collect(linspace(-5,5,100))
	new_values = Runge.(new_x)
	preds1 = getBasis(new_x,b1) * coef1
	preds2 = getBasis(new_x,b2) * coef2
	err1 = preds1 - new_values
	err2 = preds2 - new_values
	plot1 = plot(Runge, -5, 5, title = "Runge")
	plot2 = plot([ err1, err2], labels =["equally spaced nodes", "concentrated nodes"])
	global plot_total4b = plot(plot1, plot2)
end

function q5()
	function f(x)
		return(abs(x)^0.5)
	end
	plot1 = plot(f, title = "True Function")
	b1 = BSpline(13,3,-1,1)
	my_knots = collect(linspace(-1, 0.1, 5))
	push!(my_knots, 0,0,0)
	append!(my_knots, collect(linspace(0.1,1,5)))
	b2 = BSpline(sort(my_knots),3)
	nevals = 65
	nodes = collect(linspace(-1,1,nevals))
  	coef1 = getBasis(nodes, b1) \ f.(nodes)
  	coef2 = getBasis(nodes, b2) \ f.(nodes)
	new_x = collect(linspace(-1,1,100))
	new_values = f.(new_x)
	preds1 = getBasis(new_x,b1) * coef1
	preds2 = getBasis(new_x,b2) * coef2
	err1 = preds1 - new_values
	err2 = preds2 - new_values
	plot2 = plot(x = new_x, [preds1, preds2], title = "Cubic Spline Approximation", label = ["equally spaced", "stacked knots"])
	plot3 = plot([ err1, err2], labels =["equally spaced", "stacked knots"], title = "Error")
	global plot_total5 = plot(plot1, plot2, plot3)
end

		# function to run all questions
function runall()
	info("Running all questions of HWFunApprox:")
	info("Question 1:")
	q1(15)
	display(plot_total1)
	savefig("Plot1.png")
	println("Plot for question 1 added to the current directory as a .png file.")
	info("Question 2:")
	q2(15)
	display(plot_total2)
	savefig("Plot2.png")
	println("Plot for question 2 added to the current directory as a .png file.")
	info("Question 3:")
	q3()
	display(plot_total3)
	savefig("Plot3_1.png")
	display(plot_eqsp)
	savefig("Plot3_2.png")
	println("Plots for question 3 added to the current directory as .png files.")
	info("Question 4.A:")
	q4a()
	display(plot_total4a)
	savefig("Plot4a.png")
	println("Plot for question 4.A added to the current directory as a .png file.")
	info("Question 4.B:")
	q4b()
	display(plot_total4b)
	savefig("Plot4b.png")
	println("Plot for question 4.B added to the current directory as a .png file.")
	info("Question 5:")
	q5()
	display(plot_total5)
	savefig("Plot5.png")
	println("Plot for question 5 added to the current directory as a .png file.")
end
end
