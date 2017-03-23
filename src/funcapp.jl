

module funcapp

using FastGaussQuadrature
using ApproxFun
using Gadfly
using Plots
using DataFrames
using Base.Test
using ApproXD
using Interpolations

# Functions to approximate
function truef1(x)
    x + 2.0.*x.^2.0 - exp(-x)
end

function truef2(x)
    (1.0 .+ 25.0.*x.^2.0).^(-1.0)
end

function truef3(x)
    abs(x).^.5
end

# Useful functions for Chebychev approximation
ChebyT(x,deg) = cos(acos(x)*deg)
unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]

##############################################################################
################################# QUESTION 1 #################################
##############################################################################

	# CHEBYSCHEV APPROXIMATION OF A FUNCTION f
	function q1(n)

		# Definition of the paramateres of the interpolation
		J = n-1
		a = -3
		b = 3

		# Defintion of the Chebyshv interpolation nodes
		z = gausschebyshev(n)
		x = 0.5*(a+b)+0.5*(b-a)*z[:1] # Adapting the nodes to the [-3; 3] space
		y = truef1.(x)

		# Evaluating the Chebyschev polynomial using the nodes z
		phi = Matrix{Float64}(n, J+1)
		for i =1:n, j =1:(J+1)
			#phi[i,j] = cos((n-i+0.5)*(j-1)*pi/n)
			phi[i,j]=cos(acos(z[:1][i])*(j-1))
		end
		# Solve for the approximation coefficients
		c = phi\y

		###

		# Using the approximation to predict n_new points in [-3;3]
		n_new = 100
		x_new = linspace(a,b,n_new)
		z_new = unitmap(x_new,a,b) # Adapting the new points to the [-1;1] space
		phi_new = Matrix{Float64}(n_new, J+1)
		for i =1:n_new, j =1:(J+1)
			phi_new[i,j] = ChebyT(z_new[i], j-1)
		end
		y_predict = phi_new*c
		y_true = truef1.(x_new)
		y_dev = abs(y_predict - y_true)

		# Testing the accuracy of the aproximation
        @testset "Q1: Testing Accuracy of the Approximation:" begin
        @test maximum(y_dev)<1e-9
        end

		# Printing the Results
		plotq1 = Gadfly.plot(layer(x=x_new, y=y_predict, Geom.point, Theme(default_color=colorant"blue", default_point_size=1pt, highlight_width=0.05pt)), layer(truef1,a,b, Geom.line, Theme(default_color=colorant"purple", line_width=4pt)), Guide.manual_color_key("", ["Approximated Points", "True Function"], ["blue", "purple"]), Guide.Title("Q1: True and approx function"))
    plotq1_bis = Gadfly.plot(x=x_new, y=y_dev, Geom.line, Guide.Title("Q1 : Deviation"), Guide.ylabel("Deviation(in absolute terms)"))
		display(hstack(plotq1,plotq1_bis))

	end

	##############################################################################
################################# QUESTION 2 #################################
##############################################################################

	# APPROXIMATION OF A FUNCTION f USING THE APPROXFUN PACKAGE
	function q2(n)

		# Definition of the paramateres of the interpolation
		a = -3
		b = 3

		# Approximating f using Chebyshev nodes
		f_approx = Fun(truef1, Chebyshev(-3..3))

		# Using our approximator to predict n_new points
		n_new = 100
		x_new = linspace(a,b,n_new)
		y_predict=f_approx.(x_new)
		y_dev = abs(truef1.(x_new) - y_predict)

		# Testing the accuracy of the aproximation
        @testset "Q2: Testing Accuracy of the Approximation using ApproxFun:" begin
        @test maximum(y_dev)<1e-9
        end

		plotq2 = Gadfly.plot(layer(x=x_new, y=y_predict, Geom.point, Theme(default_color=colorant"blue", default_point_size=1pt, highlight_width=0.05pt)), layer(truef1,a,b, Geom.line, Theme(default_color=colorant"purple", line_width=4pt)), Guide.manual_color_key("", ["Approximated Points", "True Function"], ["blue", "purple"]), Guide.Title("Q2: True and approx using ApproxFun"))
    plotq2_bis = Gadfly.plot(x=x_new, y=y_dev, Geom.line, Guide.Title("Q2: Deviation"), Guide.ylabel("Deviation (in absolute terms)"))
		display(hstack(plotq2,plotq2_bis))


	end


##############################################################################
################################# QUESTION 3 #################################
##############################################################################

	# PLOTTING THE FIRST 9 CHEBYSCHEV BASIS POLYNOMIALS
	function q3()

		J = 9
		a = -3
		b = 3

		# Defining each the Chebyschev polynomials (of degree j)
		poly_cheb = (x, j) -> ChebyT(unitmap(x,a,b),j-1)
		# Plotting the results
 		plotq3 = Plots.plot([x -> poly_cheb(x, j) for j in 1:J], a, b, layout = (3,3), legend=false, title = ["Polynomial of degree $j" for j in 0:J-1]', titlefont = font(7))
		display(plotq3)

	end

##############################################################################
################################# QUESTION 4 #################################
##############################################################################

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

################################# Question 4a #################################

	# POLYNOMIAL INTERPOLATION : UNIFORM VS CHEBYSHEV NODES
	function q4a(deg=(5,9,15),lb=-1.0,ub=1.0, typenodes="cheby")

		# Function to approximate (Runge's function)
		f = (x) -> 1./(1+25.*x.^2)
		# Parameters
		n_new = 100
		x_new = linspace(lb,ub,n_new)
		colors = [colorant"yellow", colorant"orange", colorant"red"]

		global plot1 = Gadfly.plot(layer(f,lb,ub, Geom.line, Theme(default_color=colorant"darkblue", line_width=2pt)), Guide.YLabel("f(x)"), Guide.Title("Q4a : Runge Approximation using "* typenodes*" nodes"))

		# Approximating the function with various degree of polynomials (input "deg")
		for j = 1:length(deg)

			n = deg[j]+1
			if typenodes == "uniform"
				nodes = linspace(lb, ub, n)
			end
			if typenodes == "cheby"
				nodes = gausschebyshev(n)[:1]
			end

			# Using the function ChebyType to get the approximating function
			approx= ChebyType(nodes, deg[j],lb, ub, f)

			# Using our approximating function to make predictions on x_new points
			dict_result = predict(approx,x_new)
			y_predict =  dict_result["preds"]
			push!(plot1,layer(x=x_new, y=y_predict,Geom.line, Theme(default_color=colors[j]), order=2))

		end

		push!(plot1, Guide.manual_color_key("Legend", ["Degree 5", "Degree 9", "Degree 15"], ["yellow", "orange" , "red"]))
		# Plotting the final results
		return plot1

	end

  ################################# Question 4b #################################
	# SPLINE INTERPOLATION : KNOTS PLACEMENT
	function q4b()

		# Parameters of the interpolation
		lb = -5
		ub = 5
		nknots = 13
		deg = 3 # Cubic spline

		# Points used for computing the apprixomation function
		n = 65 # So that we have more observations than parameters to estimate
		x = collect(linspace(lb,ub,n))
		y = truef2.(x)
		# Points used for making predictions using the approximation function
		n_new = 1000
		x_new = collect(linspace(lb,ub,n_new))
		y_new = truef2.(x_new)

		## Version 1 B-Spline Approximation : equally spaced knots
		b1 = ApproXD.BSpline(nknots,deg,lb,ub) # Getting the basis function
		phi1 = ApproXD.getBasis(x,b1) #Evaluating the basis function at the points x
		c1 = phi1\y
		knots1 = collect(linspace(lb,ub,nknots)) # knots1 = nknots equally spaced knots

		## Version 2 B-Spline Aproximation : concentrated knots
		knots2 = vcat(collect(linspace(lb, lb/4,3)), collect(linspace(lb/4+eps(),ub/4-eps(),7)), collect(linspace(ub/4,ub,3))) # we construct concentrated towards 0 knots manually using linspace with different n
		b2 = ApproXD.BSpline(knots2,deg)
		phi2 = ApproXD.getBasis(x,b2)
		c2 = phi2\y


		# Making predictions using n_new points
		error_predict1 = ApproXD.getBasis(x_new,b1)*c1 - y_new
		error_predict2 = ApproXD.getBasis(x_new,b2)*c2 - y_new
		#display(error_predict2)

		# Plotting the results
		plot1 = Gadfly.plot(truef2,-5,5, Theme(line_width = 3pt, major_label_font_size = 10pt, major_label_color = colorant"white"), Guide.Title("Q4b: Runge's Function"))
    plot2 = Gadfly.plot(layer(x=x_new, y=error_predict1, Theme(line_width = 2.5pt, default_color=colorant"blue")), Geom.line, Guide.Title("Q4b : Error in Runge's function"))
		plot2 = push!(plot2, layer(x=x_new, y=error_predict2, Theme(line_width = 2.5pt, default_color=colorant"green")), Geom.line)
		plot2 = push!(plot2, layer(xintercept=knots1, Theme(default_color=colorant"blue", line_width=0.5pt)), Geom.vline)
		plot2 = push!(plot2, layer(xintercept=knots2, Theme(default_color=colorant"green", line_width=0.5pt)), Geom.vline)
    display(vstack(plot1,plot2))


	end

  ##############################################################################
  ################################# QUESTION 5 #################################
  ##############################################################################

	# SPLINES AND KINKS
	function q5(k)

		# Parameters
		a,b = (-1,1)
		nknots = 13
		deg = 3

		# Function for myknots with k number of knots set at 0
		function myknots(k)
		    if isodd(k)==true
		        vcat(-1, -.8, collect(linspace(-.5, -.1, 5-(k+1)/2)), zeros(k), collect(linspace(.1, .5, 5-(k+1)/2)), .8, 1)
		    else
		        vcat(-1, -.8, collect(linspace(-.5, -.1, 5-k/2)), zeros(k), collect(linspace(.1, .5, 4-k/2)), .8, 1)
		    end
		end

		# Getting B Splines
		bs1 = ApproXD.BSpline(nknots,deg,a,b)
		bs2 = ApproXD.BSpline(myknots(k), deg)
		nevals = 65

		# Getting coefficients c
		points = collect(linspace(a,b,nevals))
		b1 = full(ApproXD.getBasis(points,bs1))
		b2 = full(ApproXD.getBasis(points,bs2))
		c1 = b1 \ (truef3.(points))
		c2 = b2 \ (truef3.(points))

		# Approximate functions with a 1000 equally spaced elements array
		nodes = collect(linspace(a,b,1000))
		truef = truef3.(nodes)
		f1 = getBasis(nodes,bs1) * c1
		f2 = getBasis(nodes,bs2) * c2

		# Difference with the true function
		Diff5_1 = truef - f1
		Diff5_2 = truef - f2

		figure5_a = Gadfly.plot(x=nodes, y=truef, Geom.line, Guide.Title("Q5: True function with a kink"))
    figure5_b = Gadfly.plot(layer(x=nodes, y=f1, Geom.line, Theme(default_color=colorant"deepskyblue")), layer(x=nodes, y=f2, Geom.line, Theme(default_color=colorant"blue")), Guide.manual_color_key("Legend", ["Uniform", "Centered towards 0"], ["deepskyblue", "blue"]), Guide.Title("Q5: Approx with uniform and centrerd towards 0 knots distribution"))
    figure5_c = Gadfly.plot(layer(x=nodes, y=Diff5_1, Geom.line, Theme(default_color=colorant"deepskyblue")), layer(x=nodes, y=Diff5_2, Geom.line, Theme(default_color=colorant"blue")), Guide.manual_color_key("Legend", ["Uniform", "Centered towards 0"], ["deepskyblue", "blue"]), Guide.Title("Q5: Approx with uniform and centrerd towards 0 knots distribution"))

		display(figure5_a)
		display(figure5_b)
		display(figure5_c)

	end


		# function to run all questions
	function runall()
		println("running all questions of HW-funcapprox:")
		q1(15)
		q2(15)
		q3()
    plot1 = q4a((5,9,15), -1 , 1, "uniform")
		plot2 = q4a((5,9,15),-1 ,1 , "cheby")
		display(plot1)
		display(plot2)
		q4b()
		q5(1)
		println("Q5: We could not set more than one knot to 0 otherwise the inverse-matrix operator in the computation of the coefficients would yield an Array of NaN elements")
	end


end
