using FastGaussQuadrature
using Gadfly
using Base.Test
using ApproxFun
using ApproXD
using Plots

module funcapp

	# use chebyshev to interpolate this:
	function q1(n_cheb)
		#####QUESTION 1#####
		# Approximate f(x) = x + 2x^2 - exp(-x) for x in -3:3
		# n = # interpolation points
		# Approximation of degree deg=n-1
		# n = 15 Chebyshev interpolation nodes
		# Define function
		f(x::Vector{Float64}) = x + 2.*x.^2 - exp(-x)
		f(x::LinSpace{Float64}) = x + 2.*x.^2 - exp(-x)
		# Degree of approximation
		deg = n_cheb-1
		# Bounds of domain
		a = -3
		b = 3
		# Chebyshev interpolation nodes
		nodes = gausschebyshev(n_cheb)
		# Map z into x
		x = (1/2)*(a+b) + (1/2)*(b-a)*nodes[:1]
		# Get function values at x, y=f(x)
		y = f(x)
		# Evaluate the Chebyshev basis matrix z
		z = 2.*(x-a)./(b-a) - 1
		#phi = [cos((n-i+0.5)*(j-1)*pi/n) for i = 1:n_cheb, j = 1:n_cheb]
		phi = [cos(acos(z[i])j) for i = 1:n_cheb, j = 0:deg]
		transpose(phi)*phi
		# Compute approximation coeffs c by matrix inversion
		# Use the backslash operator \ to achieve this. type ?\ to know more.
		c = \(phi,y)
		# Predict
		n_new = 100
		nodes_new = linspace(a,b,n_new)
		z_new = 2.*(nodes_new-a)./(b-a) - 1
		phi_new = [cos(acos(z_new[i])*j) for i = 1:n_new, j = 0:deg]
		f_hat = phi_new*c
		y = f(nodes_new)
		epsi = y - f_hat
		# Graphs
		graph1_1 = Gadfly.plot(layer(x=nodes_new, y=y, Geom.line),
		                layer(x=nodes_new,y=f_hat, Geom.point, Theme(default_point_size=1pt)))
		graph1_2 = Gadfly.plot(x=nodes_new, y=epsi, Geom.point)
		graph1 = hstack(graph1_1,graph1_2)
		savefig("question1.png")
		display(graph1)
		# Automated test
		t1 = @test maximum(epsi)<1e-9
		display(t1)
	end

	function q2(n)
	    #####QUESTION 2#####
	    f(x) = x + 2*x^2 - exp(-x)
	    f(x::Vector{Float64}) = x + 2.*x.^2 - exp(-x)
	    f(x::LinSpace{Float64}) = x + 2.*x.^2 - exp(-x)
	    a = -3
	    b = 3
	    g = Fun(f, Chebyshev(a..b))
	    n_new = 100
	    nodes_new2 = linspace(a,b,n_new)
	    g_hat = g.(nodes_new2)
	    y2 = f.(nodes_new2)
	    epsi2 = abs(y2 - g_hat)
	    # Graphs
	    graph2_1 = Gadfly.plot(layer(x=nodes_new2, y=y2, Geom.line),
	                        layer(x=nodes_new2,y=g_hat, Geom.point, Theme(default_point_size=1pt)))
	    graph2_2 = Gadfly.plot(x=nodes_new2, y=epsi2, Geom.point)
	    graph2 = hstack(graph2_1, graph2_2)
			savefig("question2.png")
	    display(graph2)
	    # Automated test
	    t2 = @test maximum(epsi2)<1e-9
	    display(t2)
	end


	# plot the first 9 basis Chebyshev Polynomial Basis Fnctions
	function q3()
		phi3(z,j) = cos(acos(z)j)
		graph3 = Plots.plot([z -> phi3(z,j) for j in 1:9], -1, 1, layout=(3,3), label=transpose(["Basis num $j" for j in 1:9]))
		savefig("question3.png")
		display(graph3)
	end

	function q4a()
		function Cheby_4a(n)
		t(y)= 1/(1+25*(y^2))

		  # Degree of approximation
		  deg = n-1
		  # Bounds of domain
		  a = -5
		  b = 5
		  # Chebyshev interpolation nodes
		  nodes = gausschebyshev(n)
		  # Map z into x
		  x = (1/2)*(a+b) + (1/2)*(b-a)*nodes[:1]
		  # Get function values at x, y=f(x)
		  y=zeros(n)
		  for i in 1:n
		  y[i] = t(x[i])
		end
		  # Evaluate the Chebyshev basis matrix z
		  z = 2*(x-a)/(b-a) - 1
		  #phi = [cos((n-i+0.5)*(j-1)*pi/n) for i = 1:n_cheb, j = 1:n_cheb]
		  phi = [cos(acos(z[i])j) for i = 1:n, j = 0:deg]
		  transpose(phi)*phi
		  # Compute approximation coeffs c by matrix inversion
		  # Use the backslash operator \ to achieve this. type ?\ to know more.
		  c = \(phi,y)
		  g_hat= phi*c
		  errors= y - g_hat
		return (g_hat,x)
		end

		Cheby_4a(6)
		Cheby_4a(10)
		Cheby_4a(16)
		b(y)=1/(1+25*(y^2))
		graph_1= Plots.plot(b,-5,5)
		Plots.plot!(Cheby_4a(6)[2],Cheby_4a(6)[1],label="5 degrees")
		Plots.plot!(Cheby_4a(10)[2],Cheby_4a(10)[1],label="9 degrees")
		Plots.plot!(Cheby_4a(16)[2],Cheby_4a(16)[1],label="15 degrees")

		function Uniform_4a(n)
		t(y)= 1/(1+25*(y^2))
		  # Degree of approximation
		  deg = n-1
		  # Bounds of domain
		  a = -5
		  b = 5
		#Uniform interpolation nodes
		  nodes = linspace(a,b,n)
		  # Map z into x
		  x = zeros(n)
		   for i in 1:n
		     x[i]=nodes[i]
		   end
		  # Get function values at x, y=f(x)
		  y=zeros(n)
		  for i in 1:n
		  y[i] = t(x[i])
		end
		  # Evaluate the Chebyshev basis matrix z
		  z = 2*(x-a)/(b-a) - 1
		  #phi = [cos((n-i+0.5)*(j-1)*pi/n) for i = 1:n_cheb, j = 1:n_cheb]
		  phi = [cos(acos(z[i])j) for i = 1:n, j = 0:deg]
		  transpose(phi)*phi
		  # Compute approximation coeffs c by matrix inversion
		  # Use the backslash operator \ to achieve this. type ?\ to know more.
		  c = \(phi,y)
		  g_hat= phi*c
		  errors= y - g_hat
		  return (g_hat,x)
		end
		Uniform_4a(6)
		Uniform_4a(10)
		Uniform_4a(16)
		b1(y)=1/(1+25*(y^2))
		graph_2= Plots.plot(b1,-5,5)
		Plots.plot!(Uniform_4a(6)[2],Uniform_4a(6)[1],label="5 degrees")
		Plots.plot!(Uniform_4a(10)[2],Uniform_4a(10)[1],label="9 degrees")
		Plots.plot!(Uniform_4a(16)[2],Uniform_4a(16)[1],label="15 degrees")

		graph_final=Plots.plot(graph_1,graph_2,layout=2)
		savefig("question4a.png")
		display(graph_final)
	end


	function q4b()
		#First version==>Random Knots
		g(x)= 1/(1+25*(x^2))
		bs=BSpline(13,3,-5,5)
		#15 coefficients because 13(p)-3(k)-1
		B=full(getBasis(collect(linspace(-5,5,15)),bs))
		p=linspace(-5,5,15)
		v=zeros(15)
		for i in 1:15
		  v[i]=g(p[i])
		end
		B_1= inv(B)
		coeff=B_1*v

		B_new=full(getBasis(collect(linspace(-5,5,65)),bs))
		y=linspace(-5,5,65)
		v_1=zeros(65)
		for i in 1:65
		  v_1[i]=g(y[i])
		end
		approx_1= B_new*coeff
		z=v_1-approx_1

		#Second version==> we choose the Knots
		knots=linspace(-1,1,13)
		my_knots=zeros(13)
		for i in 1:13
		  my_knots[i]=knots[i]
		end
		bs2=BSpline(my_knots,3)
		#15 coefficients because 13(p)-3(k)-1
		B_2=full(getBasis(collect(linspace(-1,1,15)),bs2))
		h=linspace(-1,1,15)

		m=zeros(15)
		for i in 1:15
		  m[i]=g(h[i])
		end
		B_2inv= inv(B_2)
		coeff=B_2inv*m

		B_new2=full(getBasis(collect(linspace(-1,1,65)),bs2))
		y2=linspace(-1,1,65)
		v_2=zeros(65)
		for i in 1:65
		  v_2[i]=g(y2[i])
		end
		approx_2= B_new2*coeff
		z2=v_2-approx_2

		#Plots and Panels
		graph4b_1 = Plots.plot(g,-5,5, label="Runge's function")
		graph4b_2 = Plots.plot(y,z,ylim=(-0.5,0.5))
		Plots.plot!(y2,z2, ylim=(-0.5,0.5))
		graph4 = Plots.plot(graph4b_1, graph4b_2, layout=2)
		display(graph4)
		savefig("question4b.png")
	end

	function q5()
		j(x)=abs(x)^(1/2)
		bs=BSpline(13,3,-1,1)
		B=full(getBasis(collect(linspace(-1,1,15)),bs))
		p=linspace(-1,1,15)
		v=zeros(15)
		for i in 1:15
		  v[i]=j(p[i])
		end
		B_1= inv(B)
		coeff=B_1*v

		B_new=full(getBasis(collect(linspace(-1,1,65)),bs))
		y=linspace(-1,1,65)
		v_1=zeros(65)
		for i in 1:65
		  v_1[i]=j(y[i])
		end
		approx_1= B_new*coeff
		z=v_1-approx_1


		#Second version
		knots_neg=linspace(-1,0,6)
		knots_pos=linspace(0,1,6)

		my_knots=zeros(13)
		for i in 1:6
		  my_knots[i]=knots_neg[i]
		  my_knots[6]=0
		  my_knots[i+7]=knots_pos[i]
		  end

		bs2=BSpline(my_knots,3)

		B_2=full(getBasis(collect(linspace(-1,1,15)),bs2))

		h=linspace(-1,1,15)

		m=zeros(15)
		for i in 1:15
		  m[i]=j(h[i])
		end
		B_2inv = inv(B_2)
		coeff=B_2inv*m

		B_new2=full(getBasis(collect(linspace(-1,1,65)),bs2))
		y2=linspace(-1,1,65)
		v_2=zeros(65)
		for i in 1:65
		  v_2[i]=j(y2[i])
		end
		approx_2= B_new2*coeff
		z2=v_2-approx_2

		#Plots and Panels
		graph5_1 = Plots.plot(j,-1,1)
		graph5_2 = Plots.plot(y,approx_1,xlim=(-1,1))
		Plots.plot!(y2,approx_2,xlim=(-1,1))
		graph5_3 = Plots.plot(y,z,ylim=(-2,2))
		Plots.plot!(y2,z2, ylim=(-0.05,0.05))
		graph5 = Plots.plot(graph5_1, graph5_2, graph5_3, layout=3)
		savefig("question5.png")
		display(graph5)
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
