

module funcapp

	using Plots
	using FastGaussQuadrature
	using ApproxFun
	using ApproXD
	using FactCheck

	f(x::Float64)=x+2*x^2-exp(-x)
	f(x::Vector{Float64})=x+2*x.^2-exp(-x)

	r(x::Float64)=(1+25*x^2)^(-1)
	r(x::Vector{Float64})=(1+25*x.^2).^(-1)

	g(x::Float64)=(abs(x))^(0.5)
	g(x::Vector{Float64})=(abs(x)).^(0.5)

	# use chebyshev to interpolate this:
	function q1(n=15)
		deg=n-1
		nodes = gausschebyshev(n)[1]
		a=-3
		b=3
		denormalized_nodes=0.5*(b+a)+0.5*(b-a)*nodes

		T=zeros(n,deg)
		T[:,1]=1
		T[:,2]=nodes

		for i=2:deg-1
		  for j=1:n
		    T[j,i+1]=2*nodes[j]*T[j,i]-T[j,i-1]
		  end
		end

		y=f(denormalized_nodes)
		c=T\y
		approx_f=vec(T*c)

		return denormalized_nodes, approx_f

	end

	function q1_graph(approx_f)
		p=plot(f,label="True function",layout=2)
		plot!(p[2],approx_f,label="Approximation")
	end

	function q2()
		f=x->x+2*x.^2-exp(-x)
		x  = Fun(f,-3..3)
		space(x)

		approx_f=f(x)

		return approx_f

	end

	function q2_graph(approx_f)
		p=plot(f,label="True function",layout=2)
		plot!(p[2],approx_f,label="Approximation")
	end


	# plot the first 9 basis Chebyshev Polynomial Basisi Fnctions
	function q3()
		nodes = gausschebyshev(100)[1]

		T=zeros(100,100)
		T[:,1]=1
		T[:,2]=nodes
		deg=99

		for i=2:deg
		  for j=1:100
		    T[j,i+1]=2*nodes[j]*T[j,i]-T[j,i-1]
		  end
		end

		plot(T[:,1:9],layout=(3,3),grid=false,label=transpose(["Basis num $i" for i in 1:9]),linewidth=2,linecolor=:black)
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

	function q4a(k,a,b)
		Cheby_pol=Matrix(3,2)
		Predictions=Matrix(4,2)

		for i=1:length(k)
		  nodes = gausschebyshev(k[i]+1)[1]
		  chebyshev_nodes=0.5*(b+a)+0.5*(b-a)*nodes
		  uniform_nodes=collect(linspace(-5,5,k[i]+1))

		  Cheby_pol[i,1]=funcapp.ChebyType(chebyshev_nodes,k[i],a,b,funcapp.r)
		  Cheby_pol[i,2]=funcapp.ChebyType(uniform_nodes,k[i],a,b,funcapp.r)

		  Predictions[i,1]=funcapp.predict(Cheby_pol[i,1],chebyshev_nodes)
		  Predictions[i,2]=funcapp.predict(Cheby_pol[i,2],uniform_nodes)

		end

		uniform_nodes=collect(linspace(-5,5,16))

		p=plot(r(uniform_nodes),label="True function",layout=2)
		plot!(p[1],Predictions[1,1]["preds"],label="Approximation with degree 5")
		plot!(p[1],Predictions[2,1]["preds"],label="Approximation with degree 9")
		plot!(p[1],Predictions[3,1]["preds"],label="Approximation with degree 15")

		plot!(p[2],r(uniform_nodes),label="True function")
		plot!(p[2],Predictions[1,2]["preds"],label="Approximation with degree 5")
		plot!(p[2],Predictions[2,2]["preds"],label="Approximation with degree 9")
		plot!(p[2],Predictions[3,2]["preds"],label="Approximation with degree 15")

	end

	function q4b(a,b)

		nb_knots=13
		deg=5
		bs = BSpline(nb_knots,3,a,b)
		x=linspace(a,b,deg*nb_knots)
		B = full(getBasis(collect(x),bs))
		c=B\g(collect(x))
		approx_g=vec(B*c)

		knots1=[-5 -3 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 3 5]
		knots1=vec(knots1)
		bs1 = BSpline(knots1,3)
		B1 = full(getBasis(collect(x),bs1))
		c=B1\g(collect(x))
		approx_g1=vec(B1*c)

		p=plot(g(collect(x)),layout=2,label="Original function")
		plot!(p[2],g(collect(x))-approx_g,color=:green, label="Error with uniform nodes")
		plot!(p[2],g(collect(x))-approx_g1,color=:blue, label="Error with concentrated nodes")

	end

	function q5()

		nb_knots=13
		deg=5

		a=-1
		b=1
		y=linspace(a,b,65)

		bs_unif = BSpline(nb_knots,3,a,b)
		B_unif = full(getBasis(collect(y),bs_unif))

		c_unif=B_unif\g(collect(y))


		approx_g_unif=vec(B_unif*c_unif)

		knots_mult1=vec([-1 -0.8 -0.6 -0.4 -0.2 -0.01 0 0.01 0.2 0.4 0.6 0.8 1])
		#knots_mult2=vec([-1 -0.75 -0.5 -0.25 -0.02 -0.01 0 0.01 0.02 0.25 0.50 0.75 1])

		bs_mult= BSpline(knots_mult1,3)
		B_mult = full(getBasis(collect(y),bs_mult))

		c_mult=B_mult\g(collect(y))


		approx_g_mult=vec(B_mult*c_mult)

		p=plot(y,g(collect(y)),layout=3,label="True function")
		plot!(p[2],approx_g_unif, label="Uniform approximation")
		plot!(p[3],approx_g_mult, label="Approximation with multiplicity of nodes")




	end


		# function to run all questions
	function runall()
		println("running all questions of HW-funcapprox:")
		approx1=q1(100)
		println("In question 1, the obtained nodes when n=100 are:")
		println(approx1[1])
		println("and the corresponding approximation:")
		println(approx1[2])
		@fact abs(maximum(approx1[2]-f(approx1[1])))--> less_than(1e-9)
		plot=q1_graph(approx1[2])
		display(plot)

		approx2=q2()
		println("In question 2,the corresponding approximation with the ApproxFun package is:")
		println(approx2)
		plot2=q2_graph(approx2)
		display(plot2)

		println("The first nine basis polynomials are:")
		plot3=q3()
		display(plot3)


		plot4a=q4a((5,9,15),-5,5)
		display(plot4a)
		plot4b=q4b(-5,5)
		display(plot4b)

		plot5=q5()
		display(plot5)
	end


end
