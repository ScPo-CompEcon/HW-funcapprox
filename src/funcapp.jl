	

module funcapp

##SETUP
    
    #packages required
    using Plots
    using FastGaussQuadrature
    using Base.Test
    using ApproxFun
    using ApproXD

    #options
    pyplot()

    #define some other functions
    function ChebyBasis(ni,nj,nodes) #calculates the Chebyshev basis matrix
    V = Array(Float64, ni, nj)
    V[:,1] = 1
    V[:,2] = nodes
        for j in 3:nj
            for i in 1:ni
                V[i,j] = 2*V[i,2]*V[i,j-1]-V[i,j-2]
            end
        end
        return V
    end

    unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1     #[a,b] -> [-1,1]

    #function to get coefficients
    function get_coeffs(lb::Float64, ub::Float64, nodes::Vector, f::Function)
        n = length(nodes)
        y = f.(5.*nodes)
        V = ChebyBasis(n,n,nodes)
        c = V \ y
        return c
    end

    #function to give true and predicted y-values
    function predict(n_new::Int, nj::Int, coeffs::Vector, lb::Float64, ub::Float64, f::Function)
        x_new = linspace(lb, ub, n_new)
        true_y = map(f, x_new)
        z_new = map((x) -> unitmap(x,lb,ub), x_new)
        cheby_basis = ChebyBasis(n_new,nj,z_new)
        cheby_y = cheby_basis * coeffs
        return cheby_y, true_y
    end
    
    r(x) = 1/(1+25*x^2)
    f(x) = x + 2*x^2 -exp(-x)
    t(x) = sqrt(abs(x))

##BEGIN QUESTIONS

	# use chebyshev to interpolate this:
	function q1(n)
		# function to approximate
        eps = 1e-9 # degree of approx
        #define Chebychev nodes on [-1,1]
        nodes = gausschebyshev(n)[1]
        # map nodes onto [-3,3]
        x_nodes = map((x) -> 3*x, nodes)
        #true y values
        y = map(f, x_nodes)
        # basis matrix
        V = ChebyBasis(n,n,nodes)
        # solve to find c
        c = V \ y
        #predict values
        cheby_y, true_y = predict(100,n,c,-3.0,3.0,f)
        #test max deviation < 1e-9
        @test maxabs(true_y - cheby_y) < eps
        #plots
        x_new = linspace(-3,3,100)
        plot1 = plot(x_new, [true_y, cheby_y], 
                    label = ["f(x)" "f_hat(x)"], linestyle = [:dot :solid], linecolor = ["black" "red"])
        plot2 = plot(x_new, true_y - cheby_y, label = "f(x) - f_hat(x)", yformatter = :scientific)
        return plot1, plot2
	end

	function q2(n)
		fc = Fun(f, Interval(-3,3))
        n_new = 100
        x_new = linspace(-3, 3, n_new)
        true_y = map(f, x_new)
        cheby_y = map(fc, x_new)
        plot1 = plot(x_new, [true_y, cheby_y], label = ["f(x)" "fc(x)"])
        plot2 = plot(x_new, true_y - cheby_y, yformatter = :scientific, label = "f(x) - fc(x)")
        return plot1, plot2
	end


	# plot the first 9 basis Chebyshev Polynomial Basis Fnctions
	function q3()
        T0(x) = 1
        T1(x) = x
        T2(x) = 2*x*T1(x) - T0(x)
        T3(x) = 2*x*T2(x) - T1(x)
        T4(x) = 2*x*T3(x) - T2(x)
        T5(x) = 2*x*T4(x) - T3(x)
        T6(x) = 2*x*T5(x) - T4(x)
        T7(x) = 2*x*T6(x) - T5(x)
        T8(x) = 2*x*T7(x) - T6(x)
        plot0 = plot(T0, labels = "T0")
        plot1 = plot(T1, labels = "T1")
        plot2 = plot(T2, labels = "T2")
        plot3 = plot(T3, labels = "T3")
        plot4 = plot(T4, labels = "T4")
        plot5 = plot(T5, labels = "T5")
        plot6 = plot(T6, labels = "T6")
        plot7 = plot(T7, labels = "T7")
        plot8 = plot(T8, labels = "T8")
        return plot0,plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8
	end


	function q4a()
        unodes = Array(Vector,3)
        i = 1
        for n in [6, 10, 16]
            unodes[i] = linspace(-1,1,n)
            i = i+1
        end
        chnodes = Array(Vector,3)
        i = 1
        for n in [6, 10, 16]
            chnodes[i] = gausschebyshev(n)[1]
            i = i+1
        end
        c = get_coeffs(-5.0, 5.0, unodes[1],r)
        
        y = Array{Vector}(3,2)
        j=1
        for nodes in [unodes, chnodes]
            i = 1
            for n in [6, 10, 16]
                y[i,j] = predict(100, n, get_coeffs(-5.0,5.0,nodes[i],r), -5.0, 5.0,r)[1]
                i = i+1
            end
            j=j+1
        end
        x_new = linspace(-5,5,100)
        true_y = r.(x_new)
        plot1 = plot(x_new, [true_y, y[1,1], y[2,1], y[3,1]],
        labels = ["Runge's function" "5 deg approx" "9 deg aprox." "15 deg approx."], title = "Uniformly spaced nodes")
        plot2 = plot(x_new, [true_y, y[1,2], y[2,2], y[3,2]],
        labels = ["Runge's function" "5 deg approx" "9 deg aprox." "15 deg approx."], title = "Chebyshev nodes")
        return plot1, plot2
    end

	function q4b()
        b = BSpline(13, 3, -5.0, 5.0)
        Basis = full(getBasis(collect(linspace(-5,5,65)),b))
        y = r.(linspace(-5,5,65))
        c = Basis \ y
        x_new = linspace(-5,5,65)
        true_y = r.(x_new)
        Basis1 = full(getBasis(collect(x_new),b))
        Bs_y1 = Basis1*c

        knots = vcat( collect(linspace(-5,-1.75,3)), collect(linspace(-1,1,7)), collect(linspace(1.75,5,3)) )
        b2 = BSpline(knots, 3)
        Basis2 = full(getBasis(collect(x_new),b2))
        y = r.(x_new)
        c2 = Basis2 \ y
        Basis3 = full(getBasis(collect(x_new),b2))
        Bs_y2 = Basis3*c2
        plot1 = plot(x_new, [true_y, Bs_y1, Bs_y2], 
            linecolor = [:auto "purple" "red"], line = [:solid :dot :dot], 
            labels = ["Runge's function" "Uniform knots" "Concentrated at zero"])
        plot2 = scatter(knots, zeros(15), labels = "My knot positions")
        plot2 = plot!(x_new, [true_y - Bs_y1, true_y-Bs_y2], labels = ["Uniform knots" "My knots"])
        return plot1, plot2
		
	end

	function q5()
        true_x = linspace(-1,1,65)
        true_y = map(t,true_x)
        plot0 = plot(true_x,true_y)
        b = BSpline(13, 3, -1.0, 1.0)
        Basis = full(getBasis(collect(linspace(-1,1,65)),b))
        y = t.(linspace(-1,1,65))
        c = Basis \ y
        x_new = linspace(-1,1,65)
        true_y = t.(x_new)
        Basis1 = full(getBasis(collect(x_new),b))
        Bs_y1 = Basis1*c

        #knots = vcat( collect(linspace(-1,0,7)), 0, collect(linspace(0,1,7)) )
        knots = vcat( collect(linspace(-1,-.25,3)), collect(linspace(-.1,.1,7)), collect(linspace(.25,1,3)) )
        b2 = BSpline(knots, 3)
        Basis2 = full(getBasis(collect(x_new),b2))
        y = t.(x_new)
        c2 = Basis2 \ y
        Basis3 = full(getBasis(collect(x_new),b2))
        Bs_y2 = Basis3*c2
        plot1 = plot(x_new, [Bs_y1, Bs_y2], 
            linecolor = [:auto "purple" "red"], line = [:solid :dot :dot], 
            labels = ["Uniform knots" "Concentrated at zero"])
        plot2 = scatter(knots, zeros(15), labels = "My knot positions")
        plot2 = plot!(x_new, [true_y - Bs_y1, true_y-Bs_y2], labels = ["Uniform knots" "My knots"])
        return plot0, plot1, plot2
		
	end


		# function to run all questions
	function runall()
		println("running all questions of HW-funcapprox:")
		p1a, p1b = q1(15)
		p2a, p2b = q2(15)
		p3a, p3b, p3c, p3d, p3e, p3f, p3g, p3h, p3i = q3()
		p4a_a, p4a_b = q4a()
		p4b_a, p4b_b = q4b()
		p5a, p5b, p5c = q5()
        l = @layout [ grid(1,2); grid(1,2); grid(1,9);
        grid(1,2); grid(1,2); grid(1,3) ]
        return plot(p1a, p1b, p2a, p2b, p3a, p3b, p3c, p3d, p3e, p3f, p3g, p3h, p3i,
                p4a_a, p4a_b, p4b_a, p4b_b, p5a, p5b, p5c, layout = l, yformatter = :scientific)
	end


end

