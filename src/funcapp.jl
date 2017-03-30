## Sorry for the delay. I still don’t know how to “functionize” this, so here’s just the code for each question ##


using Plots
using FastGaussQuadrature
using Gadfly
using ApproxFun ## this prints LoadError "syntax: numeric constant 1. cannot be 			##implicitely multiplied because it ends with".""
using ApproXD

Q1/

a=-3
b=3
n=15
J = n-1 #df approximation
f(x)=x+2x^2-exp(-x)

nodes = gausschebyshev(n) #compute n cheby nodes

nodes[1] #collect cheby nodes
x=(1/2)*(a+b) + (1/2)*(b-a)*nodes[1] #map z into x
y=f.(x)
z = 2.*(x-a)./(b-a) - 1 #for practicality later
#F=map(f(x), x) #get function value at x

Gadfly.plot(x=x, y=y, Geom.line) #plot f(x)

#PHI=zeros(15,15)

#for i in 1:size(PHI)[1], j in 1:size(PHI)[2]
  #  PHI[i,j] = cos((n-z[i]+0.5)*(j-1)*pi/n)  ## why does not this work ?
#end


phi = [cos(acos(z[i])j) for i = 1:15, j = 0:J] #building basis function matrix 


c = fill(.0, n)
c= \(phi,y) #vector of coefficients
Gadfly.plot(x=y, y=c, Geom.line)  #plot coeeficients for fun

n_new = 100 #predict 100 new nodes
nodes_new = linspace(a,b,n_new)
z_new = 2.*(nodes_new-a)./(b-a) - 1
phi_new = [cos(acos(z_new[i])*j) for i = 1:n_new, j = 0:J]

#for i in 1:size(PHI)[1], j in 1:size(PHI)[2]
#phi_new= cos((n-z_new[i]+0.5)*(j-1)*pi/n)
#end

f_hat = phi_new*c
y = f.(nodes_new)

diff = y - f_hat
plot(diff)
graph1 = Gadfly.plot(layer(x=nodes_new, y=y, Geom.line),layer(x=nodes_new,y=f_hat, Geom.point), layer(x=nodes_new, y=diff, Geom.line))


Q2/
## this suddenly stopped working as my computer seems not to be able to pre-compile ApproxFun anymore.. 

  f(x) = x + 2*x^2 - exp(-x)
a = -3
b = 3
FUN = Fun(f, Chebyshev(a..b))
n_new = 100
nodes_new2 = linspace(a,b,n_new)
FUN_hat = FUN.(nodes_new2)

 n_new = 100
x_new = linspace(a,b,n_new)
y_predict=f_approx.(x_new)
y_dev = abs(truef1.(x_new) - y_predict)

Q3/

J = 9
a = -3
b = 3

poly_cheb = (x, j) -> ChebyT(unitmap(x,a,b),j-1)

Plots.plot([x -> poly_cheb(x, j) for j in 1:J], a, b, layout = (3,3))

Q4/

I couldn't do question 4.

Q5/

function myknots(k)
if isodd(k)==true
vcat(-1, -.8, collect(linspace(-.5, -.1, 5-(k+1)/2)), zeros(k), collect(linspace(.1, .5, 5-(k+1)/2)), .8, 1)
else
vcat(-1, -.8, collect(linspace(-.5, -.1, 5-k/2)), zeros(k), collect(linspace(.1, .5, 4-k/2)), .8, 1)
end
end

a=-1
b=1
p = 13
k = 3
eval=65
g(x)= abs(x)^(0.5)
Plots.plot(g)

bs_uni = BSpline(p,k,a,b) #get splines
bs_kink= BSpline(myknots(k), k)

B_uni =full(getBasis(collect(linspace(a,b,eval)),bs_uni)) #get basis functions
B_kink =full(getBasis(collect(linspace(a,b,eval)),bs_kink))

x= collect(linspace(a,b,eval)) #get coefficients
y=g.(x)
c_uni= \(B_uni,y)
c_kink= \(B_kink,y)

new_points = collect(linspace(a,b,1000))
g_new_points = g.(new_points) #true function
g_hat_uni = getBasis(new_points,bs_uni) * c_uni
g_hat_kink = getBasis(new_points,bs_kink) * c_kink

Plots.plot(new_points,g_new_points) #plot true function
Plots.plot(new_points,g_hat_uni) #plot uni approx
Plots.plot(new_points,g_hat_kink) #plot kink approx

graph2 = Gadfly.plot(layer(x=new_points, y=g_new_points, Geom.line),layer(x=new_points,y=g_hat_uni, Geom.point), layer(x=new_points, y=g_hat_kink, Geom.point))
