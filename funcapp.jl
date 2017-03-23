module HW_funapprox

using FastGaussQuadrature

ff(x::Float64)=x+2*x^2-exp(-x)
ff(x::Vector{Float64})=x+2*x.^2-exp(-x)

function question1(n)
  nodes = gausschebyshev(15)
  deg=n-1
  a=-3
  b=3
  normalized_nodes=0.5*(b+a)+0.5*(b-a)*nodes[1]
end

end
