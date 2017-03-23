


include("src/funcapp.jl")
funcapp.runall()

using FactCheck

println("Here are the graphs for each question:")
approx1=funcapp.q1(100)
@fact abs(maximum(approx1[2]-funcapp.f(approx1[1])))--> less_than(1e-9)
plot=funcapp.q1_graph(approx1[2])
display(plot)

approx2=funcapp.q2()
plot2=funcapp.q2_graph(approx2)
display(plot2)

plot3=funcapp.q3()
display(plot3)

plot4a=funcapp.q4a((5,9,15),-5,5)
display(plot4a)
plot4b=funcapp.q4b(-5,5)
display(plot4b)

plot5=funcapp.q5()
display(plot5)
funcapp.q5()
