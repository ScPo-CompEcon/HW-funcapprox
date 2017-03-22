# not very useful tests, just training

module IntTest
	using funcapp
	using Base.Test
	using Plots
		@testset "Do the questions return graphs?" begin
		@test typeof(funcapp.q1(15)) == Plots.Plot{Plots.PyPlotBackend}
		@test typeof(funcapp.q2(15)) == Plots.Plot{Plots.PyPlotBackend}
		@test typeof(funcapp.q3()) == Plots.Plot{Plots.PyPlotBackend}
		@test typeof(funcapp.q4a()) == Plots.Plot{Plots.PyPlotBackend}
		@test typeof(funcapp.q4b()) == Plots.Plot{Plots.PyPlotBackend}
		@test typeof(funcapp.q5()) == Plots.Plot{Plots.PyPlotBackend}
	end
end
