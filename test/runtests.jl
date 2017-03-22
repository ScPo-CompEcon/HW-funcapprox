# not very useful tests, just training

module IntTest
	using funcapp
	using Base.Test
	using Plots
		@testset "Do the questions return graphs?" begin
		@test typeof(funcapp.q(1)) == Plots.Plot{Plots.PyPlotBackend}
		@test typeof(funcapp.q(2)) == Plots.Plot{Plots.PyPlotBackend}
		@test typeof(funcapp.q(3)) == Plots.Plot{Plots.PyPlotBackend}
		@test typeof(funcapp.q(4)) == Plots.Plot{Plots.PyPlotBackend}
		@test typeof(funcapp.q(5)) == Plots.Plot{Plots.PyPlotBackend}
	end
end 