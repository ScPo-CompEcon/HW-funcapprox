using Base.Test

@testset "basics" begin
	n = 15
	for i in 1:length(funcapp.q1(n)[:error])
		@test funcapp.q1(n)[:error][i] < 1e-6
	end
end
