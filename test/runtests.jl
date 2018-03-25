using Base.Test

@testset "basics" begin
	@test funcapp.q1(15)[:error] .< 1e-9
end
