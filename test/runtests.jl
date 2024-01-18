using TropicalImplicitization
using Test

@testset "TropicalImplicitization.jl" begin
    A = [1 1 1 1 1 1 1 1; 0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1]
    cone_list, weight_list = get_trop_A_disc(A)
    Delta = get_polytope_from_cycle(cone_list, weight_list)
end
