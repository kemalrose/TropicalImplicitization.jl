using TropicalImplicitization
using Test

@testset "TropicalImplicitization.jl" begin
    A = [1 1 1 1 1 1 1 1; 0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1]
    cone_list, weight_list = get_trop_A_disc(A)
    Delta = get_polytope_from_cycle(cone_list, weight_list)

    A = [1 1 1 1 1 1; 0 1 0 1 0 1; 0 0 1 1 0 0; 0 0 0 0 1 1]
    cone_list, weight_list = TropicalImplicitization.get_trop_con_Var(A)
    TropicalImplicitization.extract_polar_degrees(cone_list, weight_list)
    TropicalImplicitization.compute_polar_degrees(A)
end
