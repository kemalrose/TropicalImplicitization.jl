# Authors: Kemal Rose, Bernd Sturmfels and Simon Telen
# Date: Juli 14, 2023
# Short description: This code accompanies the paper
# "Tropical Implicitization Revisited" by Kemal Rose, Bernd Sturmfels and Simon Telen


module TropicalImplicitization

using Oscar, LinearAlgebra, Combinatorics



include("functions.jl")
export get_trop_A_disc, get_polytope_from_cycle, compute_A_discriminant, get_dual_fan_of_chow_polytope, get_tropical_cycle, get_vandermonde_matrix, sample, get_chow_fan, get_aux_data, getVertex, extract_polar_degrees, get_trop_con_Var, compute_polar_degrees
include("experimental_functions.jl")
export get_tropical_conormal_variety, getAk, get_cones_and_mults_from_lin_space_and_proj_with_lineality

end
