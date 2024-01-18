# Authors: Kemal Rose, Bernd Sturmfels and Simon Telen
# Date: Juli 14, 2023
# Short description: This code accompanies the paper
# "Tropical Implicitization Revisited" by Kemal Rose, Bernd Sturmfels and Simon Telen


module TropicalImplicitization

using Oscar, LinearAlgebra

export get_trop_A_disc, get_polytope_from_cycle, compute_A_discriminant, get_dual_fan_of_chow_polytope, get_tropical_cycle, get_vandermonde_matrix, sample, get_chow_fan

include("functions.jl")

end
