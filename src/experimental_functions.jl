







struct simplicial_inters_data
    n
    d
    n_cones
    RR
    Fσ
end


function get_tropical_conormal_variety(A, k = 1)
    m = size(A)[1]-1
    n = size(A)[2]-1
    
    Ak = getAk(A, k)
    S = matrix_space(ZZ, size(Ak)...)
    B = nullspace(S(Ak))[2]
    B = Array(B)

    linear_space_mat =  [B  zeros(Int64, n+1, m+1);  zeros(Int64, m+1, size(B,2)) LinearAlgebra.I[1:m+1, 1:m+1]]
    projection_mat = [zeros(Int64, n+1,n+1) -transpose(A); LinearAlgebra.I[1:n+1, 1:n+1] transpose(A)]
    projection_mat = matrix_space(ZZ, size(projection_mat)...)(projection_mat)
    get_cones_and_mults_from_lin_space_and_proj_with_lineality(linear_space_mat, projection_mat)
end



function compute_polar_degrees(A, k)
    cone_list, weight_list = get_tropical_conormal_variety(A, k)
    pdeg_list = extract_polar_degrees(cone_list, weight_list)

    EDdeg = sum(pdeg_list)
    first_ind_nonzero = findfirst(x->x!=0, pdeg_list)
    last_ind_nonzero = findfirst(x->x!=0, reverse(pdeg_list))
    deg_X = pdeg_list[first_ind_nonzero]
    deg_X_dual = reverse(pdeg_list)[last_ind_nonzero]


    println("The toric variety is of degree $deg_X.")
    println("The generic Euclidean distance degree of order $k is $EDdeg.")
    println("The dual variety of order $k is of degree $deg_X_dual and of codimension $last_ind_nonzero.")
    println("The polar degrees of order $k are $pdeg_list.")
end






function getAk(A, k=1)
    d, n = size(A)
    simplex = Oscar.simplex(d-1,k)
    latticePoints = Oscar.lattice_points(simplex)
    Ak = transpose(ones(Int64, n))
    for latticePoint in latticePoints
        newVector = transpose(ones(Int64, n))
        for indexNewVector in 1:n
            for rowVectorIndex in 2:d
                newVector[indexNewVector] *= (A[rowVectorIndex, indexNewVector]) ^ (latticePoint[rowVectorIndex-1].d)
            end
        end
        Ak = vcat(Ak, newVector)
    end
    return Ak[2:end,:]
end






function get_inters_multiplicities_linear_spaces(data, indices_linear_spaces, cone_list, weight_list, w)
    n = data.nß
    d = data.d
    n_cones_from_var = data.n_cones
    n_cones_from_lin = length(indices_linear_spaces)

    does_intersect_mat = fill(false, n_cones_from_var, n_cones_from_lin)
    lattice_mult_mat = zeros(Int64, n_cones_from_var, n_cones_from_lin)

    Ssq = matrix_space(QQ, n, n)
    lattice_mtx = zero(Ssq)
    E = Matrix{Int64}(identity_matrix(ZZ, n))


    for i in 1:n_cones_from_var
        Rσ = data.RR[i]
        lattice_mtx = Ssq(hcat(Array(Rσ), zeros(ZZ, n, n-d)))

        for j in 1:n_cones_from_lin

            indices_linear_space = indices_linear_spaces[j]

            for row_index in 1:n 
                for col_index in 1:n-d
                    lattice_mtx[row_index, d+col_index] = -E[row_index, indices_linear_space[col_index]]
                end
            end

            lattice_mult = Int64.(abs(det(lattice_mtx)) * weight_list[i])
            lattice_mult_mat[i, j] = lattice_mult

            if lattice_mult != 0
                x = Array(Oscar.solve(lattice_mtx, w; side =:right))
                if all(x[d+1:n] .≥ 0)
                    vector = data.RR[i]*x[1:d]
                    #@assert vector in cone_list[i]
                    #inner_prods_with_defining_inequalities = -data.Fσ[i]*(vector)
                    #if all(inner_prods_with_defining_inequalities.≥0)
                    if vector in cone_list[i]
                        does_intersect_mat[i, j] = 1
                    end
                end
            end
        end
    end
    total = sum(does_intersect_mat.*lattice_mult_mat)
    return does_intersect_mat, lattice_mult_mat, total
end




