

    
using Oscar, LinearAlgebra
using TropicalImplicitization
using Combinatorics



    
# Computes the weighted linear projection of a tropicalized linear space. This is a tropical cycle
# -------------  Input:
# linear_space_mat     a matrix defining a linear space
# projection_mat       a matrix defining a linear projection
# -------------  Output:
# cone_list             a list of polyhedral cones. Need not form a fan.
# weight_list           a list of integer multiplicities. The balancing condition is satisfied.
function get_cones_and_mults_from_lin_space_and_proj_with_lineality(linear_space_mat, projection_mat)
    matroid = Oscar.Polymake.matroid.Matroid(VECTORS = Matrix{Int64}(linear_space_mat))
    ΣB = Oscar.Polymake.tropical.matroid_fan{min}(matroid)
    rayMat = ΣB.RAYS
    rayMat = convert(Array{Int64,2},rayMat)
    rayMat = matrix_space(QQ,size(rayMat)...)(rayMat)
    lineality = ΣB.LINEALITY_SPACE
    lineality = convert(Array{Int64,2},lineality)
    lineality = matrix_space(QQ,size(lineality)...)(lineality)
    lineality = vcat(lineality, -lineality)[:,2:end]
    one_vec = lineality[1:1,:].=1
    lineality = vcat(lineality, one_vec, -one_vec)

    #ones_vec = transpose(ones(eltype(lineality), size(lineality, 2)))
    #lineality = vcat(lineality, ones_vec)
    maxpols = ΣB.MAXIMAL_POLYTOPES
    cone_list = []
    for j in 1:size(maxpols,1)
        R = rayMat[findall(ℓ->ℓ,collect(maxpols[j,:]))[1:end-1],2:end]
        M = vcat(R, lineality)
        σ = positive_hull(M)
        push!(cone_list,σ)
    end
    d = maximum(dim.(cone_list))
    weight_list = ones(Int64, length(cone_list))
    project_tropical_variety(cone_list, weight_list, projection_mat)
end



function get_trop_con_Var(A)
    d, n = size(A)
    S = matrix_space(ZZ, d, n)
    B = nullspace(S(A))[2]
    B = Array(B)
    linear_space_mat =  [B  zeros(Int64, n, d);  zeros(Int64, d, n-d) LinearAlgebra.I[1:d, 1:d]]
    projection_mat = [zeros(Int64, n,n) -transpose(A); LinearAlgebra.I[1:n, 1:n] transpose(A)]
    projection_mat = matrix_space(ZZ, size(projection_mat)...)(projection_mat)
    get_cones_and_mults_from_lin_space_and_proj_with_lineality(linear_space_mat, projection_mat)
end






function extract_polar_degrees(cone_list, weight_list)
    data = get_simplicial_inters_data(cone_list, weight_list)
    n = data.n
    d = data.d
    pdeg_list = Array{Int64}([])

    for polar_index in 1:d-2
        pdim1 = polar_index
        pdim2 = d - pdim1 
        inds1 = collect(combinations(1:d,pdim1))
        inds2 = collect(combinations(d+1:n,pdim2))

        indices_linear_spaces = []
        for ind_list1 in inds1 
            for ind_list2 in inds2  
                current_inds = []
                append!(current_inds, ind_list1)
                append!(current_inds, ind_list2)
                push!(indices_linear_spaces, current_inds)
            end
        end

        w = matrix_space(QQ, n, 1)(rand(-100000:100000, n))

        does_intersect_mat, lattice_mult_mat, total = get_inters_multiplicities(data, indices_linear_spaces, cone_list, weight_list, w)
        push!(pdeg_list, total)
    end
    pdeg_list
end



struct simplicial_inters_data
    n
    d
    n_cones
    RR
    Fσ
end


function get_simplicial_inters_data(cone_list, weight_list)

    n = ambient_dim(cone_list[1])
    d = dim(cone_list[1])
    @assert 2*d == n

    RR = []
    Fσ = []
    n_cones = length(cone_list)

    for j in 1:length(cone_list)
        sigma = cone_list[j]
    
        R = get_cone_rays(sigma)[:,1:end-1]
    
        SS, UU, VV = snf_with_transform(R)
    
        Rσ = inv(UU)[:,1:d]
        RR = push!(RR,Rσ)


        fcts = collect(Oscar.facets(sigma))
        fctmtx = convert(Array{Int64,2},numerator.(transpose(hcat([Oscar.normal_vector(f) for f ∈ fcts]...))))
        T = matrix_space(QQ,size(fctmtx)...)
        Fσ = push!(Fσ,T(fctmtx))
    end

    simplicial_inters_data(n,d,n_cones,RR,Fσ)
end



function get_inters_multiplicities(data, indices_linear_spaces, cone_list, weight_list, w)
    n = data.n
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
        lattice_mtx = Ssq(hcat(Array(Rσ),zeros(ZZ, n, d)))

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
                x = Array(Oscar.solve(lattice_mtx, w))
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







function get_aux_data(cone_list, weight_list)
    n = ambient_dim(cone_list[1])
    d = dim(cone_list[1])

    Ssq = matrix_space(QQ, n, n)

    E = Matrix{Int64}(identity_matrix(ZZ, n))

    Γ = fill(zero(Ssq),length(cone_list),n)
    dets = fill(zero(QQFieldElem),length(cone_list),n)
    RR = []
    Fσ = []
    for j in 1:length(cone_list)
        sigma = cone_list[j]

        R = get_cone_rays(sigma)[:,1:end-1]

        SS, UU, VV = snf_with_transform(R)

        Rσ = inv(UU)[:,1:d]

        for i = 1:n
            mtx = Ssq(hcat(Array(Rσ),-E[:,i]))
            dmtx = det(mtx)
            if dmtx !=0
                Γ[j,i] = inv(mtx)
                dets[j,i] = abs(dmtx)*weight_list[j]
            end
        end
        fcts = collect(Oscar.facets(sigma))

        fctmtx = convert(Array{Int64,2},numerator.(transpose(hcat([Oscar.normal_vector(f) for f ∈ fcts]...))))
        T = matrix_space(QQ,size(fctmtx)...)
        RR = push!(RR,Rσ)
        Fσ = push!(Fσ,T(fctmtx))
    end

    DTA = tropical_sampling_data(n, d, Γ, dets, RR, Fσ)
    DTA
end





function compute_polar_degrees(A)
    cone_list, weight_list = get_trop_con_Var(A)
    pdeg_list = extract_polar_degrees(cone_list, weight_list)

    EDdeg = sum(pdeg_list)
    first_ind_nonzero = findfirst(x->x!=0, pdeg_list)
    last_ind_nonzero = findfirst(x->x!=0, reverse(pdeg_list))
    deg_X = pdeg_list[first_ind_nonzero]
    deg_X_dual = reverse(pdeg_list)[last_ind_nonzero]


    println("The degree of the toric variety is $deg_X and its Euclidean distance degree is $EDdeg.")
    println("The degree of the dual variety is $deg_X_dual and its codimension is $last_ind_nonzero.")
    println("The polar degrees are $pdeg_list.")
end






A = [1 1 1 1 1 1; 0 1 0 1 0 1; 0 0 1 1 0 0; 0 0 0 0 1 1]
cone_list, weight_list = get_trop_con_Var(A)
dim.(cone_list)

extract_polar_degrees(cone_list, weight_list)

compute_polar_degrees(A)

