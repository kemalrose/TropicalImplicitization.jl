



using Oscar, LinearAlgebra
using TropicalImplicitization
using Combinatorics



A = [1 1 1 1 1 1; 0 1 0 1 0 1; 0 0 1 1 0 0; 0 0 0 0 1 1]
cone_list, weight_list = get_trop_con_Var(A)
Delta = get_polytope_from_cycle(cone_list, weight_list)


cone_list, weight_list = get_trop_con_Var(A)
dim.(cone_list)


function get_trop_con_Var(A)
    d, n = size(A)
    S = matrix_space(ZZ, d, n)
    B = nullspace(S(A))[2]
    B = Array(B)
    linear_space_mat =  [B  zeros(Int64, n, d);  zeros(Int64, d, n-d) LinearAlgebra.I[1:d, 1:d]]
    projection_mat = [zeros(Int64, n,n) transpose(A); LinearAlgebra.I[1:n, 1:n] transpose(A)]
    projection_mat = matrix_space(ZZ, size(projection_mat)...)(projection_mat)
    get_cones_and_mults_from_lin_space_and_proj(linear_space_mat, projection_mat)
end



function extract_polar_degrees(cone_list, weight_list)

    pdeg_list = []
        
    simplicial_inters_data = get_simplicial_inters_data(cone_list, weight_list)
    n = simplicial_inters_data.n
    d = simplicial_inters_data.d


    for polar_index in 1:d-1
        pdeg = 0
        dim1 = polar_index+1
        dim2 = d - dim1
        inds1 = collect(combinations(1:d,dim1-1))
        inds2 = collect(combinations(d+1:n,dim2-1))

        indices_linear_spaces = []
        for ind_list1 in inds1 
            for ind_list2 in inds2  
                current_inds = []
                append!(current_inds, ind_list1)
                append!(current_inds, ind_list2)
                push!(indices_linear_spaces, current_inds)
            end
        end



        w = rand(...)
        for inds in indices_linear_spaces
            E_inds = E[:, inds]
            σ += get_inters_multiplicity(simplicial_inters_data, E_inds, w)
        end



        push!(pdeg_list, pdeg)
    end
    σ_list
end

struct get_simplicial_inters_data
    n
    d
    RR
end

function get_simplicial_inters_data(cone_list, weight_list)

    n = ambient_dim(cone_list[1])
    d = dim(cone_list[1])
    @assert 2*d == n



    RR = []

    for j in 1:length(cone_list)
        sigma = cone_list[j]
    
        R = get_cone_rays(sigma)[:,1:end-1]
    
        SS, UU, VV = snf_with_transform(R)
    
        Rσ = inv(UU)[:,1:d]
        push!(RR, Rσ)




        fcts = collect(Oscar.facets(sigma))

        fctmtx = convert(Array{Int64,2},numerator.(transpose(hcat([Oscar.normal_vector(f) for f ∈ fcts]...))))
        T = matrix_space(QQ,size(fctmtx)...)
        RR = push!(RR,Rσ)
        Fσ = push!(Fσ,T(fctmtx))
    end



    get_simplicial_inters_data()
end




function get_inters_multiplicity(simplicial_inters_data, E_inds, w)
    n = simplicial_inters_data.n
    d = simplicial_inters_data.d

    does_intersect_list = [false for i in 1:simplicial_inters_data.nr_cnes]
    lattice_determinant_list = [0 for i in 1:simplicial_inters_data.nr_cnes]


    Ssq = matrix_space(QQ, n, n)

    lattice_mtx = zero(Ssq)

    for j in 1:simplicial_inters_data.nr_cnes

    end

    return does_intersect_list, lattice_determinant_list 
end

function get_inters_multiplicity(Fσ, lattice_mtx, E_inds, w, d)
    # σ is a d-dimensional cone, w a vector in the real numbers.
    # This function computes whether w+σ intersects the simplicial cone spanned by the columns of E_inds and the lattice multiplicity of the intersection.
    # Here Fσ is a matrix whose rows are the ray generators of the dual of σ modulo lineality.
    # lattice_mtx is an nxn matrix whose first dim(σ)


    
    does_intersect, lattice_determinant 
end








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
 
end

DTA = tropical_sampling_data(n, d, Γ, dets, RR, Fσ)
DTA

for i = 1:n
    for j = 1:size(data.Γ,1)
        if data.dets[j, i] != 0
            x = Array(data.Γ[j,i] * w_mat)
            if x[end]≥0
                vec = -data.Fσ[j]*(data.RR[j]*x[1:size(data.RR[j],2)])
                if all(vec.≥0)
                    is_contained[i, j] = 1
                    if sum(vec.==0) > 0 || x[end] == 0
                        is_in_face[i, j] = 1
                    end
                end
            end
        end
    end
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












cnes, weights = get_tropical_Conormal(A, 1)





function get_tropical_Conormal(A, k)
    Ak = A
    d, n = size(Ak)
    S = matrix_space(ZZ, d, n)
    B = nullspace(S(Ak))[2]
    B = Array(B)
    linear_space_mat =  [B  zeros(Int64, n, d);  zeros(Int64, d, n-d) LinearAlgebra.I[1:d, 1:d]]
    projection_mat = [LinearAlgebra.I[1:n, 1:n] transpose(A)]
    projection_mat = matrix_space(ZZ, size(projection_mat)...)(projection_mat)
    get_cones_and_mults_from_lin_space_and_proj(linear_space_mat, projection_mat)
end



# Computes the weighted linear projection of a tropicalized linear space. This is a tropical cycle
# -------------  Input:
# linear_space_mat     a matrix defining a linear space
# projection_mat       a matrix defining a linear projection
# -------------  Output:
# cone_list             a list of polyhedral cones. Need not form a fan.
# weight_list           a list of integer multiplicities. The balancing condition is satisfied.
function get_cones_and_mults_from_lin_space_and_proj(linear_space_mat, projection_mat)
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
