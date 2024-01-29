# Authors: Kemal Rose, Bernd Sturmfels and Simon Telen
# Date: Juli 14, 2023
# Short description: This code accompanies the paper
# "Tropical Implicitization Revisited" by Kemal Rose, Bernd Sturmfels and Simon Telen
# arXiv: ...


# Computes the stable intersection of inner normal fans of polytopes.
# -------------  Input:
# pol_list     a list of polytopes
# -------------  Output:
# cone_list    a List of polyhedral cones
# weight_list  a List of multiplicities

function stable_inters(pol_list::Vector{})
    cone_list = []
    weight_list = []
    @assert length(unique(ambient_dim.(pol_list))) == 1
    n = ambient_dim(pol_list[1])
    d = n - length(pol_list)
    MKsum = sum(pol_list)
    Sigma = normal_fan(MKsum)
    d_cones = cones(Sigma, d)
    if d == 0
        d_cones = [positive_hull(zeros(QQ,n))]
    end
    vert_lists = vertices.(pol_list)
    for sigma in d_cones
        cone_rays = get_cone_rays(sigma)
        inter_pt = sum([cone_rays[:,i] for i in 1:size(cone_rays, 2)])
        opt_vals = [minimum([dot(inter_pt, vert) for vert in vert_lists[i]]) for i in 1:length(vert_lists)]
        revealed_face_verts = [ filter(vert->dot(vert,inter_pt)==opt_vals[i], vert_lists[i]) for i in 1:length(vert_lists)]
        cone_weight = general_mixed_vol(revealed_face_verts)
        if cone_weight > 0
            push!(cone_list, sigma)
            push!(weight_list, cone_weight)
        end
    end
    cone_list, weight_list
end

# Computes generators of a polyhedral cone
# -------------  Input:
# sigma        a polyhedral cone
# -------------  Output:
# cone_rays    an Oscar matrix, the column vectors generate sigma.

function get_cone_rays(sigma)
    n = ambient_dim(sigma)
    cone_rays = hcat(rays_modulo_lineality(sigma)[1]...,lineality_space(sigma)...,-lineality_space(sigma)...,zeros(ZZ, n))
    for i in 1:ncols(cone_rays)
        cone_rays[:,i] = cone_rays[:,i].*lcm(denominator.(cone_rays[:,i]))
    end
    cone_rays =  numerator.(cone_rays)
    cone_rays = Oscar.matrix(cone_rays)
    cone_rays
end


# Computes the mixed volume of polytopes in an appropriate lattice.
# -------------  Input:
# vert_lists        a list of k lists. Each list comprises the vertices of a polytope in Q^n
# -------------  Output:
# vol         the mixed volume of the polytopes, computed in the lattice defined by their affine hull.
function general_mixed_vol(vert_lists)
    n = length(vert_lists[1][1])
    vert_lists = [[vert-vert_list[1] for vert in vert_list] for vert_list in vert_lists]
    face_matrix = hcat([ hcat(verts...) for verts in vert_lists]...)
    face_matrix_Oscar = Oscar.matrix(ZZ.(face_matrix))
    D, U, V = snf_with_transform(face_matrix_Oscar)
    d = rank(D)
    projection_matrix = U[1:d,:]
    revealed_face_verts_projected = [[projection_matrix*collect(vert) for vert in verts] for verts in vert_lists]
    revealed_faces = convex_hull.(revealed_face_verts_projected)
    vol = QQ(Polymake.polytope.mixed_volume( [ pol.pm_polytope for pol in revealed_faces ]...  ))
    @assert denominator(vol) == 1
    vol = numerator(vol)
    vol
end

# Computes the lattice index of a lattice after linear projection.
# -------------  Input:
# sigma     a polyhedral cone in R^n.
# A         a k × n matrix.
# -------------  Output:
# mlt       the lattice index of the lattice A * (span(sigma) ∩ Z^n) in the lattice  (A * span(sigma)) ∩ Z^k.
function get_mult(sigma, A)
    @assert(typeof(A) <: ZZMatrix)
    B = get_cone_rays(sigma)
    BasisB = get_lattice_basis(B)
    S, U, V = snf_with_transform(A * BasisB)
    d = rank(S)
    mlt = prod([S[i,i] for i in 1:d])
    mlt
end

# Computes a lattice basis of a linear space
# -------------  Input:
# A         an integer n × k matrix of rank d.
# -------------  Output:
# Mat       an integer matrix such that (A * R^k) ∩ Z^n = Mat * Z^d
function get_lattice_basis(A)
    @assert(typeof(A) <: ZZMatrix)
    A_Oscar = A
    D, U, V = snf_with_transform(A_Oscar)
    d = rank(D)
    Mat = inv(U)[:,1:d]
    Mat
end

# Computes a linear projection of a tropical cycle
# -------------  Input:
# cone_list             a list of polyhedral cones. Need not form a fan.
# weight_list           a list of integer multiplicities. The balancing condition is satisfied.
# projection_mat        an integer matrix
# -------------  Output:
# resulting_cones       a list of polyhedral cones. Their union is the image of the tropical cycle.
# resulting_weights     a list of integer multiplicities. They define the multiplicities of the projected cycle according to the tropical pushforward formula.
function project_tropical_variety(cone_list, weight_list, projection_mat)
    @assert(typeof(projection_mat) <: ZZMatrix)
    resulting_cones = []
    resulting_weights = []
    for i in 1:length(cone_list)
        sigma = cone_list[i]
        weight = weight_list[i]
        new_rays = projection_mat * get_cone_rays(sigma)
        new_cone = positive_hull(transpose(new_rays))
        new_weight = weight * get_mult(sigma, projection_mat)
        push!(resulting_cones, new_cone)
        push!(resulting_weights, new_weight)
    end
    d = maximum(dim.(resulting_cones))
    inds = findall(sigma-> dim(sigma) == d, resulting_cones)
    resulting_cones = resulting_cones[inds]
    resulting_weights = resulting_weights[inds]
    resulting_cones, resulting_weights
end

# From a weighted codimension one fan compute the data that is needed to build a vertex oracle.
# -------------  Input:
# cone_list             a list of polyhedral cones. Need not form a fan.
# weight_list           a list of integer multiplicities. The balancing condition is satisfied.
# -------------  Output:
# DTA       an object of type "tropical_sampling_data".

function get_aux_data(cone_list, weight_list)
    n = ambient_dim(cone_list[1])
    d = dim(cone_list[1])

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

# An object comprises the data that is needed to build a vertex oracle from a codimension one tropical cycle.
# -------------  Fields:
#     n          ambient dimension of the cycle
#     d          dimension of the tropical cycle
#     dets       lattice multiplicities for the intersections of the cone with a ray e_i
#     RR         a ist of lattice bases for all lattices generated by cones
#     Fσ         a list of defining facet inequalities for each cone
#     Γ          inverses of lattice bases appearing in RR when topped up with base vectors e_i
struct tropical_sampling_data
    n
    d
    Γ
    dets
    RR
    Fσ
end

# From a weight vector w and a tropical cycle
# list the incidence of intersections w+pos(e_i) with cones in the cycle.
# -------------  Input:
# w                 an element of R^n
# data              an object of type "tropical_sampling_data" representing a tropical cycle
# -------------  Output:
# is_contained    incidence Matrix of intersections w+pos(e_i) ∩ sigma. Here sigma is a maximal cone in the tropical cycle, 1 ≤ i ≤ n.
# is_in_face      incidence Matrix of intersections w+pos(e_i) ∩ boundary(sigma). Here sigma is a maximal cone in the tropical cycle, 1 ≤ i ≤ n.
# is_interiour    incidence Matrix of intersections w+pos(e_i) ∩ interior(sigma). Here sigma is a maximal cone in the tropical cycle, 1 ≤ i ≤ n.
function cone_containments(w, data)

    w_mat = matrix_space(QQ,data.n,1)(w)
    n = data.n
    is_contained = zeros(Int64,n, size(data.Γ, 1))
    is_in_face = zeros(Int64, n, size(data.Γ, 1))
    is_interiour = zeros(Int64, n, size(data.Γ, 1))

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
    is_contained, is_in_face, is_interiour
end






#Computes the tropicalization of the image of a generic map of Laurent polynomials with prescribed Newton polytopes.
# -------------  Input:
# newton_pols             a list of k polytopes in R^n.
# -------------  Output:
# resulting_cones       a list of polyhedral cones.
# resulting_weights     a list of integer multiplicities.
"""
get_tropical_cycle(newton_pols)\\
Computes the tropicalization of the image of a generic map of Laurent polynomials with prescribed Newton polytopes.
# Examples
```jldoctest
julia> newton_pols = [convex_hull([[2], [4]]), convex_hull([[0], [8]])];
julia> cone_list, weight_list = get_tropical_cycle(newton_pols)
(Any[Polyhedral cone in ambient dimension 2, Polyhedral cone in ambient dimension 2, Polyhedral cone in ambient dimension 2, Polyhedral cone in ambient dimension 2], Any[2, 2, 8, 4])
```
"""
function get_tropical_cycle(newton_pols::Vector)
    m = ambient_dim(newton_pols[1])

    n = length(newton_pols)

    vert_list = [collect(Array.(vertices(pol))) for pol in newton_pols]
    new_vert_list = [push!([ vcat(vert, zeros(ZZ,n)) for vert in verts  ]) for verts in vert_list]
    basis_vecs = []
    for i in 1:length(vert_list)
        ei = zeros(ZZ, m+n)
        ei[m+i] = 1
        push!(basis_vecs, ei)
        push!(new_vert_list[i],ei)
    end
    A = transpose(hcat(basis_vecs...))
    A = matrix_space(ZZ, size(A)...)(A)
    pol_list = convex_hull.(new_vert_list)
    cone_list, weight_list = stable_inters(pol_list)
    resulting_cones, resulting_weights = project_tropical_variety(cone_list, weight_list, A)
    resulting_cones, resulting_weights
end

# From a list of vertices and a vertex oracle for a polytope reconstruct the polytope
function get_Polytope(V_0, getV)
    V_old, F_old = [], []
    V_new  = copy(V_0)
    F_new = []
    P = convex_hull(transpose(hcat(V_0...)))

    while size(V_old) != size(V_new)
        V_old = copy(V_new)
        F_old = copy(F_new)


        V_new, F_new, P = update(V_old, F_old, getV)
    end
    return P
end

function update(V_old, F_old, getV)
    V_new = copy(V_old)
    P = convex_hull(transpose(hcat(V_old...)))
    F_new  = [normal_vector(f) for f in facets(P)]
    keep_facets = []

    for a in setdiff(F_new, F_old)
        vert = getV(a)
        if !in(vert, V_new)

            V_new = push!(V_new, vert)
        else

            index = findfirst(ℓ -> ℓ == a, F_new)
            keep_facets = push!(keep_facets, index)
        end
    end
    copy(V_new), [copy(F_old); F_new[keep_facets]], P
end


# From a vertex oracle compute a polytope
# -------------  Input:
# data              an object of type "tropical_sampling_data" representing a tropical cycle of codimension one
# -------------  Output:
# monomials         the vertices of a polytope
function newton_pol(data::tropical_sampling_data)
    n = data.n
    V_ambient = sampleRandom(data, 2*n)

    function getV(w)
        vtx = getVertex(-w, data)
        vtx
    end
    V_0 = V_ambient
    P = get_Polytope(V_0, getV);
end

# From a vertex oracle sample random vertices of a polytope
# -------------  Input:
# data              an object of type "tropical_sampling_data" representing a tropical cycle of codimension one
# -------------  Output:
# monomials         random vertices of a polytope
function sampleRandom(sampling_data, nsamples)
    n = sampling_data.n
    monomials = []
    for i = 1:nsamples
        w = rand(-100000:100000,n)
        monomial = getVertex(w, sampling_data)
        monomials = unique!(push!(monomials,monomial))
    end
    monomials
end


# Construct a vertex oracle for a polytope from a tropical codimension one cycle
# -------------  Input:
# w                 a generic weight vector in R^n
# data              an object of type "tropical_sampling_data" representing a tropical cycle of codimension one
# -------------  Output:
# monomials         the extremal vertex of P in direction w
function getVertex(w, data)
    n = data.n

    monomial = zeros(QQFieldElem,n)
    is_contained, is_in_face = cone_containments(w, data)
    is_generic = sum(is_in_face) == 0

    if !(is_generic)

        w_new = 5000 * w + rand(-100:100, n)
        is_contained_new, is_in_face_new = cone_containments(w_new, data)
        while sum((is_contained_new - is_contained).>0) != 0
            w_new = 20000 * w + rand(-1000:1000, n)
            is_contained_new, is_in_face_new = cone_containments(w_new, data)
        end
        is_contained = is_contained_new
    end

    for i = 1:n
        for j = 1:size(data.Γ,1)
            if data.dets[j,i] != 0
                if is_contained[i, j] == 1
                    monomial[i] += data.dets[j,i]
                end
            end
        end
    end
    monomial
end


# Construct a Newton polytope from a tropical codimension one cycle
# -------------  Input:
# cone_list             a list of polyhedral cones. Need not form a fan.
# weight_list           a list of integer multiplicities. The balancing condition is satisfied.
# -------------  Output:
# pol                   a Polytope with weighted normal fan equal to the tropical cycle
"""
get_polytope_from_cycle(newton_pols)\\
From a tropical hypersurface with trivial valuation computes the, unique up to translation, Newton Polytope that is dual.
# Examples
```jldoctest
julia> A = [1 1 1 1 1 1; 2 3 5 7 11 13; 13 8 5 3 2 1];
julia> cone_list, weight_list = get_trop_A_disc(A);
julia> Delta = get_polytope_from_cycle(cone_list, weight_list)
Polyhedron in ambient dimension 6
```
"""
function get_polytope_from_cycle(cone_list, weight_list)
    n = ambient_dim(cone_list[1])
    @assert maximum(dim.(cone_list)) == n-1
    sampling_data = get_aux_data(cone_list, weight_list)
    pol = newton_pol(sampling_data)
    pol
end




# Computes the tropicalization of an A discriminant
# -------------  Input:
# A                     a matrix
# -------------  Output:
# cone_list             a list of polyhedral cones. Need not form a fan.
# weight_list           a list of integer multiplicities. The balancing condition is satisfied.
"""
get_trop_A_disc(A)\\
Computes the tropicalization of an A-discriminant.
# Examples
```jldoctest
julia> A = [1 1 1 1 1 1; 2 3 5 7 11 13; 13 8 5 3 2 1];

julia> cone_list, weight_list = get_trop_A_disc(A)
(Any[Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6, Polyhedral cone in ambient dimension 6], Any[1, 15, 2, 6, 2, 2, 1, 2, 1, 9, 1, 1, 1, 1, 1])
```
"""
function get_trop_A_disc(A)
    d, n = size(A)
    S = matrix_space(ZZ, d, n)
    B = nullspace(S(A))[2]
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
    matroid = Polymake.matroid.Matroid(VECTORS = Matrix{Int64}(linear_space_mat))
    ΣB = Polymake.tropical.matroid_fan{min}(matroid)
    rayMat = ΣB.RAYS
    rayMat = convert(Array{Int64,2},rayMat)
    rayMat = matrix_space(QQ,size(rayMat)...)(rayMat)
    lineality = ΣB.LINEALITY_SPACE
    lineality = convert(Array{Int64,2},lineality)
    lineality = matrix_space(QQ,size(lineality)...)(lineality)
    lineality = vcat(lineality, -lineality)[:,2:end]
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


# Given a tropical cycle, forms a codimension one tropical cycle in one more variable. This is done by homogenizing the fan and forming the stable sum with the negative of a tropical linear space.
# -------------  Input:
# cone_list             a list of polyhedral cones. Need not form a fan.
# weight_list           a list of integer multiplicities. The balancing condition is satisfied.
# -------------  Output:
# resulting_cone_list      a list of polyhedral cones. Need not form a fan.
# resulting_weight_list    a list of integer multiplicities. The balancing condition is satisfied.
function get_chow_fan(cone_list, weight_list)
    cone_list_proj, weight_list_proj = homogenize_tropical_cycle(cone_list,weight_list)
    d = maximum(dim.(cone_list))
    n = ambient_dim(cone_list[1])
    cone_list_lin, weight_list_lin = get_negative_projective_linear_cycle(n,n-d-1)
    resulting_cone_list, resulting_weight_list = stable_sum(cone_list_proj, weight_list_proj, cone_list_lin, weight_list_lin)
    return resulting_cone_list, resulting_weight_list
end

function get_negative_projective_linear_cycle(n,d)
    rayList = LinearAlgebra.I[1:n+1, 1:n+1]
    lineality = [transpose(ones(Int64,n+1)); -transpose(ones(Int64,n+1))]

    maxpols = [collect(st) for st in subsets(Set(1:n+1), d)]
    cone_list = []
    for pol in maxpols
        cone = positive_hull([-rayList[pol,:]; lineality])
        push!(cone_list, cone)
    end
    weight_list = ones(Int64, length(cone_list))
    cone_list, weight_list
end


function homogenize_tropical_cycle(cone_list, weight_list)
    n = ambient_dim(cone_list[1])
    #lin_space = hcat(ones(Int64,n+1), -ones(Int64,n+1))
    #lin_space = matrix_space(ZZ,size(line_space)...)(lin_space)
    resulting_cone_list = []
    resulting_weight_list = weight_list
    for sigma in cone_list
        ray_list = get_cone_rays(sigma)
        zero_vec = zero(matrix_space(ZZ, 1 , size(ray_list,2)))
        ray_list = vcat(zero_vec, ray_list)
        sigma_new = positive_hull(transpose(ray_list))
        push!(resulting_cone_list, sigma_new)
    end
    resulting_cone_list, resulting_weight_list
end


# Given two tropical cycles, forms their stable sum.
# -------------  Input:
# cone_list1             a list of polyhedral cones. Need not form a fan.
# weight_list1           a list of integer multiplicities. The balancing condition is satisfied.
# cone_list2             a list of polyhedral cones. Need not form a fan.
# weight_list2           a list of integer multiplicities. The balancing condition is satisfied.
# -------------  Output:
# cone_list      a list of polyhedral cones. Need not form a fan.
# weight_list    a list of integer multiplicities. The balancing condition is satisfied.
function stable_sum(cone_list1, weight_list1, cone_list2, weight_list2)
    cone_list = []
    weight_list = []
    for i in eachindex(cone_list1)
        for j in eachindex(cone_list2)
            sigma1 = cone_list1[i]
            sigma2 = cone_list2[j]
            weight1 = weight_list1[i]
            weight2 = weight_list2[j]
            A1 = get_cone_rays(sigma1)
            A2 = get_cone_rays(sigma2)

            new_weight = get_index_stable_sum(A1, A2)
            if new_weight != 0
                new_cone = positive_hull(transpose([A1 A2]))
                push!(cone_list, new_cone)
                push!(weight_list, new_weight * weight1 * weight2)
            end
        end
    end
    cone_list, weight_list
end

function get_index_stable_sum(A1, A2)
    M1 = get_lattice_basis(A1)
    M2 = get_lattice_basis(A2)
    D, U, V = snf_with_transform([M1 M2])
    abs(det(  (U * [M1 M2])[1:minimum(size(D)),1:minimum(size(D))]  ))
end



# Computes the coefficients and monomials of an A-discriminant and its Newton polytope over the field Fld. Is experimental, possibly slow.
# -------------  Input:
# A               a matrix
# Delta           the Newton polytope of the A discriminant belonging to A
# Fld             a field, e.g Fld = QQ, GF(101),...
# -------------  Output:
# mons            the lattice points of Delta
# coeffs          the coefficients of the A-discriminant
"""
compute_A_discriminant(A)\\
Computes the coefficients and monomials of an A-discriminant and its Newton polytope over the field Fld. Is experimental, possibly slow.
# Examples
```jldoctest
julia> A = [1 1 1 1 1; 0 1 0 4 1; 0 0 1 1 4];
julia> cone_list, weight_list = get_trop_A_disc(A);
julia> Delta = get_polytope_from_cycle(cone_list, weight_list);
julia> Fld = QQ;
julia> mons, coeffs = compute_A_discriminant(A, Delta, Fld);
```
"""
function compute_A_discriminant(A, Delta, Fld)
    A_Oscar = matrix_space(ZZ, size(A)...)(A)
    (d,n) = size(A)
    B = Array(nullspace(A_Oscar)[2])
    println("Collect all lattice points.")
    mons = collect(lattice_points(Delta))
    pts = []
    n_samples = round(length(mons)*1.2)
    println("Sample $n_samples points from Discriminant")
    for j in 1:n_samples
        v = ZZRingElem.(rand(-50:50,n-d))
        λ = ZZRingElem.(rand(-50:50,d))
        pt = horn_param(A, B, v, λ)

        push!(pts,pt)
    end
    println("Construct Vandermonde matrix.")
    Mat_Space = Oscar.matrix_space(Fld, length(pts), length(mons))
    Vand = zero(Mat_Space)
    for i in 1:length(pts)
      Vand[i, :] = [prod(pts[i].^convert(Vector{ZZRingElem}, mon)) for mon in mons]
    end
    println("Compute coefficients of the Discriminant.")
    null_dim, coeffs = nullspace(Vand)
    mons, coeffs
end

# Samples points from the A discriminant
# -------------  Input:
# A               a matrix of size d times n of full rank.
# B               Gale dual to A
# v               a vector of size n-d
# λ               a vector of size d
# -------------  Output:
# vec            an element of the A discriminant
function horn_param(A, B, v, λ)
    n = size(A,2)
    φλ = [prod(λ.^A[:,i]) for i in 1:n]
    vec = φλ.*(B*v)
    vec
end


# Samples points from a rational variety
# -------------  Input:
# f               a list of polynomials
# n_samples       the number of samples
# -------------  Output:
# P               the list of samples
function sample(f, n_samples)
    n_vars = length(gens(parent(f[1])))
    P = []
    for i in 1:n_samples
        pt = rand(-1000:1000, n_vars)
        one_sample = [evaluate(pol, pt) for pol in f]
        push!(P,one_sample)
    end
    P
end

# Computes a Vandermonde matrix from a list of monomials and from
# -------------  Input:
# B               a matrix of monomials
# P               a matrix of samples
# -------------  Output:
# V               a Vandermonde matrix
function get_vandermonde_matrix(B,P)
    pts = P
    mons = B
    T = eltype(pts[1])
    V = zeros(T, (length(pts), length(mons)))
    for i in 1:length(pts)
        V[i, :] = [prod(pts[i].^mon) for mon in mons]
        #V[i, :] = V[i, :]/norm(V[i, :])
    end
    V = matrix_space(QQ, size(V)...)(V)
    V
end


