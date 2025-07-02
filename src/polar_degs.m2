
p_ring = QQ[x, y, z]
I = ideal(x^4-x*y^2*z+z^4)






needsPackage "Resultants"
needsPackage "Quasidegrees"



A = matrix{{1, 1, 1, 1}, {0, 2, 3, 4}}


k = numgens source A
p_ring = QQ[x_1..x_k]
I = toricIdeal(A, p_ring)


J = conormalVariety I;
origvars = (gens ring J)_{0..k-1}
dualvars = (gens ring J)_{k..2*k-1}

eliminate(J, origvars)
eliminate(J, dualvars)


u = {7, 2, 1}

eliminate(I, origvars)

critIdeal =  ideal (dualvars + origvars - u)+ J
