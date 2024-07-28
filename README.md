# SFS-Tajima-s-D

Tajima’s D is  a statistical test that compares the number of segregating sites to the average number of nucleotide differences [1]. It is used to test for neutrality in a population by determining if the observed DNA sequence variation is consistent with the neutral theory of evolution [5]. This can help detect different selection models. 

Tajima’s D = (π - θ_W) / sqrt(V_D)
where:
π is the average number of pairwise differences (estimator of θ based on π)
θ_W is Watterson's estimator of θ
VD is the variance of (π - θ_W)

The code calculates these components as follows:
a) π (theta_pi):
π = Σ_i 2i(n-i)SFS[i-1] / (n(n-1))

b) θ_W (theta_w):
θ_W = S / a1
where S is the number of segregating sites, and a1 = Σ_i^(n-1) 1/i

c) Variance V_D:
VD = e_1S + e_2S(S-1)
where e_1 and e_2 are complex functions of n derived in Tajima's original paper [1]. The implementation uses JavaScript to perform these calculations in the browser.

# Please also cite:
1)	Tajima, F. 1989. Statistical method for testing the neutral   mutation hypothesis by DNA polymorphism. Genetics 123:585-595.
