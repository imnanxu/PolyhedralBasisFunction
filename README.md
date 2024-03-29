# PolyhedralBasisFunction
Software for computing the real-valued basis function of polyhedral groups [1]. If use any part of the code or solutions, please cite this paper [1] in your work.

[1] Xu, Nan*, and Peter C. Doerschuk. "Computation of real-valued basis functions which transform as irreducible representations of the polyhedral groups." SIAM Journal on Scientific Computing 2021 43:6, A3657-A3676, [doi:10.1137/20M1318183](https://epubs.siam.org/doi/abs/10.1137/20M1318183). ([arXiv](https://arxiv.org/abs/1701.01348v2))


## I. Software for computing the real irrep matrices and real basis functions
Please see the MATHEMATICA notebook file "Main.nb" for the tutorial of calling the relavent package to generate real basis function for each polyhedral group. Three MATHEMATICA software packages (i.e., RealIrrepBasisT.m, RealIrrepBasisO.m, and RealIrrepBasisI.m) were developed for computing the real irrep matrices as well as the spherical harmonics coefficients `c_{p,l,n,j,m}` (Eq. 6.1, [1]) which define the real basis functions in terms of spherical harmonics for the three polyhedral groups. The following computations can be performed:
 ### 1. The complex and the equivalent real irrep matrix
     $\Gamma_c^p$ for complex, and $\Gamma_r^p$ for real (Eq. 5.1, [1])
 ### 2. Non-orthogonized coefficients matrix 
     $\hat{D}_{l,m}^p$  (Eq. 6.10, [1])
 ### 3. Spherical harmonics coefficients matrix
     $\hat{H}_l^p=((\hat{H}_{l,1}^p)^T, ..., (\hat{H}_{l,Npl}^p)^T)^T$  (Eq. 6.11, [1])

The final independent real basis functions can be obtained by multiplying each row of $\hat{H}_l^p$ by the spherical harmonics vector (i.e., `Table[SphericalHarmonicY[l,m,\[Theta],\[Phi]],{m,-l,l}]` in MATHEMATICA). All programs were tested in MATHEMATICA v12.2. 

## II. Numerical solutions
Note that the solution of real irrep matrices and coefficients are not unique as described in [1]. One solution for each group is included ("sol_T.zip" for the tetrahedral group, "sol_O.zip" for the octahedral group, and "sol_I.zip" for the icosahedral group). The following three files are included in each "_*.zip" folder.
 ### 1. Real irrep matrices: "RealIrreps_*.txt"
     A 2-dim matrix with the 1st row {a, b} and 2nd row {c, d} has the form of {{a,b},{c,d}} in all these files.
  Note: For the tetrahedral group, only p=1 and p=4 irrep matrices are real valued.
  
 ### 2. Spherical harmonics coeffcients matrix `\hat{H}_l^p` for 0<=l<=100: "BasisCoeff<ins> </ins>*.txt"
     File format: a line of '$l$' value
                  a line of '$p$' value
                  a line of coefficient matrix '$\hat{H}_l^p$'
  Note: The `d_p*n+j`th row in `\hat{H}_l^p` includes the coefficients `{c_{p,l,n,j,m}}_{m=-l}^{m=l}` (Eq. 6.1, [1]) for `n=1,...,Npl`. 

 ### 3. Real basis functions `F_l^p` (Eq. 6.11, [1]) for 0<=l<=100 at a randomly selected `(\[Theta], \[Phi])`: "RealBasisTest_*.txt"
     File format: a line of '$l$' value
                  a line of '$p$' value
                  a line of '{$\theta$, $\phi$}     $F_l^p(\theta, \phi)$'
  Note: For the tetrahedral group, only p=1 and p=4 irrpes give real basis functions. For testing the completeness of the subspace determined by l degree, complex irrep matrices, coefficients and basis functions for p=2 and 3 irreps of the tetrahedral group were also included in above three files.
 
## III. Obtain real basis functions in MATLAB
MATLAB functions are also developed to read the coefficients file (i.e., read_coefMat.m for "BasisCoeff_*.txt") and then to compute the real basis functions (i.e., demonstrate_get_Fplnj.m and get_Fplnj.m).
