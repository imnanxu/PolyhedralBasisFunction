# PolyhedralBasisFunction
Software for computing the real-valued basis function of polyhedral groups [1]. If use any part of the code, please cite this paper [1] in your work.

[1] Xu, Nan, and Peter C. Doerschuk. "Computation of real-valued basis functions which transform as irreducible representations of the polyhedral groups." arXiv preprint arXiv:1701.01348 (2017).

I. Software for computing the real irrep matrices and real basis functions:
Three MATHEMATICA software packages (i.e., RealIrrepBasisT.m, RealIrrepBasisO.m, and RealIrrepBasisI.m) were developed for computing the real irrep matrices as well as the spherical harmonics coefficients c_{p,l,n,j,m} (Eq. 6.1) which define the real basis functions in terms of spherical harmonics for the three polyhedral groups. The following computations can be performed:
 1. The complex and the equivalent real irrep matrix, \Gamma_c^p and \Gamma_r^p (Eq. 5.1, [1])
 2. Non-orthogonized coefficients matrix \hat{D}_{l,m}^p  (Eq. 6.10, [1])
 3. Spherical harmonics coefficients matrix, \hat{H}_l^p=((\hat{H}_{l,1}^p)^T, ..., (\hat{H}_{l,Npl}^p)^T)^T  (Eq. 6.11, [1])

The final real basis functions can be obtained by multiplying each row of \boldsymbol{\hat\calH}_{l}^{p} by the spherical harmonics vector (i.e., Table[SphericalHarmonicY[l,m,\[Theta],\[Phi]],{m,-l,l}] in MATHEMATICA). Please see the notebook file Main.nb for the tutorial of calling these packages to generate real basis function for each polyhedral group. 

II. Numerical solution:
Note that the solution of real irrep matrices and coefficients are not unique as described in [1]. One solution for each group is included:
 1. Real irrep matrices: "RealIrreps_T.txt", "RealIrreps_O.txt", and "RealIrreps_I.txt" for groups T, O and I, respectively.
(*A 2-dim matrix with the 1st row {a, b} and 2nd row {c, d} has the form of {{a,b},{c,d}} in these files.*)
 2. Spherical harmonics coeffcients matrix \hat{H}_l^p for l=0,...,100: "BasisCoeff_T.txt" (for T), "BasisCoeff_O.txt" (for O), "BasisCoeff_I.txt" (for I).
(*File format: a line of 'l' value, a line of 'p' value, and then a line of coefficient matrix '\hat{H}_l^p'.*)
 3.Real basis functions {F}_{l}^p(\[Theta], \[Phi]) at randomly selected (\[Theta], \[Phi])'s for l=0,...,100: "RealBasisTest_T.txt" (for T), "RealBasisTest_O.txt" (for O), and "RealBasisTest_I.txt" (for I).
(*File format: a line of 'l' value, a line of 'p' value, and then a line of '{(\[Theta], \[Phi]}   {F}_{l}^p(\[Theta], \[Phi])'.*)

III. Obtain real basis functions in MATLAB:
MATLAB functions are also developed to read the coefficients file (i.e., read_coefMat.m for "BasisCoeff_*.txt") and then to compute the real basis functions (i.e., demonstrate_get_Fplnj.m and get_Fplnj.m). 
