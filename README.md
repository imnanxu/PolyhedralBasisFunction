# PolyhedralBasisFunction
Matlab and Mathematica Software for computing the complete basis function for the article [1]. If use any part of the code, please cite this paper [1] in your work.

[1] Xu, Nan, and Peter C. Doerschuk. "Computation of real basis functions for the 3-D rotational polyhedral groups T, O, and I." arXiv preprint arXiv:1701.01348 (2017).

I. Software in Mathematica has been written to compute the real irrep matrices and then the spherical harmonics coeffcients ${\hat\calH}_{l,n}^p$ in [1, Eq. 6.10] which define the symmetrical basis in terms of spherical harmonics for these three groups.
1. For icosahedral group I: In sequence,
 (a) Run LiuPingChen_toUnitaryReal_1012.m.
 (b) Run OthoBasisFunctionRealirre_1012.m.
 (c) Run BasisRealFunctionCoeffMatrix[l_,m_,p_] which computes the coeffcients for c_{p=p_,l=l_,n,j,m=m_}.

2. For octahedral group O: In sequence,
 (a) Run OctRealIrredBasisMat 081115.wl.
 (b) Run BasisRealFunctionCoeffMatrixO[l_,m_,p_] which computes the coefficients for c_{p=p_,l=l_,n,jmm=m_}.

3. For tetrahedral group T: In sequence,
 (a) Run TecRealIrredBasisMat_081115.wl
 (b) Run BasisRealFunctionCoeffMatrixT[l_,m_,p_] which computes the coefficients for c_{p=p_,l=l_,n,j,m=m_}. 
Note, only the p = 1 and p = 4 irreps are real valued and lead to real-valued functions.

II. Numerical solution:
Note that the solution of real irrep matrices and coefficients are not unique as described in [1]. One solution for each group is included:
1. Real irrep matrices: "RealIrreps_T.txt", "RealIrreps_O.txt", and "RealIrreps_I.txt" for groups T, O and I, respectively.
 (*A 2-dim matrix with the 1st row {a, b} and 2nd row {c, d} is demonstrated in the form of {{a,b},{c,d}} in these three files.*)
2. Spherical harmonics coeffcients ${\hat\calH}_{l,n}^p$ in [1, Eq. 6.10] that construct a real basis: "BasisCoeff_T.txt" (for T), "BasisCoeff_O.txt" (for O), "BasisCoeff_I.txt" (for I).
(*The coefficient matrices {{\hat\calH}_{l,n=1}^p, ..., {\hat\calH}_{l,n=N_{p;l}}^p} are listed after a line of l value and a line of p value.*)

 
III. Get the value for basis functions for a polyhedral group in Matlab:
1. Example: The coefficients for the Icosahedral group generated from above with 0<=l<=45: IcosahedralRealBasisFunctionCoeff.txt
2. The Matlab that reads the coefficients .txt file: read_coefficients.m
3. The Matlab that computes the basis function values F_{p,l,n,j}(\theta;\phi) in [1, Eq. 9]: demonstrate_get_Fplnj.m and get_Fplnj.m. 


