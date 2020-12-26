# PolyhedralBasisFunction
Matlab and Mathematica Software for computing the complete basis function for the article [1]. If use any part of the code, please cite this paper [1] in your work.

[1] Xu, Nan, and Peter C. Doerschuk. "Computation of real basis functions for the 3-D rotational polyhedral groups T, O, and I." arXiv preprint arXiv:1701.01348 (2017).

Software in Mathematica has been written to compute the coeffcients c_{p,l,n,j,m} in [1, Eq. 9] which define the symmetrical basis in terms of spherical harmonics.
1. For icosahedral symmetry: In sequence,
 (a) Run LiuPingChen_toUnitaryReal_1012.m.
 (b) Run OthoBasisFunctionRealirre_1012.m.
 (c) Run BasisRealFunctionCoeffMatrix[l_,m_,p_] which computes the coeffcients for c_{p=p_,l=l_,n,j,m=m_}.

2. For octahedral symmetry: In sequence,
 (a) Run OctRealIrredBasisMat 081115.wl.
 (b) Run BasisRealFunctionCoeffMatrixO[l_,m_,p_] which computes the coefficients for c_{p=p_,l=l_,n,jmm=m_}.

3. For tetrahedral symmetry: In sequence,
 (a) Run TecRealIrredBasisMat_081115.wl
 (b) Run BasisRealFunctionCoeffMatrixT[l_,m_,p_] which computes the coefficients for c_{p=p_,l=l_,n,j,m=m_}. 
Note, only the p = 1 and p = 4 irreps are real valued and lead to real-valued functions.
 
Get the value for basis functions for a polyhedral group:
1. Example: The coefficients for the Icosahedral group generated from above with 0<=l<=45: IcosahedralRealBasisFunctionCoeff.txt
2. The Matlab that reads the coefficients .txt file: read_coefficients.m
3. The Matlab that computes the basis function values F_{p,l,n,j}(\theta;\phi) in [1, Eq. 8]: demonstrate get Fplnj.m and get Fplnj.m. 
