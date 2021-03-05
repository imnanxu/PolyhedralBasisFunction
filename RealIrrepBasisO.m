(* ::Package:: *)

(*Written by Nan Xu Aug 11, 2015*)
(*Modified by Nan Xu Mar 02, 2021*)
(*ClearAll["Global`*"]*)
BeginPackage["RealIrrepBasisO`"]

  \[CapitalGamma]o0::usage="\[CapitalGamma]o0[p,g] returns the complex irrep matrix of the octahedral group element g (g=1,2,...,24) of the p'th (p=1,2,3,4,5) irrep.";

  \[CapitalGamma]o::usage="\[CapitalGamma]o[p,g] returns the equivalent real irrep matrix of the octahedral group element g (g=1,2,...,24) of the p'th (p=1,2,3,4,5) irrep.";

  BasisRealFunctionCoeffMatrixO::usage="BasisRealFunctionCoeffMatrixO[l,m,p] computes the non-orthogonalized coefficents matrix of the octahedral group that associates with the p'th irrep (p=1,2,3,4,5), the spherical harmonics degree l (l>=0) and order m (-l<=m<=l).";
  
  BasisRealFunctionOthoCoeffMatrixO::usage="BasisRealFunctionOthoCoeffMatrixO[l,p] computes the Gram\[Dash]Schmidt orthogonalized coefficents matrix of the octahedral group that associates with the p'th irrep (p=1,2,3,4,5), the spherical harmonics degree l (l>=0).";

  Begin["`Private`"]


\[Epsilon]=10^(-14);
\[Phi]=Pi/2;
\[Gamma]=Cos[\[Phi]]+I*Sin[\[Phi]];
SimilarTransMatrix[\[Gamma]_,n_]:=(J=SparseArray[{},{n,n}]; For[i=1,i<=n,i++,J[[i,n-i+1]]=IdentityMatrix[n][[i,i]]]; aGreek=SparseArray[{},{2n,2n}]; aGreek[[1;;n,1;;n]]=Re[\[Gamma]]*J; aGreek[[n+1;;2n,n+1;;2n]]=-Re[\[Gamma]]*J; aGreek[[1;;n,n+1;;2n]]=Im[\[Gamma]]*J; aGreek[[n+1;;2n,1;;n]]=Im[\[Gamma]]*J; Sall=Eigenvectors[aGreek]; St=SparseArray[{},{n,n}]; For[k=1,k<=n,k++,St[[k,1;;n]]=Sall[[k,1;;n]]+I*Sall[[k,n+1;;2n]]]; S=Normal[Transpose[St]])
MatrixForm[S3=SimilarTransMatrix[\[Gamma],2]/Sqrt[2]]
MatrixForm[S3=Orthogonalize[SimilarTransMatrix[\[Gamma],2]]](*same as the last one*)
(*MatrixForm[S3=SimilarTransMatrix[\[Gamma],3]/Sqrt[2]]
MatrixForm[S3=Orthogonalize[SimilarTransMatrix[\[Gamma],3]]](*same as the last one*)*)


Nx=6;
NgO=24;
dimO={1,1,2,3,3};
xvec=RandomReal[{0,1},{Nx,3}];
x=xvec[[2, All]];

Oired1={{{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}};
Oired2={{{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{-1}}, {{-1}}, {{-1}}, {{-1}}, {{-1}}, {{-1}}, {{-1}}, {{-1}}, {{-1}}, {{-1}}, {{-1}}, {{-1}}};
KeyImaginary=0.866025403784438*I;
Oired3={{{ 1, 0 }, 
{ 0, 1 }}, 
{{ 1, 0 }, 
{ 0, 1 }}, 
{{ 1, 0 }, 
{ 0, 1 }}, 
{{ 1, 0 }, 
{ 0, 1 }}, 
{{-0.50+KeyImaginary, 0 }, 
{ 0, -0.50-KeyImaginary }}, 
{{-0.50+KeyImaginary, 0 }, 
{ 0, -0.50-KeyImaginary }}, 
{{-0.50+KeyImaginary, 0 }, 
{ 0, -0.50-KeyImaginary }}, 
{{-0.50+KeyImaginary, 0 }, 
{ 0, -0.50-KeyImaginary }}, 
{{-0.50-KeyImaginary, 0 }, 
{ 0, -0.50+KeyImaginary }}, 
{{-0.50-KeyImaginary, 0 }, 
{ 0, -0.50+KeyImaginary }}, 
{{-0.50-KeyImaginary, 0 }, 
{ 0, -0.50+KeyImaginary }}, 
{{-0.50-KeyImaginary, 0 }, 
{ 0, -0.50+KeyImaginary }}, 
{{ 0, 1 }, 
{ 1, 0 }}, 
{{ 0, 1 }, 
{ 1, 0 }}, 
{{ 0, 1 }, 
{ 1, 0 }}, 
{{ 0, 1 }, 
{ 1, 0 }}, 
{{ 0, -0.50-KeyImaginary }, 
{-0.50+KeyImaginary, 0 }}, 
{{ 0, -0.50-KeyImaginary }, 
{-0.50+KeyImaginary, 0 }}, 
{{ 0, -0.50-KeyImaginary }, 
{-0.50+KeyImaginary, 0 }}, 
{{ 0, -0.50-KeyImaginary }, 
{-0.50+KeyImaginary, 0 }}, 
{{ 0, -0.50+KeyImaginary }, 
{-0.50-KeyImaginary, 0 }}, 
{{ 0, -0.50+KeyImaginary }, 
{-0.50-KeyImaginary, 0 }}, 
{{ 0, -0.50+KeyImaginary }, 
{-0.50-KeyImaginary, 0 }}, 
{{ 0, -0.50+KeyImaginary }, 
{-0.50-KeyImaginary, 0 }}};

Oired4={{{ 1, 0, 0 }, 
{ 0, 1, 0 }, 
{ 0, 0, 1 }}, 
{{ 1, 0, 0 }, 
{ 0, -1, 0 }, 
{ 0, 0, -1 }}, 
{{-1, 0, 0 }, 
{ 0, 1, 0 }, 
{ 0, 0, -1 }}, 
{{-1, 0, 0 }, 
{ 0, -1, 0 }, 
{ 0, 0, 1 }}, 
{{ 0, 0, 1 }, 
{ 1, 0, 0 }, 
{ 0, 1, 0 }}, 
{{ 0, 0, -1 }, 
{ 1, 0, 0 }, 
{ 0, -1, 0 }}, 
{{ 0, 0, -1 }, 
{-1, 0, 0 }, 
{ 0, 1, 0 }}, 
{{ 0, 0, 1 }, 
{-1, 0, 0 }, 
{ 0, -1, 0 }}, 
{{ 0, 1, 0 }, 
{ 0, 0, 1 }, 
{ 1, 0, 0 }}, 
{{ 0, -1, 0 }, 
{ 0, 0, -1 }, 
{ 1, 0, 0 }}, 
{{ 0, 1, 0 }, 
{ 0, 0, -1 }, 
{-1, 0, 0 }}, 
{{ 0, -1, 0 }, 
{ 0, 0, 1 }, 
{-1, 0, 0 }}, 
{{ 1, 0, 0 }, 
{ 0, 0, -1 }, 
{ 0, -1, 0 }}, 
{{ 1, 0, 0 }, 
{ 0, 0, 1 }, 
{ 0, 1, 0 }}, 
{{-1, 0, 0 }, 
{ 0, 0, 1 }, 
{ 0, -1, 0 }}, 
{{-1, 0, 0 }, 
{ 0, 0, -1 }, 
{ 0, 1, 0 }}, 
{{ 0, 0, 1 }, 
{ 0, -1, 0 }, 
{-1, 0, 0 }}, 
{{ 0, 0, -1 }, 
{ 0, 1, 0 }, 
{-1, 0, 0 }}, 
{{ 0, 0, -1 }, 
{ 0, -1, 0 }, 
{ 1, 0, 0 }}, 
{{ 0, 0, 1 }, 
{ 0, 1, 0 }, 
{ 1, 0, 0 }}, 
{{ 0, 1, 0 }, 
{-1, 0, 0 }, 
{ 0, 0, -1 }}, 
{{ 0, -1, 0 }, 
{-1, 0, 0 }, 
{ 0, 0, 1 }}, 
{{ 0, 1, 0 }, 
{ 1, 0, 0 }, 
{ 0, 0, 1 }}, 
{{ 0, -1, 0 }, 
{ 1, 0, 0 }, 
{ 0, 0, -1 }}};

Oired5={{{ 1, 0, 0 }, 
{ 0, 1, 0 }, 
{ 0, 0, 1 }}, 
{{ 1, 0, 0 }, 
{ 0, -1, 0 }, 
{ 0, 0, -1 }}, 
{{-1, 0, 0 }, 
{ 0, 1, 0 }, 
{ 0, 0, -1 }}, 
{{-1, 0, 0 }, 
{ 0, -1, 0 }, 
{ 0, 0, 1 }}, 
{{ 0, 0, 1 }, 
{ 1, 0, 0 }, 
{ 0, 1, 0 }}, 
{{ 0, 0, -1 }, 
{ 1, 0, 0 }, 
{ 0, -1, 0 }}, 
{{ 0, 0, -1 }, 
{-1, 0, 0 }, 
{ 0, 1, 0 }}, 
{{ 0, 0, 1 }, 
{-1, 0, 0 }, 
{ 0, -1, 0 }}, 
{{ 0, 1, 0 }, 
{ 0, 0, 1 }, 
{ 1, 0, 0 }}, 
{{ 0, -1, 0 }, 
{ 0, 0, -1 }, 
{ 1, 0, 0 }}, 
{{ 0, 1, 0 }, 
{ 0, 0, -1 }, 
{-1, 0, 0 }}, 
{{ 0, -1, 0 }, 
{ 0, 0, 1 }, 
{-1, 0, 0 }}, 
{{-1, 0, 0 }, 
{ 0, 0, 1 }, 
{ 0, 1, 0 }}, 
{{-1, 0, 0 }, 
{ 0, 0, -1 }, 
{ 0, -1, 0 }}, 
{{ 1, 0, 0 }, 
{ 0, 0, -1 }, 
{ 0, 1, 0 }}, 
{{ 1, 0, 0 }, 
{ 0, 0, 1 }, 
{ 0, -1, 0 }}, 
{{ 0, 0, -1 }, 
{ 0, 1, 0 }, 
{ 1, 0, 0 }}, 
{{ 0, 0, 1 }, 
{ 0, -1, 0 }, 
{ 1, 0, 0 }}, 
{{ 0, 0, 1 }, 
{ 0, 1, 0 }, 
{-1, 0, 0 }}, 
{{ 0, 0, -1 }, 
{ 0, -1, 0 }, 
{-1, 0, 0 }}, 
{{ 0, -1, 0 }, 
{ 1, 0, 0 }, 
{ 0, 0, 1 }}, 
{{ 0, 1, 0 }, 
{ 1, 0, 0 }, 
{ 0, 0, -1 }}, 
{{ 0, -1, 0 }, 
{-1, 0, 0 }, 
{ 0, 0, -1 }}, 
{{ 0, 1, 0 }, 
{-1, 0, 0 }, 
{ 0, 0, 1 }}};


\[Epsilon]=10^(-14);
Oired3Real=SparseArray[{},{NgO,2,2}];
For[k=1,k<=NgO,k++,Oired3Real[[k, All, All]]=Inverse[S3] . Oired3[[k, All, All]] . S3]
MatrixForm[Oired1Real=Chop[Normal[Oired1]*1.0,\[Epsilon]]];
MatrixForm[Oired2Real=Chop[Normal[Oired2]*1.0,\[Epsilon]]];
MatrixForm[Oired3Real=Chop[Normal[Oired3Real]*1.0,\[Epsilon]]];
MatrixForm[Oired4Real=Chop[Normal[Oired4]*1.0,\[Epsilon]]];
MatrixForm[Oired5Real=Chop[Normal[Oired5]*1.0,\[Epsilon]]];
\[CapitalGamma]o[p_,g_]:=(Piecewise[{{{{1}},p==1},{Normal[Oired2Real[[g]]],p==2},{Normal[Oired3Real[[g]]],p==3},{Normal[Oired4Real[[g,All,All]]],p==4},{Normal[Oired5Real[[g,All,All]]],p==5}}])
\[CapitalGamma]o0[p_,g_]:=(Piecewise[{{{{1}},p==1},{Normal[Oired2[[g]]],p==2},{Normal[Oired3[[g]]],p==3},{Normal[Oired4[[g,All,All]]],p==4},{Normal[Oired5[[g,All,All]]],p==5}}])


(*CHECK for Unitary
For[p=1, p<=5, p++, For[g=1, g<=24, g++, If[ConjugateTranspose[\[CapitalGamma]o[p,g]] . \[CapitalGamma]o[p,g]== IdentityMatrix[dimO[[p]]], Print[{p, g}]]]]*)


(*CHECK for Group Homomorphism*)
IsThere[M_,g1_,g2_]:=(\[Epsilon]=10^(-14); R2=Chop[M[[g1, All, All]] . M[[g2, All, All]],\[Epsilon]]; n=0; flg=0; For[k=1,k<=NgO,k++, test=M[[k, All, All]]; If[test==R2,flg=flg+1; n=k,flg=flg]]; n)
g1=15; g2=17;
IsThere[Oired3Real,g1,g2]
IsThere[Oired5Real,g1,g2]


(*Use O(432) p=5 as the rotation matrices*)
Ro=SparseArray[{},{NgO,3,3}];
For[k=1,k<=NgO, k++, Ro[[k, All, All]]=Oired5[[k, All, All]]; 
(*Print[{k, MatrixForm[Ro[[k, All, All]]], MatrixForm[Oired4[[k, All, All]]], MatrixForm[Oired5[[k, All, All]]]}]*)
]


(*Transformation Matrix To Euler Angle Converter*)
Ritb[\[Alpha]_,\[Beta]_,\[Gamma]_]:={{Cos[\[Gamma]]Cos[\[Beta]]Cos[\[Alpha]]-Sin[\[Gamma]]Sin[\[Alpha]], -(Cos[\[Gamma]]Cos[\[Beta]]Sin[\[Alpha]]+Sin[\[Gamma]]Cos[\[Alpha]]), Cos[\[Gamma]]Sin[\[Beta]]},
				{Sin[\[Gamma]]Cos[\[Beta]]Cos[\[Alpha]]+Cos[\[Gamma]]Sin[\[Alpha]],-Sin[\[Gamma]]Cos[\[Beta]]Sin[\[Alpha]]+Cos[\[Gamma]]Cos[\[Alpha]], Sin[\[Gamma]]Sin[\[Beta]]},
				{-Sin[\[Beta]]Cos[\[Alpha]], Sin[\[Beta]]Sin[\[Alpha]], Cos[\[Beta]]}}
MatrixToEulerAngleITB[M_]:=(Angles=SparseArray[{},{3,1}];
b1=ArcCos[M[[3,3]]]; 
a1=Piecewise[{{ArcTan[M[[1,1]],M[[2,1]]], b1==0}, {0, b1==Pi}, {ArcTan[M[[3,1]],-M[[3,2]]], b1!=0||Pi}}];
c1=Piecewise[{{0, b1==0}, {ArcTan[M[[2,2]],-M[[1,2]]], b1==Pi}, {ArcTan[M[[1,3]],M[[2,3]]], b1!=0||Pi}}];
Angles[[1]]=a1;Angles[[2]]=b1;Angles[[3]]=c1; Angles)
Rotest=SparseArray[{},{NgO,3,3}];
\[Alpha]o=SparseArray[{},{NgO,1}];
\[Beta]o=SparseArray[{},{NgO,1}];
\[Gamma]o=SparseArray[{},{NgO,1}];
For[g=1, g<=NgO, g++, M=Inverse[Ro[[g,All,All]]]; a=Chop[MatrixToEulerAngleITB[M][[1]],\[Epsilon]]; b=Chop[MatrixToEulerAngleITB[M][[2]],\[Epsilon]]; c=Chop[MatrixToEulerAngleITB[M][[3]],\[Epsilon]]; \[Alpha]o[[g]]=a; \[Beta]o[[g]]=b; \[Gamma]o[[g]]=c; Rotest[[g,All,All]]=Chop[Ritb[a,b,c],\[Epsilon]]]


(*For[g=1, g<=NgO, g++, Print[{g, MatrixForm[Ro[[g, All, All]]], MatrixForm[Inverse[Ro[[g, All, All]]]]}]]*)


(*Transformation Matrix To Euler Angle Converter*)
Ritb[\[Alpha]_,\[Beta]_,\[Gamma]_]:={{Cos[\[Gamma]]Cos[\[Beta]]Cos[\[Alpha]]-Sin[\[Gamma]]Sin[\[Alpha]], -(Cos[\[Gamma]]Cos[\[Beta]]Sin[\[Alpha]]+Sin[\[Gamma]]Cos[\[Alpha]]), Cos[\[Gamma]]Sin[\[Beta]]},
				{Sin[\[Gamma]]Cos[\[Beta]]Cos[\[Alpha]]+Cos[\[Gamma]]Sin[\[Alpha]],-Sin[\[Gamma]]Cos[\[Beta]]Sin[\[Alpha]]+Cos[\[Gamma]]Cos[\[Alpha]], Sin[\[Gamma]]Sin[\[Beta]]},
				{-Sin[\[Beta]]Cos[\[Alpha]], Sin[\[Beta]]Sin[\[Alpha]], Cos[\[Beta]]}}
MatrixToEulerAngleITB[M_]:=(Angles=SparseArray[{},{3,1}];
b1=ArcCos[M[[3,3]]]; 
a1=Piecewise[{{ArcTan[M[[1,1]],M[[2,1]]], b1==0}, {ArcTan[M[[2,2]],M[[1,2]]], b1==Pi}, {ArcTan[M[[3,1]],-M[[3,2]]], b1!=0||Pi}}];
c1=Piecewise[{{0, b1==0}, {0, b1==Pi}, {ArcTan[M[[1,3]],M[[2,3]]], b1!=0||Pi}}];
Angles[[1]]=a1;Angles[[2]]=b1;Angles[[3]]=c1; Angles)
Rotest=SparseArray[{},{NgO,3,3}];
\[Alpha]o=SparseArray[{},{NgO,1}];
\[Beta]o=SparseArray[{},{NgO,1}];
\[Gamma]o=SparseArray[{},{NgO,1}];
For[g=1, g<=NgO, g++, M=Inverse[Ro[[g,All,All]]]; a=Chop[MatrixToEulerAngleITB[M][[1]],\[Epsilon]]; b=Chop[MatrixToEulerAngleITB[M][[2]],\[Epsilon]]; c=Chop[MatrixToEulerAngleITB[M][[3]],\[Epsilon]]; \[Alpha]o[[g]]=a; \[Beta]o[[g]]=b; \[Gamma]o[[g]]=c; Rotest[[g,All,All]]=Chop[Ritb[a,b,c],\[Epsilon]]]


\[Alpha]oM=SparseArray[{},{NgO,1}];
\[Beta]oM=SparseArray[{},{NgO,1}];
\[Gamma]oM=SparseArray[{},{NgO,1}];
For[g=1, g<=NgO, g++, M=Ro[[g,All,All]]; a=EulerAngles[M][[1]]; b=EulerAngles[M][[2]]; c=EulerAngles[M][[3]]; \[Alpha]oM[[g]]=a; \[Beta]oM[[g]]=b; \[Gamma]oM[[g]]=c; Rotest[[g,All,All]]=Chop[Ritb[c,b,a],\[Epsilon]]]
(*
\[Alpha]oM
\[Beta]oM
\[Gamma]oM
Rotest==Ro
*)


(*Test for spherical harmonics rotation property*)
(*
For[l=0,l<=4, l++, Print[l]; 
For[g=1, g<=NgO, g++,
Ry=Y[Inverse[Ro[[g,All,All]]] . x,l,Table[m,{m, -l,l}]];
y=Y[x,l,Table[m,{m, -l,l}]];
MatrixForm[MWigD=Table[WignerD[{l, mp, m},-\[Alpha]oM[[g]], -\[Beta]oM[[g]], -\[Gamma]oM[[g]]],{m,-l,l},{mp,-l,l}]]
MatrixForm[MWigD=Table[WignerD[{l, mp, m},\[Gamma]oM[[g]], \[Beta]oM[[g]], \[Alpha]oM[[g]]],{m,-l,l},{mp,-l,l}]]
MatrixForm[MWigD=Table[Dwig[l,mp,m,\[Alpha]o[[g]],\[Beta]o[[g]],\[Gamma]o[[g]]],{m,-l,l},{mp,-l,l}]]
If[Chop[MWigD . y-Ry,\[Epsilon]]!=ConstantArray[0, 2l+1], Print[{l, MWigD . y, Ry}]]]]
*)


(*TEST for Multiplication Table*)
(*
M4=Oired4Real; M3=Oired3Real; M2=Oired2Real; Mr=Ro;
For[g1=1, g1<= NgO, g1++, 
	For[g2=1,g2<=NgO,g2++, 
		g3=IsThere[Oired5Real,g1,g2];
		R4=Chop[M4[[g1, All, All]] . M4[[g2, All, All]],\[Epsilon]];
		R3=Chop[M3[[g1, All, All]] . M3[[g2, All, All]],\[Epsilon]];
		R2=Chop[M2[[g1, All, All]] . M2[[g2, All, All]],\[Epsilon]];
		Mrc=Chop[Mr[[g1, All, All]] . Mr[[g2, All, All]],\[Epsilon]];
		If[R4!= M4[[g3, All, All]], Print[{4, g1, g2, g3}] ]
		If[R3!= M3[[g3, All, All]], Print[{3, g1, g2, g3}] ]
		If[R2!=M2[[g3, All, All]], Print[{2, g1, g2, g3}] ]
		If[Mrc!=Mr[[g3, All, All]], Print[{"R Mat", g1, g2, g3}] ]
]]
*)


(*TEST for Rotation Matrix Set and 3D Irep Mat*)
(*
For[g1=1, g1<= NgO, g1++, 
	R=Ro[[g1, All, All]]; 
	Yes1=0; Yes2=0;
	For[g2=1,g2<=NgO,g2++, 
		If[R==Oired4[[g2, All, All]], n1=g1; Yes1=1, If[Yes1==0, n1=0]]
		If[R==Oired5[[g2, All, All]], n2=g2; Yes2=1, If[Yes2==0, n2=0]]
	]; 
	Print[{g1, n1, n2}];
]*)


(*Spherical Harmonics*)
Clear[d]
\[Epsilon]=10^(-14);
Y[x_,l_,m_]:=( r=Sqrt[x[[1]]^2+x[[2]]^2+x[[3]]^2];  \[Theta]=ArcCos[x[[3]]/r];\[Phi]=ArcTan[x[[1]],x[[2]]] ;SphericalHarmonicY[l,m,\[Theta],\[Phi]])
(*D coefficiets*)
powerFunc=Unevaluated[#1^#2]/.HoldPattern[0^0]:>1&;
d[l_,mp_,m_,\[Beta]_]:=Sqrt[(l+m)!*(l-m)!*(l+mp)!*(l-mp)!]*Sum[(-1)^k/((l-mp-k)!*(l+m-k)!*(k+mp-m)!*k!)*powerFunc[(Cos[\[Beta]/2]),(2l+m-mp-2k)]*powerFunc[(-Sin[\[Beta]/2]),(mp-m+2k)],{k,Max[0,(m-mp)],Min[(l-mp),(l+m)]}]
Dwig[l_,mp_,m_,\[Alpha]_,\[Beta]_,\[Gamma]_]:=Exp[-I*mp*\[Alpha]]*d[l,mp,m,\[Beta]]*Exp[-I*m*\[Gamma]]
(*Group Projection Operator*)
(*\[CapitalDelta]o[l_,mp_,m_,p_,i_,j_]:=dimO[[p]]/NgO*Sum[(Conjugate[\[CapitalGamma]o[p,g][[i,j]]])*Dwig[l, m, mp, \[Alpha]o[[g]], \[Beta]o[[g]], \[Gamma]o[[g]]], {g,NgO}]*)
\[CapitalDelta]o[l_,mp_,m_,p_,i_,j_]:=dimO[[p]]/NgO*Sum[(Conjugate[\[CapitalGamma]o[p,g][[i,j]]])*WignerD[{l, mp, m}, \[Alpha]o[[g]], \[Beta]o[[g]], \[Gamma]o[[g]]],{g,NgO}] 
\[CapitalDelta]o[l_,mp_,m_,p_,i_,j_]:=dimO[[p]]/NgO*Sum[(Conjugate[\[CapitalGamma]o[p,g][[i,j]]])*WignerD[{l, mp, m}, -\[Alpha]oM[[g]], -\[Beta]oM[[g]], -\[Gamma]oM[[g]]],{g,NgO}] 
\[CapitalDelta]o[l_,mp_,m_,p_,i_,j_]:=dimO[[p]]/NgO*Sum[(Conjugate[\[CapitalGamma]o[p,g][[i,j]]])*Dwig[l, mp, m, \[Alpha]oM[[g]], \[Beta]oM[[g]], \[Gamma]oM[[g]]],{g,NgO}] 

(*MatrixForm[MWigD=Table[WignerD[{l, mp, m},-\[Alpha]oM[[g]], -\[Beta]oM[[g]], -\[Gamma]oM[[g]]],{m,-l,l},{mp,-l,l}]];*)
ProjOperatorO[l_,m_,p_,i_,j_,x_]:=Chop[Table[\[CapitalDelta]o[l,mp,m,p,i,j],{mp,-l,l}] . Y[x,l,Table[k,{k,-l,l}]],\[Epsilon]]
ProjOperatorTestO[l_,m_,p_,i_,j_,x_]:=dimO[[p]]/NgO*Chop[Sum[Conjugate[\[CapitalGamma]o[p,g][[i,j]]]*Y[Ro[[g,All,All]] . x,l,m],{g,NgO}],\[Epsilon]]
LengthProjOperatorO[l_,m_,p_,n_]:=Chop[Sum[Abs[\[CapitalDelta]o[l,mp,m,p,n,n]]^2,{mp, -l,l}],\[Epsilon]]

BasisFunctionMatrixO[l_,m_,p_]:=(Matrix=SparseArray[{},{dimO[[p]],2l+1}]; n=0;
	For[k=1,k<=dimO[[p]],k++, 
		(*cn=LengthProjOperatorO[l,m,p,k];*)
		ChoiceN=Table[\[CapitalDelta]o[l,mp,m,p,k,k],{mp, -l,l}];
		If[ChoiceN!=ConstantArray[0, 2l+1] ,n=k; Break[]]
	];
	If[n!=0, 
		Matrix[[n,All]]=Chop[Table[\[CapitalDelta]o[l,mp,m,p,n,n],{mp,-l,l}],\[Epsilon]];
		For[np=1,np<=dimO[[p]],np++, 
			If[np!=n,
				Matrix[[np,All]]=Chop[Table[\[CapitalDelta]o[l,mp,m,p,np,n],{mp,-l,l}],\[Epsilon]],
				Continue[]
			]
		],
	(*cn=0*)];
	Normal[Matrix])
BasisRealFunctionCoeffMatrixO[l_,m_,p_]:=(Mat=BasisFunctionMatrixO[l,m,p];
		afunc[mp_]:=Mat[[All,mp+l+1]]; lh=0; rh=0; 
		For[mp=-l, mp<=-1, mp++, lh=lh+afunc[mp]; rh=rh+(-1)^Abs[mp]*Conjugate[afunc[-mp]];];
		angle=0;
		For[i=1,i<=Length[rh],i++,If[lh[[i]]==0,Continue[],Break[]];];
		If[i<= Length[rh], angle=Log[rh[[i]]/lh[[i]]]/2];
	MatNew=Exp[angle]*Mat)
BasisRealFunctionOthoCoeffMatrixO[l_,p_]:=(Mat0={};
	(*WriteString["stdout", "l=",l,", p=", p, ", m={"];*)
	For[m=-l, m<=l, m++,
		Mattest=BasisRealFunctionCoeffMatrixO[l, m, p];
		(*WriteString["stdout", m, " "];*)
		(*Print[{m,timeCur[[1]]}]*)
		Mat0=Join[Mat0, Mattest];
	]; (*WriteString["stdout", "}\n"];*)
	(*Mat0orth=Select[Orthogonalize[N[Mat0,precisionN]],#!=ConstantArray[0,2*(l)+1]&]*)
	Mat0orth=Select[Orthogonalize[Chop[Mat0*1.0,\[Epsilon]]],#!=ConstantArray[0,2*(l)+1]&]
	)
RealBasisFunctionsO[l_,p_,\[Theta]_,\[Phi]_]:=(
	CoeffMat=BasisRealFunctionOthoCoeffMatrixO[l,p];
	If[Dimensions[CoeffMat]!={0},
		y=Table[SphericalHarmonicY[l,m,\[Theta],\[Phi]],{m,-l,l}];
		realBasis=Chop[CoeffMat . y],
		realBasis={}
	]
)


  End[]
  EndPackage[]
