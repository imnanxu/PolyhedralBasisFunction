(* ::Package:: *)

(*Written by Nan Xu Aug 11, 2015*)
(*Modified by Nan Xu Feb 28, 2021*)
(*ClearAll["Global`*"]*)
BeginPackage["RealIrrepBasisT`"]

  \[CapitalGamma]t0::usage="\[CapitalGamma]t0[p,g] returns the complex irrep matrix of the tetrahedral group element g (g=1,2,...,12) of the p'th (p=1,2,3,4) irrep.";

  \[CapitalGamma]t::usage="\[CapitalGamma]t[p,g] returns the equivalent real irrep matrix of the tetrahedral group element g (g=1,2,...,12) of the p'th (p=1,2,3,4) irrep.";

  BasisRealFunctionCoeffMatrixT::usage="BasisRealFunctionCoeffMatrixT[l,m,p] computes the non-orthogonalized coefficents matrix of the tetrahedral group that associates with the p'th irrep (p=1,2,3,4), the spherical harmonics degree l (l>=0) and order m (-l<=m<=l).";
  
  BasisRealFunctionOthoCoeffMatrixT::usage="BasisRealFunctionOthoCoeffMatrixT[l,p] computes the Gram\[Dash]Schmidt orthogonalized coefficents matrix of the tetrahedral group that associates with the p'th irrep (p=1,2,3,4), the spherical harmonics degree l (l>=0).";
		
 (* RealBasisFunctionsT::usage="RealBasisFunctionsT[l,p,\[Theta],\[Phi]] computes the real basis functions of the tetrahedral group that associates with the p'th irrep (p=1,2,3,4), the spherical harmonics degree l (l>=0) at the angle of (\[Theta],\[Phi]) where 0<=\[Theta]<=Pi and 0<=\[Phi]<=2Pi.";*)

  Begin["`Private`"]



(*Similarity Matrix Sp*)
\[Epsilon]=10^(-14);
\[Phi]=Pi/2;
\[Gamma]=Cos[\[Phi]]+I*Sin[\[Phi]];
SimilarTransMatrix[\[Gamma]_,n_]:=(J=SparseArray[{},{n,n}]; For[i=1,i<=n,i++,J[[i,n-i+1]]=IdentityMatrix[n][[i,i]]]; aGreek=SparseArray[{},{2n,2n}]; aGreek[[1;;n,1;;n]]=Re[\[Gamma]]*J; aGreek[[n+1;;2n,n+1;;2n]]=-Re[\[Gamma]]*J; aGreek[[1;;n,n+1;;2n]]=Im[\[Gamma]]*J; aGreek[[n+1;;2n,1;;n]]=Im[\[Gamma]]*J; Sall=Eigenvectors[aGreek]; St=SparseArray[{},{n,n}]; For[k=1,k<=n,k++,St[[k,1;;n]]=Sall[[k,1;;n]]+I*Sall[[k,n+1;;2n]]]; S=Normal[Transpose[St]])
S2=SimilarTransMatrix[\[Gamma],1]/Sqrt[2];
S2=Orthogonalize[SimilarTransMatrix[\[Gamma],1]];(*same as the last one*)
S3=SimilarTransMatrix[\[Gamma],3]/Sqrt[2];
S3=Orthogonalize[SimilarTransMatrix[\[Gamma],3]];(*same as the last one*)


(*Complex Irreps from the literature*)
NgT=12;
dimT={1,1,1,3};
Tired1={{{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}};
Tired2={{{1}}, {{1}}, {{1}}, {{1}}, {{-0.5+0.866025403784439*I}}, {{-0.5+0.866025403784439*I}}, {{-0.5+0.866025403784439*I}}, {{-0.5+0.866025403784439*I}}, {{-0.5-0.866025403784438*I}}, {{-0.5-0.866025403784438*I}}, {{-0.5-0.866025403784438*I}}, {{-0.5-0.866025403784438*I}}};
Tired3={{{1}}, {{1}}, {{1}}, {{1}}, {{-0.5-0.866025403784438*I}}, {{-0.5-0.866025403784438*I}}, {{-0.5-0.866025403784438*I}}, {{-0.5-0.866025403784438*I}}, {{-0.5+0.866025403784439*I}}, {{-0.5+0.866025403784439*I}}, {{-0.5+0.866025403784439*I}}, {{-0.5+0.866025403784439*I}}};
Tired4={{{ 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }},
{{ 1, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 }},
{{-1, 0, 0 }, { 0, -1, 0 }, { 0, 0, 1 }},
{{-1, 0, 0 }, { 0, 1, 0 }, { 0, 0, -1 }},
{{ 0, 0, 1 }, { 1, 0, 0 }, { 0, 1, 0 }},
{{ 0, 0, -1 }, { 1, 0, 0 }, { 0, -1, 0 }},
{{ 0, 0, 1 }, {-1, 0, 0 }, { 0, -1, 0 }},
{{ 0, 0, -1 }, {-1, 0, 0 }, { 0, 1, 0 }},
{{ 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 0 }},
{{ 0, -1, 0 }, { 0, 0, -1 }, { 1, 0, 0 }},
{{ 0, -1, 0 }, { 0, 0, 1 }, {-1, 0, 0 }}, 
{{ 0, 1, 0 }, { 0, 0, -1 }, {-1, 0, 0 }}};
(*For[k=1, k<= NgT, k++, Print[MatrixForm[Tired1[[k,All, All]]]]]
For[k=1, k<= NgT, k++, Print[MatrixForm[Tired2[[k,All, All]]]]]
For[k=1, k<= NgT, k++, Print[MatrixForm[Tired3[[k,All, All]]]]]
For[k=1, k<= NgT, k++, Print[{k, MatrixForm[Tired4[[k,All, All]]]}]]*)


(*Transform to Real irreps*)
Tired2Real=SparseArray[{},{NgT,1}];
For[k=1,k<=NgT,k++,Tired2Real[[k]]=Inverse[S2] . Tired2[[k, All, All]] . S2]
Tired3Real=SparseArray[{},{NgT,1}];
For[k=1,k<=NgT,k++,Tired3Real[[k]]=Inverse[S2] . Tired3[[k, All, All]] . S2]
Tired1Real=Normal[Tired1]*1.0;
Tired2Real=Normal[Tired2Real]*1.0;
Tired3Real=Normal[Tired3Real]*1.0;
Tired4Real=Normal[Tired4]*1.0;
\[CapitalGamma]t0[p_,g_]:=(Piecewise[{{{{1}},p==1},{Normal[Tired2[[g, All, All]]],p==2},{Normal[Tired3[[g, All, All]]],p==3},{Normal[Tired4[[g, All,All]]],p==4}}])
\[CapitalGamma]t[p_,g_]:=(Piecewise[{{{{1}},p==1},{Normal[Tired2Real[[g, All, All]]],p==2},{Normal[Tired3Real[[g, All, All]]],p==3},{Normal[Tired4Real[[g, All,All]]],p==4}}])


(*Test for Real irrep's homomorphism*)
(*
IsThere[M_,g1_,g2_]:=(\[Epsilon]=10^(-14); R2=Chop[M[[g1, All, All]] . M[[g2, All, All]],\[Epsilon]]; n=0; flg=0; For[k=1,k<=NgT,k++, test=M[[k, All, All]];If[test==R2,flg=flg+1;n=k,flg=flg]]; n)
g1=3; g2=11;
IsThere[Tired3Real,g1,g2]
*)


(*Use O(432) p=5 as the rotation matrices*)
Rt=SparseArray[{},{NgT,3,3}];
For[k=1,k<=NgT, k++, Rt[[k, All, All]]=Tired4[[k, All, All]]; 
	(*Print[{k, MatrixForm[Rt[[k, All, All]]], MatrixForm[Tired4[[k, All, All]]]}]*)
	]

Ritb[\[Alpha]_,\[Beta]_,\[Gamma]_]:={{Cos[\[Gamma]]Cos[\[Beta]]Cos[\[Alpha]]-Sin[\[Gamma]]Sin[\[Alpha]], -(Cos[\[Gamma]]Cos[\[Beta]]Sin[\[Alpha]]+Sin[\[Gamma]]Cos[\[Alpha]]), Cos[\[Gamma]]Sin[\[Beta]]},
				{Sin[\[Gamma]]Cos[\[Beta]]Cos[\[Alpha]]+Cos[\[Gamma]]Sin[\[Alpha]],-Sin[\[Gamma]]Cos[\[Beta]]Sin[\[Alpha]]+Cos[\[Gamma]]Cos[\[Alpha]], Sin[\[Gamma]]Sin[\[Beta]]},
				{-Sin[\[Beta]]Cos[\[Alpha]], Sin[\[Beta]]Sin[\[Alpha]], Cos[\[Beta]]}}
MatrixToEulerAngleITB[M_]:=(Angles=SparseArray[{},{3,1}];
b1=ArcCos[M[[3,3]]]; 
a1=Piecewise[{{ArcTan[M[[1,1]],M[[2,1]]], b1==0}, {0, b1==Pi}, {ArcTan[-M[[3,1]],M[[3,2]]], b1!=0||Pi}}];
c1=Piecewise[{{0, b1==0}, {ArcTan[M[[2,2]],-M[[1,2]]], b1==Pi}, {ArcTan[M[[1,3]],M[[2,3]]], b1!=0||Pi}}];
Angles[[1]]=a1;Angles[[2]]=b1;Angles[[3]]=c1; Angles)
Rttest=SparseArray[{},{NgT,3,3}];
\[Alpha]t=SparseArray[{},{NgT,1}];
\[Beta]t=SparseArray[{},{NgT,1}];
\[Gamma]t=SparseArray[{},{NgT,1}];
 \[Epsilon]=10^(-14);
For[g=1, g<=NgT, g++,M=Inverse[Rt[[g,All,All]]]; a=Chop[MatrixToEulerAngleITB[M][[1]],\[Epsilon]]; b=Chop[MatrixToEulerAngleITB[M][[2]],\[Epsilon]]; c=Chop[MatrixToEulerAngleITB[M][[3]],\[Epsilon]]; \[Alpha]t[[g]]=a; \[Beta]t[[g]]=b; \[Gamma]t[[g]]=c;Rttest[[g,All,All]]=Ritb[a,b,c]]


(*TEST for Multiplication Table*)
(*
M3=Tired3Real; M2=Tired2Real; Mr=Rt;
For[g1=1, g1<= NgT, g1++, 
	For[g2=1,g2<=NgT,g2++, 
		g3=IsThere[Tired4Real,g1,g2];
		R3=Chop[M3[[g1, All, All]] . M3[[g2, All, All]],\[Epsilon]];
		R2=Chop[M2[[g1, All, All]] . M2[[g2, All, All]],\[Epsilon]];
		Mrc=Chop[Mr[[g1, All, All]] . Mr[[g2, All, All]],\[Epsilon]];
		If[R3!= M3[[g3, All, All]], Print[{3, g1, g2, g3}] ]
		If[R2!=M2[[g3, All, All]], Print[{2, g1, g2, g3}] ]
		If[Mrc!=Mr[[g3, All, All]], Print[{"rotation", g1, g2, g3}] ]
]]
*)


(*Test for spherical harmonics rotation property*)
(*
For[l=0,l<=45, l++, Print[l]; 
For[g=1, g<=NgT, g++,
Ry=Y[Inverse[Rt[[g,All,All]]].x,l,Table[m,{m, -l,l}]];
y=Y[x,l,Table[m,{m, -l,l}]];
MatrixForm[MWigD=Table[WignerD[{l, mp, m}, \[Alpha]t[[g]], \[Beta]t[[g]], \[Gamma]t[[g]]],{m,-l,l},{mp,-l,l}]];
If[Chop[MWigD.y-Ry,\[Epsilon]]!=ConstantArray[0, 2l+1], Print[{l, MWigD.y, Ry}]]]]*)


(*Real Basis Functions for Spherical Harmonics*)
Clear[d]
Off[General::munfl]; Off[General::stop];
Y[x_,l_,m_]:=( r=Sqrt[x[[1]]^2+x[[2]]^2+x[[3]]^2];  \[Theta]=ArcCos[x[[3]]/r];\[Phi]=ArcTan[x[[1]],x[[2]]]; SphericalHarmonicY[l,m,\[Theta],\[Phi]])
(*D coefficiets*)
powerFunc=Unevaluated[#1^#2]/.HoldPattern[0^0]:>1&;
d[l_,mp_,m_,\[Beta]_]:=Sqrt[(l+m)!*(l-m)!*(l+mp)!*(l-mp)!]*Sum[(-1)^k/((l-mp-k)!*(l+m-k)!*(k+mp-m)!*k!)*powerFunc[(Cos[\[Beta]/2]),(2l+m-mp-2k)]*powerFunc[(-Sin[\[Beta]/2]),(mp-m+2k)],{k,Max[0,(m-mp)],Min[(l-mp),(l+m)]}]
Dwig[l_,mp_,m_,\[Alpha]_,\[Beta]_,\[Gamma]_]:=Exp[-I*mp*\[Alpha]]*d[l,mp,m,\[Beta]]*Exp[-I*m*\[Gamma]]
(*Group Projection Operator*)
(*\[CapitalDelta]t[l_,mp_,m_,p_,i_,j_]:=dimT[[p]]/NgT*Sum[(Conjugate[\[CapitalGamma]t[p,g][[i,j]]])*Dwig[l, mp, m, \[Alpha]t[[g]], \[Beta]t[[g]], \[Gamma]t[[g]]], {g,NgT}]*)
\[CapitalDelta]t[l_,mp_,m_,p_,i_,j_]:=dimT[[p]]/NgT*Sum[(Conjugate[\[CapitalGamma]t[p,g][[i,j]]])*WignerD[{l, mp, m}, \[Alpha]t[[g]], \[Beta]t[[g]], \[Gamma]t[[g]]],{g,NgT}] 
ProjOperatorT[l_,m_,p_,i_,j_,x_]:=Chop[Table[\[CapitalDelta]t[l,mp,m,p,i,j],{mp,-l,l}] . Y[x,l,Table[k,{k,-l,l}]],\[Epsilon]]
ProjOperatorTestT[l_,m_,p_,i_,j_,x_]:=dimT[[p]]/NgT*Chop[Sum[Conjugate[\[CapitalGamma]t[p,g][[i,j]]]*Y[Rt[[g,All,All]] . x,l,m],{g,NgT}],\[Epsilon]]
LengthProjOperatorT[l_,m_,p_,n_]:=Chop[Sum[Abs[\[CapitalDelta]t[l,mp,m,p,n,n]]^2,{mp, -l,l}],\[Epsilon]]

BasisFunctionMatrixT[l_,m_,p_]:=(Matrix=SparseArray[{},{dimT[[p]],2l+1}]; n=0;
	For[k=1,k<=dimT[[p]],k++, 
		(*cn=LengthProjOperatorO[l,m,p,k];*)
		ChoiceN=Table[\[CapitalDelta]t[l,mp,m,p,k,k],{mp, -l,l}];
		If[ChoiceN!=ConstantArray[0, 2l+1] ,n=k; Break[]]
	];
	If[n!=0, 
		Matrix[[n,All]]=Table[\[CapitalDelta]t[l,mp,m,p,n,n],{mp,-l,l}];
		For[np=1,np<=dimT[[p]],np++, 
			If[np!=n,
				Matrix[[np,All]]=Table[\[CapitalDelta]t[l,mp,m,p,np,n],{mp,-l,l}],
				Continue[]
			]
		],
	(*cn=0*)];
	Normal[Matrix])
BasisRealFunctionCoeffMatrixT[l_,m_,p_]:=(Mat=BasisFunctionMatrixT[l,m,p];
		afunc[mp_]:=Mat[[All,mp+l+1]]; lh=0; rh=0; 
		For[mp=-l, mp<=-1, mp++, lh=lh+afunc[mp]; rh=rh+(-1)^Abs[mp]*Conjugate[afunc[-mp]*1.0];];
		angle=0;
		For[i=1,i<=Length[rh],i++,If[lh[[i]]==0,Continue[],Break[]];];
		If[i<= Length[rh], angle=Log[rh[[i]]/lh[[i]]]/2];
	MatNew=Exp[angle]*Mat)
BasisRealFunctionOthoCoeffMatrixT[l_,p_]:=(Mat0={};
	(*WriteString["stdout", "l=",l,", p=", p, ", m={"];*)
	For[m=-l, m<=l, m++,
		Mattest=BasisRealFunctionCoeffMatrixT[l, m, p];
		(*WriteString["stdout", m, " "];*)
		(*Print[{m,timeCur[[1]]}]*)
		Mat0=Join[Mat0, Mattest];
	]; (*WriteString["stdout", "}\n"];*)
	Mat0orth=Select[Orthogonalize[Chop[Mat0*1.0,\[Epsilon]]],#!=ConstantArray[0,2*(l)+1]&])
(*RealBasisFunctionsT[l_,p_,\[Theta]_,\[Phi]_]:=(
	CoeffMat=BasisRealFunctionOthoCoeffMatrixT[l,p];
	If[Dimensions[CoeffMat]!={0},
		y=Table[SphericalHarmonicY[l,m,\[Theta],\[Phi]],{m,-l,l}];
		realBasis=Chop[CoeffMat . y, \[Epsilon]],
		realBasis={}
	]
)*)


  End[]

  EndPackage[]
