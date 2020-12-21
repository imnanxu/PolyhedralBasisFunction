(* ::Package:: *)

(*Written by Nan Xu Oct 12, 2015*)
(*Parameters Set-Up*)
(* <<KtoI_irredi.m *)
l=5;
Nx=6;
Ng=60;
dim={1,3,3,4,5};
xvec=RandomReal[{0,1},{Nx,3}];
x=xvec[[2,All]];


Ryb=SparseArray[{},{60,3,3}];
S={{Cos[2*Pi/5],-Sin[2*Pi/5],0},{Sin[2*Pi/5],Cos[2*Pi/5],0},{0,0,1}};
U={{1/Sqrt[5],0,2/Sqrt[5]},{0,1,0},{-2/Sqrt[5],0,1/Sqrt[5]}};
P={{-1,0,0},{0,1,0},{0,0,-1}};
T=U.S.Inverse[U];
Ryb[[1,All,All]]={{1,0,0},{0,1,0},{0,0,1}};
Ryb[[2, All,All]]=S;
Ryb[[3,All,All]]=S.Ryb[[2, All,All]];
Ryb[[4,All,All]]=S.Ryb[[3, All,All]];
Ryb[[5,All,All]]=S.Ryb[[4, All,All]];
Ryb[[6,All,All]]=S.T;
Ryb[[7,All,All]]=T.Ryb[[6,All,All]];
Ryb[[8,All,All]]=T.Ryb[[7,All,All]];
Ryb[[9,All,All]]=Inverse[T].Ryb[[6,All,All]];
Ryb[[10,All,All]]=Inverse[T].Ryb[[9,All,All]];
Ryb[[11,All,All]]=S.Ryb[[6,All,All]];
Ryb[[12,All,All]]=S.Ryb[[7,All,All]];
Ryb[[13,All,All]]=S.Ryb[[8,All,All]];
Ryb[[14,All,All]]=S.Ryb[[9,All,All]];
Ryb[[15,All,All]]=S.Ryb[[10,All,All]];
Ryb[[16,All,All]]=S.Ryb[[11,All,All]];
Ryb[[17,All,All]]=S.Ryb[[12,All,All]];
Ryb[[18,All,All]]=S.Ryb[[13,All,All]];
Ryb[[19,All,All]]=S.Ryb[[14,All,All]];
Ryb[[20,All,All]]=S.Ryb[[15,All,All]];
Ryb[[21,All,All]]=Inverse[S].Ryb[[6,All,All]];
Ryb[[22,All,All]]=Inverse[S].Ryb[[7,All,All]];
Ryb[[23,All,All]]=Inverse[S].Ryb[[8,All,All]];
Ryb[[24,All,All]]=Inverse[S].Ryb[[9,All,All]];
Ryb[[25,All,All]]=Inverse[S].Ryb[[10,All,All]];
Ryb[[26,All,All]]=Inverse[S].Ryb[[21,All,All]];
Ryb[[27,All,All]]=Inverse[S].Ryb[[22,All,All]];
Ryb[[28,All,All]]=Inverse[S].Ryb[[23,All,All]];
Ryb[[29,All,All]]=Inverse[S].Ryb[[24,All,All]];
Ryb[[30,All,All]]=Inverse[S].Ryb[[25,All,All]];
For[i=1,i<=30,i++, Ryb[[i+30,All,All]]=P.Ryb[[i,All,All]];];


(*from, /home/yara/doerschuk/doerschu/yibinsoftware/EM3.matlab*)
Rpeter=SparseArray[{},{60,3,3}];
Rpeter[[1, All,All]]={{1,0,0},{0,1,0},{0,0,1}};
Rpeter[[2, All,All]]={{0.3090169943749474,-0.951056516295153,0},{0.951056516295153,0.3090169943749474,0},{0,0,1}};
Rpeter[[3,All,All]]={{-0.809016994374947,-0.5877852522924732,0},{0.5877852522924732,-0.809016994374947,0},{0,0,1}};
Rpeter[[4,All,All]]={{0.3090169943749474,0.951056516295153,0},{-0.951056516295153,0.3090169943749474,0},{0,0,1}};
Rpeter[[5,All,All]]={{-0.809016994374947,0.5877852522924732,0},{-0.5877852522924732,-0.809016994374947,0},{0,0,1}};
Rpeter[[6,All,All]]={{-0.1381966011250104,-0.4253254041760199,0.894427190999916},{0.951056516295154,-0.3090169943749474,0},{0.2763932022500211,0.85065080835204,0.447213595499958}};
Rpeter[[7,All,All]]={{-0.4472135954999578,0,0.894427190999916},{0,-1,0},{0.894427190999916,0,0.447213595499958}};
Rpeter[[8,All,All]]={{-0.1381966011250104,0.4253254041760199,0.894427190999916},{-0.951056516295153,-0.3090169943749474,0},{0.2763932022500212,-0.85065080835204,0.447213595499958}};
Rpeter[[9,All,All]]={{0.3618033988749896,-0.2628655560595667,0.894427190999916},{0.5877852522924732,0.809016994374947,0},{-0.723606797749979,0.5257311121191335,0.447213595499958}};
Rpeter[[10,All,All]]={{0.3618033988749896,0.2628655560595667,0.894427190999916},{-0.587785252292473,0.809016994374947,0},{-0.723606797749979,-0.5257311121191334,0.447213595499958}};
Rpeter[[11,All,All]]={{-0.947213595499958,0.1624598481164531,0.2763932022500211},{0.1624598481164533,-0.5,0.85065080835204},{0.2763932022500211,0.85065080835204,0.447213595499958}};
Rpeter[[12,All,All]]={{-0.1381966011250105,0.951056516295153,0.2763932022500211},{-0.4253254041760198,-0.3090169943749474,0.85065080835204},{0.894427190999916,0,0.447213595499958}};
Rpeter[[13,All,All]]={{0.861803398874989,0.4253254041760199,0.2763932022500212},{-0.4253254041760198,0.3090169943749474,0.85065080835204},{0.2763932022500212,-0.85065080835204,0.447213595499958}};
Rpeter[[14,All,All]]={{-0.447213595499958,-0.85065080835204,0.276393202250021},{0.5257311121191338,0,0.85065080835204},{-0.723606797749979,0.5257311121191335,0.447213595499958}};
Rpeter[[15,All,All]]={{0.6708203932499368,-0.6881909602355866,0.2763932022500212},{0.1624598481164532,0.5,0.85065080835204},{-0.723606797749979,-0.5257311121191334,0.447213595499958}};
Rpeter[[16,All,All]]={{-0.4472135954999582,0.5257311121191335,-0.7236067977499788},{-0.85065080835204,0,0.5257311121191336},{0.2763932022500211,0.85065080835204,0.447213595499958}};
Rpeter[[17,All,All]]={{0.3618033988749893,0.5877852522924732,-0.723606797749979},{-0.2628655560595667,0.809016994374947,0.5257311121191338},{0.894427190999916,0,0.447213595499958}};
Rpeter[[18,All,All]]={{0.6708203932499368,-0.1624598481164531,-0.7236067977499792},{0.6881909602355867,0.5,0.5257311121191339},{0.2763932022500212,-0.85065080835204,0.447213595499958}};
Rpeter[[19,All,All]]={{-0.6381966011250108,-0.2628655560595669,-0.723606797749979},{-0.2628655560595667,-0.809016994374947,0.5257311121191338},{-0.723606797749979,0.5257311121191335,0.447213595499958}};
Rpeter[[20,All,All]]={{0.052786404500042,-0.6881909602355866,-0.7236067977499792},{0.6881909602355867,-0.4999999999999999,0.5257311121191339},{-0.723606797749979,-0.5257311121191334,0.447213595499958}};
Rpeter[[21,All,All]]={{0.86180339887499,-0.4253254041760198,0.2763932022500211},{0.4253254041760199,0.3090169943749474,-0.85065080835204},{0.2763932022500211,0.85065080835204,0.447213595499958}};
Rpeter[[22,All,All]]={{-0.1381966011250104,-0.951056516295153,0.276393202250021},{0.4253254041760198,-0.3090169943749474,-0.85065080835204},{0.894427190999916,0,0.447213595499958}};
Rpeter[[23,All,All]]={{-0.947213595499958,-0.1624598481164532,0.2763932022500212},{-0.1624598481164533,-0.5,-0.85065080835204},{0.2763932022500212,-0.85065080835204,0.447213595499958}};
Rpeter[[24,All,All]]={{0.6708203932499371,0.6881909602355867,0.2763932022500211},{-0.1624598481164532,0.5,-0.85065080835204},{-0.723606797749979,0.5257311121191335,0.447213595499958}};
Rpeter[[25,All,All]]={{-0.4472135954999578,0.85065080835204,0.2763932022500212},{-0.5257311121191338,0,-0.85065080835204},{-0.723606797749979,-0.5257311121191334,0.447213595499958}};
Rpeter[[26,All,All]]={{0.670820393249937,0.1624598481164531,-0.7236067977499788},{-0.6881909602355869,0.4999999999999999,-0.5257311121191336},{0.2763932022500211,0.85065080835204,0.447213595499958}};
Rpeter[[27,All,All]]={{0.3618033988749894,-0.5877852522924732,-0.723606797749979},{0.2628655560595667,0.809016994374947,-0.5257311121191338},{0.894427190999916,0,0.447213595499958}};
Rpeter[[28,All,All]]={{-0.447213595499958,-0.5257311121191335,-0.7236067977499792},{0.85065080835204,0,-0.5257311121191339},{0.2763932022500212,-0.85065080835204,0.447213595499958}};
Rpeter[[29,All,All]]={{0.05278640450004211,0.6881909602355867,-0.723606797749979},{-0.6881909602355869,-0.5,-0.5257311121191338},{-0.723606797749979,0.5257311121191335,0.447213595499958}};
Rpeter[[30,All,All]]={{-0.6381966011250108,0.2628655560595669,-0.7236067977499792},{0.2628655560595666,-0.809016994374947,-0.5257311121191339},{-0.723606797749979,-0.5257311121191334,0.447213595499958}};
Rpeter[[31,All,All]]={{-1,0,0},{0,1,0},{0,0,-1}};
Rpeter[[32,All,All]]={{-0.3090169943749474,0.951056516295153,0},{0.951056516295153,0.3090169943749474,0},{0,0,-1}};
Rpeter[[33,All,All]]={{0.809016994374947,0.5877852522924732,0},{0.5877852522924732,-0.809016994374947,0},{0,0,-1}};
Rpeter[[34,All,All]]={{-0.3090169943749474,-0.951056516295153,0},{-0.951056516295153,0.3090169943749474,0},{0,0,-1}};
Rpeter[[35,All,All]]={{0.809016994374947,-0.5877852522924732,0},{-0.5877852522924732,-0.809016994374947,0},{0,0,-1}};
Rpeter[[36,All,All]]={{0.1381966011250104,0.4253254041760199,-0.894427190999916},{0.951056516295154,-0.3090169943749474,0},{-0.2763932022500211,-0.85065080835204,-0.447213595499958}};
Rpeter[[37,All,All]]={{0.4472135954999579,0,-0.894427190999916},{0,-1,0},{-0.894427190999916,0,-0.447213595499958}};
Rpeter[[38,All,All]]={{0.1381966011250104,-0.4253254041760199,-0.894427190999916},{-0.951056516295153,-0.3090169943749474,0},{-0.2763932022500212,0.85065080835204,-0.447213595499958}};
Rpeter[[39,All,All]]={{-0.3618033988749896,0.2628655560595667,-0.894427190999916},{0.5877852522924732,0.809016994374947,0},{0.723606797749979,-0.5257311121191335,-0.447213595499958}};
Rpeter[[40,All,All]]={{-0.3618033988749895,-0.2628655560595667,-0.894427190999916},{-0.587785252292473,0.809016994374947,0},{0.723606797749979,0.5257311121191335,-0.447213595499958}};
Rpeter[[41,All,All]]={{0.947213595499958,-0.1624598481164531,-0.2763932022500211},{0.1624598481164533,-0.5,0.85065080835204},{-0.2763932022500211,-0.85065080835204,-0.447213595499958}};
Rpeter[[42,All,All]]={{0.1381966011250105,-0.951056516295153,-0.2763932022500211},{-0.4253254041760198,-0.3090169943749474,0.85065080835204},{-0.894427190999916,0,-0.447213595499958}};
Rpeter[[43,All,All]]={{-0.861803398874989,-0.4253254041760199,-0.2763932022500212},{-0.4253254041760198,0.3090169943749474,0.85065080835204},{-0.2763932022500212,0.85065080835204,-0.447213595499958}};
Rpeter[[44,All,All]]={{0.447213595499958,0.85065080835204,-0.276393202250021},{0.5257311121191338,0,0.85065080835204},{0.723606797749979,-0.5257311121191335,-0.447213595499958}};
Rpeter[[45,All,All]]={{-0.6708203932499368,0.6881909602355866,-0.2763932022500212},{0.1624598481164532,0.5,0.85065080835204},{0.723606797749979,0.5257311121191335,-0.447213595499958}};
Rpeter[[46,All,All]]={{0.4472135954999582,-0.5257311121191335,0.7236067977499788},{-0.85065080835204,0,0.5257311121191336},{-0.2763932022500211,-0.85065080835204,-0.447213595499958}};
Rpeter[[47,All,All]]={{-0.3618033988749893,-0.5877852522924732,0.723606797749979},{-0.2628655560595667,0.809016994374947,0.5257311121191338},{-0.894427190999916,0,-0.447213595499958}};
Rpeter[[48,All,All]]={{-0.6708203932499368,0.1624598481164531,0.7236067977499792},{0.6881909602355867,0.5,0.5257311121191339},{-0.2763932022500212,0.85065080835204,-0.447213595499958}};
Rpeter[[49,All,All]]={{0.6381966011250108,0.2628655560595669,0.723606797749979},{-0.2628655560595667,-0.809016994374947,0.5257311121191338},{0.723606797749979,-0.5257311121191335,-0.447213595499958}};
Rpeter[[50,All,All]]={{-0.052786404500042,0.6881909602355866,0.7236067977499792},{0.6881909602355867,-0.4999999999999999,0.5257311121191339},{0.723606797749979,0.5257311121191335,-0.447213595499958}};
Rpeter[[51,All,All]]={{-0.86180339887499,0.4253254041760198,-0.2763932022500211},{0.4253254041760199,0.3090169943749474,-0.85065080835204},{-0.2763932022500211,-0.85065080835204,-0.447213595499958}};
Rpeter[[52,All,All]]={{0.1381966011250104,0.951056516295153,-0.276393202250021},{0.4253254041760198,-0.3090169943749474,-0.85065080835204},{-0.894427190999916,0,-0.447213595499958}};
Rpeter[[53,All,All]]={{0.947213595499958,0.1624598481164532,-0.2763932022500212},{-0.1624598481164533,-0.5,-0.85065080835204},{-0.2763932022500212,0.85065080835204,-0.447213595499958}};
Rpeter[[54,All,All]]={{-0.6708203932499371,-0.6881909602355867,-0.2763932022500211},{-0.1624598481164532,0.5,-0.85065080835204},{0.723606797749979,-0.5257311121191335,-0.447213595499958}};
Rpeter[[55,All,All]]={{0.4472135954999578,-0.85065080835204,-0.2763932022500212},{-0.5257311121191338,0,-0.85065080835204},{0.723606797749979,0.5257311121191335,-0.447213595499958}};
Rpeter[[56,All,All]]={{-0.670820393249937,-0.1624598481164531,0.7236067977499788},{-0.6881909602355869,0.4999999999999999,-0.5257311121191336},{-0.2763932022500211,-0.85065080835204,-0.447213595499958}};
Rpeter[[57,All,All]]={{-0.3618033988749894,0.5877852522924732,0.723606797749979},{0.2628655560595667,0.809016994374947,-0.5257311121191338},{-0.894427190999916,0,-0.447213595499958}};
Rpeter[[58,All,All]]={{0.447213595499958,0.5257311121191335,0.7236067977499792},{0.85065080835204,0,-0.5257311121191339},{-0.2763932022500212,0.85065080835204,-0.447213595499958}};
Rpeter[[59,All,All]]={{-0.05278640450004211,-0.6881909602355867,0.723606797749979},{-0.6881909602355869,-0.5,-0.5257311121191338},{0.723606797749979,-0.5257311121191335,-0.447213595499958}};
Rpeter[[60,All,All]]={{0.6381966011250108,-0.2628655560595669,0.7236067977499792},{0.2628655560595666,-0.809016994374947,-0.5257311121191339},{0.723606797749979,0.5257311121191335,-0.447213595499958}};
permutation4p2={1, 2, 5, 9, 17, 10, 27, 13, 21, 18, 24, 15, 26, 3, 4, 48, 45, 56, 54, 49, 60, 36, 52, 42, 38, 14, 16, 47, 40, 46, 55, 41, 53, 20, 29, 6, 12, 57, 39, 8, 22, 44, 58, 28, 25, 11, 31, 59, 33, 30, 19, 43, 35, 34, 37, 23, 7, 50, 32, 51};
permutation4p2N={1, 2, 4, 9, 17, 10, 27, 13, 21, 18, 24, 15, 26, 3, 5, 48, 45, 56, 54, 49, 60, 36, 52, 42, 38, 14, 16, 47, 40, 46, 55, 41, 53, 20, 29, 6, 12, 57, 39, 8, 22, 44, 58, 28, 25, 11, 31, 59, 33, 30, 19, 43, 34, 35, 37, 23, 7, 50, 32, 51};
permutation4p3={1, 4, 3, 36, 52, 38, 42, 49, 60, 54, 56, 48, 45, 2, 5, 24, 18, 15, 26, 21, 13, 10, 27, 17, 9, 46, 55, 22, 8, 25, 28, 20, 29, 53, 41, 40, 47, 12, 6, 39, 57, 16, 14, 44, 58, 50, 31, 11, 32, 43, 51, 19, 33, 35, 7, 59, 37, 23, 34, 30};
Rzdpermute=Rpeter[[permutation4p2N,All,All]];


(*See the difference between Rpeter and Ryb For[i=1,i<=60,i++, Print[{i,MatrixForm[1.0*Normal[Ryb[[i,All,All]]]],MatrixForm[Rpeter[[i,All,All]]]}]];*)


Rzdpermute=Rpeter[[permutation4p2N,All,All]];
Rybpermute=Ryb[[permutation4p2,All,All]];
(*For[i=1,i<=60,i++, Print[{i,MatrixForm[Chop[1.0*Normal[Rzdpermute[[i,All,All]]]-1.0*Rybpermute[[i,All,All]]]]}]];*)
For[i=1,i<=60,i++, Print[{i,MatrixForm[Chop[1.0*Normal[Rzdpermute[[i,All,All]]]-1.0*Normal[Inverse[Rzdpermute[[i,All,All]]]]]]}]];


(*Transformation Matrix To Euler Angle Converter*)
Rrose[\[Alpha]_,\[Beta]_,\[Gamma]_]:={{Cos[\[Alpha]]Cos[\[Beta]]Cos[\[Gamma]]-Sin[\[Alpha]]Sin[\[Gamma]], Sin[\[Alpha]]Cos[\[Beta]]Cos[\[Gamma]]+Cos[\[Alpha]]Sin[\[Gamma]], -Sin[\[Beta]]Cos[\[Gamma]]},{-Cos[\[Alpha]]Cos[\[Beta]]Sin[\[Gamma]]-Sin[\[Alpha]]Cos[\[Gamma]],-Sin[\[Alpha]]Cos[\[Beta]]Sin[\[Gamma]]+Cos[\[Alpha]]Cos[\[Gamma]], Sin[\[Beta]]Sin[\[Gamma]]},{Cos[\[Alpha]]Sin[\[Beta]], Sin[\[Alpha]]Sin[\[Beta]], Cos[\[Beta]]}}
MatrixToEulerAngle[M_]:=(Angles=SparseArray[{},{3,1}];
b1=ArcCos[M[[3,3]]]; 
a1=Piecewise[{{ArcTan[M[[1,1]],-M[[2,1]]], b1==0}, {0, b1==Pi}, {ArcTan[M[[3,1]]/Sin[b1],M[[3,2]]/Sin[b1]], b1!=0||Pi}}];
c1=Piecewise[{{0, b1==0}, {ArcTan[-M[[1,1]],M[[2,1]]], b1==Pi}, {ArcTan[M[[1,3]]/(-Sin[b1]),M[[2,3]]/Sin[b1]], b1!=0||Pi}}];
Angles[[1]]=a1;Angles[[2]]=b1;Angles[[3]]=c1; Angles)
Rzdtest=SparseArray[{},{60,3,3}];
\[Alpha]zd=SparseArray[{},{60,1}];
\[Beta]zd=SparseArray[{},{60,1}];
\[Gamma]zd=SparseArray[{},{60,1}];
 \[Epsilon]=10^(-14);
For[g=1, g<=Ng, g++,M=Inverse[Rzdpermute[[g,All,All]]]; a=Chop[MatrixToEulerAngle[M][[1]],\[Epsilon]]; b=Chop[MatrixToEulerAngle[M][[2]],\[Epsilon]]; c=Chop[MatrixToEulerAngle[M][[3]],\[Epsilon]]; \[Alpha]zd[[g]]=a; \[Beta]zd[[g]]=b;  \[Gamma]zd[[g]]=c; Rzdtest[[g,All,All]]=Chop[Rrose[a,b,c],\[Epsilon]]]
(*For[g=1, g<=Ng, g++,M=Inverse[Rzdpermute[[g,All,All]]]; a=Chop[MatrixToEulerAngle[M][[1]],\[Epsilon]]; b=Chop[MatrixToEulerAngle[M][[2]],\[Epsilon]]; c=Chop[MatrixToEulerAngle[M][[3]],\[Epsilon]]; If[a>= 0, \[Alpha]zd[[g]]=a, \[Alpha]zd[[g]]=a+2Pi]; \[Beta]zd[[g]]=b;  If[c>= 0, \[Gamma]zd[[g]]=c, \[Gamma]zd[[g]]=c+2Pi]; Rzdtest[[g,All,All]]=Chop[Rrose[a,b,c],\[Epsilon]]]*)
MatrixForm[{\[Alpha]zd,\[Beta]zd,\[Gamma]zd}]


g=2;
M=Inverse[Rzdpermute[[g,All,All]]];
MatrixToEulerAngle[M]
M=Rzdpermute[[g,All,All]];
MatrixToEulerAngle[M]


\[Alpha]yb=SparseArray[{},{60,1}];
\[Beta]yb=SparseArray[{},{60,1}];
\[Gamma]yb=SparseArray[{},{60,1}];
For[g=1, g<=Ng, g++,M=Rybpermute[[g,All,All]]; a=Chop[EulerAngles[M][[1]],\[Epsilon]]; b=Chop[EulerAngles[M][[2]],\[Epsilon]]; c=Chop[EulerAngles[M][[3]],\[Epsilon]]; \[Alpha]yb[[g]]=a; \[Beta]yb[[g]]=b;  \[Gamma]yb[[g]]=c; Rzdtest[[g,All,All]]=Chop[Rrose[a,b,c],\[Epsilon]]]
MatrixForm[{\[Alpha]yb,\[Beta]yb,\[Gamma]yb}*1.0]


(*Test: Rzdpermute instead of Inverse[Rzdpermute] share the same multiplication table with LiuPing2 and Liuping3
M2=LiuPing2; M3=LiuPing3;
For[g1=1, g1<= Ng, g1++, 
	For[g2=1,g2<=Ng,g2++, 
		g3=IsThere[Rzdpermute,g1,g2];
		R2=M2[[g1, All, All]].M2[[g2, All, All]];
		R3=M3[[g1, All, All]].M3[[g2, All, All]];	
        If[R3!= M3[[g3, All, All]], Print[{3, g1, g2, g3}]];
		If[R2!=M2[[g3, All, All]], Print[{2, g1, g2, g3}]];	
]]*)


(*Test: Rotated Spherical Harmonics= corresponding WigD times unrotated ones*)
For[l=1,l<=4, l++, Print[l]; 
For[g=1, g<=Ng, g++,
Ry=Y[Inverse[Rzdpermute[[g,All,All]]].x,l,Table[m,{m, -l,l}]];
y=Y[x,l,Table[m,{m, -l,l}]];
MatrixForm[MWigD=Table[WignerD[{l, mp, m}, -\[Alpha]zd[[g]], -\[Beta]zd[[g]], -\[Gamma]zd[[g]]],{m,-l,l},{mp,-l,l}]];
(*MatrixForm[MWigD=Table[WignerD[{l, mp, m}, \[Alpha]zd[[g]], \[Beta]zd[[g]], \[Gamma]zd[[g]]],{m,-l,l},{mp,-l,l}]];*)
(*MatrixForm[MWigD=Table[WignerD[{l, mp, m}, -\[Gamma]zd[[g]], -\[Beta]zd[[g]], -\[Alpha]zd[[g]]],{m,-l,l},{mp,-l,l}]];*)
If[Chop[MWigD.y-Ry,\[Epsilon]]!=ConstantArray[0, 2l+1], Print[{l, MWigD.y, Ry}]]]]


(*Spherical Harmonics*)
Clear[d]
Y[x_,l_,m_]:=( r=Sqrt[x[[1]]^2+x[[2]]^2+x[[3]]^2];  \[Theta]=ArcCos[x[[3]]/r];\[Phi]=ArcTan[x[[1]],x[[2]]] ;SphericalHarmonicY[l,m,\[Theta],\[Phi]])
(*D coefficiets*)
powerFunc=Unevaluated[#1^#2]/.HoldPattern[0^0]:>1&;
d[l_,mp_,m_,\[Beta]_]:=Sqrt[(l+m)!*(l-m)!*(l+mp)!*(l-mp)!]*Sum[(-1)^k/((l-mp-k)!*(l+m-k)!*(k+mp-m)!*k!)*powerFunc[(Cos[\[Beta]/2]),(2l+m-mp-2k)]*powerFunc[(-Sin[\[Beta]/2]),(mp-m+2k)],{k,Max[0,(m-mp)],Min[(l-mp),(l+m)]}]
Dwig[l_,mp_,m_,\[Alpha]_,\[Beta]_,\[Gamma]_]:=Exp[-I*mp*\[Alpha]]*d[l,mp,m,\[Beta]]*Exp[-I*m*\[Gamma]]
(*Group Projection Operator*)
\[CapitalDelta][l_,mp_,m_,p_,i_,j_]:=dim[[p]]/Ng*Sum[ (Conjugate[\[CapitalGamma]r[p,g][[i,j]]])*WignerD[{l, mp, m}, -\[Alpha]zd[[g]], -\[Beta]zd[[g]],-\[Gamma]zd[[g]]],{g,Ng}] 
(*\[CapitalDelta][l_,mp_,m_,p_,i_,j_]:=dim[[p]]/Ng*Sum[ (Conjugate[\[CapitalGamma][p,g][[i,j]]])*WignerD[{l, mp, m}, \[Gamma]zd[[g]], \[Beta]zd[[g]], \[Alpha]zd[[g]]],{g,Ng}] *)
(*\[CapitalDelta][l_,mp_,m_,p_,i_,j_]:=dim[[p]]/Ng*Sum[ (\[CapitalGamma]r[p,g][[i,j]])*Dwig[l,mp,m,\[Alpha]zd[[g]],\[Beta]zd[[g]],\[Gamma]zd[[g]]],{g,Ng}] Tested finally-Oct 12, 2014 !!!*)
ProjOperator[l_,m_,p_,i_,j_,x_]:=Chop[Table[\[CapitalDelta][l,mp,m,p,i,j],{mp,-l,l}].Y[x,l,Table[k,{k,-l,l}]],\[Epsilon]]
ProjOperatorTest[l_,m_,p_,i_,j_,x_]:=dim[[p]]/Ng*Chop[Sum[Conjugate[\[CapitalGamma]r[p,g][[i,j]]]*Y[Inverse[Rzdpermute[[g,All,All]]].x,l,m],{g,Ng}],\[Epsilon]]
LengthProjOperator[l_,m_,p_,n_]:=Chop[Sum[Abs[\[CapitalDelta][l,mp,m,p,n,n]]^2,{mp, -l,l}],\[Epsilon]]
BasisFunctionMatrix[l_,m_,p]:=(\[Epsilon]=10^(-14);Matrix=SparseArray[{},{dim[[p]],2l+1}]; n=0;
	For[k=1,k<=dim[[p]],k++, 
		cn=LengthProjOperator[l,m,p,k];
		If[cn!=0,n=k; Break[]]
	];
	If[n!=0, 
		Matrix[[n,All]]=Chop[Table[\[CapitalDelta][l,mp,m,p,n,n],{mp,-l,l}]/cn,\[Epsilon]];
		For[np=1,np<=dim[[p]],np++, 
			If[np!=n,
				Matrix[[np,All]]=Chop[Table[\[CapitalDelta][l,mp,m,p,np,n],{mp,-l,l}]/cn,\[Epsilon]],
				Continue[]
			]
		],
	cn=0];
	Normal[Matrix])
BasisRealFunctionCoeffMatrix[l_,p_,m_]:=(Mat=BasisFunctionMatrix[l,m,p];
		afunc[mp_]:=Mat[[All,mp+l+1]]; lh=0; rh=0; 
		For[mp=-l, mp<=-1, mp++, lh=lh+afunc[mp]; rh=rh+(-1)^Abs[mp]*Conjugate[afunc[-mp]];];
		angle=0;
		For[i=1,i<=Length[rh],i++,If[lh[[i]]==0,Continue[],Break[]];];
		If[i<= Length[rh], angle=Log[rh[[i]]/lh[[i]]]/2];
	MatNew=Exp[angle]*Mat)



x=xvec[[1,All]]
l=6; p=1; 
Y=Table[Y[x,l,m],{m,-l,l}]
Dmat=Table[BasisFunctionMatrix[l,m,p],{m,-l,l}]



