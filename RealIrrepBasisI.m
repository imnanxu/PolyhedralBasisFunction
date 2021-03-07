(* ::Package:: *)

(*Written by Nan Xu Oct 12, 2015*)
(*Modified by Nan Xu Mar 02, 2021*)
(*ClearAll["Global`*"]*)
BeginPackage["RealIrrepBasisI`"]

  \[CapitalGamma]0::usage="\[CapitalGamma]0[p,g] returns the complex irrep matrix of the icosahedral group element g (g=1,2,...,60) of the p'th (p=1,2,3,4,5) irrep.";

  \[CapitalGamma]r::usage="\[CapitalGamma]r[p,g] returns the equivalent real irrep matrix of the icosahedral group element g (g=1,2,...,60) of the p'th (p=1,2,3,4,5) irrep.";

  BasisRealFunctionCoeffMatrixI::usage="BasisRealFunctionCoeffMatrixI[l,m,p] computes the non-orthogonalized coefficents matrix of the icosahedral group that associates with the p'th irrep (p=1,2,3,4,5), the spherical harmonics degree l (l>=0) and order m (-l<=m<=l).";
  
  BasisRealFunctionOthoCoeffMatrixI::usage="BasisRealFunctionOthoCoeffMatrixI[l,p] computes the Gram\[Dash]Schmidt orthogonalized coefficents matrix of the icosahedral group that associates with the p'th irrep (p=1,2,3,4,5), the spherical harmonics degree l (l>=0).";

  Begin["`Private`"]



Ng=60;
dim={1,3,3,4,5};


(*%Fa Liu, Jia-Lun Ping, Jin-Quan Chen, "Application of the eigenfunction method to the icosahedral group", J. Math. Phys 31(5):1065-1075, May 1990.*)
(*parameters setup*)
z=Exp[Sqrt[-1]*2*Pi/5];
P=(Sqrt[5]-1)/2;
Q=-(Sqrt[5]+1)/2;
A=1/Sqrt[5];
B=-A*Q;
myD=-Sqrt[2]*A;
myE=A*P;
F=-A^2;
G=B^2;
H=A^2*P^2;
L=Sqrt[6]*A^2;
M=-2*A^2*P;
myN=-2*A^2*Q;
a=1/z;
b=-B*z^2;
d=1/z^2;
e=-myE*z;
f=myD*z^2;
g=myD/z;
h=myE*z^2;
j=B/z;
k=-A/z;
m=A*z^2;
n=G/z;
p=H*z^2;
q=myN/z;
r=myN/z^2;
s=M*z^2;
t=M*z;
u=L/z;
v=L/z^2;
w=H/z;
y=G*z^2;
nrows=1+9+9+16+25;
ncols=30;
part1=SparseArray[{},{nrows,ncols}];
part2=SparseArray[{},{nrows,ncols}];

part1[[1, All]]=1;

part1[[2, All]]={
1, a, Conjugate[a], b, Conjugate[b], Conjugate[b], b, Conjugate[b], b, b, Conjugate[b], b, Conjugate[b], d, Conjugate[d], e, Conjugate[e], Conjugate[e], e, Conjugate[e], e, e, Conjugate[e], e, Conjugate[e], -Conjugate[j], -j, -h, -Conjugate[h], -Conjugate[h]};

part1[[3, All]]={
0, 0, 0, f, myD, Conjugate[f], myD, g, Conjugate[g], g, f, Conjugate[f], Conjugate[g], 0, 0, -g, -Conjugate[f], -Conjugate[g], -f, -f, -Conjugate[f], -Conjugate[g], -myD, -myD, -g, f, Conjugate[g], -myD, -Conjugate[f], -Conjugate[g]};

part1[[4, All]]={
0, 0, 0, h, h, Conjugate[h], Conjugate[h], myE, myE, -e, -e, -Conjugate[e], -Conjugate[e], 0, 0, -b, -b, -Conjugate[b], -Conjugate[b], B, B, Conjugate[j], Conjugate[j], j, j, Conjugate[h], Conjugate[h], -Conjugate[b], -Conjugate[b], j};

part1[[5, All]]={
0, 0, 0, myD, Conjugate[f], myD, f, g, Conjugate[g], Conjugate[f], Conjugate[g], g, f, 0, 0, -f, -Conjugate[g], -Conjugate[f], -g, -f, -Conjugate[f], -myD, -g, -Conjugate[g], -myD, g, Conjugate[f], -f, -myD, -f};

part1[[6, All]]={
1, 1, 1, A, A, A, A, A, A, A, A, A, A, 1, 1, -A, -A, -A, -A, -A, -A, -A, -A, -A, -A, A, A, -A, -A, -A};

part1[[7, All]]={
0, 0, 0, myD, f, myD, Conjugate[f], Conjugate[g], g, f, g, Conjugate[g], Conjugate[f], 0, 0, -Conjugate[f], -g, -f, -Conjugate[g], -Conjugate[f], -f, -myD, -Conjugate[g], -g, -myD, Conjugate[g], f, -Conjugate[f], -myD, -Conjugate[f]};

part1[[8, All]]={
0, 0, 0, Conjugate[h], Conjugate[h], h, h, myE, myE, -Conjugate[e], -Conjugate[e], -e, -e, 0, 0, -Conjugate[b], -Conjugate[b], -b, -b, B, B, j, j, Conjugate[j], Conjugate[j], h, h, -b, -b, Conjugate[j]};

part1[[9, All]]={
0, 0, 0, Conjugate[f], myD, f, myD, Conjugate[g], g, Conjugate[g], Conjugate[f], f, g, 0, 0, -Conjugate[g], -f, -g, -Conjugate[f], -Conjugate[f], -f, -g, -myD, -myD, -Conjugate[g], Conjugate[f], g, -myD, -f, -g};

part1[[10, All]]={
1, Conjugate[a], a, Conjugate[b], b, b, Conjugate[b], b, Conjugate[b], Conjugate[b], b, Conjugate[b], b, Conjugate[d], d, Conjugate[e], e, e, Conjugate[e], e, Conjugate[e], Conjugate[e], e, Conjugate[e], e, -j, -Conjugate[j], -Conjugate[h], -h, -h};

part1[[11, All]]={
1, d, Conjugate[d], Conjugate[e], e, e, Conjugate[e], e, Conjugate[e], Conjugate[e], e, Conjugate[e], e, Conjugate[a], a, b, Conjugate[b], Conjugate[b], b, Conjugate[b], b, b, Conjugate[b], b, Conjugate[b], -h, -Conjugate[h], -j, -Conjugate[j], -Conjugate[j]};

part1[[12, All]]={
0, 0, 0, -g, -myD, -Conjugate[g], -myD, -Conjugate[f], -f, -Conjugate[f], -g, -Conjugate[g], -f, 0, 0, Conjugate[f], Conjugate[g], f, g, g, Conjugate[g], f, myD, myD, Conjugate[f], -g, -f, myD, Conjugate[g], f};

part1[[13, All]]={
0, 0, 0, j, j, Conjugate[j], Conjugate[j], B, B, -b, -b, -Conjugate[b], -Conjugate[b], 0, 0, -Conjugate[e], -Conjugate[e], -e, -e, myE, myE, h, h, Conjugate[h], Conjugate[h], Conjugate[j], Conjugate[j], -e, -e, Conjugate[h]};

part1[[14, All]]={
0, 0, 0, -myD, -Conjugate[g], -myD, -g, -Conjugate[f], -f, -Conjugate[g], -f, -Conjugate[f], -g, 0, 0, g, f, Conjugate[g], Conjugate[f], g, Conjugate[g], myD, Conjugate[f], f, myD, -Conjugate[f], -Conjugate[g], g, myD, g};

part1[[15, All]]={
1, 1, 1, -A, -A, -A, -A, -A, -A, -A, -A, -A, -A, 1, 1, A, A, A, A, A, A, A, A, A, A, -A, -A, A, A, A};

part1[[16, All]]={
0, 0, 0, -myD, -g, -myD, -Conjugate[g], -f, -Conjugate[f], -g, -Conjugate[f], -f, -Conjugate[g], 0, 0, Conjugate[g], Conjugate[f], g, f, Conjugate[g], g, myD, f, Conjugate[f], myD, -f, -g, Conjugate[g], myD, Conjugate[g]};

part1[[17, All]]={
0, 0, 0, Conjugate[j], Conjugate[j], j, j, B, B, -Conjugate[b], -Conjugate[b], -b, -b, 0, 0, -e, -e, -Conjugate[e], -Conjugate[e], myE, myE, Conjugate[h], Conjugate[h], h, h, j, j, -Conjugate[e], -Conjugate[e], h};

part1[[18, All]]={
0, 0, 0, -Conjugate[g], -myD, -g, -myD, -f, -Conjugate[f], -f, -Conjugate[g], -g, -Conjugate[f], 0, 0, f, g, Conjugate[f], Conjugate[g], Conjugate[g], g, Conjugate[f], myD, myD, f, -Conjugate[g], -Conjugate[f], myD, g, Conjugate[f]};

part1[[19, All]]={
1, Conjugate[d], d, e, Conjugate[e], Conjugate[e], e, Conjugate[e], e, e, Conjugate[e], e, Conjugate[e], a, Conjugate[a], Conjugate[b], b, b, Conjugate[b], b, Conjugate[b], Conjugate[b], b, Conjugate[b], b, -Conjugate[h], -h, -Conjugate[j], -j, -j};

part1[[20, All]]={
1, d, Conjugate[d], k, Conjugate[k], Conjugate[k], k, Conjugate[k], k, k, Conjugate[k], k, Conjugate[k], Conjugate[a], a, m, Conjugate[m], Conjugate[m], m, Conjugate[m], m, m, Conjugate[m], m, Conjugate[m], -m, -Conjugate[m], -k, -Conjugate[k], -Conjugate[k]};

part1[[21, All]]={
0, 0, 0, Conjugate[e], -Conjugate[h], e, -h, -h, -Conjugate[h], e, -myE, -myE, Conjugate[e], 0, 0, -B, b, -B, Conjugate[b], -Conjugate[j], -j, b, -j, -Conjugate[j], Conjugate[b], -Conjugate[h], -myE, b, -Conjugate[j], -j};

part1[[22, All]]={
0, 0, 0, j, -b, Conjugate[j], -Conjugate[b], j, Conjugate[j], B, -Conjugate[b], -b, B, 0, 0, -e, myE, -Conjugate[e], myE, h, Conjugate[h], h, -e, -Conjugate[e], Conjugate[h], B, j, Conjugate[h], -e, myE};

part1[[23, All]]={
0, 0, 0, k, k, Conjugate[k], Conjugate[k], -A, -A, -m, -m, -Conjugate[m], -Conjugate[m], 0, 0, -k, -k, -Conjugate[k], -Conjugate[k], A, A, m, m, Conjugate[m], Conjugate[m], Conjugate[k], Conjugate[k], -Conjugate[k], -Conjugate[k], Conjugate[m]};

part1[[24, All]]={
0, 0, 0, -h, e, -Conjugate[h], Conjugate[e], -h, -Conjugate[h], -myE, Conjugate[e], e, -myE, 0, 0, Conjugate[b], -B, b, -B, -Conjugate[j], -j, -Conjugate[j], Conjugate[b], b, -j, -myE, -h, -j, Conjugate[b], -B};

part1[[25, All]]={
1, a, Conjugate[a], m, Conjugate[m], Conjugate[m], m, Conjugate[m], m, m, Conjugate[m], m, Conjugate[m], d, Conjugate[d], Conjugate[k], k, k, Conjugate[k], k, Conjugate[k], Conjugate[k], k, Conjugate[k], k, -Conjugate[k], -k, -m, -Conjugate[m], -Conjugate[m]};

part1[[26, All]]={
0, 0, 0, m, m, Conjugate[m], Conjugate[m], A, A, -Conjugate[k], -Conjugate[k], -k, -k, 0, 0, -m, -m, -Conjugate[m], -Conjugate[m], -A, -A, Conjugate[k], Conjugate[k], k, k, Conjugate[m], Conjugate[m], -Conjugate[m], -Conjugate[m], k};

part1[[27, All]]={
0, 0, 0, -b, j, -Conjugate[b], Conjugate[j], Conjugate[j], j, -Conjugate[b], B, B, -b, 0, 0, myE, -e, myE, -Conjugate[e], Conjugate[h], h, -e, h, Conjugate[h], -Conjugate[e], j, B, -e, Conjugate[h], h};

part1[[28, All]]={
0, 0, 0, -Conjugate[b], Conjugate[j], -b, j, j, Conjugate[j], -b, B, B, -Conjugate[b], 0, 0, myE, -Conjugate[e], myE, -e, h, Conjugate[h], -Conjugate[e], Conjugate[h], h, -e, Conjugate[j], B, -Conjugate[e], h, Conjugate[h]};

part1[[29, All]]={
0, 0, 0, Conjugate[m], Conjugate[m], m, m, A, A, -k, -k, -Conjugate[k], -Conjugate[k], 0, 0, -Conjugate[m], -Conjugate[m], -m, -m, -A, -A, k, k, Conjugate[k], Conjugate[k], m, m, -m, -m, Conjugate[k]};

part1[[30, All]]={
1, Conjugate[a], a, Conjugate[m], m, m, Conjugate[m], m, Conjugate[m], Conjugate[m], m, Conjugate[m], m, Conjugate[d], d, k, Conjugate[k], Conjugate[k], k, Conjugate[k], k, k, Conjugate[k], k, Conjugate[k], -k, -Conjugate[k], -Conjugate[m], -m, -m};

part1[[31, All]]={
0, 0, 0, -Conjugate[h], Conjugate[e], -h, e, -Conjugate[h], -h, -myE, e, Conjugate[e], -myE, 0, 0, b, -B, Conjugate[b], -B, -j, -Conjugate[j], -j, b, Conjugate[b], -Conjugate[j], -myE, -Conjugate[h], -Conjugate[j], b, -B};

part1[[32, All]]={
0, 0, 0, Conjugate[k], Conjugate[k], k, k, -A, -A, -Conjugate[m], -Conjugate[m], -m, -m, 0, 0, -Conjugate[k], -Conjugate[k], -k, -k, A, A, Conjugate[m], Conjugate[m], m, m, k, k, -k, -k, m};

part1[[33, All]]={
0, 0, 0, Conjugate[j], -Conjugate[b], j, -b, Conjugate[j], j, B, -b, -Conjugate[b], B, 0, 0, -Conjugate[e], myE, -e, myE, Conjugate[h], h, Conjugate[h], -Conjugate[e], -e, h, B, Conjugate[j], h, -Conjugate[e], myE};

part1[[34, All]]={
0, 0, 0, e, -h, Conjugate[e], -Conjugate[h], -Conjugate[h], -h, Conjugate[e], -myE, -myE, e, 0, 0, -B, Conjugate[b], -B, b, -j, -Conjugate[j], Conjugate[b], -Conjugate[j], -j, b, -h, -myE, Conjugate[b], -j, -Conjugate[j]};

part1[[35, All]]={
1, Conjugate[d], d, Conjugate[k], k, k, Conjugate[k], k, Conjugate[k], Conjugate[k], k, Conjugate[k], k, a, Conjugate[a], Conjugate[m], m, m, Conjugate[m], m, Conjugate[m], Conjugate[m], m, Conjugate[m], m, -Conjugate[m], -m, -Conjugate[k], -k, -k};

part1[[36, All]]={
1, d, Conjugate[d], n, Conjugate[n], Conjugate[n], n, Conjugate[n], n, n, Conjugate[n], n, Conjugate[n], Conjugate[a], a, p, Conjugate[p], Conjugate[p], p, Conjugate[p], p, p, Conjugate[p], p, Conjugate[p], y, Conjugate[y], w, Conjugate[w], Conjugate[w]};

part1[[37, All]]={
0, 0, 0, q, r, Conjugate[q], Conjugate[r], Conjugate[r], r, Conjugate[q], myN, myN, q, 0, 0, M, s, M, Conjugate[s], t, Conjugate[t], s, Conjugate[t], t, Conjugate[s], r, myN, s, t, Conjugate[t]};

part1[[38, All]]={
0, 0, 0, u, L, Conjugate[u], L, v, Conjugate[v], v, u, Conjugate[u], Conjugate[v], 0, 0, v, Conjugate[u], Conjugate[v], u, u, Conjugate[u], Conjugate[v], L, L, v, u, Conjugate[v], L, Conjugate[u], Conjugate[v]};

part1[[39, All]]={
0, 0, 0, Conjugate[t], s, t, Conjugate[s], Conjugate[t], t, M, Conjugate[s], s, M, 0, 0, Conjugate[q], myN, q, myN, Conjugate[r], r, Conjugate[r], Conjugate[q], q, r, M, Conjugate[t], r, Conjugate[q], myN};

part1[[40, All]]={
0, 0, 0, w, w, Conjugate[w], Conjugate[w], H, H, p, p, Conjugate[p], Conjugate[p], 0, 0, n, n, Conjugate[n], Conjugate[n], G, G, y, y, Conjugate[y], Conjugate[y], Conjugate[w], Conjugate[w], Conjugate[n], Conjugate[n], Conjugate[y]};

part1[[41, All]]={
0, 0, 0, Conjugate[r], Conjugate[q], r, q, Conjugate[r], r, myN, q, Conjugate[q], myN, 0, 0, Conjugate[s], M, s, M, t, Conjugate[t], t, Conjugate[s], s, Conjugate[t], myN, Conjugate[r], Conjugate[t], Conjugate[s], M};

part1[[42, All]]={
1, a, Conjugate[a], p, Conjugate[p], Conjugate[p], p, Conjugate[p], p, p, Conjugate[p], p, Conjugate[p], d, Conjugate[d], Conjugate[n], n, n, Conjugate[n], n, Conjugate[n], Conjugate[n], n, Conjugate[n], n, Conjugate[w], w, y, Conjugate[y], Conjugate[y]};

part1[[43, All]]={
0, 0, 0, -Conjugate[v], -L, -v, -L, -u, -Conjugate[u], -u, -Conjugate[v], -v, -Conjugate[u], 0, 0, -u, -v, -Conjugate[u], -Conjugate[v], -Conjugate[v], -v, -Conjugate[u], -L, -L, -u, -Conjugate[v], -Conjugate[u], -L, -v, -Conjugate[u]};

part1[[44, All]]={
0, 0, 0, y, y, Conjugate[y], Conjugate[y], G, G, Conjugate[n], Conjugate[n], n, n, 0, 0, p, p, Conjugate[p], Conjugate[p], H, H, Conjugate[w], Conjugate[w], w, w, Conjugate[y], Conjugate[y], Conjugate[p], Conjugate[p], w};

part1[[45, All]]={
0, 0, 0, s, Conjugate[t], Conjugate[s], t, t, Conjugate[t], Conjugate[s], M, M, s, 0, 0, myN, Conjugate[q], myN, q, r, Conjugate[r], Conjugate[q], Conjugate[r], r, q, Conjugate[t], M, Conjugate[q], r, Conjugate[r]};

part1[[46, All]]={
0, 0, 0, L, Conjugate[u], L, u, v, Conjugate[v], Conjugate[u], Conjugate[v], v, u, 0, 0, u, Conjugate[v], Conjugate[u], v, u, Conjugate[u], L, v, Conjugate[v], L, v, Conjugate[u], u, L, u};

part1[[47, All]]={
0, 0, 0, -L, -v, -L, -Conjugate[v], -u, -Conjugate[u], -v, -Conjugate[u], -u, -Conjugate[v], 0, 0, -Conjugate[v], -Conjugate[u], -v, -u, -Conjugate[v], -v, -L, -u, -Conjugate[u], -L, -u, -v, -Conjugate[v], -L, -Conjugate[v]};

part1[[48, All]]={
1, 1, 1, F, F, F, F, F, F, F, F, F, F, 1, 1, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F};

part1[[49, All]]={
0, 0, 0, -L, -Conjugate[v], -L, -v, -Conjugate[u], -u, -Conjugate[v], -u, -Conjugate[u], -v, 0, 0, -v, -u, -Conjugate[v], -Conjugate[u], -v, -Conjugate[v], -L, -Conjugate[u], -u, -L, -Conjugate[u], -Conjugate[v], -v, -L, -v};

part1[[50, All]]={
0, 0, 0, L, u, L, Conjugate[u], Conjugate[v], v, u, v, Conjugate[v], Conjugate[u], 0, 0, Conjugate[u], v, u, Conjugate[v], Conjugate[u], u, L, Conjugate[v], v, L, Conjugate[v], u, Conjugate[u], L, Conjugate[u]};

part1[[51, All]]={
0, 0, 0, Conjugate[s], t, s, Conjugate[t], Conjugate[t], t, s, M, M, Conjugate[s], 0, 0, myN, q, myN, Conjugate[q], Conjugate[r], r, q, r, Conjugate[r], Conjugate[q], t, M, q, Conjugate[r], r};

part1[[52, All]]={
0, 0, 0, Conjugate[y], Conjugate[y], y, y, G, G, n, n, Conjugate[n], Conjugate[n], 0, 0, Conjugate[p], Conjugate[p], p, p, H, H, w, w, Conjugate[w], Conjugate[w], y, y, p, p, Conjugate[w]};

part1[[53, All]]={
0, 0, 0, -v, -L, -Conjugate[v], -L, -Conjugate[u], -u, -Conjugate[u], -v, -Conjugate[v], -u, 0, 0, -Conjugate[u], -Conjugate[v], -u, -v, -v, -Conjugate[v], -u, -L, -L, -Conjugate[u], -v, -u, -L, -Conjugate[v], -u};

part1[[54, All]]={
1, Conjugate[a], a, Conjugate[p], p, p, Conjugate[p], p, Conjugate[p], Conjugate[p], p, Conjugate[p], p, Conjugate[d], d, n, Conjugate[n], Conjugate[n], n, Conjugate[n], n, n, Conjugate[n], n, Conjugate[n], w, Conjugate[w], Conjugate[y], y, y};

part1[[55, All]]={
0, 0, 0, r, q, Conjugate[r], Conjugate[q], r, Conjugate[r], myN, Conjugate[q], q, myN, 0, 0, s, M, Conjugate[s], M, Conjugate[t], t, Conjugate[t], s, Conjugate[s], t, myN, r, t, s, M};

part1[[56, All]]={
0, 0, 0, Conjugate[w], Conjugate[w], w, w, H, H, Conjugate[p], Conjugate[p], p, p, 0, 0, Conjugate[n], Conjugate[n], n, n, G, G, Conjugate[y], Conjugate[y], y, y, w, w, n, n, y};

part1[[57, All]]={
0, 0, 0, t, Conjugate[s], Conjugate[t], s, t, Conjugate[t], M, s, Conjugate[s], M, 0, 0, q, myN, Conjugate[q], myN, r, Conjugate[r], r, q, Conjugate[q], Conjugate[r], M, t, Conjugate[r], q, myN};

part1[[58, All]]={
0, 0, 0, Conjugate[u], L, u, L, Conjugate[v], v, Conjugate[v], Conjugate[u], u, v, 0, 0, Conjugate[v], u, v, Conjugate[u], Conjugate[u], u, v, L, L, Conjugate[v], Conjugate[u], v, L, u, v};

part1[[59, All]]={
0, 0, 0, Conjugate[q], Conjugate[r], q, r, r, Conjugate[r], q, myN, myN, Conjugate[q], 0, 0, M, Conjugate[s], M, s, Conjugate[t], t, Conjugate[s], t, Conjugate[t], s, Conjugate[r], myN, Conjugate[s], Conjugate[t], t};

part1[[60, All]]={
1, Conjugate[d], d, Conjugate[n], n, n, Conjugate[n], n, Conjugate[n], Conjugate[n], n, Conjugate[n], n, a, Conjugate[a], Conjugate[p], p, p, Conjugate[p], p, Conjugate[p], Conjugate[p], p, Conjugate[p], p, Conjugate[y], y, Conjugate[w], w, w};

part2[[1, All]]=1;

part2[[2, All]]={
-h, -h, -Conjugate[h], -Conjugate[j], -j, -Conjugate[j], -j, -Conjugate[h], -h, -j, -Conjugate[j], -Conjugate[h], -h, -Conjugate[j], -j, -B, 0, -myE, 0, -B, -B, -myE, 0, 0, -myE, -B, -B, -myE, 0, -myE};

part2[[3, All]]={
-Conjugate[f], -Conjugate[g], -g, Conjugate[f], f, Conjugate[g], myD, -myD, -f, g, myD, -f, -g, g, Conjugate[f], Conjugate[g], 0, -f, 0, Conjugate[f], f, -g, 0, 0, -myD, g, myD, -Conjugate[f], 0, -Conjugate[g]};

part2[[4, All]]={
j, B, B, myE, myE, -e, -e, -b, -b, -Conjugate[e], -Conjugate[e], Conjugate[j], Conjugate[j], h, h, h, -1, j, -d, -e, -Conjugate[e], -Conjugate[b], -Conjugate[a], -Conjugate[d], B, Conjugate[h], myE, Conjugate[j], -a, -b};

part2[[5, All]]={
-g, -Conjugate[g], -g, Conjugate[f], f, myD, g, -Conjugate[f], -myD, myD, Conjugate[g], -Conjugate[g], -Conjugate[f], f, Conjugate[g], g, 0, -Conjugate[f], 0, f, Conjugate[f], -Conjugate[g], 0, 0, -myD, Conjugate[g], myD, -f, 0, -g};

part2[[6, All]]={
-A, -A, -A, A, A, A, A, -A, -A, A, A, -A, -A, A, A, A, -1, -A, -1, A, A, -A, -1, -1, -A, A, A, -A, -1, -A};

part2[[7, All]]={
-Conjugate[g], -g, -Conjugate[g], f, Conjugate[f], myD, Conjugate[g], -f, -myD, myD, g, -g, -f, Conjugate[f], g, Conjugate[g], 0, -f, 0, Conjugate[f], f, -g, 0, 0, -myD, g, myD, -Conjugate[f], 0, -Conjugate[g]};

part2[[8, All]]={
Conjugate[j], B, B, myE, myE, -Conjugate[e], -Conjugate[e], -Conjugate[b], -Conjugate[b], -e, -e, j, j, Conjugate[h], Conjugate[h], Conjugate[h], -1, Conjugate[j], -Conjugate[d], -Conjugate[e], -e, -b, -a, -d, B, h, myE, j, -Conjugate[a], -Conjugate[b]};

part2[[9, All]]={
-f, -g, -Conjugate[g], f, Conjugate[f], g, myD, -myD, -Conjugate[f], Conjugate[g], myD, -Conjugate[f], -Conjugate[g], Conjugate[g], f, g, 0, -Conjugate[f], 0, f, Conjugate[f], -Conjugate[g], 0, 0, -myD, Conjugate[g], myD, -f, 0, -g};

part2[[10, All]]={
-Conjugate[h], -Conjugate[h], -h, -j, -Conjugate[j], -j, -Conjugate[j], -h, -Conjugate[h], -Conjugate[j], -j, -h, -Conjugate[h], -j, -Conjugate[j], -B, 0, -myE, 0, -B, -B, -myE, 0, 0, -myE, -B, -B, -myE, 0, -myE};

part2[[11, All]]={
-j, -j, -Conjugate[j], -h, -Conjugate[h], -h, -Conjugate[h], -Conjugate[j], -j, -Conjugate[h], -h, -Conjugate[j], -j, -h, -Conjugate[h], -myE, 0, -B, 0, -myE, -myE, -B, 0, 0, -B, -myE, -myE, -B, 0, -B};

part2[[12, All]]={
Conjugate[g], f, Conjugate[f], -Conjugate[g], -g, -f, -myD, myD, g, -Conjugate[f], -myD, g, Conjugate[f], -Conjugate[f], -Conjugate[g], -f, 0, g, 0, -Conjugate[g], -g, Conjugate[f], 0, 0, myD, -Conjugate[f], -myD, Conjugate[g], 0, f};

part2[[13, All]]={
Conjugate[h], myE, myE, B, B, -b, -b, -Conjugate[e], -Conjugate[e], -Conjugate[b], -Conjugate[b], h, h, j, j, j, -1, Conjugate[h], -Conjugate[a], -b, -Conjugate[b], -e, -Conjugate[d], -a, myE, Conjugate[j], B, h, -d, -Conjugate[e]};

part2[[14, All]]={
Conjugate[f], f, Conjugate[f], -Conjugate[g], -g, -myD, -Conjugate[f], Conjugate[g], myD, -myD, -f, f, Conjugate[g], -g, -f, -Conjugate[f], 0, Conjugate[g], 0, -g, -Conjugate[g], f, 0, 0, myD, -f, -myD, g, 0, Conjugate[f]};

part2[[15, All]]={
A, A, A, -A, -A, -A, -A, A, A, -A, -A, A, A, -A, -A, -A, -1, A, -1, -A, -A, A, -1, -1, A, -A, -A, A, -1, A
};

part2[[16, All]]={
f, Conjugate[f], f, -g, -Conjugate[g], -myD, -f, g, myD, -myD, -Conjugate[f], Conjugate[f], g, -Conjugate[g], -Conjugate[f], -f, 0, g, 0, -Conjugate[g], -g, Conjugate[f], 0, 0, myD, -Conjugate[f], -myD, Conjugate[g], 0, f};

part2[[17, All]]={
h, myE, myE, B, B, -Conjugate[b], -Conjugate[b], -e, -e, -b, -b, Conjugate[h], Conjugate[h], Conjugate[j], Conjugate[j], Conjugate[j], -1, h, -a, -Conjugate[b], -b, -Conjugate[e], -d, -Conjugate[a], myE, j, B, Conjugate[h], -Conjugate[d], -e};

part2[[18, All]]={
g, Conjugate[f], f, -g, -Conjugate[g], -Conjugate[f], -myD, myD, Conjugate[g], -f, -myD, Conjugate[g], f, -f, -g, -Conjugate[f], 0, Conjugate[g], 0, -g, -Conjugate[g], f, 0, 0, myD, -f, -myD, g, 0, Conjugate[f]};

part2[[19, All]]={
-Conjugate[j], -Conjugate[j], -j, -Conjugate[h], -h, -Conjugate[h], -h, -j, -Conjugate[j], -h, -Conjugate[h], -j, -Conjugate[j], -Conjugate[h], -h, -myE, 0, -B, 0, -myE, -myE, -B, 0, 0, -B, -myE, -myE, -B, 0, -B};

part2[[20, All]]={
-k, -k, -Conjugate[k], -m, -Conjugate[m], -m, -Conjugate[m], -Conjugate[k], -k, -Conjugate[m], -m, -Conjugate[k], -k, -m, -Conjugate[m], -A, 0, A, 0, -A, -A, A, 0, 0, A, -A, -A, A, 0, A};

part2[[21, All]]={
-B, Conjugate[b], b, Conjugate[e], e, -h, Conjugate[e], Conjugate[b], -j, -Conjugate[h], e, -B, -Conjugate[j], -myE, -h, e, 0, b, 0, -Conjugate[h], -h, -j, 0, 0, -B, Conjugate[e], -myE, Conjugate[b], 0, -Conjugate[j]};

part2[[22, All]]={
h, -e, -Conjugate[e], -Conjugate[b], -b, -b, Conjugate[j], h, -Conjugate[e], -Conjugate[b], j, Conjugate[h], myE, Conjugate[j], B, -Conjugate[b], 0, -e, 0, j, Conjugate[j], h, 0, 0, myE, -b, B, -Conjugate[e], 0, Conjugate[h]};

part2[[23, All]]={
Conjugate[m], A, A, -A, -A, -m, -m, -k, -k, -Conjugate[m], -Conjugate[m], m, m, k, k, k, -1, Conjugate[m], -Conjugate[a], -m, -Conjugate[m], -Conjugate[k], -Conjugate[d], -a, A, Conjugate[k], -A, m, -d, -k};

part2[[24, All]]={
-Conjugate[j], Conjugate[b], b, Conjugate[e], e, e, -Conjugate[h], -Conjugate[j], b, Conjugate[e], -h, -j, -B, -Conjugate[h], -myE, Conjugate[e], 0, Conjugate[b], 0, -h, -Conjugate[h], -Conjugate[j], 0, 0, -B, e, -myE, b, 0, -j};

part2[[25, All]]={
-m, -m, -Conjugate[m], -Conjugate[k], -k, -Conjugate[k], -k, -Conjugate[m], -m, -k, -Conjugate[k], -Conjugate[m], -m, -Conjugate[k], -k, A, 0, -A, 0, A, A, -A, 0, 0, -A, A, A, -A, 0, -A};

part2[[26, All]]={
k, -A, -A, A, A, -Conjugate[k], -Conjugate[k], -m, -m, -k, -k, Conjugate[k], Conjugate[k], m, m, m, -1, k, -d, -Conjugate[k], -k, -Conjugate[m], -Conjugate[a], -Conjugate[d], -A, Conjugate[m], A, Conjugate[k], -a, -m};

part2[[27, All]]={
myE, -Conjugate[e], -e, -b, -Conjugate[b], Conjugate[j], -b, -Conjugate[e], h, j, -Conjugate[b], myE, Conjugate[h], B, Conjugate[j], -Conjugate[b], 0, -e, 0, j, Conjugate[j], h, 0, 0, myE, -b, B, -Conjugate[e], 0, Conjugate[h]};

part2[[28, All]]={
myE, -e, -Conjugate[e], -Conjugate[b], -b, j, -Conjugate[b], -e, Conjugate[h], Conjugate[j], -b, myE, h, B, j, -b, 0, -Conjugate[e], 0, Conjugate[j], j, Conjugate[h], 0, 0, myE, -Conjugate[b], B, -e, 0, h};

part2[[29, All]]={
Conjugate[k], -A, -A, A, A, -k, -k, -Conjugate[m], -Conjugate[m], -Conjugate[k], -Conjugate[k], k, k, Conjugate[m], Conjugate[m], Conjugate[m], -1, Conjugate[k], -Conjugate[d], -k, -Conjugate[k], -m, -a, -d, -A, m, A, k, -Conjugate[a], -Conjugate[m]};

part2[[30, All]]={
-Conjugate[m], -Conjugate[m], -m, -k, -Conjugate[k], -k, -Conjugate[k], -m, -Conjugate[m], -Conjugate[k], -k, -m, -Conjugate[m], -k, -Conjugate[k], A, 0, -A, 0, A, A, -A, 0, 0, -A, A, A, -A, 0, -A};

part2[[31, All]]={
-j, b, Conjugate[b], e, Conjugate[e], Conjugate[e], -h, -j, Conjugate[b], e, -Conjugate[h], -Conjugate[j], -B, -h, -myE, e, 0, b, 0, -Conjugate[h], -h, -j, 0, 0, -B, Conjugate[e], -myE, Conjugate[b], 0, -Conjugate[j]};

part2[[32, All]]={
m, A, A, -A, -A, -Conjugate[m], -Conjugate[m], -Conjugate[k], -Conjugate[k], -m, -m, Conjugate[m], Conjugate[m], Conjugate[k], Conjugate[k], Conjugate[k], -1, m, -a, -Conjugate[m], -m, -k, -d, -Conjugate[a], A, k, -A, Conjugate[m], -Conjugate[d], -Conjugate[k]};

part2[[33, All]]={
Conjugate[h], -Conjugate[e], -e, -b, -Conjugate[b], -Conjugate[b], j, Conjugate[h], -e, -b, Conjugate[j], h, myE, j, B, -b, 0, -Conjugate[e], 0, Conjugate[j], j, Conjugate[h], 0, 0, myE, -Conjugate[b], B, -e, 0, h};

part2[[34, All]]={
-B, b, Conjugate[b], e, Conjugate[e], -Conjugate[h], e, b, -Conjugate[j], -h, Conjugate[e], -B, -j, -myE, -Conjugate[h], Conjugate[e], 0, Conjugate[b], 0, -h, -Conjugate[h], -Conjugate[j], 0, 0, -B, e, -myE, b, 0, -j};

part2[[35, All]]={
-Conjugate[k], -Conjugate[k], -k, -Conjugate[m], -m, -Conjugate[m], -m, -k, -Conjugate[k], -m, -Conjugate[m], -k, -Conjugate[k], -Conjugate[m], -m, -A, 0, A, 0, -A, -A, A, 0, 0, A, -A, -A, A, 0, A};

part2[[36, All]]={
w, w, Conjugate[w], y, Conjugate[y], y, Conjugate[y], Conjugate[w], w, Conjugate[y], y, Conjugate[w], w, y, Conjugate[y], G, 0, H, 0, G, G, H, 0, 0, H, G, G, H, 0, H};

part2[[37, All]]={
M, Conjugate[s], s, q, Conjugate[q], Conjugate[r], q, Conjugate[s], Conjugate[t], r, Conjugate[q], M, t, myN, Conjugate[r], Conjugate[q], 0, s, 0, r, Conjugate[r], Conjugate[t], 0, 0, M, q, myN, Conjugate[s], 0, t};

part2[[38, All]]={
Conjugate[u], Conjugate[v], v, Conjugate[u], u, Conjugate[v], L, L, u, v, L, u, v, v, Conjugate[u], Conjugate[v], 0, u, 0, Conjugate[u], u, v, 0, 0, L, v, L, Conjugate[u], 0, Conjugate[v]};

part2[[39, All]]={
Conjugate[r], Conjugate[q], q, Conjugate[s], s, s, t, Conjugate[r], q, Conjugate[s], Conjugate[t], r, myN, t, M, Conjugate[s], 0, Conjugate[q], 0, Conjugate[t], t, Conjugate[r], 0, 0, myN, s, M, q, 0, r};

part2[[40, All]]={
Conjugate[y], G, G, H, H, p, p, n, n, Conjugate[p], Conjugate[p], y, y, w, w, w, 1, Conjugate[y], Conjugate[a], p, Conjugate[p], Conjugate[n], Conjugate[d], a, G, Conjugate[w], H, y, d, n};

part2[[41, All]]={
t, Conjugate[s], s, q, Conjugate[q], Conjugate[q], r, t, s, q, Conjugate[r], Conjugate[t], M, r, myN, q, 0, Conjugate[s], 0, Conjugate[r], r, t, 0, 0, M, Conjugate[q], myN, s, 0, Conjugate[t]};

part2[[42, All]]={
y, y, Conjugate[y], Conjugate[w], w, Conjugate[w], w, Conjugate[y], y, w, Conjugate[w], Conjugate[y], y, Conjugate[w], w, H, 0, G, 0, H, H, G, 0, 0, G, H, H, G, 0, G};

part2[[43, All]]={
-v, -Conjugate[u], -u, -v, -Conjugate[v], -Conjugate[u], -L, -L, -Conjugate[v], -u, -L, -Conjugate[v], -u, -u, -v, -Conjugate[u], 0, -Conjugate[v], 0, -v, -Conjugate[v], -u, 0, 0, -L, -u, -L, -v, 0, -Conjugate[u]};

part2[[44, All]]={
w, H, H, G, G, Conjugate[n], Conjugate[n], p, p, n, n, Conjugate[w], Conjugate[w], y, y, y, 1, w, d, Conjugate[n], n, Conjugate[p], Conjugate[a], Conjugate[d], H, Conjugate[y], G, Conjugate[w], a, p};

part2[[45, All]]={
myN, q, Conjugate[q], s, Conjugate[s], t, s, q, Conjugate[r], Conjugate[t], Conjugate[s], myN, r, M, t, Conjugate[s], 0, Conjugate[q], 0, Conjugate[t], t, Conjugate[r], 0, 0, myN, s, M, q, 0, r};

part2[[46, All]]={
v, Conjugate[v], v, Conjugate[u], u, L, v, Conjugate[u], L, L, Conjugate[v], Conjugate[v], Conjugate[u], u, Conjugate[v], v, 0, Conjugate[u], 0, u, Conjugate[u], Conjugate[v], 0, 0, L, Conjugate[v], L, u, 0, v};

part2[[47, All]]={
-u, -Conjugate[u], -u, -v, -Conjugate[v], -L, -u, -v, -L, -L, -Conjugate[u], -Conjugate[u], -v, -Conjugate[v], -Conjugate[u], -u, 0, -v, 0, -Conjugate[v], -v, -Conjugate[u], 0, 0, -L, -Conjugate[u], -L, -Conjugate[v], 0, -u};

part2[[48, All]]={
F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, 1, F, 1, F, F, F, 1, 1, F, F, F, F, 1, F};

part2[[49, All]]={
-Conjugate[u], -u, -Conjugate[u], -Conjugate[v], -v, -L, -Conjugate[u], -Conjugate[v], -L, -L, -u, -u, -Conjugate[v], -v, -u, -Conjugate[u], 0, -Conjugate[v], 0, -v, -Conjugate[v], -u, 0, 0, -L, -u, -L, -v, 0, -Conjugate[u]};

part2[[50, All]]={
Conjugate[v], v, Conjugate[v], u, Conjugate[u], L, Conjugate[v], u, L, L, v, v, u, Conjugate[u], v, Conjugate[v], 0, u, 0, Conjugate[u], u, v, 0, 0, L, v, L, Conjugate[u], 0, Conjugate[v]};

part2[[51, All]]={
myN, Conjugate[q], q, Conjugate[s], s, Conjugate[t], Conjugate[s], Conjugate[q], r, t, s, myN, Conjugate[r], M, Conjugate[t], s, 0, q, 0, t, Conjugate[t], r, 0, 0, myN, Conjugate[s], M, Conjugate[q], 0, Conjugate[r]};

part2[[52, All]]={
Conjugate[w], H, H, G, G, n, n, Conjugate[p], Conjugate[p], Conjugate[n], Conjugate[n], w, w, Conjugate[y], Conjugate[y], Conjugate[y], 1, Conjugate[w], Conjugate[d], n, Conjugate[n], p, a, d, H, y, G, w, Conjugate[a], Conjugate[p]};

part2[[53, All]]={
-Conjugate[v], -u, -Conjugate[u], -Conjugate[v], -v, -u, -L, -L, -v, -Conjugate[u], -L, -v, -Conjugate[u], -Conjugate[u], -Conjugate[v], -u, 0, -v, 0, -Conjugate[v], -v, -Conjugate[u], 0, 0, -L, -Conjugate[u], -L, -Conjugate[v], 0, -u};

part2[[54, All]]={
Conjugate[y], Conjugate[y], y, w, Conjugate[w], w, Conjugate[w], y, Conjugate[y], Conjugate[w], w, y, Conjugate[y], w, Conjugate[w], H, 0, G, 0, H, H, G, 0, 0, G, H, H, G, 0, G};

part2[[55, All]]={
Conjugate[t], s, Conjugate[s], Conjugate[q], q, q, Conjugate[r], Conjugate[t], Conjugate[s], Conjugate[q], r, t, M, Conjugate[r], myN, Conjugate[q], 0, s, 0, r, Conjugate[r], Conjugate[t], 0, 0, M, q, myN, Conjugate[s], 0, t};

part2[[56, All]]={
y, G, G, H, H, Conjugate[p], Conjugate[p], Conjugate[n], Conjugate[n], p, p, Conjugate[y], Conjugate[y], Conjugate[w], Conjugate[w], Conjugate[w], 1, y, a, Conjugate[p], p, n, d, Conjugate[a], G, w, H, Conjugate[y], Conjugate[d], Conjugate[n]};

part2[[57, All]]={
r, q, Conjugate[q], s, Conjugate[s], Conjugate[s], Conjugate[t], r, Conjugate[q], s, t, Conjugate[r], myN, Conjugate[t], M, s, 0, q, 0, t, Conjugate[t], r, 0, 0, myN, Conjugate[s], M, Conjugate[q], 0, Conjugate[r]};

part2[[58, All]]={
u, v, Conjugate[v], u, Conjugate[u], v, L, L, Conjugate[u], Conjugate[v], L, Conjugate[u], Conjugate[v], Conjugate[v], u, v, 0, Conjugate[u], 0, u, Conjugate[u], Conjugate[v], 0, 0, L, Conjugate[v], L, u, 0, v};

part2[[59, All]]={
M, s, Conjugate[s], Conjugate[q], q, r, Conjugate[q], s, t, Conjugate[r], q, M, Conjugate[t], myN, r, q, 0, Conjugate[s], 0, Conjugate[r], r, t, 0, 0, M, Conjugate[q], myN, s, 0, Conjugate[t]};

part2[[60, All]]={
Conjugate[w], Conjugate[w], w, Conjugate[y], y, Conjugate[y], y, w, Conjugate[w], y, Conjugate[y], w, Conjugate[w], Conjugate[y], y, G, 0, H, 0, G, G, H, 0, 0, H, G, G, H, 0, H};
parts1and2=SparseArray[{},{nrows,2*ncols}];
parts1and2[[1;;nrows,1;;ncols]]=part1;
parts1and2[[1;;nrows,ncols+1;;2*ncols]]=part2;


LiuPingMat2=parts1and2[[2;;10, 1;;2*ncols]];
LiuPingMat3=parts1and2[[11;;19, 1;;2*ncols]];(*//MatrixForm*)
LiuPingMat4=parts1and2[[20;;35, 1;;2*ncols]];(*//MatrixForm*)
LiuPingMat5=parts1and2[[36;;60, 1;;2*ncols]];(*//MatrixForm*)
LiuPing2=SparseArray[{},{60,3,3}];
LiuPing3=SparseArray[{},{60,3,3}];
LiuPing4=SparseArray[{},{60,4,4}];
LiuPing5=SparseArray[{},{60,5,5}];
For[k=1,k<=60,k++, LiuPing2[[k,All,All]]=Transpose[ArrayReshape[LiuPingMat2[[All,k]], {3, 3}]]; LiuPing3[[k,All,All]]=Transpose[ArrayReshape[LiuPingMat3[[All,k]], {3, 3}]]; LiuPing4[[k,All,All]]=Transpose[ArrayReshape[LiuPingMat4[[All,k]], {4, 4}]]; LiuPing5[[k,All,All]]=Transpose[ArrayReshape[LiuPingMat5[[All,k]], {5, 5}]]; ]
MatrixForm[LiuPing2=Normal[LiuPing2]];
MatrixForm[LiuPing3=Normal[LiuPing3]];
MatrixForm[LiuPing4=Normal[LiuPing4]];
MatrixForm[LiuPing5=Normal[LiuPing5]];
(*LiuPinglist2=ArrayReshape[LiuPingMat2, {3, 3, 60}];
LiuPinglist3=ArrayReshape[LiuPingMat3, {3, 3, 60}];
LiuPinglist4=ArrayReshape[LiuPingMat4, {4, 4, 60}];
LiuPinglist5=ArrayReshape[LiuPingMat5, {5, 5, 60}];*)



\[Epsilon]=10^(-14);
\[Phi]=0;
\[Gamma]=Cos[\[Phi]]+I*Sin[\[Phi]];
SimilarTransMatrix[\[Gamma]_,n_]:=(J=SparseArray[{},{n,n}]; For[i=1,i<=n,i++,J[[i,n-i+1]]=IdentityMatrix[n][[i,i]]]; aGreek=SparseArray[{},{2n,2n}]; aGreek[[1;;n,1;;n]]=Re[\[Gamma]]*J; aGreek[[n+1;;2n,n+1;;2n]]=-Re[\[Gamma]]*J; aGreek[[1;;n,n+1;;2n]]=Im[\[Gamma]]*J; aGreek[[n+1;;2n,1;;n]]=Im[\[Gamma]]*J; Sall=Eigenvectors[aGreek]; St=SparseArray[{},{n,n}]; For[k=1,k<=n,k++,St[[k,1;;n]]=Sall[[k,1;;n]]+I*Sall[[k,n+1;;2n]]]; S=Normal[Transpose[St]])
(*SimilarTransMatrix[\[Gamma]_,n_]:=(J=IdentityMatrix[n]; aGreek=SparseArray[{},{2n,2n}]; aGreek[[1;;n,1;;n]]=Re[\[Gamma]]*J; aGreek[[n+1;;2n,n+1;;2n]]=-Re[\[Gamma]]*J; aGreek[[1;;n,n+1;;2n]]=Im[\[Gamma]]*J; aGreek[[n+1;;2n,1;;n]]=Im[\[Gamma]]*J; Sall=Eigenvectors[aGreek]; St=SparseArray[{},{n,n}]; For[k=1,k<=n,k++,St[[k,1;;n]]=Sall[[k,1;;n]]+I*Sall[[k,n+1;;2n]]]; S=Normal[Transpose[St]])*)
MatrixForm[S2=SimilarTransMatrix[\[Gamma],3]/Sqrt[2]];
MatrixForm[S2=Orthogonalize[SimilarTransMatrix[\[Gamma],3]]];(*2D Similiarty Matrix*)
MatrixForm[S3=SimilarTransMatrix[\[Gamma],3]/Sqrt[2]];
MatrixForm[S3=Orthogonalize[SimilarTransMatrix[\[Gamma],3]]];(*3D Similiarty Matrix*)
MatrixForm[S4=SimilarTransMatrix[\[Gamma],4]/Sqrt[2]];
MatrixForm[S4=Orthogonalize[SimilarTransMatrix[\[Gamma],4]]];(*4D Similiarty Matrix*)
MatrixForm[S5=SimilarTransMatrix[\[Gamma],5]/Sqrt[2]];
MatrixForm[S5=Orthogonalize[SimilarTransMatrix[\[Gamma],5]]];(*5D Similiarty Matrix*)


LiuPing2Real=SparseArray[{},{60,3,3}]; 
For[k=1,k<=60,k++,LiuPing2Real[[k,1;;3,1;;3]]=Inverse[S2] . LiuPing2[[k,1;;3,1;;3]] . S2]
MatrixForm[LiuPing2Real=Chop[Normal[LiuPing2Real]*1.0],\[Epsilon]];
LiuPing3Real=SparseArray[{},{60,3,3}];
For[k=1,k<=60,k++,LiuPing3Real[[k,1;;3,1;;3]]=Inverse[S3] . LiuPing3[[k,1;;3,1;;3]] . S3]
MatrixForm[LiuPing3Real=Chop[Normal[LiuPing3Real]*1.0],\[Epsilon]];
LiuPing4Real=SparseArray[{},{60,4,4}];
For[k=1,k<=60,k++,LiuPing4Real[[k,1;;4,1;;4]]=Inverse[S4] . LiuPing4[[k,1;;4,1;;4]] . S4]
MatrixForm[LiuPing4Real=Chop[Normal[LiuPing4Real]*1.0],\[Epsilon]];
LiuPing5Real=SparseArray[{},{60,5,5}];
For[k=1,k<=60,k++,LiuPing5Real[[k,1;;5,1;;5]]=Inverse[S5] . LiuPing5[[k,1;;5,1;;5]] . S5]
MatrixForm[LiuPing5Real=Chop[Normal[LiuPing5Real]*1.0],\[Epsilon]];
\[CapitalGamma]0[p_,g_]:=(Piecewise[{{{{1}},p==1},{Normal[LiuPing2[[g,All,All]]],p==2},{Normal[LiuPing3[[g,All,All]]],p==3},{Normal[LiuPing4[[g,All,All]]],p==4},{Normal[LiuPing5[[g,All,All]]],p==5}}])
\[CapitalGamma]r[p_,g_]:=(Piecewise[{{{{1}},p==1},{Normal[LiuPing2Real[[g,All,All]]],p==2},{Normal[LiuPing3Real[[g,All,All]]],p==3},{Normal[LiuPing4Real[[g,All,All]]],p==4},{Normal[LiuPing5Real[[g,All,All]]],p==5}}])


(*Test Homomorphism*)
(*
IsThere[M_,g1_,g2_]:=(\[Epsilon]=10^(-14); R2=Chop[M[[g1,All,All]].M[[g2,All,All]],\[Epsilon]]; flg=0; For[k=1,k<=60,k++, test=M[[k,All,All]];If[test==R2,flg=flg+1;n=k,flg=flg]]; flg n)
IsThereInv[M_,g1_]:=(\[Epsilon]=10^(-14); R2=Chop[Inverse[M[[g1,All,All]]],\[Epsilon]]; n=0; For[k=1,k<=60,k++, test=M[[k,All,All]];If[test==R2,n=k;Break[],flg=flg]]; Print[{g1, MatrixForm[Inverse[R2]],n,MatrixForm[test]}])
g1=59; g2=48;
IsThere[LiuPing3Real,g1,g2]
For[g1=45,g1<=60, g1=g1+1, IsThereInv[LiuPing3Real,g1]]
*)


(*Written by Nan Xu Oct 12, 2015*)
(*Parameters Set-Up*)
(* <<KtoI_irredi.m *)
Ryb=SparseArray[{},{60,3,3}];
S={{Cos[2*Pi/5],-Sin[2*Pi/5],0},{Sin[2*Pi/5],Cos[2*Pi/5],0},{0,0,1}};
U={{1/Sqrt[5],0,2/Sqrt[5]},{0,1,0},{-2/Sqrt[5],0,1/Sqrt[5]}};
P={{-1,0,0},{0,1,0},{0,0,-1}};
T=U . S . Inverse[U];
Ryb[[1,All,All]]={{1,0,0},{0,1,0},{0,0,1}};
Ryb[[2, All,All]]=S;
Ryb[[3,All,All]]=S . Ryb[[2, All,All]];
Ryb[[4,All,All]]=S . Ryb[[3, All,All]];
Ryb[[5,All,All]]=S . Ryb[[4, All,All]];
Ryb[[6,All,All]]=S . T;
Ryb[[7,All,All]]=T . Ryb[[6,All,All]];
Ryb[[8,All,All]]=T . Ryb[[7,All,All]];
Ryb[[9,All,All]]=Inverse[T] . Ryb[[6,All,All]];
Ryb[[10,All,All]]=Inverse[T] . Ryb[[9,All,All]];
Ryb[[11,All,All]]=S . Ryb[[6,All,All]];
Ryb[[12,All,All]]=S . Ryb[[7,All,All]];
Ryb[[13,All,All]]=S . Ryb[[8,All,All]];
Ryb[[14,All,All]]=S . Ryb[[9,All,All]];
Ryb[[15,All,All]]=S . Ryb[[10,All,All]];
Ryb[[16,All,All]]=S . Ryb[[11,All,All]];
Ryb[[17,All,All]]=S . Ryb[[12,All,All]];
Ryb[[18,All,All]]=S . Ryb[[13,All,All]];
Ryb[[19,All,All]]=S . Ryb[[14,All,All]];
Ryb[[20,All,All]]=S . Ryb[[15,All,All]];
Ryb[[21,All,All]]=Inverse[S] . Ryb[[6,All,All]];
Ryb[[22,All,All]]=Inverse[S] . Ryb[[7,All,All]];
Ryb[[23,All,All]]=Inverse[S] . Ryb[[8,All,All]];
Ryb[[24,All,All]]=Inverse[S] . Ryb[[9,All,All]];
Ryb[[25,All,All]]=Inverse[S] . Ryb[[10,All,All]];
Ryb[[26,All,All]]=Inverse[S] . Ryb[[21,All,All]];
Ryb[[27,All,All]]=Inverse[S] . Ryb[[22,All,All]];
Ryb[[28,All,All]]=Inverse[S] . Ryb[[23,All,All]];
Ryb[[29,All,All]]=Inverse[S] . Ryb[[24,All,All]];
Ryb[[30,All,All]]=Inverse[S] . Ryb[[25,All,All]];
For[i=1,i<=30,i++, Ryb[[i+30,All,All]]=P . Ryb[[i,All,All]];];


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


(*See the difference between Rpeter and Ryb *)
(*For[i=1,i<=60,i++, Print[{i,MatrixForm[1.0*Normal[Ryb[[i,All,All]]]],MatrixForm[Rpeter[[i,All,All]]]}]];*)


Rzdpermute=Rpeter[[permutation4p2N,All,All]];
Rybpermute=Ryb[[permutation4p2,All,All]];
(*For[i=1,i<=60,i++, Print[{i,MatrixForm[Chop[1.0*Normal[Rzdpermute[[i,All,All]]]-1.0*Rybpermute[[i,All,All]]]]}]];*)
(*For[i=1,i<=60,i++, Print[{i,MatrixForm[Chop[1.0*Normal[Rzdpermute[[i,All,All]]]-1.0*Normal[Inverse[Rzdpermute[[i,All,All]]]]]]}]];*)


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
MatrixForm[{\[Alpha]zd,\[Beta]zd,\[Gamma]zd}];
(*
\[Alpha]yb=SparseArray[{},{60,1}];
\[Beta]yb=SparseArray[{},{60,1}];
\[Gamma]yb=SparseArray[{},{60,1}];
For[g=1, g<=Ng, g++,M=Rybpermute[[g,All,All]]; a=Chop[EulerAngles[M][[1]]]; b=Chop[EulerAngles[M][[2]]]; c=Chop[EulerAngles[M][[3]]]; \[Alpha]yb[[g]]=a; \[Beta]yb[[g]]=b;  \[Gamma]yb[[g]]=c; Rzdtest[[g,All,All]]=Chop[Rrose[a,b,c]]]
MatrixForm[Chop[{\[Alpha]yb,\[Beta]yb,\[Gamma]yb}*1.0]]
*)


(*Test: Rzdpermute instead of Inverse[Rzdpermute] share the same multiplication table with LiuPing2 and Liuping3*)
(*
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
(*
Nx=6; xvec=RandomReal[{0,1},{Nx,3}]; x=xvec[[2,All]];
For[l=1,l<=4, l++, Print[l]; 
For[g=1, g<=Ng, g++,
Ry=Y[Inverse[Rzdpermute[[g,All,All]]] . x,l,Table[m,{m, -l,l}]];
y=Y[x,l,Table[m,{m, -l,l}]];
MatrixForm[MWigD=Table[WignerD[{l, mp, m}, -\[Alpha]zd[[g]], -\[Beta]zd[[g]], -\[Gamma]zd[[g]]],{m,-l,l},{mp,-l,l}]];
(*MatrixForm[MWigD=Table[WignerD[{l, mp, m}, \[Alpha]zd[[g]], \[Beta]zd[[g]], \[Gamma]zd[[g]]],{m,-l,l},{mp,-l,l}]];*)
(*MatrixForm[MWigD=Table[WignerD[{l, mp, m}, -\[Gamma]zd[[g]], -\[Beta]zd[[g]], -\[Alpha]zd[[g]]],{m,-l,l},{mp,-l,l}]];*)
If[Chop[MWigD . y-Ry,\[Epsilon]]!=ConstantArray[0, 2l+1], Print[{l, MWigD . y, Ry}]]]]
*)


(*Real basis functions*)
Clear[d]
Off[General::munfl]; Off[General::stop];
\[Epsilon]=10^(-14);
Y[x_,l_,m_]:=( r=Sqrt[x[[1]]^2+x[[2]]^2+x[[3]]^2];  \[Theta]=ArcCos[x[[3]]/r];\[Phi]=ArcTan[x[[1]],x[[2]]] ;SphericalHarmonicY[l,m,\[Theta],\[Phi]])
(*D coefficiets*)
powerFunc=Unevaluated[#1^#2]/.HoldPattern[0^0]:>1&;
d[l_,mp_,m_,\[Beta]_]:=Sqrt[(l+m)!*(l-m)!*(l+mp)!*(l-mp)!]*Sum[(-1)^k/((l-mp-k)!*(l+m-k)!*(k+mp-m)!*k!)*powerFunc[(Cos[\[Beta]/2]),(2l+m-mp-2k)]*powerFunc[(-Sin[\[Beta]/2]),(mp-m+2k)],{k,Max[0,(m-mp)],Min[(l-mp),(l+m)]}]
Dwig[l_,mp_,m_,\[Alpha]_,\[Beta]_,\[Gamma]_]:=Exp[-I*mp*\[Alpha]]*d[l,mp,m,\[Beta]]*Exp[-I*m*\[Gamma]]
(*Group Projection Operator*)
(*\[CapitalDelta][l_,mp_,m_,p_,i_,j_]:=dim[[p]]/Ng*Sum[ (Conjugate[\[CapitalGamma]r[p,g][[i,j]]])*WignerD[{l, mp, m}, -\[Alpha]zd[[g]], -\[Beta]zd[[g]],-\[Gamma]zd[[g]]],{g,Ng}]*)
(*\[CapitalDelta][l_,mp_,m_,p_,i_,j_]:=dim[[p]]/Ng*Sum[ (Conjugate[\[CapitalGamma][p,g][[i,j]]])*WignerD[{l, mp, m}, \[Gamma]zd[[g]], \[Beta]zd[[g]], \[Alpha]zd[[g]]],{g,Ng}] *)
\[CapitalDelta][l_,mp_,m_,p_,i_,j_]:=dim[[p]]/Ng*Sum[ (\[CapitalGamma]r[p,g][[i,j]])*Dwig[l,mp,m,\[Alpha]zd[[g]],\[Beta]zd[[g]],\[Gamma]zd[[g]]],{g,Ng}] (*Tested finally-Oct 12, 2014 !!!*)
ProjOperator[l_,m_,p_,i_,j_,x_]:=Chop[Table[\[CapitalDelta][l,mp,m,p,i,j],{mp,-l,l}] . Y[x,l,Table[k,{k,-l,l}]],\[Epsilon]]
ProjOperatorTest[l_,m_,p_,i_,j_,x_]:=dim[[p]]/Ng*Chop[Sum[Conjugate[\[CapitalGamma]r[p,g][[i,j]]]*Y[Inverse[Rzdpermute[[g,All,All]]] . x,l,m],{g,Ng}],\[Epsilon]]
LengthProjOperator[l_,m_,p_,n_]:=Chop[Sum[Abs[\[CapitalDelta][l,mp,m,p,n,n]]^2,{mp, -l,l}],\[Epsilon]]
BasisFunctionMatrix[l_,m_,p_]:=(\[Epsilon]=10^(-14);Matrix=SparseArray[{},{dim[[p]],2l+1}]; n=0;
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
BasisRealFunctionCoeffMatrixI[l_,m_,p_]:=(Mat=BasisFunctionMatrix[l,m,p];
		afunc[mp_]:=Mat[[All,mp+l+1]]; lh=0; rh=0; 
		For[mp=-l, mp<=-1, mp++, lh=lh+afunc[mp]; rh=rh+(-1)^Abs[mp]*Conjugate[afunc[-mp]];];
		angle=0;
		For[i=1,i<=Length[rh],i++,If[lh[[i]]==0,Continue[],Break[]];];
		If[i<= Length[rh], angle=Log[rh[[i]]/lh[[i]]]/2];
	MatNew=Exp[angle]*Mat)

(*Below were added on Mar 02, 2021*)
BasisRealFunctionOthoCoeffMatrixI[l_,p_]:=(Mat1={};
	For[m=-l, m<=l, m++,
		Mattest=BasisRealFunctionCoeffMatrixI[l, m, p];
		AppendTo[Mat1,ArrayReshape[Mattest,{dim[[p]]*(2l+1)}]];
	];	
	Mat1oth=Select[Chop[Orthogonalize[Chop[Mat1]]],#!=ConstantArray[0,dim[[p]]*(2l+1)]&];
	Mat1oth=ArrayReshape[Mat1oth,{dim[[p]],2l+1}]
	)


(*Test 1: generate real basis functions for the random angle*)
(*
lmin=0; lmax=5; 
For[l=lmin, l<=lmax, l++, 
  For[p=1, p<=5, p++,
	\[Theta]1=RandomReal[{0,Pi}];\[Phi]1=RandomReal[{0,2*Pi}]; 
	tcur=AbsoluteTiming[test1=RealBasisFunctionsI[l,p,\[Theta]1,\[Phi]1]];
	If[Dimensions[test1]!= {0},
		outLine1=StringForm["(l,p)=(``,``): (``,``)-``s>>``\n",l,p,TextString[\[Theta]1], TextString[\[Phi]1],TextString[tcur[[1]]],TextString[test1]];
		WriteString["stdout", outLine1];
	]
  ]
]
*)


  End[]
  EndPackage[]
