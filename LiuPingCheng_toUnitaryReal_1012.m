(* ::Package:: *)

(*%Fa Liu, Jia-Lun Ping, Jin-Quan Chen, "Application of the eigenfunction method to the icosahedral group", J. Math. Phys 31(5):1065-1075, May 1990.*)
(*Written by Nan Xu Oct 12, 2015*)
ClearAll["Global`*"]
Ng=60;
dim={1,3,3,4,5};


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
\[Phi]=0; \[Gamma]=Cos[\[Phi]]+I*Sin[\[Phi]];
n=4;


J=SparseArray[{},{n,n}]; For[i=1,i<=n,i++,J[[i,n-i+1]]=IdentityMatrix[n][[i,i]]]; 
MatrixForm[J]
aGreek=SparseArray[{},{2n,2n}]; aGreek[[1;;n,1;;n]]=Re[\[Gamma]]*J; aGreek[[n+1;;2n,n+1;;2n]]=-Re[\[Gamma]]*J; aGreek[[1;;n,n+1;;2n]]=Im[\[Gamma]]*J; aGreek[[n+1;;2n,1;;n]]=Im[\[Gamma]]*J; 
Sall=Eigenvectors[aGreek]; 
St=SparseArray[{},{n,n}]; For[k=1,k<=n,k++,St[[k,1;;n]]=Sall[[k,1;;n]]+I*Sall[[k,n+1;;2n]]]; 
S=Normal[Transpose[St]]
MatrixForm[Orthogonalize[S]]


n=4; J=SparseArray[{},{n,n}]; For[i=1,i<=n,i++,J[[i,n-i+1]]=IdentityMatrix[n][[i,i]]]; J=Normal[J]
k=28;
MatrixForm[A=Conjugate[LiuPing4[[k,1;;n,1;;n]]]];
MatrixForm[B=J.LiuPing4[[k,1;;n,1;;n]].J];
MatrixForm[Chop[(A-B)*1.0]]


k=28;
MatrixForm[A=Conjugate[LiuPing3[[k,1;;3,1;;3]]]]
MatrixForm[B=J.LiuPing3[[k,1;;3,1;;3]].J]
MatrixForm[Chop[(A-B)*1.0]]


(*Test*)
sum2=0; For[k=1,k<=60,k++,sum2=sum2+LiuPing2[[k,1;;3,1;;3]].Transpose[LiuPing2[[k,1;;3,1;;3]]]];
sum3=0; For[k=1,k<=60,k++,sum3=sum3+LiuPing3[[k,1;;3,1;;3]].Transpose[LiuPing3[[k,1;;3,1;;3]]]];
sum4=0; For[k=1,k<=60,k++,sum4=sum4+LiuPing4[[k,1;;4,1;;4]].J.Transpose[LiuPing4[[k,1;;4,1;;4]]]];
sum5=0; For[k=1,k<=60,k++,sum5=sum5+LiuPing5[[k,1;;5,1;;5]].Transpose[LiuPing5[[k,1;;5,1;;5]]]];
Simplify[Chop[sum2*1.0]/Ng]
Simplify[Chop[sum3*1.0]/Ng]
Simplify[Chop[sum4*1.0]/Ng]
Simplify[Chop[sum5*1.0]/Ng]





\[Epsilon]=10^(-14);
\[Phi]=0;
\[Gamma]=Cos[\[Phi]]+I*Sin[\[Phi]];
SimilarTransMatrix[\[Gamma]_,n_]:=(J=SparseArray[{},{n,n}]; For[i=1,i<=n,i++,J[[i,n-i+1]]=IdentityMatrix[n][[i,i]]]; aGreek=SparseArray[{},{2n,2n}]; aGreek[[1;;n,1;;n]]=Re[\[Gamma]]*J; aGreek[[n+1;;2n,n+1;;2n]]=-Re[\[Gamma]]*J; aGreek[[1;;n,n+1;;2n]]=Im[\[Gamma]]*J; aGreek[[n+1;;2n,1;;n]]=Im[\[Gamma]]*J; Sall=Eigenvectors[aGreek]; St=SparseArray[{},{n,n}]; For[k=1,k<=n,k++,St[[k,1;;n]]=Sall[[k,1;;n]]+I*Sall[[k,n+1;;2n]]]; S=Normal[Transpose[St]])
(*SimilarTransMatrix[\[Gamma]_,n_]:=(J=IdentityMatrix[n]; aGreek=SparseArray[{},{2n,2n}]; aGreek[[1;;n,1;;n]]=Re[\[Gamma]]*J; aGreek[[n+1;;2n,n+1;;2n]]=-Re[\[Gamma]]*J; aGreek[[1;;n,n+1;;2n]]=Im[\[Gamma]]*J; aGreek[[n+1;;2n,1;;n]]=Im[\[Gamma]]*J; Sall=Eigenvectors[aGreek]; St=SparseArray[{},{n,n}]; For[k=1,k<=n,k++,St[[k,1;;n]]=Sall[[k,1;;n]]+I*Sall[[k,n+1;;2n]]]; S=Normal[Transpose[St]])*)
MatrixForm[S2=SimilarTransMatrix[\[Gamma],3]/Sqrt[2]]
MatrixForm[S2=Orthogonalize[SimilarTransMatrix[\[Gamma],3]]](*same as the last one*)
MatrixForm[S3=SimilarTransMatrix[\[Gamma],3]/Sqrt[2]]
MatrixForm[S3=Orthogonalize[SimilarTransMatrix[\[Gamma],3]]](*same as the last one*)
MatrixForm[S4=SimilarTransMatrix[\[Gamma],4]/Sqrt[2]]
MatrixForm[S4=Orthogonalize[SimilarTransMatrix[\[Gamma],4]]](*same as the last one*)
MatrixForm[S5=SimilarTransMatrix[\[Gamma],5]/Sqrt[2]]
MatrixForm[S5=Orthogonalize[SimilarTransMatrix[\[Gamma],5]]](*same as the last one*)



LiuPing2Real=SparseArray[{},{60,3,3}];
For[k=1,k<=60,k++,LiuPing2Real[[k,1;;3,1;;3]]=Inverse[S2].LiuPing2[[k,1;;3,1;;3]].S2]
MatrixForm[LiuPing2Real=Chop[Normal[LiuPing2Real]*1.0,\[Epsilon]]]
LiuPing3Real=SparseArray[{},{60,3,3}];
For[k=1,k<=60,k++,LiuPing3Real[[k,1;;3,1;;3]]=Inverse[S3].LiuPing3[[k,1;;3,1;;3]].S3]
MatrixForm[LiuPing3Real=Chop[Normal[LiuPing3Real]*1.0,\[Epsilon]]]
LiuPing4Real=SparseArray[{},{60,4,4}];
For[k=1,k<=60,k++,LiuPing4Real[[k,1;;4,1;;4]]=Inverse[S4].LiuPing4[[k,1;;4,1;;4]].S4]
MatrixForm[LiuPing4Real=Chop[Normal[LiuPing4Real]*1.0,\[Epsilon]]]
LiuPing5Real=SparseArray[{},{60,5,5}];
For[k=1,k<=60,k++,LiuPing5Real[[k,1;;5,1;;5]]=Inverse[S5].LiuPing5[[k,1;;5,1;;5]].S5]
MatrixForm[LiuPing5Real=Chop[Normal[LiuPing5Real]*1.0,\[Epsilon]]]
\[CapitalGamma]r[p_,g_]:=(Piecewise[{{{{1}},p==1},{Normal[LiuPing2Real[[g,All,All]]],p==2},{Normal[LiuPing3Real[[g,All,All]]],p==3},{Normal[LiuPing4Real[[g,All,All]]],p==4},{Normal[LiuPing5Real[[g,All,All]]],p==5}}])


\[CapitalGamma][p_,g_]:=(Piecewise[{{{{1}},p==1},{Normal[LiuPing2[[g,All,All]]],p==2},{Normal[LiuPing3[[g,All,All]]],p==3},{Normal[LiuPing4[[g,All,All]]],p==4},{Normal[LiuPing5[[g,All,All]]],p==5}}])


p=2; For[g=1, g<=60, g++, M=Det[\[CapitalGamma][p,g]]; Print[{p, M}]] 


(*p=2;q=1; j=1; k=2;
Clear[g]
For[p=1,p<=5,p++,
	For[j=1,j<= dim[[p]],j++,
		For[k=1,k<= dim[[p]],k++,
			For[q=1,q<=5,q++,
				For[s=1,s<= dim[[q]],s++,
					For[t=1,t<= dim[[q]],t++,
						Cur=Simplify[Sum[Conjugate[\[CapitalGamma][p,g][[j,k]]]*\[CapitalGamma][q,g][[s,t]],{g,Ng}]]*1.0; 
						If[Cur!=0, 
							Print[{"p", p, j, k, "q", q, s, t,  Cur}]
						]
					]
				]
			]
		]
	]
]*)


IsThere[M_,g1_,g2_]:=(\[Epsilon]=10^(-14); R2=Chop[M[[g1,All,All]].M[[g2,All,All]],\[Epsilon]]; flg=0; For[k=1,k<=60,k++, test=M[[k,All,All]];If[test==R2,flg=flg+1;n=k,flg=flg]]; flg n)
g1=59; g2=48;
IsThere[LiuPing3Real,g1,g2]


IsThereInv[M_,g1_]:=(\[Epsilon]=10^(-14); R2=Chop[Inverse[M[[g1,All,All]]],\[Epsilon]]; n=0; For[k=1,k<=60,k++, test=M[[k,All,All]];If[test==R2,n=k;Break[],flg=flg]]; Print[{g1, MatrixForm[Inverse[R2]],n,MatrixForm[test]}])
g1=59;
IsThereInv[LiuPing3Real,g1]


61-47


For[g1=45,g1<=60, g1=g1+1, IsThereInv[LiuPing3Real,g1]]






