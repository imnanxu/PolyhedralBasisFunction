function Fplnj=get_Fplnj(p,l,n,j,thetavalues,phivalues,tilde_b)
%function Fplnj=get_Fplnj(p,l,n,j,thetavalues,phivalues,tilde_b)
%Peter C. Doerschuk, 2018-01-14
%icosahedral symmetry
%p=which irrep, possible values 1,2,3,4,5
%l=spherical harmonic order, possible values 0,1,2,...
%n=which vector of length d_p, possible values 1,2,3,...
%j=which component of the vector, possible values 1,2,...,d_p
%thetavalues and phivalues are vectors of angle values
%whether a harmonic with indices p,l,n,j exists is not tested

assert(~isinteger(p));
assert(~isinteger(l));
assert(~isinteger(n));
assert(~isinteger(j));
assert(all(isreal(thetavalues)));
assert(all(isreal(phivalues)));

dp=[1 3 3 4 5]; %d_p=dimension of the p th irrep for icosahedral symmetry

assert(1<=p && p<=length(dp));
assert(0<=l);
assert(1<=n);
assert(1<=j && j<=dp(p));

costhetavalues=cos(thetavalues);
nj=(n-1)*dp(p) + j; %n and j should not occur after this statement

%compute Y_{l,m} for one value of l and all values of m (-l<=m<=+l)
%Matlab legendre function provides only has 0<=m<=+l not -l<=m<=+l
%Matlab legendre function writes P_n^m while PCD writes P_{l,m}:
%Matlab m = PCD m = order
%Matlab n = PCD l = degree
%P: row index is m+1, column index is index to costhetavalues
P=legendre(l,costhetavalues);
assert(ndims(P)==2);
assert(size(P,1)==l+1);
assert(size(P,2)==length(thetavalues));

%expimphi: row index is m+1, column index is index to phivalues
expimphi=exp(sqrt(-1).*[0:l]'*phivalues'); %[0:l] are the values of m
assert(ndims(expimphi)==2);
assert(size(expimphi,1)==l+1);
assert(size(expimphi,2)==length(phivalues));

Ftable=factorial(0:2*l).'; %table of factorial values, result is column vector
tmp=(2*l+1)/(4*pi); %scalar
normalizer=sqrt(tmp*Ftable(l-[0:l]+1)./Ftable(l+[0:l]+1)); %[0:l] are the values of m, result is column vector

%Y_{l,m}(\theta,\phi) for fixed l, varying m\in\{0,\dots,l\} , (\theta,\phi)=[thetavalues(:) phivalues(:)]
%Y: row index is m+1, column index is index to [costhetavalues(:) phivalues(:)]
Y=P.*expimphi;
Y=bsxfun(@times,Y,normalizer);
assert(ndims(Y)==2);
assert(size(Y,1)==l+1);
assert(size(Y,2)==length(phivalues));
assert(size(Y,2)==length(thetavalues));

%Can compute linear combinations of Y's for all angle pairs by doing weighted sums of the rows.

cpl_allnj_allm=tilde_b{l+1,p}; %matrix of c_{p,l,(n,j)=row,m=col}
cplnj_allm=cpl_allnj_allm(nj,:); %row vector of c_{p,l,(n,j),m=col} where m\in\{-l,\dots,+l\}
assert(size(cplnj_allm,1)==1);
assert(size(cplnj_allm,2)==2*l+1);

%evaluate the linear combination
%F_{p,l,n,j}(\theta,\phi) for fixed indices but all [thetavalues(:) phivalues(:)], result is row vector
%add the terms for m^\prime=0,...,l for all values of [thetavalues(:) phivalues(:)]
Fplnj=cplnj_allm(l+1:end)*Y; %row vector times matrix, result is row vector with length(thetavalues) columns
%add the terms for m^\prime=-l,...,-1 (note reversed order of cplnj_allm elements)
signfromm=1-2*mod(1:l,2); %answer is empty when l=0
assert(size(signfromm,1)==1);
assert(size(signfromm,2)==l);
Fplnj=Fplnj+(signfromm.*cplnj_allm(l:-1:1))*conj(Y(2:end,:)); %row vector times matrix, result is row vector with length(thetavalues) columns
assert(size(Fplnj,1)==1);
assert(size(Fplnj,2)==length(thetavalues));

Fplnj=real(Fplnj);
