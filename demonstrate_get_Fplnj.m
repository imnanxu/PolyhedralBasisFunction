%Peter C. Doerschuk, 2018-01-14
%This harmonic is shown in Figure 1(c)
p=2;
l=10;
n=1;
j=1;
thetavalues=[0.1 0.2 0.3 0.4]'*pi;
phivalues=[0.15 0.25 0.35 0.45]'*pi;

filename='IcosahedralRealBasisFunctionCoeff.txt';
lmax=l;
[tilde_b, space_check, ll, Ipl]=read_coefficients(filename, lmax);
F=get_Fplnj(p,l,n,j,thetavalues,phivalues,tilde_b);
fprintf(1,'       theta          phi        Fplnj for p %d l %d n %d j %d\n',p,l,n,j);
fprintf(1,'%12.5g %12.5g %12.5g\n',[thetavalues phivalues F']');
