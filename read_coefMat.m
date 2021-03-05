function [tildeb, space_check, ls, Ips]=read_coefMat(b_fn, lmax, nIp)
%Input:
% b_fn: the file name of tildeb data file. The file should have the format
% below:
%     l1    //First l value
%     Ip1   //1st Ip value which has nonzero basis function coefficients under l1
%     {{}}   //basis coefficients matrix under l1, Ip1
%     l1    //First l value
%     Ip2   //2nd Ip value which has nonzero basis function coefficients under l1
%     {{}}  //basis coefficients matrix under l1, Ip2
%     ...
%     l2    //Second l value
%     Ip1'  //1st Ip value which has nonzero basis function coefficients under l2
%     {{},{},{},...}    //basis coefficients matrix under l2, Ip1'
%     
% lmax: max l

%Output: 
% tildeb: a cell array of size (lmax+1) X Ip, where each cell entry
% includes a matrix of basis coefficient matrix. Row of the 
% matrix is the vector of basis coefficient for m=-l,...,l.
%
% space_check: a (lmax+1) X 1 vector to test if the total number of basis
% function under certain l span the whole space: (0, yes; 1 no.)
%
% ls: a vector of list of all ls 
% Ips: a vector of list of all Ips
%
% Test:
% b_fn='BasisCoeff_T.txt'; lmax=30; nIp=4;
% Nan Xu
% 10/16/2014
% modified by Nan Xu on 03/02/2021

tildeb=cell(lmax+1, nIp); %lmax+1-->l, nIp-->p, nmax-->given l, each Ip has n basis, 2l
space_check=zeros(lmax+1,1);
[fid,fopenmsg]=fopen(b_fn,'r');
if fid == -1
  error(['read_coefficients: fid ' num2str(fid) ' fopenmsg ' fopenmsg ' opening ' b_fn]);
end

% TotalLength=2*sum(0:lmax)+lmax+1;
% ll=zeros(TotalLength,1); Ipl=zeros(TotalLength,1); 
tline = fgetl(fid); rd_l1=str2num(tline); 
tline = fgetl(fid); rd_Ip=str2num(tline);  %read the l and Ip
if isempty(rd_Ip) || isempty(rd_Ip)
    error('read_coefficients: file opening reading error\n');
end
% if (rd_l1~=0 || rd_Ip~=1)
%     error('read_coefficients: l 0 rd_l %d Ip 1 rd_Ip %d\n',rd_l1,rd_Ip);
% end
ll_ct=0;
while ~feof(fid)   
    ll_ct=ll_ct+1;
    ls(ll_ct)=rd_l1;
    Ips(ll_ct)=rd_Ip;
    n=0;  
    
    tline = fgetl(fid); tildeb_curv=eval(tline); Nlp=size(tildeb_curv,2);
    tildeb_curm=zeros(Nlp,2*rd_l1+1);
    for j=1:Nlp
        tildeb_curm(j,:)=cell2mat(tildeb_curv{j});
    end
    
    tildeb{rd_l1+1,rd_Ip}=tildeb_curm;
    n=n+Nlp;
    if  ~feof(fid)
        tline = fgetl(fid); 
        rd_l1=str2num(tline); %tline is the new l
        tline = fgetl(fid);      
        rd_Ip=str2num(tline); %tline is the new Ip
    end
    
    if n~=2*rd_l1+1
        space_check(rd_l1+1)=1;
    end
end
fclose(fid);
