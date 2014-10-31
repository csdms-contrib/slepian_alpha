function [lon,lat,colnr,rownr]=ind2cor(indices,c11,cmn,m,n)
% [lon,lat,colnr,rownr]=IND2COR(indices,c11,cmn,m,n)
%
% Given 'c11' and 'cmn', transforms the 'indices' into a matrix
% to PIXEL CENTER  longitude and latitude coordinates for indices from 
% a matrix with 'm' rows and 'n' columns
%
% c11 and cmn are PIXEL CENTERED coordinates of center of first and last element
% of the matrix, in (x-y).  Indices may be a matrix or a vector.
%
% Test on:
%
% c11=[90 -10]; cmn=[100 -30]; m=4 ; n=6;
% [lon,lat]=ind2cor(1:m*n,c11,cmn,m,n)
%
% Last modified by fjsimons-at-alum.mit.edu

lonint=(cmn(1)-c11(1))/(n-1);
latint=(cmn(2)-c11(2))/(m-1);

colnr=ceil(indices/m);
rownr=indices-(colnr-1)*m;

lon=c11(1)+(colnr-1)*lonint;
lat=c11(2)+(rownr-1)*latint;

% Reinvented this Dec. 4th 2001. 
%  lon=C11(1)+floor(reginds/M)/N*(CMN(1)-C11(1));
%  lat=C11(2)-mod(reginds,M)/M*(C11(2)-CMN(2));



  
