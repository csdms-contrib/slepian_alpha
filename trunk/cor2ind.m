function [ind,colnr,rownr]=cor2ind(lon,lat,c11,cmn,m,n)
% [ind,colnr,rownr]=cor2ind(lon,lat,c11,cmn,m,n)
%
% Transforms longitudes and latitudes given an 
% mXn array with c11 and cmn the BOUNDARY
% coordinates of the first and last elements of the matrix
% into a running index and column and row number.
%
% EXAMPLE:
%
% cor2ind('demo');
%
% See also IND2COR
%
% Written by fjsimons-at-mit.edu, October 11th 2000

if ~isstr(lon)
  lonspan=(cmn(1)-c11(1));
  latspan=(cmn(2)-c11(2));
  
  % m and n are the elements of the matrix underlying the grid;
  % the grid has n+1 and m+1 elements
  lonint=lonspan/n;
  latint=latspan/m;
  
  colnr=ceil((lon-c11(1))/lonint);
  rownr=ceil((lat-c11(2))/latint);
  
  ind=rownr+(colnr-1)*m;
elseif strmatch(lon,'demo')
  clf
  c11=[10 8]; cmn=[38 -42]; n=32; m=27;
  % Make grid with these GRID BOUNDARIES
  long=linspace(c11(1),cmn(1),n); lont=indeks(diff(long),1);
  latg=linspace(cmn(2),c11(2),m); lant=indeks(diff(latg),1);
  [lonlon,latlat]=meshgrid(long,latg);
  
  lon=c11(1)+rand*[cmn(1)-c11(1)];
  lat=cmn(2)+rand*[c11(2)-cmn(2)];
  [ind,colnr,rownr]=cor2ind(lon,lat,c11,cmn,m-1,n-1)
  mat=repmat(NaN,m-1,n-1); mat(rownr,colnr)=0;
  imagef(c11+[lont -lant]/2,cmn+[-lont lant]/2,mat);
  hold on ;  plot(lon,lat,'+'); 
  fridplot(lonlon,latlat,'Color','w'); hold off
end


