function [x,y]=mollweide(lon,lat,lonref)
% [x,y]=MOLLWEIDE(lon,lat,lonref)
%
% Mollweide projection of spherical coordinates
%
% INPUT:
%
% lon,lat     longitude and latitude arrays [radians]
% lonref      the projection reference [default: pi]
%
% OUTPUT:
%
% x,y        the coordinates after Mollweide projection
%
% Last modified by fjsimons-at-alum.mit.edu , 05/26/2021

defval('lonref',pi)

% Mapping to the right range
i=(lonref-lon)>pi;
lon(i)=lon(i)+2*pi;

i=(lon-lonref)>pi;
lon(i)=lon(i)-2*pi;

% The projection
theta=lat;
dtheta=1;

% I do believe this routine originated with Hrafknell Karason
while max(abs(dtheta))>1e-8,
  dtheta=-(theta+sin(theta)-pi*sin(lat))./(1+cos(theta));
  theta=theta+dtheta;
end
theta=theta/2;

% The final output
x=2*sqrt(2)*(lon-lonref).*cos(theta)/pi;
y=sqrt(2)*sin(theta);
