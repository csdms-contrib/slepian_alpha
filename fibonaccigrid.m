function [lon,lat]=fibonaccigrid(N)
% [lon,lat]=FIBONACCIGRID(N)
%
% Construct a Fibonacci grid, see Álvaro González, "Measurement of
% areas on a sphere using Fibonacci and latitude–longitude lattices",
% Mathematical Geosciences (2010) 42: 49–64, 10.1007/s11004-009-9257-x
%
% INPUT:
%
% N      The number of points will be 2N+1 [default: 1000]
%
% Last modified by Cornelis Slobbe, 07/12/2012 
% Last modified by fjsimons-at-alum.mit.edu, 01/26/2023

defval('N',1000);
defval('WhichAngle','ComplGoldenAngle')

% Golden ratio
Phi=(1+sqrt(5))/2;

lat=asin((2*[-N:1:N])/(2*N+1))'*(180/pi);

if strcmp(WhichAngle,'ComplGoldenAngle')
    lon=mod([-N:1:N],Phi)'*(360/Phi);
elseif strcmp(WhichAngle,'GoldenAngle')
    lon=-mod([-N:1:N],Phi^2)'*(360/Phi^2);
end

lon(lon<0)  =lon(lon<0)  +360;
lon(lon>360)=lon(lon>360)-360;

% axesm globe
% plotm(lat,lon,'.')
% plot(lon,lat,'.')

% Prepare output
vars={lon,lat};
varargout=vars(1:nargout);
