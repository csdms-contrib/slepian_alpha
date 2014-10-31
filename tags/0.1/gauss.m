function z=gauss(x,s)
% z=GAUSS(x,s)
%
% Retuns a gaussian defined on 'x' with standard deviation 's'.
%
% Example:
%
% x=-10:0.1:10; 
% for index=1:4 ; plot(x,gauss(x,index)) ; hold on ; end ; grid
% lin2col(gca)

% Written by FJS, August 6th 1998

z=exp(-x.^2/(2*s^2))/(s*sqrt(2*pi));