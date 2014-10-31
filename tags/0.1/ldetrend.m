function [y1,a,b]=ldetrend(y0)
% [y1,a,b]=ldetrend(y0)
%
% Remove linear trend from data by finding the 
% best fitting coeficients a,b, that describe
% y=ax+b
%
% Last modified by fjsimons-at-alum.mit.edu, 15.11.2004

y0=y0(:);
m=length(y0);
x=1:m;

% Make Jacobian
A=[x(:) repmat(1,m,1)];

% Invert it
ab=inv(A'*A)*A'*y0;

a=ab(1);
b=ab(2);

% Calculate and remove linear trend
y1=y0-a*x(:)-b;

