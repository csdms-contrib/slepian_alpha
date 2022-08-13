function [x,y]=randcirc(xm,ym,r,dr,N)
% [x,y]=randcirc(xm,ym,r,dr,N)
%
% Gives the coordinates of a random blobby circle.
%
% INPUT:
%
% xm           Horizontal position of the center [default: 0]
% ym           Vertical position of the center [default: 0]
% r            Radius [default: 1]
% dr           Size of random perturbations around the radius [default: 0.1]
% N            Number of random spike points [default: 10]
%
% OUPUT:
%
% [x,y]        The coordinates
%
% EXAMPLE:
%
% [x,y]=randcirc(1,3,2,0.3,128);
% plot(x,y); axis image; grid on
%
% SEE ALSO: 
%
% BLOB
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)
% Last modified by fjsimons-at-alum.mit.edu, 08/04/2022

defval('xm',0)
defval('ym',0)
defval('r',1)
defval('dr',0.1)
defval('N',10)

nr=100;

r=repmat(r,N,1);
r=r+dr*(rand(N,1)-0.5)*2;
t=linspace(0,2*pi*N/(N+1),N)';

r=[r;r;r];
t=[t-2*pi;t;t+2*pi];

tt=linspace(0,2*pi,nr);

if verLessThan('matlab','6')
  rr=interp1(t,r,tt,'cubic');
else
  % rr=interp1(t,r,tt,'v5cubic');
  rr=interp1(t,r,tt,'pchip');
end

% Ouput
x=xm+rr.*cos(tt);
y=ym+rr.*sin(tt);

