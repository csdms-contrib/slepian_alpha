function phint=dphpatch(th,thR,th0,ph0)
% phint=DPHPATCH(th,thR,th0,ph0)
%
% Calculates the integration domain of a spherical patch.
% If th=th0=pi/2 then correctly gives range(phint)=2*thR.
%
% INPUT:
%
% th         Colatitude where you want it evaluated, in radians
% thR        Radius of the cap, in radians
% th0        Colatitude of the cap center, in radians
% ph0        Longitude of the cap center, in radians
%
% OUTPUT:
%
% phint      Longitudinal integration interval, in radians
%
% Last modified by fjsimons-at-alum.mit.edu, May 12th, 2004

dph=real(acos((cos(thR)-cos(th)*cos(th0))./(sin(th)*sin(th0))));
phint=[ph0-dph(:) ph0+dph(:)];

