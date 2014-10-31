function [lmcosi,C,V]=boxcap(TH,L,ms)
% [lmcosi,C,V]=boxcap(TH,L,ms)
%
% OBSOLETE.
%
% Returns the spherical harmonic coefficients of a cylindrical 
% boxcar, adequately normalized to unity on the unit sphere.
%
% INPUT:
%
% TH       Colatitudinal radius of the cap, in degrees
% L        Bandwidth
% ms       0 Coefficients normalized to unit power [default]
%          1 Coefficients normalized according to Mark Simons
%
% OUTPUT:
%
% lmcosi   Spherical harmonic coefficient array, normalized to unity
% C        Just the zonal coefficients vector, normalized to unity
% V        The "eigenvalue" - in this case the energy leakage parameter
%
% See also SSHCAP
%
% Last modified by fjsimons-at-alum.mit.edu, 05.05.2005

defval('TH',20)
defval('L',8)
defval('ms',0)

% Evaluate Legendre polynomials at one point only
x=cos(TH*pi/180);
[dems,dels,mzero,lmcosi]=addmon(L);
lmcosi(1,3)=sqrt(4*pi);

% Simons, Solomon and Hager Eq. 49
C(1)=sqrt(4*pi);
Plm=repmat(NaN,1,addmup(L));
for l=0:L+2
  Pl(l+1)=rindeks(legendre(l,x),1);
end

% Simons, Solomon and Hager Eq. 50
fax=Pl(1)-Pl(2);
for l=1:L
  C(l+1)=sqrt(4*pi/(2*l+1))*(Pl(l)-Pl(l+2))/fax;
  lmcosi(mzero(l+1),3)=C(l+1);
end

C=C(:);

if ms==0
  % Normalization to unity on the unit sphere
  C=C/sqrt(C'*C);
  
  % Check that all is right, and determine the leakage parameter
  [ngl1,ngl2,com,V]=orthocheck(C,[],TH*pi/180,0);
else
    V=NaN;
    disp('Normalization according to Mark Simons, GJI 1996')
end

lmcosi(mzero,3)=C;

