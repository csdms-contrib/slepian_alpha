function lmrlim=rsh2cpx(lmcosi)
% lmrlim=RSH2CPX(lmcosi)
%
% Transforms the coefficients of a field expanded in REAL spherical
% harmonics into an expansion of COMPLEX spherical harmonics.
% Incomplete, only used in PLM2ROT.
%
% INPUT:
%
% lmcosi        Standard matrix listing l, m, Ccos and Csin
%
% OUTPUT:
%
% lmrlim        Matrix listing coefficients for m=-l:l as
%               l abs(m) Creal Cimag
%
% SEE ALSO:
%
% UMMP, ULMMP
%
% Last modified by fjsimons-at-alum.mit.edu, 19.05.2006

[l,m,mz]=addmon(max(lmcosi(:,1)));

lmrlim=lmcosi;

% Divide everything by sqrt(2)...
lmrlim(:,3:4)=lmcosi(:,3:4)/sqrt(2);
% ... except the m=0 component
lmrlim(mz,3:4)=lmcosi(mz,3:4);




