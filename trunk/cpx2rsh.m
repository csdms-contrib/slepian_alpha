function lmcosi=cpx2rsh(lmrlim)
% lmcosi=CPX2RSH(lmrlim)
%
% Transforms the coefficients of a field expanded in COMPLEX spherical
% harmonics into an expansion of REAL spherical harmonics. 
% Incomplete, only used in PLM2ROT.
%
% INPUT:
%
% lmrlim        Matrix listing coefficients for m=-l:l as
%               l abs(m) Creal Cimag
%
% OUTPUT:
%
% lmcosi        Standard matrix listing l, m, Ccos and Csin
%
% SEE ALSO:
%
% UMMP, ULMMP
%
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012

% Find out where the m=0 coefficients sit in this matrix
[l,m,mz]=addmon(max(lmrlim(:,1)));

lmcosi=lmrlim;

% Multiply everything by sqrt(2)...
lmcosi(:,3:4)=lmrlim(:,3:4)*sqrt(2);

% ... except the m=0 component
lmcosi(mz,3:4)=lmrlim(mz,3:4);
