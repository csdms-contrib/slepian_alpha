function lmcosi=coef2lmcosi(coef,onorout)
% lmcosi=coef2lmcosi(coef,onorout)
%
% Transform coefficients given as a vector "coef" into the lmcosi format.
% Requires the coef vector to be "complete" as in having (L+1)^2 entries.
%
% INPUT:
%
% coef      the spherical-harmonic coefficients ordered in either addmon or
%           addmout format.
% onorout   if coef in addmout format, set this to 1. Otherwise (addmon) 
%           leave out or set it to zero
% 
% OUTPUT:
%
% lmcosi    the coefficients ordered in lmcosi format (l m cos sin)
%
% See also lmcosi2coef, coef2blmclm, blmclm2coef
%
% Last modified by plattner-at-alumni.ethz.ch, 12/09/2014

defval('onorout',0)

coef=coef(:);
Lmax=sqrt(length(coef))-1;

if onorout
    coef=out2on(coef,Lmax);
end    

[demsz,delsz,mz,lmc,mzin]=addmon(Lmax);

lmcosi=[delsz demsz reshape(insert(coef,0,mzin),2,length(demsz))'];


