function coef=lmcosi2coef(lmcosi,onorout)
% coef=lmcosi2coef(lmcosi,onorout)
%
% Transform coefficients given as lmcosi into a single vector of 
% coefficients in the format of your choosing (addmon or addmout)
%
% INPUT:
%
% lmcosi    the coefficients ordered in lmcosi format (l m cos sin)
% 
% onorout   if you would like your coef vector in addmon format, you can
%           omit this variable or set it to zero. If you want your vector
%           in addmout, set this to 1.
% 
% OUTPUT:
%
% coef      the spherical-harmonic coefficients ordered in either addmon or
%           addmout format.
%
% See also coef2lmcosi, coef2blmclm, blmclm2coef
%
% Last modified by plattner-at-alumni.ethz.ch, 12/09/2014

defval('onorout',0)

Lmax=lmcosi(end,1);

[demsz,delsz,mz,lmc,mzin,mzo,bigm,bigl,rinm,ronm,demin]=addmon(Lmax); 
cosivec=reshape(lmcosi(:,3:4),1,2*length(demsz));
coef=cosivec(mzo);

if onorout
    coef=coef(rinm);
end
coef=coef(:);
