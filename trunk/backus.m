function [ap1,th,th0,ap2]=backus(L,N,norma)
% [ap1,th,th0,ap2]=BACKUS(L,N,norma)
%
% Backus/Robin asymptotic approximation to the Legendre polynomials of
% large degree L, which works for 0<theta<pi, for ZONAL polynomials. 
%
% INPUT:
%
% L          Degree of the Legendre polynomial
% N          Number of points in ]0,pi[
% norma      'sch' Schmidt-normalized polynomials 
%            'fnr' Fully-normalized real
%            'fnc' Fully-normalized complex
%
% OUTPUT:
%
% ap1        The approximation due to DT (B.86)
% th         The abscissa
% th0        The validity range is [th0 pi-th0]
% ap2        A better approximation due to DT (B.87)
%
% See Dahlen and Tromp, Appendix B p 855.
%
% EXAMPLE:
%
% libbrecht('demo4')
%
% See also DAHLEN, HILBXLM.
%
% By fjsimons-at-alum.mit.edu, Feb 11th, 2004

defval('norma','sch')
defval('N',1000);

% DT Eq. B. 80
th0=asin(1/(sqrt(L*(L+1))));
th=linspace(0+eps,pi-eps,N);

% Extra normalization factors, when m is always 0
switch norma
 case 'sch'
  fac=sqrt(4*pi/(2*L+1));
 case 'fnr'
  fac=sqrt(4*pi);
 case 'fnc'
  fac=1;
 otherwise
  error('Specify valid normalization')
end

% The approximation DT Eq. B. 86
arg=[(L+1/2)*th-pi/4];
ap1=fac/pi./sqrt(sin(th)).*(cos(arg));
ap2=fac/pi./sqrt(sin(th)).*(cos(arg)+1/8/(L+1/2).*cot(th).*sin(arg));

