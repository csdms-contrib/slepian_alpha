function [fth,theta]=pl2th(cofs,nth,nrm)
% [fth,theta]=pl2th(cofs,nth,nrm)
%
% Inverse Legendre transform, undone by TH2PL
%
% Construct spatial functions from zonal coefficients
% (given starting with L=0). Resulting spherical harmonic
% is normalized over the unit sphere.
%
% f(\theta)=\sum_{l=0}^{L}f_lP_l(\cos(theta))
%
% INPUT:
%
% cofs      Vector/Matrix with zonal spherical harmonics coefficients
% nth       Number of colatitudes on interval [0 pi] [default: 720]
%           If nth is a vector, calculates it at specific points
% nrm       Normalization factor for Ylm [default: 4pi], may want 1
%
% OUTPUT:
%
% fth       The colatitudinal function, given by the weighted sum of ZSH
% theta     The data grid
%
% EXAMPLE:
% 
% [fth,theta]=pl2th([0 0 0 1 0 0 0]',128);
% [C1,rnk]=th2pl(fth,10,'gl');
% [C2,rnk]=th2pl(fth,30,'im');
% [C3,rnk]=th2pl(fth,60,'im');
% [C4,rnk]=th2pl(fth,128,'im');
% [C5,rnk]=th2pl(fth,250,'im');
% [C6,rnk]=th2pl(fth,30,'simpson');
%% Now watch out: for this:
% [C7,rnk]=th2pl(fth,7,'gl');
% [C8,rnk]=th2pl(fth,120,'simpson');
%
% See also TH2PL, PLM2XYZ, XYZ2PLM. Watch the factors, e.g. SDWCAPT2.
%
% Last modified by fjsimons-at-alum.mit.edu, 04/02/2009

defval('nth',720);
defval('nrm',4*pi);
if length(nth)==1 % If this is a number, in other words
  % Make sure nth is more than the Nyquist sampling
  if nth < (size(cofs,1)+1)
    warning('Sample finer to avoid aliasing')
  end
  % Define data grid
  theta=linspace(0,pi,nth);
  as=1; % Linearly spaced
else
  theta=nth(:)';
  nth=length(theta);
  as=0; % Not evenly spaced
end

% Construct matrix with spherical harmonics
% Here put modification to preload the Legendre functions
% Get Legendre function values if equally spaced on 0->pi
L=size(cofs,1)-1;
fnpl=sprintf('%s/LSS0-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'LEGENDRE'),L,nth);

if exist(fnpl,'file')==2 & as==1
  eval(sprintf('load %s',fnpl))
else  
  Pl=repmat(NaN,nth,L+1);
  for lin=0:L
    % This is SDW (2006) equation (3.2) without the sqrt(4pi)
    % BUT IT DOES INCLUDE sqrt(2-dom) !
    Pl(:,lin+1)=...
	rindeks(legendre(lin,cos(theta),'sch')*sqrt(2*lin+1),1)';
  end
  % For equally spaced data, save the polynomials in the original
  % normalization  
  if as==1
    eval(sprintf('save %s Pl',fnpl)) 
  end
end
% Watch the normalization in case it's changed
% This now is SDW (2006) equation (3.2) with the sqrt(4pi) if nrm=1
% BUT IT ALSO INCLUDES sqrt(2-dom) !
Pl=Pl/sqrt(4*pi)*sqrt(nrm);

% Expand to colatitudinal functions
fth=Pl*cofs;




