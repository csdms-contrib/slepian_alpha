function [C,rnk]=th2pl(fth,L,method,intv,normfac) 
% [C,rnk]=TH2PL(fth,L,method,intv,normfac)
%
% Forward Legendre transform, undone by PL2TH.
%
% Calculates ZONAL spherical harmonic coefficients (given starting with
% L=0) from an input spatial function. Normalization is to the area of the
% unit sphere. We use 'schmidt'*sqrt(2l+1) to get spherical harmonics
% normalized to 4\pi. So the reconstruction integral has a factor of
% (4-2*(m==0)) in it reflecting the normalization of the polynomials.
%
% f_l=\int_{0}^{\pi}f(\theta)P_l(\cos(theta))\sin(\theta)\,d\theta
%
% INPUT:
%
% fth       The colatitudinal function, starting at 0 and ending on pi 
%           This could be a matrix, down the columns
% L         Maximum degree of expansion [default is Nyquist or 256]
% method    'gl' by Gauss-Legendre integration (not recommended)
%           'im' by inversion [default]
%           'simpson' method (not recommended)
% inv       Interval (default is all: [-1 1])
% normfac   Leave alone for the Legendre transform corresponding to PL2TH
%           but may set to any other value if you just want to integrate
% 
% OUTPUT:
%
% C         Vector/Matrix with zonal spherical harmonics coefficients
% rnk       Rank of the SVD
%
% EXAMPLE: Compare some integration rules:
%
% th=linspace(0,pi,100);
% th2pl(sin(th),0,'gl',[],1)
% th2pl(sin(th),0,'im')
% th2pl(sin(th),0,'simpson')
% gausslegendre([0 pi],inline('sin(x)'),15)
% simpson(th,sin(th))
% trapeze(th,sin(th))
%
% What you need to do is the integral from 0 to pi of the function
% multiplied by the Legendre polynomial. This can be done by Gaussian
% integration. On the other hand, you can invert for the coefficients by
% writing the expansion as a matrix equation. This is more accurate.
%
% See also PL2TH, PLM2XYZ, XYZ2PLM. Watch the factors, e.g. SDWCAPT2.
% 
% Last modified by fjsimons-at-alum.mit.edu, 14.07.2005

if length(fth)==prod(size(fth))
  fth=fth(:);
end

nth=size(fth,1);
theta=linspace(0,pi,nth);

% Decide on the Nyquist frequency
defval('L',min(nth-1,255));
defval('method','im')
defval('intv',[-1 1])

% Only zonal terms, unassociated polynomials
m=0; % And normalization in line with PL2TH
defval('normfac',(4-2*(m==0)));

if L>(nth-1)
  warning('Function undersampled. Aliasing will occur.')
end

% Define evaluation points
switch method
 case 'gl'
  disp('This method may lead to highly inaccurate results')
  % Highest degree of integrand will always be 2*L
  [w,x]=gausslegendrecof(2*L,[],intv);
  % This is cheating and why the inversion method is best
  % We're coming up with f evaluated at the GL points
  fthx=interp1(cos(theta),fth,x,'spline');
  as=0; % Not equally spaced
 case 'im'
  % Where to evaluate the Legendre polynomials
  if ~all(intv==[-1 1])
    error('Choose GL method')
  end
  x=cos(theta);
  as=1; % Equally spaced
 case 'simpson'
  % Where to evaluate the Legendre polynomials
  warning('This method may lead to highly inaccurate results')
  if ~all(intv==[-1 1])
    error('Choose GL method')
  end
  x=cos(theta);
  as=1; % Equally spaced
 otherwise
  error('Specify valid method.')
end

% Get zonal Legendre function values
fnpl=sprintf('%s/LSS0-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'LEGENDRE'),L,length(x));

%disp(sprintf('Lnyq= %i ; expansion out to degree L= %i',nth-1,L))
  
if exist(fnpl,'file')==2 & as==1
  eval(sprintf('load %s',fnpl))
else  
  % Evaluate Legendre polynomials at (Gauss-Legendre) points
  Pl=repmat(NaN,length(x),L+1);
  if L>200
    h=waitbar(0,'Evaluating zonal Legendre polynomials');
  end
  for l=0:L
    Pl(:,l+1)=rindeks(legendre(l,x(:)','sch')*sqrt(2*l+1),1)';
    if L>200
      waitbar((l+1)/(L+1),h)
    end
  end
  if L>200
    delete(h)
  end
  if as==1
    eval(sprintf('save %s Pl',fnpl)) 
end
end

% Calculate expansion
switch method
 case 'gl'
  C=(fthx'*diag(w)*Pl)'*(intv(2)-intv(1))/2/normfac;
  rnk=[]; 
 case 'im'
  [C,merr,mcov,chi2,L2err,rnk,dw]=datafit(Pl,fth);
 case 'simpson'
  C=simpson(theta,repmat(fth(:).*sin(theta(:)),1,size(Pl,2)).*Pl)/normfac;
  rnk=[]; 
end

