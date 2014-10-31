function [C,rnk]=th2plm(fth,L,m)
% [C,rnk]=TH2PLM(fth,L,m)
%
% Forward spherical harmonic transform for a single m.
% Normalization to unity over unit sphere.
%
% INPUT:
%
% fth          Colatitudinal functions on equally spaced grid.
%              matrix of size [M,N], N functions evaluated at M points
% L            Maximum degree of expansion [default is Nyquist or 256]
% m            The only order that is involved
%
% OUTPUT:
%
% C            Matrix with spherical harmonic coefficients
%
% This requires a bit of work to get the normalization right, in line
% with TH2PL. Just live with it for now.
%
% Last modified by fjsimons-at-alum.mit.edu, 06/19/2008

if length(fth)==prod(size(fth))
  fth=fth(:);
end

nth=size(fth,1);
theta=linspace(0,pi,nth); % Always equally spaced
as=1;

% Calculate Nyquist degree
defval('L',min(nth-1,255));

% Choice of algorithm
if L>255
  disp('Using Libbrecht algorithm')
  libb=1;
else 
  libb=0;
end

if L>(nth-1)
  warning('Function undersampled. Aliasing will occur.')
end

% Get Legendre function values
fnpl=sprintf('%s/LSSM-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'LEGENDRE'),L,length(theta));

%disp(sprintf('Lnyq= %i ; expansion out to degree L= %i',nth-1,L))

waitmax=100;

if exist(fnpl,'file')==2 & as==1
  eval(sprintf('load %s',fnpl))
else  
  % Evaluate Legendre polynomials at selected points
  Plm=repmat(NaN,length(theta),addmup(L));
  if L>waitmax
    h=waitbar(0,'Evaluating all Legendre polynomials');
  end
  in1=0;
  in2=1;
  % Always start from the beginning
  for l=0:L
    if libb==0
      Plm(:,in1+1:in2)=(legendre(l,cos(theta(:)'),'sch')*sqrt(2*l+1))';
    else
      Plm(:,in1+1:in2)=(libbrecht(l,cos(theta(:)'),'sch')*sqrt(2*l+1))';
    end
    in1=in2;
    in2=in1+l+2;
    if L>waitmax
      waitbar((l+1)/(L+1),h)
    end
  end
  if L>waitmax
    delete(h)
  end
  if as==1
    eval(sprintf('save %s/LSSM-%i-%i Plm ',...
		 fullfile(getenv('IFILES'),'LEGENDRE'),L,length(theta)))
  end
end

[dems,dels,mz,lmcosi]=addmon(L);

% Locate the harmonics of order m
Pm=Plm(:,mz([m:L]+1)+m)/sqrt(4*pi);

% Retrieve the coefficients
[C,merr,mcov,chi2,L2err,rnk,Dw]=datafit(Pm,fth);
