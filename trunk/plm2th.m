function varargout=plm2th(cofs,nth,m,nrm) 
% [r,nlon,lat]=plm2th(cofs,nth,m,nrm)
%
% Inverse spherical harmonic transform for a single m.
%
% Like PLM2XYZ, except only returns the th-dependent part,
% thus leaving out the cos(m*phi) or sin(m*phi) part.
%
% INPUT:
%
% cofs          Matrix or vector listing cosine or sine coefficients
% nth           Number of colatitudes on interval [0 pi] [default: 720]
%               If nth is a vector, calculates it at specific points
% m             The only order that is involved, m>=0
% nrm           Normalization factor for Ylm [default: 4pi], may want 1
% 
% OUTPUT:
%
% r             The colatitudinal function
% nlon,lat      The grid info, in degrees
% 
% EXAMPLE:
%
%
% Last modified by fjsimons-at-alum.mit.edu, 05/17/2011
%
% Could make this even faster, of course, but not now.

defval('nth',720);
defval('nrm',4*pi);
defval('m',1)

m=abs(m);

if length(nth)==1
  % Make sure nth is more than the Nyquist sampling
  if nth < (size(cofs,1)+1)
    warning('Sample finer to avoid aliasing')
  end
  % Define data grid
  theta=linspace(0,pi,nth);
  as=1; % Equally spaced
else
  theta=nth(:)';
  nth=length(theta);
  as=0; % Not equally spaced
end

% Assumes the listing is complete
L=size(cofs,1)+m-1;
lmin=0;

% Fills in the matrix properly... m added is latest fix
[dems,dels,mzero,lmcosi]=addmon(L,m);
cosi=repmat(lmcosi(:,3),1,size(cofs,2));
for lrnk=1:size(cofs,2)
  cosi(mzero(m+1:end)+m,lrnk)=cofs(:,lrnk);
end
cofs=cosi;

if L>255
%  disp('Using Libbrecht algorithm')
  libb=1;
else 
  libb=0;
end

nlon=2*nth-1;

% Initialize grid
r=repmat(0,nth,size(cofs,2)); 

% Get Legendre function values if equally spaced on 0->pi
fnpl=sprintf('%s/LSSM-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'LEGENDRE'),L,nth);

if exist(fnpl,'file')==2 & as==1
  %eval(sprintf('load %s',fnpl))
  load(fnpl)
else  
  % Evaluate Legendre polynomials at selected points
  Plm=repmat(NaN,nth,addmup(L));
  if L>200
    h=waitbar(0,'Evaluating all Legendre polynomials');
  end
  in1=0;
  in2=1;
  % Always start from the beginning
  for l=0:L
    % This is SDW (2006) equation (3.2) without the sqrt(4pi)
    % BUT IT DOES INCLUDE sqrt(2-dom) !
    if libb==0
      Plm(:,in1+1:in2)=(legendre(l,cos(theta(:)'),'sch')*sqrt(2*l+1))';
    else
      Plm(:,in1+1:in2)=(libbrecht(l,cos(theta(:)'),'sch')*sqrt(2*l+1))';
    end
    in1=in2;
    in2=in1+l+2;
    if L>200
      waitbar((l+1)/(L+1),h)
    end
  end
  if L>200
    delete(h)
  end
  if as==1
    eval(sprintf('save %s Plm ',fnpl))
  end
end

% Fix normalization if need be
% This now is SDW (2006) equation (3.2) with the sqrt(4pi) if nrm=1
% BUT IT DOES INCLUDE sqrt(2-dom) !
Plm=Plm/sqrt(4*pi)*sqrt(nrm);

% Loop over the degrees
more off
% disp(sprintf('PLM2TH Expansion from %i to %i',lmin,L))

for l=lmin:L
  % Always start from the beginning
  b=addmup(l-1)+1;
  e=addmup(l);

  % Sum over all orders and (through loop) over all degrees
  r=r+Plm(:,b:e)*cofs(b:e,:);
end

lat=90-theta*180/pi;

vars={'r','nlon','lat'};
for index=1:nargout
  varargout{index}=eval(vars{index});
end

