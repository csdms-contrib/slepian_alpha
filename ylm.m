function [Y,theta,phi,dems,dels]=ylm(l,m,theta,phi,check,tol,blox,irr)
% [Y,theta,phi,dems,dels]=YLM(l,m,theta,phi,check,tol,blox,irr)
%
% Calculates unit-normalized real spherical harmonics, 
% Dahlen and Tromp (1998), Theoretical Global Seismology, eq. (B.72).
% Here, Y00=1/sqrt(4*pi), i.e. difer(indeks(ylm(0,0),1)*sqrt(4*pi)-1)
%
% INPUT:
%
% l      degree(s) (0 <= l <= infinity) [default: random]
% m      order (-l <= m <= l)        [default: all orders -l:l]
%        l and m can be vectors, but not both at the same time
%        If l is [0 L] and m is empty, do all degrees from 0 to L
%        If l is a single number and m is empty, it does m=0:l
% theta  colatitude vector (0<=theta<=pi) [default: 181 linearly spaced; not NaN!]
% phi    longitude vector (0<=phi<=2*pi) [default: 361 linearly spaced]
%        Unless irr=1, we assume you mean a 2-D (theta,phi) grid.
%        But if irr=1, length(theta(:)) must be equal to length(phi(:)).
% check  1 optional normalization check by Gauss-Legendre quadrature
%        0 no normalization check [default]
% tol    Tolerance for optional normalization checking [default: 1e-10]
% blox   0 Standard (lm) ordering, as ADDMOUT, l=0:L, m=-l:l [default]
%        1 Block-diagonal ordering, m=-L:L, l=abs(m):L
% irr    0 Regular grid, no matter how you input theta and phi [default]
%        1 Irregular grid, input interpreted as distinct pairs of theta, phi
%        Note that irr is a variable in the financial toolbox
%
% OUTPUT:
%
% Y      The real spherical harmonics at the desired argument(s):
%           As matrix with dimensions of 
%           length(theta) x length(phi) x max(length(m),length(l)) OR
%           (L+1)^2 x (length(theta)*length(phi)) if you put in
%           a degree l=[0 L] and an order []: lists orders -l to l OR
%           max(length(m),length(l)) or (L+1)^2 x length(theta) for irr=1
% theta  The latitude(s), which you might or not have specified
% phi    The longitude(s), which you might or not have specified
% dems   The orders to which the Ylms belong
% dels   The degrees to which the Ylms belong [if input needed interpreting]
%
% NOTES:
%
% This is a recursive algorithm which uses recursion. Get it?
%
% EXAMPLES:
%
% ylm('demo1')
% ylm('demo2') % Compare with PLM2XYZ which is differently normalized/phased
% 
% SEE ALSO:
%
% LIBBRECHT, PLM, XLM, PLM2XYZ
%
% Last modified by fjsimons-at-alum.mit.edu, 06/23/2011

% Default values
defval('l',round(rand*10))

if ~isstr(l)
  defval('m',[])
  defval('theta',linspace(0,pi,181))
  defval('phi',linspace(0,2*pi,361))
  defval('check',0)
  defval('tol',1e-10)
  defval('blox',0)
  defval('irr',0)

  if blox~=0 & blox~=1
    error('Specify valid block-sorting option ''blox''')
  end

  if irr==1 & ~all(size(theta(:))==size(phi(:)))
    error('Input arrays must have the same dimensions for irregular grids')
  end

  % Make sure phi is a row vector
  phi=phi(:)';

  % If the degrees go from 0 to some L and m is empty, know what to do...
  if min(l)==0 & max(l)>0 & isempty(m)
    if irr==0
      % Here you assume a regular grid
      [PH,TH]=meshgrid(phi,theta);
    else
      PH=phi(:);
      TH=theta(:);
    end
    Y=repmat(NaN,(max(l)+1)^2,length(TH(:)));
    for thel=0:max(l)
      % Because here, too, you do the orders explicitly
      theY=reshape(ylm(thel,-thel:thel,theta,phi,check,tol,[],irr),...
		   length(TH(:)),2*thel+1)';
      Y(thel^2+1:(thel+1)^2,:)=theY;
    end
    theta=TH(:);
    phi=PH(:);
    
    [dems,dels,mz,blkm]=addmout(max(l));
    if blox==1
      Y=Y(blkm,:);
      dems=dems(blkm);
      dels=dels(blkm);
    end
    return
  end

  % ...if not, you're doing a single degree, perhaps as part of the above

  % Error handling common to PLM, XLM, YLM - note this resets defaults
  [l,m,mu,check,tol]=pxyerh(l,m,cos(theta),check,tol);

  % Straight to the calculation, check normalization on the XLM
  % Can make this faster, obviously, we're doing twice the work here
  [X,theta,dems]=xlm(l,abs(m),theta,check);
  % The second dimension is always length(theta)

  if irr==0
    % Initialize the matrix with the spherical harmonics
    Y=repmat(NaN,[length(theta) length(phi) max(length(m),length(l))]);
    
    % Make the longitudinal phase: ones, sines or cosines, sqrt(2) or not 
    % The negative m is the cosine
    P=repmat(diag(sqrt(2-(m(:)==0)))*...
	     cos(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi))),length(l),1);
    
    % Make the real spherical harmonics
    % Here too, you assume a regular grid
    if prod(size(l))==1 & prod(size(m))==1
      Y=X'*P;
    else
      for index=1:max(length(m),length(l))
	Y(:,:,index)=X(index,:)'*P(index,:);
      end
    end
  else
    % Initialize the matrix with the spherical harmonics
    Y=repmat(NaN,[length(theta) max(length(m),length(l))]);
    
    % Make the longitudinal phase: ones, sines or cosines, sqrt(2) or not 
    % The negative m is the cosine
    P=repmat(diag(sqrt(2-(m(:)==0)))*...
	     cos(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi))),length(l),1);

    % Make the real spherical harmonics
    % Here too, you assume a regular grid
    if prod(size(l))==1 & prod(size(m))==1
      Y=X.*P;
    else
      for index=1:max(length(m),length(l))
	Y(:,index)=[X(index,:).*P(index,:)]';
      end
    end
  end
elseif strcmp(l,'demo1')
  plotplm(ylm(2,-1,[],[],1),[],[],1); colormap(flipud(gray(7)))
  Y=ylm(0:10);
  Y=ylm(2,-1,[pi/2-linspace(pi/6,-pi/6,100)],[linspace(pi/2,pi,100)]);
  Y=ylm(0:10,[],[pi/2-linspace(pi/6,-pi/6,100)],[linspace(pi/2,pi,100)]);
  theta=linspace(0,pi,90); phi=linspace(0,2*pi,90);
  difer(diag(ylm(3,2,theta,phi,[],[],[],0))'-ylm(3,2,theta,phi,[],[],[],1))
  [lon2,lat2]=grcazim;
  Y=ylm([0 3],[],[90-lat2]*pi/180,lon2*pi/180,[],[],[],1);
elseif strcmp(l,'demo2')
  L=60;
  % Load a filtered geoid
  v=plmfilt(kindeks(fralmanac('EGM96','SHM'),1:4),L);
  v(1,3)=0;
  % Expand and plot one way
  [r1,lon,lat]=plm2xyz(v,1);
  clf
  subplot(211)
  imagefnan([lon(1) lat(1)],[lon(end) lat(end)],r1)
  plotcont
  % Expand and plot the other way
  v=[zeros(3,4) ; v];
  [Y,ph,th]=ylm([0 L+2],[],[],[]);
  % Fix the phase and the normalization
  EM=addmout(L+2);
  Y=[Y.*repmat((-1).^EM,1,size(Y,2))]*sqrt(4*pi);
  [dems,dels,mz,lmcosi,mzi,mzo,bigm,bigl,rinm,ronm,demin]=addmon(L+2);
  v=v(:,3:4);
  r2=reshape(v(ronm)'*Y,181,361);
  subplot(212)
  imagefnan([lon(1) lat(1)],[lon(end) lat(end)],r2)
  plotcont
  difer(r1-r2)
end
