function [X,theta,dems]=xlm(l,m,theta,xver,tol,blox)
% [X,theta,dems]=XLM(l,m,theta,xver,tol,blox)
%
% Calculates normalized (associate) Legendre functions, DT (B.58/B.60).
%
% INPUT:
%
% l      degree (0 <= l <= infinity) [default: random]
% m      order (-l <= m <= l)        [default: all orders 0<=l]
%        l and m can be vectors, but not both at the same time
% theta  argument (0 <= theta <= pi) [default: 181 linearly spaced; not NaN!]
% xver   1 optional normalization check by Gauss-Legendre quadrature
%        0 no normalization check [default]
% tol    Tolerance for optional normalization checking
% blox   0 Standard (lm) ordering, l=0:L, m=-l:l
%        1 Block-diagonal ordering, m=-L:L, l=abs(m):L
%
% OUTPUT:
%
% X      The associated normalized Legendre function at the desired argument(s):
%           as a scalar or a row vector with length(theta) columns, OR
%           as a matrix with length(m) rows and length(theta) columns, OR 
%           as a matrix with length(l) rows and length(theta) columns, OR
%           as a 3D matrix of size [length(l) size(theta), OR
%           (L+1)^2 x length(theta) if you put in
%           a degree l=[0 L] and an order []: lists orders -l to l.
% theta  The argument(s), which you might or not have specified
% dems   The orders to which the Xlms belong (to verify input or block sorting)
%
% EXAMPLES:
%
% plot(xlm([0:5],0)')
% plot(xlm(5,[])')
%
% SEE ALSO:
%
% LIBBRECHT, PLM, YLM
%
% Last modified by fjsimons-at-alum.mit.edu, 01/18/2008

% For single m, this will accept 2D theta's and return 3D results X
% See Research Notebook VI p 77ff

% Default values
defval('l',round(rand*10))
defval('m',[])
defval('theta',linspace(0,pi,181))
% Never make this one or it will reach internal recursion limit
defval('xver',0)
defval('tol',1e-10)
defval('blox',0)

if blox~=0 & blox~=1
  error('Specify valid block-sorting option ''blox''')
end

% Revert back to cos(theta)
mu=cos(theta);

switch xver
 case 0
  % If the degrees go from 0 to some L and m is empty, know what to do
  if min(l)==0 & max(l)>0 & isempty(m)
    X=repmat(NaN,(max(l)+1)^2,length(theta));
    for thel=0:max(l)
      X(thel^2+1:(thel+1)^2,:)=xlm(thel,-thel:thel,theta,xver,tol);
    end
    [dems,dels,mz,blkm]=addmout(max(l));
    if blox==1
      X=X(blkm,:);
      dems=dems(blkm);
    end
    return
  end

  dems=m;
  
  % Error handling common to PLM, XLM, YLM - note this resets defaults
  [l,m,mu,xver,tol]=pxyerh(l,m,mu,xver,tol);

  % Calculation for m>0
  if prod(size(l))==1 & prod(size(m))==1 % SINGLE L AND M
    % Note that Matlab 'sch' has the sqrt(2-(m==0)) in there so we get
    % rid of it; this option also has gotten rid of the Condon-Shortley
    % phase which we now need to put back in
    X=(-1)^m*sqrt(2*l+1)/sqrt(2-(m==0))/sqrt(4*pi)*...
      rindeks(legendre(l,mu,'sch'),abs(m)+1);
    % This straight from the rule DT B.60
    if m<0 
      X=(-1)^m*X;
    end
  elseif prod(size(l))==1 % SINGLE L MULTIPLE M
    X=sqrt(2*l+1)/sqrt(4*pi)*...
      (rindeks(legendre(l,mu,'sch'),abs(m)+1)'*...
       diag(((-1).^m)./sqrt(2-(m==0))))';
    for index=find(m<0)
      X(index,:)=(-1)^m(index)*X(index,:);
    end
  elseif prod(size(m))==1 % MUTIPLE L SINGLE M
    if min(size(mu))==1 % VECTOR ARGUMENT
      X=repmat(NaN,[length(l) length(mu) 1]);
      ini=1;
    else % MATRIX ARGUMENT
      X=repmat(NaN,[length(l) size(mu)]);
      % The following is to bypass a dimensionality convention
      X(1,:,:)=shiftdim((-1)^m/sqrt(2-(m==0))/sqrt(4*pi)*...
			legendre(0,mu,'sch'),-1);
      ini=2;
    end
    for index=ini:length(l)
      X(index,:,:)=(-1)^m*sqrt(2*l(index)+1)/sqrt(2-(m==0))/sqrt(4*pi)*...
	  rindeks(legendre(l(index),mu,'sch'),abs(m)+1);
    end
    if m<0
      X=(-1)^m*X;
    end
  else
    error('Specify valid option')
  end
 case 1
  % Check normalization... only for different 'ms
  pxynrm(l,unique(m),tol,'X');
  % Still also need to get you the answer at the desired argument
  [X,theta,dems]=xlm(l,m,theta,0);
  % So you can't return norms anymore - since you force not to
end

