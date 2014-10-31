function [P,mu,norms]=plm(l,m,mu,check,tol)
% [P,mu,norms]=PLM(l,m,mu,check,tol)
%
% Calculates (associated) Legendre functions, DT (B.48/B.56/B.71).
%
% INPUT:
%
% l      degree (0 <= l <= infinity) [default: random]
% m      order (-l <= m <= l)        [default: all orders 0<=l]
%        l and m can be vectors, but not both at the same time
% mu     argument (-1 <= mu <= 1)    [default: 181 linearly spaced]
% check  1 optional normalization check by Gauss-Legendre quadrature
%        0 no normalization check [default]
% tol    Tolerance for optional normalization checking
%
% OUTPUT:
%
% P      The associated Legendre function at the desired argument(s):
%           as a scalar or a row vector with length(mu) columns, OR
%           as a matrix with length(m) rows and length(mu) columns, OR 
%           as a matrix with length(l) rows and length(mu) columns, OR
%           as a 3D matrix of size [length(l) size(theta), OR
%           (L+1)^2 x length(theta) if you put in
%           a degree l=[0 L] and an order []: lists orders -l to l.
% mu     The argument(s), which you might or not have specified
% norms  The normalization matrix, which should be the identity matrix
%
% EXAMPLES:
%
% plot(plm([0:5],0,[],1)')
% plot(plm(5,[])')
%
% SEE ALSO:
%
% LIBBRECHT, XLM, YLM
%
% Last modified by fjsimons-at-alum.mit.edu, 10/07/2006

% Default values
defval('l',round(rand*10))
defval('m',[])
defval('mu',linspace(-1,1,181))
defval('check',0)
defval('tol',1e-10)

% Error handling common to PLM, XLM, YLM
[l,m,mu,check,tol]=pxyerh(l,m,mu,check,tol);

% Warning against using such unnormalized polynomials
if any(2*l > 21) && any(m<0)
  disp('Factorials large, negative orders will be inaccurate')
end

switch check
  case 0
   % Calculation for m>=0
   if prod(size(l))==1 & prod(size(m))==1 % SINGLE L SINGLE M
     % Note that Matlab's unnormalized functions have the (-1)^m
     % Condon-Shortley phase in there so we get rid of it
     P=(-1)^(-m)*rindeks(legendre(l,mu),abs(m)+1);
     if m < 0
       P=(-1)^(-m)*factorial(l+m)./factorial(l-m)*P;
     end
   elseif prod(size(l))==1 % SINGLE L MULTIPLE M
     P=(rindeks(legendre(l,mu),abs(m)+1)'*diag((-1).^(-m)))';
     for index=find(m<0)
       P(index,:)=(-1)^(-m(index))*...
	   factorial(l+m(index))./factorial(l-m(index))*P(index,:);
     end
   elseif prod(size(m))==1 % MUTIPLE L SINGLE M
     if min(size(mu))==1 % VECTOR ARGUMENT
       P=repmat(NaN,[length(l) length(mu) 1]);
       ini=1;
     else % MATRIX ARGUMENT
       P=repmat(NaN,[length(l) size(mu)]);
       % The following is to bypass a dimensionality convention
       P(1,:,:)=shiftdim((-1)^m*legendre(0,mu),-1);
       ini=2;
     end
     for index=ini:length(l)
       P(index,:,:)=(-1)^m*rindeks(legendre(l(index),mu),abs(m)+1);
       if m<0
	 P(index,:,:)=(-1)^m*...
	     factorial(l(index)+m)./factorial(l(index)-m)*P(index,:,:);
       end
     end
   else
     error('Specify valid option')
   end
   norms=[];
 case 1
  % Check normalization 
  norms=pxynrm(l,m,tol,'P');
  % Still also need to get you the answer at the desired argument
  P=plm(l,m,mu,0);
end
