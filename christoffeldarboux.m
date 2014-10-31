function varargout=christoffeldarboux(ks,m,L,mu,mup)
% CHRISTOFFELDARBOUX(ks,m,L,mu,mup)
% [m,L,D]=CHRISTOFFELDARBOUX(ks,m,L,mu,mup)
%
% INPUT:
% ks   1 Standard formula with (mu-mup)
%      2 Modified formula with (mu-mup)*(mu+mup)
%      3 Modified as in 2 but with only [m:2:L]
%      4 Modified as in 2 but with only [m+1:2:L]
% 
% OUTPUT:
%
% m    Order
% L    Maximum degree
% D    Difference between both approaches
%
% Several versions of the Christoffel-Darboux formula quoted by Simons,
% Dahlen and Wieczorek (2006) eqs (3.10/5.15)
%
% Last modified by fjsimons-at-alum.mit.edu, 04.05.2006

defval('ks',ceil(rand*4))
defval('m',round(rand*20))
defval('L',max(round(rand*50),m))
defval('mu',(rand-0.5)*2)
defval('mup',(rand-0.5)*2)

if length(mu)>1 | length(mup)>1
  error('Arguments must be scalar')
end

if m>L
  error('Order must be smaller than degree')
end

% The left hand side of the equation
switch ks
 case 3
  msum=[m:2:L];
 case 4
  msum=[m+1:2:L];
 otherwise
  msum=[m:L];
end

msg=sprintf('k= %i ; m = %2.2i ; L= %2.2i',ks,m,L);
disp(msg)

if ~isempty(msum)
  K1=0;
  for l=msum
    Pl=rindeks(legendre(l,[mu mup]),m+1);
    K1=K1+alm(l,m)*prod(Pl);
  end
  switch ks
   case 1
    K1=(mu-mup)*K1;
   otherwise
    K1=(mu^2-mup^2)*K1;
  end
  
  % The right hand side of the equation
  % This is the upper degree limit of the summation
  Lp=msum(end);
  % Calculate the required Legendre polynomials
  PL=Pl;
  PL1=rindeks(legendre(Lp+1,[mu mup]),m+1);
  if ks>1
    PL2=rindeks(legendre(Lp+2,[mu mup]),m+1);
    B1=clm(Lp,m)*[PL2(1)*PL(2)-PL(1)*PL2(2)];
    % If this is the evens, you would have don the odds now
    if Lp>=1 & m<Lp
      PLm1=rindeks(legendre(Lp-1,[mu mup]),m+1);
      B2=clm(Lp-1,m)*[PL1(1)*PLm1(2)-PLm1(1)*PL1(2)];
    else
      B2=0;
    end
  end
  
  switch ks
   case 1
    K2=blm(Lp,m)*[PL1(1)*PL(2)-PL(1)*PL1(2)];
   case 2
    K2=B1+B2;
   otherwise
    K2=B1;
  end
  
  % Evaluate the difference
  D=abs(K1-K2);
else
  D=NaN;
end

difer(D)

% OUTPUT:
varn={'m','L','D'};
for index=1:nargout
  varargout{index}=eval(varn{index});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=alm(l,m)
% Computes the first prefactor, eq. (5) of SD2006
fax=l-m+1:l+m;
if ~isempty(fax) & ~~prod(fax)
  A=(2*l+1)/prod(fax);
else
  A=(2*l+1)*factorial(l-m)/factorial(l+m);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B=blm(l,m)
% Computes the second prefactor, for the regular case
% This is eq. (5.15) of SDW2006
fax=l-m+2:l+m;
if ~isempty(fax) & ~~prod(fax)
  B=1/prod(fax);
else
  B=factorial(l-m+1)/factorial(l+m)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C=clm(l,m)
% Computes the third prefactor, for the special case
% This is eq. (A.13) of SD2006
fax=l-m+3:l+m;
if ~isempty(fax) & ~~prod(fax)
  C=1/prod(fax)/(2*l+3);
else
  C=factorial(l-m+2)/factorial(l+m)/(2*l+3);
end
