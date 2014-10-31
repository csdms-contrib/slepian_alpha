function [R,C,Cnorm]=legendrecof(L,method)
% [R,C,Cnorm]=legendrecof(L,method)
%
% Calculates roots and coefficients of Legendre polynomials,
% by recursion. (Not the Associated Legendre functions!)
% Coefficients ordered with maximum power first, constant last.
% Roots lie between -1 and 1.
%
% INPUT:
%
% L         Maximum degree of Legendre polynomial
% method    'jacobi' stable algorithm for roots only
%           'cofrec' unstable algorithm for roots and coefficients
%
% OUTPUT:
%
% R         Roots of Legendre polynomial of degree L
% C         Coefficients of polynomial such that P(1)=1 and power 2/(2L+1)
% Cnorm     Coefficients of polynomial such that its power equals 1
%
%
% Note that the "standard" normalization is P(1)=1 and then,
% the inner product is returns 2/(2l+1). We also return 
% polynomials whose inner product over [-1 1] equals 1.
%
% EXAMPLE:
%
% plm=legendre(10,linspace(-1,1,100));
% plot(linspace(-1,1,100),plm(1,:)); hold on
% [R,C]=legendrecof(10);
% plot(R,0,'o'); grid on
%
%
% From
% http://perso.wanadoo.fr/jean-pierre.moreau/Basic/legendre_bas.txt
% http://dip.sun.ac.za/~weideman/research/differ.html
%
% By fjsimons-at-alum.mit.edu, August 21st, 2003

defval('method','jacobi')

switch method
 case 'jacobi'
  R=legendreroots(L);
  C=0;
  Cnorm=0;
 case 'cofrec'
  % The coefficients grow too large to be represented in the machine
  if L>40
    warning('This recursion relation is not valid for L> 40')
  end
  % Initialize
  B=zeros(L+1);
  B(1,1)=1;
  B(2,1)=0;
  B(2,2)=1;
  % Recursive relation
  % Columns are increasing powers of x ; first column is constant
  if L>=2
    for l=2:L
      B(l+1,1)=-(l-1)*B(l-1,1)/l;
      for p=1:l
	B(l+1,p+1)=((2*l-1)*B(l,p)-(l-1)*B(l-1,p+1))/l;
      end
    end
    % Take the last row
    for l=0:L
      C(l+1)=B(L+1,l+1);
    end
  elseif L==1
    C=[0 1];
  elseif L==0
    C=[1];
  end
  % Output - order matters!
  C=fliplr(C);
  % We've warned you; over 40 get complex roots 
  % sometimes exceeding [-1 1]
  R=real(roots(C));
  R(R<-1)=-1;
  R(R>1)=1;
  % Normalization to give inner product of one
  Cnorm=C*sqrt((2*L+1)/2);
end

