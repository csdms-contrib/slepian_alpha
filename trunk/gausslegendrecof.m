function [w,x,N]=gausslegendrecof(l,method,intv)
% [w,x,N]=gausslegendrecof(l,method,intv)
%
% Weights w and abscissas x for Gauss-Legendre integration of a
% polynomial function f(x) of degree l. The approximation w(:)'*f(x(:))
% to the integral W(x)f(x), W(x)=1 will be exact on [-1 1] and need only
% be rescaled to the desired interval. If intv is specified,
% w(:)'*f(x(:)) will do, as we do the scaling for you.
%
% INPUT:
%
% l           Degree of polynomial in the integrand
% method      'jacobi' stable algorithm (Legendre roots only) [default]
%             'cofrec' algorithm for l<40 (roots and coefficients)
% intv        If integration interval specified, returns scaled [w,x]
%             This may be an Nx2 matrix if you want multiple intervals
%
% OUTPUT:
%
% w           Weights (vector, or matrix in case of more than one TH)
% x           Abscissas (roots of Legendre polynomial of degree N)
% N           Number of points in the integration (length(x))
%
% SEE ALSO: GAUSSLEGENDRE
% 
% The abcissa's are the roots of a Legendre polynomials defined on the
% same interval. In particular then, this routine can be used to
% integrate (products) of Legendre functions themselves, which is useful
% in the analysis of spherical harmonics. This returns the N-point
% Gauss-Legendre integration points and weights, which is accurate for
% polynomials up to degree 2*N-1. Note that these things are symmetric so
% for the N-point integration there are only ceil(N/2) unique
% weights. The weights are positive, symmetric, and should sum to 2.
% As N nodes integrate a polynomial of degree 2N-1 exactly, the number of
% nodes returned for a polynomial degree l is ceil((l+1)/2).
%
% EXAMPLE: (from Numerical Recipes [qgaus, p. 141, Chapter 4.5])
%
% [w,x]=gausslegendrecof(19);
%
% Last modified by fjsimons-at-alum.mit.edu, 03/17/2009

defval('method','jacobi')
defval('intv',[])

% Figure out what the N will need to be
N=ceil((l+1)/2);

switch method
 case 'jacobi'
  % Find roots of Legendre polynomial
  x=legendrecof(N,method);
  % Find weight values from the derivative of the Legendre polynomial
  % Calculate Legendre function at N-1; note Pl(1)=1, which is not
  % the normalization we want in order to make the inner product of the
  % spherical harmonic be 4 pi! Now, for m=0, the inner product is
  % 2/(2N+1). These belong to the standard, listed Legendre polynomials.
  Pl=rindeks(legendre(N-1,x,'sch'),1);  
  Pdiff=-N*Pl(:)./(x(:).^2-1);
 case 'cofrec'
  % Find roots of Legendre polynomial of degree N in [-1 1]
  % These are the abscissa points
  [x,C1]=legendrecof(N,method);
  % Find the value of the derivative of the Legendre function at these
  % points, by involving the lower-degree polynomial as well
  [y,C2]=legendrecof(N-1,method);
  x=sort(x);
  % This is the equation for the derivative of the Legendre polynomial
  % as found in Numerical Recipes and other texts. Wolfram has this
  % equation under Legendre-Gauss Quadrature, not under Legendre Polynomials.
  % Pdiff=N*(x.*polyval(C1,x)-polyval(C2,x))./(x.^2-1);
  % But we are noting that we are evaluating this at the roots of C1
  % so we can save one computation
  % Note we are using coefficients of Legendre Polynomials that integrate
  % to 2/(2L+1), i.e. the "classical ones"
  Pdiff=-N*polyval(C2,x)./(x.^2-1);
end

% Find weights
% Numerical Recipes, Equation (4.5.16)
w=2./(1-x.^2)./Pdiff.^2;

% Scale to interval
if ~isempty(intv) 
  if length(intv(:))==2
    a=intv(1);
    b=intv(2);
    % Rescale nodes
    x=a+(x+1)/2*(b-a);
    % Rescale weights
    w=w*(b-a)/2;
  else
    a=intv(:,1);
    b=intv(:,2);
    xg=x;wg=w;
    x=repmat(NaN,length(xg),length(a));
    w=repmat(NaN,length(wg),length(a));
    for index=1:length(a)
      x(:,index)=a(index)+(xg+1)/2*(b(index)-a(index));
      w(:,index)=wg*(b(index)-a(index))/2;
    end
  end
end

