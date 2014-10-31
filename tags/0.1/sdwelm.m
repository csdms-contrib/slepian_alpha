function varargout=sdwelm(TH,L)
% [lrnk,mrnk,lval,VV,Vsum]=SDWELM(TH,L)
%
% Computes global ranking of eigenvalues for the concentration problem
% to the SINGLE SPHERICAL CAP considering all angular orders m.
%
% INPUT:
%
% TH      Colatitudinal radius of the concentration region
% L       Bandwidth (maximum spherical harmonic degree)
% 
% OUTPUT:
% 
% lrnk    The index of the eigenvalue within its own m
% mrnk    The m of the eigenvalue
% lval    The sorted eigenvalues
% VV      The not globally sorted eigenvalue matrix (lrank-by-m)
% Vsum    Total sum of all the eigenvalues, including double counts
%
% Last modified by fjsimons-at-alum.mit.edu, 04/24/2009

defval('TH',40);
defval('L',18);

% Initialize eigenvalue matrix
VV=repmat(NaN,L+1,L+1);
for m=0:L;
  % Make nth really small since you won't use it, really
  [E,V]=sdwcap(TH,L,m,0);
  % Now compare with the EXACT grunbaum ordering
  % [E,Vg,th,C,T,V]=grunbaum(TH,L,m,nth,grd)
  % No wait, this won't make a difference as I cannot figure out how the
  % Grunbaum eigenvalues should be sorted BETWEEN them - their sorting
  % works within a single order... too bad, really.
  % All the eigenvalues; sometimes slightly negative
  VV(1:min(L-m+1,length(V)),m+1)=V(:);
  % This is only filled for possible positive degrees from 0 to L
end

% Figure out GLOBAL rank ordering for the eigenvalues, with the repeated m
[a,b]=sort(VV(:));
% Put this in a matrix
mrnk=repmat(0:L+1,L+1,1);
lrnk=repmat([1:L+1]',1,L+1);
b=b(~isnan(a)); b=flipud(b);
a=a(~isnan(a)); a=flipud(a);
% The ranked eigenvalues belong to this m
mrnk=mrnk(b);
% And they represent this number of nth eigenvalue
lrnk=lrnk(b);
lval=a;

% Figure out total sum of the eigenvalues including double counts
Vsum=VV;
Vsum(:,2:end)=Vsum(:,2:end)*2;
Vsum=sum(Vsum(~isnan(Vsum)));

% Output
varn={lrnk,mrnk,lval,VV,Vsum};
varargout=varn(1:nargout);

