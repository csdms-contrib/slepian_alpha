function R=rots(L,V,EM,gammas)
% R=ROTS(L,V,EM,gammas)
%
% Makes a pole-rotation matrix for Slepian functions on axisymmetric
% polar domains. It does this by combining non-zero orders.
%
% INPUT:
%
% L        Bandwidth of the Slepian basis
% V        Eigenvalues of the Slepian basis
% EM       List of orders in which the basis is presented (by, e.g. GALPHA)
% gammas   List of radian angles over which you want to rotate each function
%
% OUTPUT:
%
% R        A blocky rotation matrix to multiply into the G of GALPHA
%
% EXAMPLE:
%
% L=18; 
% [G,V,EM]=galpha(40,L,1,linspace(0,pi,50),linspace(0,2*pi,100),'global');
% R=rots(L,V,EM,rand(addmup(L),1)); GR=R*G; subplot(211);
% imagesc(reshape(G(2,:),50,100)); subplot(212); imagesc(reshape(GR(2,:),50,100)); 
%
% SEE ALSO: GALPHA
%
% Last modified by fjsimons-at-alum.mit.edu, 08/19/2008

defval('L',3)

if length(V)~=length(EM)
  error('V and EM must be the same length')
end

if ~exist('V') || (exist('V') & isempty(V))
  % Just to get some defaults going
  [G,V,EM]=galpha(40,L,1,0,NaN,'global');
end

defval('gammas',pi/4)

if length(gammas)==1
  gammas=repmat(gammas,addmup(L),1);
end

% Which are the abs(orders) in question here?
% EM=abs(EM(logical([1 ~~diff(abs(EM))]))); % Nope
% The previous often backfire if for some reason two m=0 or two |m|
% are in direct succession... the real discriminant is:
EM=abs(EM(logical([1 ~~diff(V)])));

% Check the dimensions are right...
difer(length(EM)-length(unique(V)))
difer(length(EM)-length(gammas))

% Now construct the grand old rotation matrix
% This growth might have to be optimized later
R=[];
for index=1:length(EM)
  if EM(index)==0
    lR=1;
  else
    gams=gammas(index);
    lR=[cos(gams) -sin(gams) ; sin(gams) cos(gams)];
  end
  R=blkdiag(R,lR);
end
