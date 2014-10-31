function [Cp,Sp]=rotcof(C,S,angl)
% [Cp,Sp]=ROTCOF(C,S,angl)
%
% Rotates real spherical harmonic coefficients over some angle
% using essentially DT (C.245): multiplication by exp(im*angle).
%
% INPUT:
%
% C      Cosine coefficients, as in LMCOSI(:,3), must start at 0
% S      Since coefficients, as in LMCOSI(:,4), must start at 0
% angl   Angle [radians]
%
% OUTPUT:
%
% Cp     Cosine coefficients of rotated field
% Sp     Sine coefficients of rotated field
%
% SEE ALSO: ROTTP, PLM2ROT
%
% Last modified by fjsimons-at-alum.mit.edu, 07/19/2010

% Figure out maximum spherical harmonic degree
L=addmup(length(C),'r');

% Initialize
Cp=C;
Sp=S;

% Calculate exponential
cangl=cos([0:L]*angl)';
sangl=sin([0:L]*angl)';

% warning('Changes. Updates consistent with ROTTP')

for l=0:L
  % Extract coefficients
  [Cl,b,e]=shcos(C,l);
  Sl=shsin(S,l);
  
  % Rotate and collect - sign flipped 03/11/2010
  Cp(b:e,1)=Cl.*cangl(1:l+1)+Sl.*sangl(1:l+1);
  Sp(b:e,1)=Sl.*cangl(1:l+1)-Cl.*sangl(1:l+1);
end

