function [dv,foldi,be]=degamini(v)
%  [dv,foldi,be]=DEGAMINI(v)
%
% Inverse operation for GAMINI
%
% INPUT
%
% v       A vector with repeated entries (not necessarily sorted)
% 
% OUTPUT:
%
% dv      The same vector with all the repeats set to 1
% foldi   A same-dimensional vector with how many repeats there were
% be      A matrix with begin and end indices into the original vector
%
% SEE ALSO:
%
% GAMINI, MATRANGES
%
% Last modified by fjsimons-at-alum.mit.edu, April 17th, 2007

v=v(:)';
indi=[1 find(~~diff(v))+1];

dv=v(indi);

if nargout > 1,
  foldi=diff([indi length(v)+1]);
  if length(v)~=sum(foldi)
    error('Incorrect')
  end
end

if nargout > 2,
  be=[indi(:) cumsum(foldi(:))];
end
