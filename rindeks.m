function x=rindeks(M,i)
% x=RINDEKS(M,i)
%
% Returns the ROW(s) i of a matrix M.
% Returns the ith 1-st dimension(s) of a matrix M.
%
% EXAMPLE:
%
% M=rand(randij(12))
% witsj=randi(size(M,1))
% rindeks(M,witsj)
%
% See also INDEKS, KINDEKS, TINDEKS, DINDEKS, SQUEEZE
%
% Last modified by fjsimons-at-alum.mit.edu, 01/30/2019

x=M(i(:),:,:);
