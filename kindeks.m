function x=kindeks(M,i)
% x=KINDEKS(M,i)
%
% Returns the COLUMN(s) i of a matrix M.
% Returns the ith 2-nd dimension(s) of a matrix M.
%
% EXAMPLE:
%
% M=rand(randij(12))
% witsj=randi(size(M,2))
% kindeks(M,witsj)
%
% See also INDEKS, RINDEKS, TINDEKS, DINDEKS, SQUEEZE.
%
% Last modified by fjsimons-at-alum.mit.edu, 01/30/2019

x=M(:,i(:),:);

