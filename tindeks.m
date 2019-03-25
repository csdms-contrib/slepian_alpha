function x=tindeks(M,i)
% x=TINDEKS(M,i)
%
% Returns the ith 3-rd dimension(s) of a 3D matrix M.
%
% EXAMPLE:
%
% M=rand([randij(12) randi(4)])
% witsj=randi(size(M,3))
% tindeks(M,witsj)
%
% See also INDEKS, RINDEKS, KINDEKS, DINDEKS, SQUEEZE
%
% Last modified by fjsimons-at-alum.mit.edu, 01/30/2019

x=M(:,:,i(:));

