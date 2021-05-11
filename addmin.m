function vec=addmin(invec)
% vec=ADDMIN(invec)
%
% Repeats the entries of a column vector the row number of times.
% Useful for real spherical harmonics. 
%
% INPUT: 
%
% invec       Input array with numbers or strings
%
% OUTPUT:
%
% vec         Output array with repeated entries
%
% Last modified by fjsimons-at-alum.mit.edu, 05/11/2021

% So this function is really just an interface to GAMINI
vec=gamini(invec(:),1:length(invec(:)));


