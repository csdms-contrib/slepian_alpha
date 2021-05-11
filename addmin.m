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

% Ensure column vector input
invec=invec(:);

% Think - maximum spherical harmonic degree
L=size(invec,1)-1;

% So this function is really just an interface to GAMINI
vec=gamini(invec,1:L+1);


