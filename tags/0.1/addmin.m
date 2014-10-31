function vec=addmin(invec)
% vec=ADDMIN(invec)
%
% Repeats the entries of an input vector INDEX times
% where INDEX is the array index, and constructs a new vector with these
% repeated array elements (these could be numbers or strings).
% Useful for real spherical harmonics. 
%
% Last modified by fjsimons-at-alum.mit.edu, Jan 27th, 2003

invec=invec(:);

L=size(invec,1)-1;

el=1:L+1;

vec=gamini(invec,el);

vec=vec(:);

