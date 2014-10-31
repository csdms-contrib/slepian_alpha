function T=tridiag(a,b,c)
% T=TRIDIAG(a,b,c)
%
% Makes a tridiagonal matrix
%
% INPUT:
%
% a    The subdiagonal
% b    The main diagonal
% c    The superdiagonal
%
% OUTPUT:
%
% T    The tridiagonal matrix
%
% Last modified by fjsimons-at-alum.mit.edu, 03/03/2009

T=diag(b)+diag(a,-1)+diag(c,1);
