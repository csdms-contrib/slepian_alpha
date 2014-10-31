function [r,J]=legendreroots(N);
% [r,J]=LEGENDREROOTS(N);
%
% The function r = legendreroots(N) computes the roots of the 
% Legendre polynomial of degree N.
% Also returns the Jacobi matrix whose eigenvalues they are.
%
% J.A.C. Weideman, S.C. Reddy 1998.
%
% fjsimons-at-alum.mit.edu is still wondering if you can use the
% Jacobi matrix to back out the coefficients of the Polynomials.

n = [1:N-1];                   %  Indices
d = n./sqrt(4*n.^2-1);         %  Create subdiagonals
J = diag(d,1)+diag(d,-1);      %  Create Jacobi matrix
r = sort(eig(sparse(J)));      %  Compute eigenvalues
