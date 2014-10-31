function [L,C,S]=spr2lcs(A)
% [L,C,S]=SPR2LCS(A)
%
% The reverse of LCS2SPR

[m,n]=size(A);
[L,C,S]=find(A);
[L,j]=sort(L);
S=S(j);
C=C(j);
L=full(cumsum(sparse(L,1,1,m,1)));
