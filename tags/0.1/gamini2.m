function g=gamini2(mat,fscals)
% g=GAMINI2(mat,fscals)
%
% Folds the elements of a MATRIX fscals(1)Xfscals(2) times
%
% See also GAMINI
%
% sum(gamini2(1:5,[1 2])~=gamini(1:5,2))
%
% Last modified by fjsimons-at-alum.mit.edu, June 8th 2000

jrep=repmat(fscals(2),size(mat,2),1)';
jgelp=zeros(1,sum(jrep));
jgelp([1 cumsum(jrep(1:end-1))+1])=1; 

irep=repmat(fscals(1),size(mat,1),1)';
igelp=zeros(1,sum(irep));
igelp([1 cumsum(irep(1:end-1))+1])=1; 

g=mat(cumsum(igelp),cumsum(jgelp));
