function addo=adrc(mat,val)
% addo=ADRC(mat,val)
%
% Adds an extra row and column to a matrix, e.g. for use wih PCOLOR or to
% use in optimization with a Lagrange multitplier.
%
% INPUT:
%
% mat      The original matrix
% val      The scalar value of the additional row and column [default: 0]
%
% OUTPUT:
%
% addo     The new matrix with the row and column added
%
% Last modified by fjsimons-at-alum.mit.edu, 09/18/2006

% Note that for Lagrange optimization, you may want addo(end,end)=0...

defval('val',0);

[m,n]=size(mat);
addo=[mat  repmat(val,m,1) ; repmat(val,1,n+1)];
