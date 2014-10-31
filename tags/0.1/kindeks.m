function x=kindeks(M,i)
% x=KINDEKS(M,i)
%
% Returns the COLUMN i of a matrix M.
% Returns the ith 2-nd dimension of a 3D matrix.
%
% See also INDEKS, RINDEKS.
%
% Last modified by fjsimons-at-alum.mit.edu, 09/28/2006

x=M(:,i,:);

