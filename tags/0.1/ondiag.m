function d=ondiag(m)
% d=ONDIAG(m)
% d=ONDIAG(matrix)
%
% Returns logical indices to the diagonal elements of a matrix of m rows
% by m columns.
%
% By fjsimons-at-alum.mit.edu, December 16th, 2003

if length(m)==1 & ~isnan(m)
  d=~~diag(ones(1,m));
else  
  if size(m,1)~=size(m,2)
    error('Only for square matrices')
  end
  d=~~diag(ones(1,size(m,1)));
end





