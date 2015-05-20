function C=cellnan(J,M,N)
% C=CELLNAN(J,M,N)
%
% Initializes a cell array with nans
%
% INPUT:
%
% J      Number of COLUMNS (if J a scalar), or ROWS and COLUMNS [J(1) J(2)]
% M,N    The number of rows and columns of each element of this cell structure;
%        if these are vectors, they will be applied to every element in
%        each of the ROWs. Check the behavior, the dimensions are
%        somewhat nonintuitively expanded. 
%
% OUTPUT:
%
% C      The initialized cell array
%
% EXAMPLE:
%
% cellnan(3,2,1)
% cellnan([1 3],2,1)
% cellnan([3 1],2,1)
% cellnan([3 2],2,1)
% cellnan(3,[1 2 3],[3 2 1])
% cellnan([3 1],[1 2 3],[3 2 1])
% cellnan([1 4],[1 2 3],[3 2 1])
% cellnan([3 4],[1 2 3],[3 2 1])
% 
% SEE ALSO:
%
% STRUCTEMPTY
%
% Last modified by fjsimons-at-alum.mit.edu, 02/05/2015

% Defaults
defval('J',3)
defval('M',4)
defval('N',5)

% Do it!
if isscalar(J)
  % Now you're making J COLUMNS
  C=cell(1,J);
else
  C=cell(J(1),J(2));
end

if isscalar(M) && isscalar(N)
  % If all of them have the same number of dimensions
  [C{:}]=deal(nan(M,N));
else
  % Now every ROW gets the same initialization and there must be one M
  % and one N for each of the rows. Later on, can extend this.
  for ind=1:J(1)
    [C{ind,:}]=deal(nan(M(ind),N(ind)));
  end
end
