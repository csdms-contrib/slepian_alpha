function a=minmax(matrix)
% a=MINMAX(matrix)
%
% Returns the limits of a matrix
% See also RANGE

a=[min(matrix(:)) max(matrix(:))];
