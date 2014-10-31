function flipt=flipflop(matrix)
% flipt=FLIPFLOP(matrix)
% flipt=FLIPFLOP(matrix)
%
% Flips and duplicates a matrix suitable for spectral analysis

% Last modified by fjsimons-at-mit.edu, Jan 11, 2000

flipt=[matrix fliplr(matrix)];
flipt=[flipt ; flipud(flipt)]; % Changed this order such that original is
% in upper left quadrant





