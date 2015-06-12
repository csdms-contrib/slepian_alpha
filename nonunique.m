function [C,IA] = nonunique(A,varargin)
% [C,AI] = NONUNIQUE(A)
% [C,AI] = NONUNIQUE(A,'sorted')
% [C,AI] = NONUNIQUE(A,'stable')
%
% A complement to MATLAB's UNIQUE function. For a given matrix or cell
% array of strings, this function finds its nonunique members, including 
% the first occurence of each nonunique element (which, under a strict
% definition of nonuniqueness, might be exluded), as well as an index
% vector to input "A" such that C = A(IA). This function mirrors UNIQUE
% in terms of sorted/stable output and the output's orientation.
%
% INPUT:
%
% A          Matrix/cell array from which to find nonunique members
% varargin   Control on how the output nonunique values are sorted
%               <empty>: Nonunique values are sorted
%               'sorted': Same as the <empty>
%               'stable': Nonunique values are not sorted and are in the
%                         same order in which they appear in "A"
%
% OUTPUT:
%
% C     Vector/cell array of the nonunique elements of "A"
% AI    Indices of "A" to nonunique values such that A(IA) = C
%
% See also ISMEMBER, SORT, UNIQUE
%
% Last modified by gleggers-at-alumni.princeton.edu, 05/21/2014

% If "A" is a matrix, vectorize it and find nonunique values of the new "A"
if ~any(size(A) == 1)
    A = A(:);
end

% Make an index vector for all the elements of "A"
indA = (1:length(A))';

% Get the unique values of "A" and their linear indices
[uA,iuA] = unique(A);

% Get indices in the overall index vector of unique values, and NOT these
% indices will be the indices of the nonunique values of "A"
inuA = indA(~ismember(indA,iuA));

% Get the strictly defined nonunique values, which excludes the first
% occurence of nonunique values (since at that point they are still unique)
nuA = A(inuA);

% Get indices of unique values also found among the nonunique values
inboth = iuA(ismember(uA,nuA));

% Combine the previous indices of nonunique values with the indices of 
% unique values which are repeated, and use this new index array to find
% the loosely defined nonunique values
IA = sort([inuA(:); inboth(:)]);
C = A(IA);

% To match the UNIQUE function, if the 'sorted' flag or no flag is given,
% then sort the nonunique values in ascending order
if isempty(varargin) || strcmp(varargin{1},'sorted')
    [C,sortI] = sort(C);
    IA = IA(sortI);

% If no 'stable' flag, then pass an error for improper input (nothing needs
% be done because the nonunique values are already not sorted)
elseif ~strcmp(varargin{1},'stable')
    error('Variable input ''%s'' was not of the proper form',varargin{1})
end