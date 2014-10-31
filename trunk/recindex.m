function [lineara,linearb]=recindex(a,b)
% [lineara,linearb]=RECINDEX(a,b)
%
% Returns the mapping of 'b' into the elements of 'a'
% such that 'b(linearb) == a(lineara)' -
% 'a' and 'b' may not both contain zeroes.
% Usually, 'a' is the bigger set and 'b' is a possible subset,
% 'a' may have some zeroes provided 'b' doesn't
% Compare with the function 'ismember' and 'intersect'
% 'a' must be a unique set
%
% Small matrices only! Needs max(length(a),length(b))^2 locations.
% Round off floats 'cos they'll get squared and then rooted
%
% Last modified by fjsimons-at-alum.mit.edu, 06/18/2007

% 'a' and 'b' may not contain zeroes.
% But if b does not contain any then it's ok
if ~prod(a) & ~prod(b)
  error('Input matrices must not both contain zero elements')
end

% Map 'a' and 'b' into a row vector
a=a(:)';
b=b(:)';

% Matrix 'a' must contain non-repeating elements
if any(size(unique(a))~=size(a))
 error('Matrix a must be a unique set') 
end
 
% Make 'a' the greater set of the two
if length(a)<length(b)
  a=[a repmat(NaN,1,length(b)-length(a))];
end

% Find the indices 'lineara' of elements in 'a' that correspond to 
% elements in 'b' and give their indices 'linearb'
helpmat=a'*b==repmat(a'.^2,1,length(b));
[lineara,j]=find(helpmat);

lineara=lineara';				
% It's ok after all if a has a number of zeros but b doesn't
if ~prod(a)
  if ~~prod(b)
    lineara=lineara(a(lineara)~=0);
%    lineara=lineara(lineara~=1);
    linearb=find(sum(helpmat,1)-sum(~a));
  end
else
  linearb=find(sum(helpmat,1));
end

% Check we've done the right thing...
if ~isempty(linearb) | ~isempty(lineara)
  if ~all(b(linearb) == a(lineara))
    error('Something wrong here')
  end
end
