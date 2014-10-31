function outmat=strcell(incell)
% outmat=STRCELL(incell)
%
% Creates character array from cell array of strings
%
% The reverse of CELLSTR.
%
% Compare to STR2CELL.
%
% Last modified by fjsimons-at-alum.mit.edu, 02/05/03

s=length(incell);
for index=1:s
  S(index)=length(incell{index});
end
MS=max(S);
for index=1:s
  outmat(index,:)=[incell{index} repmat(' ',1,MS-S(index))];
end

