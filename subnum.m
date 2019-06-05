function sbnums=subnum(m,n)
% sbnums=SUBNUM(m,n)
%
% Returns a set of figure panel subnumbers
%
% INPUT:
%
% m          A number of rows across the figure
% n          A number of columns across the figure
%
% OUTPUT:
%
% sbnums     A vector of numbers or a matrix of strings, e.g. 
%
% See also: KRIJETEM
%
% Last modified by fjsimons-at-alum.mit.edu, 06/05/2019

% Running index
indo=1:m*n;

warning off MATLAB:conversionToLogical

if max(indo)<=9
  sbnums=m*100+n*10+indo;
else
   indo=parse([str2mat(num2str(indo)) ' '],' ');
   indo=indo(logical(sum(indo(:,:)~=32,2)),:);
   sbnums=...
       [repmat([num2str(m) ','],m*n,1) repmat([num2str(n) ','],m*n,1) indo];
end

warning on MATLAB:conversionToLogical



