function sbnums=subnum(m,n)
% sbnums=SUBNUM(m,n)
%
% Returns the complete set of mxn figure
% subplot numbers. Will return a vector
% of numbers or a matrix of strings, 
% both of which can go into 'krijetem'
%
% See also: KRIJETEM
%
% Last modified by fjsimons-at-alum.mit.edu, June 3rd, 2004

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



