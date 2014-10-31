function letmat=letter(index,kees)
% letmat=LETTER(index,kees)    
%
% This last option returns upper case if 1
%
% Gives letters of the alphabet in a vector of strings
%
% Last modified by fjsimons-at-alum.mit.edu, 09/18/2012

allc=97:122;
allup=65:90;

defval('kees',0)

switch kees
 case 1
  letmat=str2mat(allup(index));
 case 0
  letmat=str2mat(allc(index));
end
