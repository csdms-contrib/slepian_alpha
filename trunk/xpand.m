function ax=xpand(axl,num)
% ax=XPAND(axl,num)
%
% Expands a pair-vector, e.g. axis limits, by a percentage of the range
% in each pair of its entries. Currently only for one pair or two pairs.
%
% INPUT:
%
% axl    The old vector (axis limits)
% num    The percentage by which x and y are inflated [default: 10]
%
% OUTPUT:
%
% ax     The new vector (axis limits)
%
% EXAMPLE:
%
% axis(xpand(axis))
% xpand([-4 -5])
%
% SEE ALSO: OPENUP
%
% Last modified by fjsimons-at-alum.mit.edu, 11/24/2010

defval('num',10)

% Twas a percentage!
num=num/100;

% Tis a vector of pairs
axl=axl(:)';
if mod(length(axl),2)
  error('This must be a vector of pairs')
end

xran=axl(2)-axl(1);
if length(axl)==4
  yran=axl(4)-axl(3);
  ax=[axl(1)-xran*num axl(2)+xran*num axl(3)-yran*num axl(4)+yran*num];
else
  ax=[axl(1)-xran*num axl(2)+xran*num];
end
