function top(hand,aha)
% TOP(hand,aha)
%
% Moves an object, child to a figure or axis, to the top of the visual pile.
%
% INPUT:
%
% hand     An object handle
% gcf      A figure or axis handle [default: gcf]
%
% SEE ALSO: BOTTOM, UISTACK
% 
% Last modified by fjsimons-at-alum.mit.edu 08/17/2020

defval('aha',gcf)

kids=getkids(aha);

try
  set(aha,'Children',[kids(kids==hand); kids(kids~=hand)]);
catch
  [~,~,ikids]=intersect(hand,kids);
  set(aha,'Children',[kids(ikids); kids(~ismember(kids,hand))]);
end

