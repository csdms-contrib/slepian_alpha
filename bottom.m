function bottom(hand,aha)
% BOTTOM(hand,aha)
%
% Moves an object, child to a figure or axis, to the BOTTOM of the visual pile.
%
% INPUT:
%
% hand     An object handle
% gcf      A figure or axis handle [default: gcf]
%
% SEE ALSO: TOP, UISTACK
% 
% Last modified by fjsimons-at-alum.mit.edu 08/17/2020

defval('aha',gcf)

kids=getkids(aha);

try
  set(aha,'Children',[kids(kids~=hand) ; kids(kids==hand)]);
catch
  [~,~,ikids]=intersect(hand,kids);
  set(aha,'Children',[kids(~ismember(kids,hand)) ; kids(ikids)]);
end


