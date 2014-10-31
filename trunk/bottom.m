function bottom(hand,aha)
% BOTTOM(hand)
% BOTTOM(hand,aha)
%
% Moves an object 'hand', child to gcf (default) or axis 'aha'
% to the bottom of the visual pile.
%
% SEE ALSO: TOP
%
% Last modified by fjsimons-at-alum.mit.edu, 03/19/2012

defval('aha',gcf)

kids=getkids(aha);

try
  set(aha,'Children',[kids(kids~=hand) ; kids(kids==hand)]);
catch
  [~,~,ikids]=intersect(hand,kids);
  set(aha,'Children',[kids(~ismember(kids,hand)) ; kids(ikids)]);
end


