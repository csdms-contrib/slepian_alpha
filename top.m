function top(hand,aha)
% TOP(hand)
% TOP(hand,aha)
%
% Moves an object 'hand', child to gcf (default) or axis 'aha'
% to the top of the visual pile.
%
% SEE ALSO: BOTTOM
% 
% Last modified by fjsimons-at-alum.mit.edu 03/19/2012

defval('aha',gcf)

kids=getkids(aha);

try
  set(aha,'Children',[kids(kids==hand); kids(kids~=hand)]);
catch
  [~,~,ikids]=intersect(hand,kids);
  set(aha,'Children',[kids(ikids); kids(~ismember(kids,hand))]);
end

