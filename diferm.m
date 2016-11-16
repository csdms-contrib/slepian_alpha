function diferm(a,b,tolex)
% DIFERM(a,b,tolex)
% 
% The mute form of DIFER, i.e.
% difer(a,[],[],NaN)
% difer(a-b,[],[],NaN)
% difer(a-b,tolex,[],NaN)
%
% Last modified by fjsimons-at-alum.mit.edu, 11/15/2016

defval('b',0)
defval('tolex',[])
difer(a-b,tolex,[],NaN)
