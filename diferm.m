function diferm(a,b,c)
% DIFERM(a,b,c)
% 
% The mute form of DIFER, i.e.
% difer(a,[],[],NaN)
% difer(a-b,[],[],NaN)
% difer(a-b,c,[],NaN)

defval('b',0)
defval('c',[])
difer(a-b,c,[],NaN)
