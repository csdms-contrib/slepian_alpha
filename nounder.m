function sin=nounder(sin,selse)
% sin=NOUNDER(sin,selse)
%
% Changes underscores to dashes (or selse) in a string
%
% Last modified by fjsimons-at-alum.mit.edu, 03/03/2016

defval('selse','-')

if length(selse)>1 && strcmp(selse,'\_')
  sin=insert(sin,92,find(abs(sin)==95));
else
  sin(find(abs(sin)==95))=selse;
end

