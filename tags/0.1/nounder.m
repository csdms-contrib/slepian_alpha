function sin=nounder(sin,selse)
% sin=NOUNDER(sin,selse)
%
% Changes underscores to dashes (or selse) in a string
%
% Last modified by fjsimons-at-alum.mit.edu, 03/18/2013

defval('selse','-')
sin(find(abs(sin)==95))=selse;

