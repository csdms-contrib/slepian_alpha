function newpos=shrink(handle,widthfactor,heightfactor)
% SHRINK(handle,widthfactor,heightfactor)
% newpos=SHRINK(handle,widthfactor,heightfactor)
%
% This function reduces the width and height of a graphics object
% by some factor. 
%
% INPUT:
%
% handle        The handle(s) to the graphics object (default: gca)
% widthfactor   The factor by which to reduce the width (default: 2)
% heightfactor  The factor by which to reduce the height (default: 2)
%
% OUTPUT:
%
% newpos        Optional output new position (lrtb)
%
% Last modified by fjsimons-at-mit.edu, 14.03.2006

defval('handle',gca)
defval('widthfactor',2)
defval('heightfactor',2)

property= 'Position';
handle=handle(:);
if prod(size(widthfactor))*prod(size(widthfactor))~=1
  error('Cannot deal with this multiple input yet')
end
for index=1:length(handle)  
  oldvalue=get(handle(index),property);
  newvalue=[oldvalue(1)+oldvalue(3)*((1-1/widthfactor)/2)  ...
	    oldvalue(2)+(oldvalue(4)*(1-1/heightfactor))/2 ... 
	    oldvalue(3)/widthfactor ...
	    oldvalue(4)/heightfactor];
  set(handle(index),property,newvalue)
  if nargout>0
    newpos(index,:)=newvalue;
  end
end
