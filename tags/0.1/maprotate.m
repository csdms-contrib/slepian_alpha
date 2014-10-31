function [dd,lola]=maprotate(d,c11cmn)
% [dd,lola]=MAPROTATE(d,c11cmn)
%
% Takes a world map and turns it around, e.g. when it crosses the
% Greenwich meridian with negative longitudes. If the map is empty, will
% at least do this for the continental outlines of the world.
%
% INPUT:
%
% d       A world map
% c11cmn  The map corners [defaulted]
%
% OUTPUT:
%
% dd      The new map
% lola    The new continental outlines
%
% Last modified by fjsimons-at-alum.mit.edu, 01/15/2010

defval('c11cmn',[-169 90 191 -90])

% Figure out where the map starts and ends
dlon=360/(size(d,2)-1);
[a,i]=min(abs([0:dlon:360]-c11cmn(3)));

% Get rid of the redundant last column
dd=d(:,1:end-1);

% Put in the redundancy again
dd=[d(:,i:end) d(:,1:i)];

% Note the Greenwich trick
yesorno=360*[c11cmn(1)<0 0];
[ax,f,lola1]=plotcont(c11cmn(1:2)+yesorno,...
		      c11cmn(3:4)+yesorno,[],-360);
delete(f)
[ax,f,lola2]=plotcont([0 c11cmn(2)],c11cmn(3:4));
delete(f)

lola=[lola1 ; NaN NaN ; lola2];
