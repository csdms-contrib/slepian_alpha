function [ax,axl,loc]=laxis(ah,xrat,yrat)
% [ax,axl,loc]=LAXIS(ah,xrat,yrat)
%
% Creates an extra axis on top of another one, used by XTRAXIS.
%
% INPUT:
%
% ah      Axis handle (by default, the current one)
% xrat    Fraction of the WIDTH of the original axis by which the new
%         axis is displaced to the LEFT wrt to the original [default 1/4]
% yrat    Fraction of the HEIGHT of the original axis by which the new
%         axis is displaced UP wrt to the original one [default 1/4]
% 
% OUTPUT:
%
% ax      New axis handle
% axl     New axis limits
% loc     And adequate location of a text annotation to the orginal axis
%
% SEE ALSO: XTRAXIS, XTRAXIS1D, XTRAXIS2D 
% 
% Last modified by fjsimons-at-alum.mit.edu, 08/02/2012

defval('ah',gca)
defval('xrat',1/4);
defval('yrat',1/4);

axes(ah)
axl=axis;
popo=getpos(ah);
ax=axes('position',...
	[popo(1)-xrat*popo(3) popo(2)+yrat*popo(4) popo(3) popo(4)]);
axis(axl)
axis equal off
loc=[axl(1) axl(3)+1/2*(axl(4)-axl(3))];
