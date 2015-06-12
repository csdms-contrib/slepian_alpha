function varargout=circ(radius,sectorwidth,origin,np,varargin)
% [H,cords,theta]=CIRC(radius,sectorwidth,origin,np,'Value','Property')
%
% Plots circles and returns handles to it.
%
% INPUT:
%
% radius          Radius, can be a vector
% sectorwidth     In radians, makes it plot only a sector [default: 2pi]
%                 If two numbers, from this to that angle.
% origin          Coordinates of origin [x0 y0] [default: 0,0]
% np              Number of points [default: 500]
% 'Value', 'Property' ... to spruce it up
%
% OUTPUT:
%
% H               Plot handles to the circle and its origin
% cords           Coordinates plotted
% theta           Angles plotted, in radians
%
% Last modified by fjsimons-at-alum.mit.edu, 12/01/2013

defval('radius',1)
defval('sectorwidth',2*pi)
defval('origin',[0 0])
defval('np',500)

if nargin<=1 
  theta=linspace(0,2*pi,np);
else
  if length(sectorwidth)==2
    theta=linspace(sectorwidth(1),sectorwidth(2),np);
  else
    theta=linspace(pi/2-sectorwidth/2,pi/2+sectorwidth/2,np);
  end
end

x=(radius'*cos(theta))';
y=(radius'*sin(theta))';

x=x+origin(1);
y=y+origin(2);
hands=plot(x,y,'Color','k');

hold on

if nargin>2
  hands=[hands ;  plot(origin(1),origin(2),'k+')];
end

cords=[x y];

varns={hands,cords,theta};
varargout=varns(1:nargout);

if nargin>=5
  set(hands,varargin{1},varargin{2})
end
