function [arro,arroperp,cend]=arrow(X,Y,U,V,ori,sca,rot,tran)
% [arro,arroperp,cend]=ARROW(X,Y,U,V,ori,sca,rot,tran)
%
% Like QUIVER, but with better arrow heads, 
% for horizontal or vertical arrows 
%
% INPUT:
%
% X,Y,U,V   As in QUIVER, default is horizontal line
% ori       Orientation, for arrowhead scaling [default: h]
% sca       Scale factor for the arrowhead [default: 1]
%           This really only works for single arrows along the axes...
% rot       Rotation angle over z, in degrees [default: 0]
% tran      Translation, optional [default: [0 0]]
%
% Last modified by fjsimons-at-alum.mit.edu, 10/15/2007

defval('X',0)
defval('Y',0)
defval('U',10)
defval('V',0)
defval('ori','h')
defval('sca',1)
defval('rot',0)
defval('tran',[0 0])
% Stupid Matlab backwards compatibility
arro=quiver('v6',X(:),Y(:),U(:),V(:));
g=get(arro(2),'XData');
h=get(arro(2),'YData');

% Make this ready for vectorial processing... right now no luck
if strcmp(ori,'h')
  set(arro(2),'YData',h/sca)
  set(arro(2),'XData',g+[1 0 1 0]*(g(2)-g(1))*(sca-1)/sca)
elseif strcmp(ori,'v')
  set(arro(2),'YData',h+[1 0 1 0]*(h(2)-h(1))*(sca-1)/sca)
  set(arro(2),'XData',g/sca)
else
  error('Specify valid option')
end

X1=get(arro(1),'XData');
X2=get(arro(2),'XData');
Y1=get(arro(1),'YData');
Y2=get(arro(2),'YData');

if rot~=0
  rots=rotz(rot*pi/180);
  rots=rots(1:2,1:2);
  XY1=rots*[X1 ; Y1];
  XY2=rots*[X2 ; Y2];
  perp=rotz(pi/2);
  arroperp=perp(1:2,1:2)*XY1;
else
  XY1=[X1 ; Y1];
  XY2=[X2 ; Y2];
end
if ~all(tran==[0 0])
  set(arro(1),'XData',XY1(1,:)+tran(1))
  set(arro(1),'YData',XY1(2,:)+tran(2))
  set(arro(2),'XData',XY2(1,:)+tran(1))  
  set(arro(2),'YData',XY2(2,:)+tran(2))
end

cend=[max(get(arro(1),'XData')) max(get(arro(1),'YData'))];
