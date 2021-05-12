function varargout=arrow(X,Y,U,V,ori,sca,rot,tran)
% [arro,arroperp,cend]=ARROW(X,Y,U,V,ori,sca,rot,tran)
%
% Like QUIVER, but with better arrow heads, for horizontal or vertical
% arrows. Severely different behaviors between MATLAB versions!
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
% OUTPUT:
%
% arro      An object, or a handle, to the arrow
% arroperp  Coordinates defining the perpendicular to the stalk
% cend      Some idea of the extent of this picture
%
% EXAMPLE:
%
% arrow('demo1')
%
% SEE ALSO:
%
% GETX, SETX
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)
% Last modified by fjsimons-at-alum.mit.edu, 06/22/2016

defval('X',0)

if ~isstr(X)
  defval('Y',0)
  defval('U',10)
  defval('V',0)
  defval('ori','h')
  defval('sca',1)
  defval('rot',0)
  defval('tran',[0 0])
  defval('arroperp',NaN);
      
  if verLessThan('matlab','7')
    % What comes out is a numeric HANDLE
    arro=quiver(X(:),Y(:),U(:),V(:));
  elseif verLessThan('matlab','9')
    warning off
    arro=quiver('v6',X(:),Y(:),U(:),V(:));
    warning on
  else
    % What comes out is an OBJECT HANDLE
    arro=quiver(X(:),Y(:),U(:),V(:));
  end

  % Rotate in the plane around the requested angle
  if rot~=0
    rots=rotz(rot*pi/180);
    rots=rots(1:2,1:2);
  else
    rots=1;
  end
  % Rotate in the plane around 90 degrees
  perp=rotz(pi/2);
  perp=perp(1:2,1:2);
  
  % Old version
  if verLessThan('matlab','9')
    % Properties of the head
    g=get(arro(2),'XData');
    h=get(arro(2),'YData');

    % Make this ready for vectorial processing... not phenomenal
    if prod(X)==0 && strcmp(ori,'h')
      set(arro(2),'YData',h/sca)
      set(arro(2),'XData',g+[1 0 1 0]*(g(2)-g(1))*(sca-1)/sca)
    elseif prod(X)==0 && strcmp(ori,'v')
      set(arro(2),'YData',h+[1 0 1 0]*(h(2)-h(1))*(sca-1)/sca)
      set(arro(2),'XData',g/sca)
    end

    % Properties of the stalk and the head
    XY1=rots*[get(arro(1),'XData') ; get(arro(1),'YData')];
    XY2=rots*[get(arro(2),'XData') ; get(arro(2),'YData')];

    % Figure out the perpendicular
    arroperp=perp*XY1;

    if ~all(tran==[0 0]) || rot~=0
      set(arro(1),'XData',XY1(1,:)+tran(1))
      set(arro(1),'YData',XY1(2,:)+tran(2))
      set(arro(2),'XData',XY2(1,:)+tran(1))  
      set(arro(2),'YData',XY2(2,:)+tran(2))
    end
    
    % Get the exent
    cend=[max(get(arro(1),'XData')) max(get(arro(1),'YData'))];
  else    
    % New version
    % Scale
    arro.AutoScaleFactor=sca;
    % Rotate the stalk
    XY1=rots*[arro.UData ; arro.VData];
    arro.UData=XY1(1,:);
    arro.VData=XY1(2,:);
    % Translate
    arro.XData=arro.XData+tran(1);
    arro.YData=arro.YData+tran(2);
    cend=[max(arro.XData) max(arro.YData)];
  end
    
  % Prepare optional output
  varns={arro,arroperp,cend};
  varargout=varns(1:nargout);

elseif strcmp(X,'demo1')
  arrow  
end
