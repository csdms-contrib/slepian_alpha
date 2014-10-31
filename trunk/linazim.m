function varargout=linazim(ang,senter,rx,ry,varargin)
% [x,y]=LINAZIM(ang,senter,rx,ry)
% [x,y]=LINAZIM(ang,senter,rx,ry,nxy)
% ph=LINAZIM(ang,senter,rx,ry)
%
% Plots lines of given azimuth at a set of points.
%
% INPUT:
%
% ang              Angles in degrees NORTH from EAST
% senter           [x(:) y(:)] location of center points
% [rx,ry]          Radii of a box in which the arrow fits
%                  If ry=[] then scaled in circle of radius rx
% nxy              Number of points calculated [Default: 100]
% 
% OUPUT:
%
% [x,y]            Coordinates of lines to be plotted
% ph               Handles to the lines plotted
%
%
% EXAMPLE:
%
% ph=linazim([0:10:179],[5 8],5,3); axis image
% ph=linazim([0:10:179],[5 8],5,[]); axis image
% ph=linazim([0:10:179],[rand(18,1) rand(18,1)],0.15,0.15); axis image
% ph=linazim([0:10:179],[rand(18,1) rand(18,1)],0.15,[]); axis image
% ph=linazim([0:10:179],[rand(18,1) rand(18,1)],0.2*rand(18,1),[]); axis image
% ph=linazim([0:10:179],[(1:18)' repmat(0,18,1)],1,[]); axis image
%
%
% See also RUMPKER
%
% Last modified by fjsimons-at-alum.mit.edu, Feb 21st, 2004.

defval('senter',[0 0])
defval('rx',10)
defval('ry',[])

ang=ang(:)*pi/180;
if size(senter,2)~=2
  error('Wrong center specification')
end
if length(ang)>1 & size(senter,1)>1
  if length(ang)~=size(senter,1)
      error('If more than one center, size must match number of azimuths')
  end
end

if length(ang)>1 & size(senter,1)==1
  senters=repmat(senter,length(ang),1);
else
  senters=senter;
end

if nargin==5
  nxy=ceil(varargin{1}/2);
else
  nxy=ceil(100/2);
end

% IF BOX
if ~isempty(ry)
  for index=1:length(ang)
    angm=ang(index);
    senter=senters(index,:);
    if abs(tan(angm))<ry/rx
      % Intersects parallels to Y-axis
      if nargout==1
	p1=[rx rx*tan(angm)];      
      else
	X=linspace(0,rx,nxy);
	Y=X*tan(angm);
      end   
    else
      % Intersects parallels to X-axis
      if nargout==1
	p1=[ry/tan(angm) ry];
      else
	Y=linspace(0,ry,nxy);
	X=Y/tan(angm);
      end
    end  
    
    if nargout==1
      % Draw line through two points, return handle
      p2=-p1;
      p1=senter+p1;
      p2=senter+p2;
      varargout{1}(index)=line([p1(1) p2(1)],[p1(2) p2(2)]);
      set(varargout{1}(index),'Color','k')
    else
      % Return points
      X=senter(1)+[-fliplr(X) X];
      Y=senter(2)+[-fliplr(Y) Y];
      varargout{1}(:,index)=X(:);
      varargout{2}(:,index)=Y(:);
    end
  end
else % IF CIRCLE
  % Then xy is not counted double as in box
if nargout==1
    nxy=2;
  else
    if nargin==5
      nxy=varargin{1};
    else
      nxy=100;
    end
  end
  % Should be the first nonzero one not just the first one  
  if length(rx)==1
    ra=linspace(-rx(1),rx(1),nxy);
    X=ra(:)*cos(ang(:)');
    Y=ra(:)*sin(ang(:)');
  elseif length(rx)==length(ang)
    ra=linspace(-mean(rx),mean(rx),nxy);
    ra=repmat(ra(:)',length(rx),1);
    % Magnifiy even more with this
    skol=1;
    skale=repmat(rx./mean(rx),1,nxy)*skol;
    ra=(ra.*skale)';
    ang=repmat(ang(:)',nxy,1);
    X=ra.*cos(ang);
    Y=ra.*sin(ang);
  else
    error('Wrong scaling specification')
  end
  X=X+repmat(senters(:,1)',nxy,1);
  Y=Y+repmat(senters(:,2)',nxy,1);
  if nargout==1
    varargout{1}=plot(X,Y,'k');
  else
    varargout{1}=X;
    varargout{2}=Y;
  end
end
