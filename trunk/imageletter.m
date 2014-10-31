function varargout=imageletter(w,n1,n2)
% v=IMAGELETTER(w,n1,n2)
%
% Makes a n1xn2 bitmap of a certain character string.
% Note that the outcome here depends on your screen resolution.
% The output is an image, if you interpret this as coordinates you may
% have to flip it upside down!
%
% Last modified by fjsimons-at-alum.mit.edu, 07/15/2009

defval('w',1)
defval('n1',65)
defval('n2',65)

% This fixes the letter conventions
switch w
 case 1; wl='X+'; rota=[0]; axl=[602   679   429   490];
 case 2; wl='Z-'; rota=[0]; axl=[602   679   429   490];
 case 3; wl='Y+'; rota=[0]; axl=[602   679   429   490];
 case 4; wl='X-'; rota=[0]; axl=[602   679   429   490];
 case 5; wl='Z+'; rota=[0]; axl=[602   679   429   490];
 case 6; wl='Y-'; rota=[0]; axl=[602   679   429   490];
end

% See if we've already got it
fnpl=fullfile(getenv('EPS'),sprintf('imageletter_%i.tif',w));

if exist(fnpl)~=2
  clf

  % Make it
  p=text(0.5,0.5,wl,'fonts',80,...
	 'horizontala','l','verticala','t',...
	 'rotation',rota);
  axis off
  fig2print(gcf,'portrait')
  
  figdisp([],w,[],1,'tiff')
end

% Load it
v=imread(fnpl);

%imagesc(v) ; keyboard
% Crop the image
% On the Dell desktop
v=v(axl(3):axl(4),axl(1):axl(2));

% Interpolate
v1i=linspace(1,size(v,1),n1);
v2i=linspace(1,size(v,2),n2);

v=interp2(v,v2i,v1i','nearest');

% Produce output
v=1-double(v)/255;

% Take a look if no output required
if nargout==0
  imagefnan(v)
end

% Produce desired output
varns={v};
varargout=varns(1:nargout);



