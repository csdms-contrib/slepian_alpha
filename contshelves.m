function varargout=contshelves(res,rotb)
% [XY,lonc,latc]=CONTSHELVES(res,rotb)
% CONTSHELVES(...) % Only makes a plot
%
% Finds the coordinates of the continental lithosphere and rotates those
% away from the poles so KERNELC doesn't choke on the calculation.
%
% INPUT:
%
% res      0 The standard, default values
%          N Splined values at N times the resolution
% rotb     0 The rotated version [default] needed by e.g. KERNELCP
%          1 The unrotated version suitable for geographic plots 
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the continental shelves, two
%          closed curves separated by a row of NaNs.
%
% Written by Jarno Saarimaki and Ciaran Beggan, 2011
% Last modified by fjsimons-at-alum.mit.edu, 07/12/2012

defval('res',0)
defval('rotb',0)

% The directory where you keep the coordinates
whereitsat=fullfile(getenv('IFILES'),'COASTS');

if res==0
  fnpl=fullfile(whereitsat,'ContShelves.mat');
else
  fnpl=fullfile(whereitsat,sprintf('ContShelves-%i.mat',res));
end

if rotb==1
  XY=load(fullfile(whereitsat,'ContShelves.txt'));
  XY(:,1)=XY(:,1)+[XY(:,1)<0]*360;
  [X,Y]=penlift(XY(:,1),XY(:,2));
  XY=[X Y];
end

if exist(fnpl,'file')==2 
  if rotb==0; load(fnpl); end
  if nargout==0
    plot(XY(:,1),XY(:,2),'k-'); axis equal; axis tight;
  else
    varns={XY,lonc,latc};
    varargout= varns(1:nargout);
  end
else
  if res == 0
    % There are two continuous regions in the data separated by a row
    % of NaNs. The other region is defined in two (by NaNs) so that the
    % original XY could be plotted as continuous lines.
    XY=load(fullfile(whereitsat,'ContShelves.txt'));
    
    % Get rid of NaN between North America and asia and leave the one
    % between North America and Antarctica
    h = find(isnan(XY(:,1)));
    XY(h(1),:) = [];
    XY = [ XY(1:h(2)-2,:) ; XY(1,:) ; XY(h(2)-1:end,:)];
    XY(end+1,:) = XY(h(2)+1,:);
    
    % Convert to Cartesian coordinates
    [X,Y,Z]=sph2cart(XY(:,1)*pi/180,XY(:,2)*pi/180,1);

    % The rotation angles
    lonc = -220; latc = -130;

    % Apply the rotation
    xyzp=[rotz(lonc*pi/180)*roty(-latc*pi/180)*[X(:) Y(:) Z(:)]']';

    % Transform back to spherical coordinates
    [phi,piminth,r]=cart2sph(xyzp(:,1),xyzp(:,2),xyzp(:,3));
    lon=phi*180/pi; lat=piminth*180/pi;

    % Output in the usual format
    XY=[lon lat];
    
    % Eyeball
    plot(XY(:,1),XY(:,2),'LineW',2,'Color','k');
    
    axis equal
    axis tight
    
    save(fnpl,'XY','lonc','latc')
  else
    [XY,lonc,latc]=contshelves(0);
    % Interpolate separately each region that is separated by NaNs
    h = find(isnan(XY(:,1)));
    XYb = [ bezier( XY(1:h(1)-1,:) , res ) ; [NaN NaN] ];
    if length(h) > 1
      for i = 2:length(h)
	XYb = [ XYb ; bezier( XY(h(i-1)+1:h(i)-1,:) , res ) ; ...
		[NaN NaN] ];
      end
    end
    XY = [ XYb ; bezier( XY(h(end)+1:end,:) , res ) ];
    save(fnpl,'XY','lonc','latc')
  end
  varns={XY,lonc,latc};
  varargout= varns(1:nargout);
end
