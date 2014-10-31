function varargout=antarctica(res)
% [XY,lonc,latc]=antarctica(res)
% ANTARCTICA(...) % Only makes a plot
%
% Finds the coordinates of Antarctica, but rotates them to an equatorial
% location so KERNELC doesn't choke on the calculation.
%
% INPUT:
%
% res      0 The standard, default values
%          N Splined values at N times the resolution
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the continent
% lonc     The amount by which you need to rotate it back over z
% latc     The amount by which you need to rotate it back over y
%
% See also PLM2ROT, GEOBOXCAP, KLMLMP2ROT
% 
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012

defval('res',0)

% The directory where you keep the coordinates
whereitsat=fullfile(getenv('IFILES'),'COASTS');

if res==0
  fnpl=fullfile(whereitsat,'Antarctica.mat');
else
  fnpl=fullfile(whereitsat,sprintf('Antarctica-%i.mat',res));
end

if exist(fnpl,'file')==2 
  load(fnpl)
  if nargout==0
    plot(XY(:,1),XY(:,2),'k-'); axis equal; grid on
    axis([-30-lonc 30-lonc -30 20])
  else
    varns={XY,lonc,latc};
    varargout=varns(1:nargout);
  end
else
  if res==0
    % First part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c11=[0 -62];
    cmn=[360 -83];
    clf
    [axlim,handl,XY]=plotcont(c11,cmn,[],0);
    delete(handl)
    % Get rid of common NaNs
    XY=XY(~isnan(XY(:,1)) & ~isnan(XY(:,2)),:);
 
    % Form closed contour
    XY=XY(~isnan(XY(:,1)) & ~isnan(XY(:,2)),:);

    % Find the geographical center and the area
    [lonc,latc,A]=rcenter([XY(:,1) XY(:,2)]);
    
    % Convert to Cartesian coordinates
    [X,Y,Z]=sph2cart(XY(:,1)*pi/180,XY(:,2)*pi/180,1);
    [Xc,Yc,Zc]=sph2cart(lonc*pi/180,latc*pi/180,1);

    % Just make it not straddle the equator; flip the sign; stay close to
    % the value calculated by RCENTER
    lonc=-45;
    
    % Apply the rotation to put it on or near the equator
    xyzp=[rotz(lonc*pi/180)*roty(-latc*pi/180)*[X(:) Y(:) Z(:)]']';
    xyzc=[rotz(lonc*pi/180)*roty(-latc*pi/180)*[Xc   Yc   Zc  ]']';
    % See LOCALIZATION and KLMLMP2ROT for the counterrotation
    
    % Transform back to spherical coordinates
    [phi,piminth,r]=cart2sph(xyzp(:,1),xyzp(:,2),xyzp(:,3));
    lon=phi*180/pi; lat=piminth*180/pi;
    [phic,piminthc]=cart2sph(xyzc(1),xyzc(2),xyzc(3));
    loncp=phic*180/pi; latcp=piminthc*180/pi;

    % Output in the usual format
    XY=[lon lat];

    % Now make sure the distances aren't huge
    [X,Y,~,p]=penlift(XY(:,1),XY(:,2),3);
    XY=[X Y];
    
    % A bit of handiwork here... inspect the bezier work also
    XY=[flipud(XY(p+2:p+3,:)) ;  XY([[(1:p(1))]' ; [p(1)+5:size(XY,1)]'],:)];
    xx=XY(:,1); yy=XY(:,2); 
    d=sqrt((xx(2:end)-xx(1:end-1)).^2+(yy(2:end)-yy(1:end-1)).^2);
    XY=XY(find(d>1e-14),:);
    XY=[XY(2:227,:) ; XY(232:end,:) ; XY(2,:)];
    
    % Form closed contour if you haven't already
    XY=XY(~isnan(XY(:,1)) & ~isnan(XY(:,2)),:);

    % Eyeball
    plot(XY(:,1),XY(:,2),'LineW',2,'Color','k');
    hold on
    plot(loncp,latcp,'bo')
    axis equal 
    axis([-30-lonc 30-lonc -30 20])
    
    % Try this also: 
    % [X,Y,Z]=sph2cart(XY(:,1)*pi/180,XY(:,2)*pi/180,1);
    % plot3(X,Y,Z)

    % Check this out %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    curvecheck(XY(:,1),XY(:,2)) 
        
    save(fnpl,'XY','lonc','latc')
  else
    [XY,lonc,latc]=antarctica(0);
    XYb=bezier(XY,res);
    XY=XYb;
    save(fnpl,'XY','lonc','latc')
  end
  varns={XY,lonc,latc};
  varargout=varns(1:nargout);
end
