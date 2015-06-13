function XY=regselect(regn,c11,cmn,xunt,res,buf,ofs)
% XY=REGSELECT(regn,c11,cmn,xunt,res,buf,ofs)
%
% Returns coordinates of a certain named specified region
% 
% INPUT:
%
% regn     The region name string, e.g. 'greenland'
% c11      The (lon,lat) coordinates of the top left corner
%          or a [lon(:)] set of coordinates for a region
% cmn      The (lon,lat) coordinates of the bottom right corner
%          or a [lat(:)] set of coordinates for a region
% xunt     The particuluar indices required for this region
% res      0 The standard, default values
%          N Splined values at N times the resolution
% buf      Distance in degrees that the region outline will be enlarged
%          by BUFFERM, not necessarily integer, possibly negative
%          [default: 0]
% ofs      Offset(s) to be considered by PLOTCONT [default: none]
%
% OUTPUT:
%
% XY       The requested coordinates
% 
% Last modified by charig-at-princeton.edu, 05/14/2015
% Last modified by fjsimons-at-alum.mit.edu, 06/13/2015

% The directory where you keep the coordinates
whereitsat=fullfile(getenv('IFILES'),'COASTS');

% Capitalize
Regn=[upper(regn(1)) regn(2:end)];

% Offsets, resolution, and buffering
defval('ofs',0)
defval('res',0)
defval('buf',0)

% Curve speed... when doing interactive mode
defval('spd',0.1)

% Revert to original name if unbuffered
if res==0 && buf==0
  fnpl=fullfile(whereitsat,sprintf('%s.mat',Regn));
elseif buf==0;
  fnpl=fullfile(whereitsat,sprintf('%s-%i.mat',Regn,res));
elseif buf~=0
  fnpl=fullfile(whereitsat,sprintf('%s-%i-%g.mat',Regn,res,buf));
end

% If you already have a file
if exist(fnpl,'file')==2
  load(fnpl)
else
  % You are about to make a file
  if res==0
    % First part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clf
    % This could be 1x2, 2x1 or 2x2 or 2x3 indeed
    if length(c11)<=3 && length(cmn)<=3
      XY=[];
      for index=1:size(c11,1)
	[~,handl,XY2]=plotcont(c11(index,:),cmn(index,:),[],ofs(index));
	% Get rid of common NaNs
	XY2=XY2(~isnan(XY2(:,1)) & ~isnan(XY2(:,2)),:);
	delete(handl)
	XY=[XY ; XY2]; clear XY2
      end
    else
      XY(:,1)=c11;
      XY(:,2)=cmn;
      c11=[min(XY(:,1)) max(XY(:,2))];
      cmn=[max(XY(:,1)) min(XY(:,2))];
    end
        
    % Get rid of common NaNs
    XY=XY(~isnan(XY(:,1)) & ~isnan(XY(:,2)),:);
    
    % If this is antarctica, we want to rotate to the equator
    if strcmp(regn,'antarctica')
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
        % Here, lonc and latc are the same values given out from antarctica
        % for the rotation back to the pole.
    end

    % Now make sure the distances aren't huge
    [X,Y,~,p]=penlift(XY(:,1),XY(:,2),3);
    XY=[X Y];

    % Some handiwork specifically for antarctica
    if strcmp(regn,'antarctica')
        XY=[flipud(XY(p+2:p+3,:)) ;  XY([[(1:p(1))]' ; [p(1)+5:size(XY,1)]'],:)];
        xx=XY(:,1); yy=XY(:,2); 
        d=sqrt((xx(2:end)-xx(1:end-1)).^2+(yy(2:end)-yy(1:end-1)).^2);
        XY=XY(find(d>1e-14),:);
        XY=[XY(2:227,:) ; XY(232:end,:) ; XY(2,:)];
    end
    
    % Experiment with cutoff -> See "check this out"
    defval('xunt',1:length(XY))
    XY=XY(xunt,:);

    % Definitely close the contour again
    XY=XY(~isnan(XY(:,1)) & ~isnan(XY(:,2)),:);
    if XY(end,:) ~= XY(1,:)
      XY=[XY ; XY(1,:)];
    end

    % For Ellesmere Island we want a three part polygon
    if strcmp(regn,'ellesmere') 
        % Undo what we just did matching the end to the start
        % and insert NaNs in order to make this a 3 parts
        XY=XY(1:end-1,:);
        XY=insert(XY,[NaN NaN NaN NaN],[24 150 209 335]); 
        XY=reshape(XY,[],2); 
    end

    % And definitely make this go clockwise
    [X2,Y2]=poly2cw(XY(:,1),XY(:,2));
    XY=[X2 Y2]; clear X2 Y2
   
    plot(XY(:,1),XY(:,2),'LineW',2,'Color','k');
    axis equal 

    % Check this out %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on ; curvecheck(XY(:,1),XY(:,2),spd); hold off
  else
    XY=bezier(eval(sprintf('%s(0)',regn)),res);
  end
  % Change the domain into a slightly larger area
  % Do we want to buffer out or buffer inside the curve?
  if buf ~= 0
    if buf > 0
      inout='out';
    else
      inout='in';
    end
    % Make some buffered coordinates and save them for later use
    disp('Buffering the coastlines... this may take a while');

    % You might look into REDUCEM to make this easier
    % Note that BUFFERM has gone through quite a few revisions
    % The cell output is no longer supported these days
    [LatB,LonB]=bufferm(XY(:,2),XY(:,1),abs(buf),inout);
    % XY=[LonB{2}+180*[LonB{2}<0] LatB{2}];
    
    % Note that, if due to BEZIER there might be a pinched-off loop in
    % the XY you will get an extra NaN and will need to watch it
    % If this shouldn't work, keep it all unchanged in other words
    try
      % You'll need a line for every possible version behavior
      % Note how POLY2CW has disappeared from BUFFERM
      if sign(buf)<0 || ~isempty(strfind(version,'2010a'))
	% Take the last bit of non-NaNs; there might have been pinched
        % off sections in the original
	LonB=LonB(indeks(find(isnan(LonB)),'end')+1:end);
	LatB=LatB(indeks(find(isnan(LatB)),'end')+1:end);
      elseif ~isempty(strfind(version,'2011a')) || ~isempty(strfind(version,'2012a'))
	LonB=LonB(1:find(isnan(LonB))-1);
	LatB=LatB(1:find(isnan(LatB))-1);
      end
    catch
      disp('BUFFERM failed to buffer as expected. Version update?')
    end
    % Periodize the right way
    XY=[LonB+360*any(LonB<0) LatB];
    
    % Definitely get rid of the NaNs again? Should be overkill at this point
    %XY=XY(~isnan(XY(:,1)) & ~isnan(XY(:,2)),:);
  end

  % Save the file
  save(fnpl,'XY')
end
