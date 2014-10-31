function varargout=sphere2cube(lon,lat,alfa,bita,gama);
% [xi,eta,fid]=SPHERE2CUBE(lon,lat,alfa,bita,gama);
% 
% Converts a set of spherical coordinates to the equivalent standard
% coordinates on the faces of the cubed sphere.
%
% INPUT:
%
% lon,lat    Same-length vectors (:) with coordinates [degrees]
% alfa       First Euler angle of wholesale tilt of all tiles [defaulted]
% beta       Second Euler angle of wholesale tilt of all tiles [defaulted]
% gama       Third Euler angle of wholesale tilt of all tiles [defaulted]
% 
% OUTPUT:
%
% xi,eta      6-dimensional cell arays with coordinates [radians] per chunk
% fid         A list of the faces to which the original requests belong
%
% EXAMPLE:
%
% sphere2cube('demo1') % Continents and hot spots
% sphere2cube('demo2') % ... with face recognition
%
% SEE ALSO: BOXCUBE
%
% Last modified by fjsimons-at-alum.mit.edu, 03/21/2010

if ~isstr(lon)
  % Specify defaults
  defval('alfa',[]);
  defval('bita',[]);
  defval('gama',[]);

  % Turn them all into vectors
  lon=lon(:);
  lat=lat(:);
  
  % Convert the requested points to Cartesian coordinates
  [X,Y,Z]=sph2cart(lon*pi/180,lat*pi/180,1);
  coordd=[X(:)' ; Y(:)' ; Z(:)'];

  % Get the rotation matrices for this construction
  [rottot,mats,legs]=cubemats(alfa,bita,gama);

  % Initialize the new coordinate matrices
  [x,y,z]=deal(zeros(length(lon),6));
  
  % In which face are the requested coordinates?
  fid=zeros(size(lon));
  
  % Loop over all faces
  for f=1:6
    % Do the inverse coordinate transform all at once
    stuff=inv(mats{f})*inv(rottot)*coordd;
    % Now distribute over the three three-dimensional vectors
    x(:,f)=stuff(1,:)';
    y(:,f)=stuff(2,:)';
    z(:,f)=stuff(3,:)';
    % Convert back to spherical coordinates
    [phi,piminth,r]=cart2sph(x(:,f),y(:,f),z(:,f));
    lon=phi*180/pi; 
    lat=piminth*180/pi;
    % Check by plotting here
    defval('xver',0)
    if xver==1
      clf
      plot(lon,lat,'.','markers',2); axis image
      set(gca,'xgrid','on','xtick',[-45 45],'ygrid','on','ytick',[-45 45])
      pause
    end
    % Restrict to the frontal x+ chunk longitudes
    xplus=[lon>=-45 & lon<=45];

    if nargout<3
      % We won't need the face number, so shrink arrays
      lon=lon(xplus);
      lat=lat(xplus);
    else
      xkeep=xplus;
    end

    % Convert from longitude and latitude to xi and eta, remember flip
    xi{f}=-lon*pi/180;
    eta{f}=atan(1./tan(pi/2-lat*pi/180)./cos(xi{f}));

    % Restrict to the frontal x+ chunk longitudinally varying latitudes
    xplus=[eta{f}>=-pi/4 & eta{f}<=pi/4];
    if nargout<3
      xi{f}=xi{f}(xplus);
      eta{f}=eta{f}(xplus);
    else
      xi{f}=xi{f}(xplus & xkeep);
      eta{f}=eta{f}(xplus & xkeep);
      % Keep a tab of where they end up
      fid=fid+[xplus & xkeep]*f;
    end
  end
  
  % Prepare output
  varns={xi,eta,fid};
  varargout=varns(1:nargout);
elseif strcmp(lon,'demo1')
  % Load and convert continents
  [a,b,XYZ]=plotcont([],[],3); 
  [phi,piminth,r]=cart2sph(XYZ(:,1),XYZ(:,2),XYZ(:,3));
  lonc=phi*180/pi; latc=piminth*180/pi;
  [xic,etac]=sphere2cube(lonc,latc);
  % Load and convert plates
  [handl,XY]=plotplates([],[],1);
  [xip,etap]=sphere2cube(XY(:,1),XY(:,2));
  % Load and convert hotspots
  hotspots=load(fullfile(getenv('IFILES'),'GEOLOGY','hotspots17'));
  [xih,etah]=sphere2cube(hotspots(:,2),hotspots(:,1));
  % Plot the dudes
  clf
  [pc,pgc]=plotonchunk(xic,etac);
  hold on
  [pp,pgp]=plotonchunk(xip,etap);
  delete(pgc{:})
  hold on
  [ph,pgh]=plotonchunk(xih,etah);
  delete(pgh{:})
  set(ph(~isnan(ph)),'MarkerS',8,'MarkerF','r','MarkerE','r','Marker','o')
  axis image
  fig2print(gcf,'portrait')
  hold off
  figdisp([],1,[],1)
elseif strcmp(lon,'demo2')
  % Load and convert continents
  [a,b,XYZ]=plotcont([],[],3); clf
  [phi,piminth,r]=cart2sph(XYZ(:,1),XYZ(:,2),XYZ(:,3));
  lonc=phi*180/pi; latc=piminth*180/pi;
  % Find the cubed-sphere coordinates and faces of the continents
  [xic,etac,fid]=sphere2cube(lonc,latc);
  % Plot while acknowledging the faces
  lon=lonc/180*pi;
  lat=latc/180*pi;
  rad=repmat(1.001,size(lat));
  [xx,yy,zz]=sph2cart(lon,lat,rad);
  % Six different colors
  cols={'r' 'y' 'k' 'b' 'm' 'g'};
  for ind=1:6
    a(ind)=plot3(xx(fid==ind),yy(fid==ind),zz(fid==ind),'.',...
		 'Color',cols{ind},'MarkerS',5);
    hold on ; axis image
  end
  hold off; axis off; view(-10,30)
  % There may be some fid=0; these should be NaN
  figdisp([],2,[],1)
end
