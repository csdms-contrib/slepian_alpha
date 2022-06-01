function varargout=plotcont(c11,cmn,res,ofs,pcol,lolax)
% [axlim,handl,XYZ,xyze]=PLOTCONT(c11,cmn,res,ofs,pcol,lolax)
% 
% Plots continental outlines (from a database) on a map
% 
% INPUT:
%
% c11    [x,y]/[lon,lat] of upper left node [defaulted]
% cmn    [x,y]/[lon,lat] of bottom right node [defaulted]
% res    0 regular coastlines (default) 
%        1 (slightly) higher-resolution coastlines and lakes, too
%        2 Global Mollweide projection centered on Pacific
%        3 Three-dimensional coordinates, no great lakes (to return)
%        4 Three-dimensional coordinates, Northern hermisphere only
%        5 High-resolution coast line of California (including Baja, on WGS84)
%        6 Simple restriction to Australia
%        7 The world's continents, opaquely as patches 
%        8 The world's coastlines from the NOAA>NESDIS>NGDC>MGGD database
%        9 Plot this on the two-dimensional cubed sphere
%        10 Global Mollweide projection centered on Greenwich
%        11 Three-dimensional flattened obstructed view centered on an axis
%           [defaulted: looking down upon Guyot Hall]
% ofs    Longitude offset, e.g. 360 degrees [only for 0,1,5,6,7]
% pcol   The patch color in case option 7 is chosen [default: grey]
% lolax  Longitude and latitude of the view axis when res==11
%
% Longitude ranges from 0 to 360; plots are centered on the PACIFIC.
%
% OUTPUT:
%
% axlim  Axis that fits snugly around the data.
% handl  Handle to the plotted line or individual patches, and for
%        option Mollweide, also the handle to the box around it
% XYZ    The actual data points plotted (2D or 3D)
% xyze   The equatorial data points plotted (2D or 3D)
% az,el  Azimuth and elevation of the VIEW axis 
%
% AUSTRALIA:
% plotcont([90 10],[180 -60])
%
% WORLD:
% plotcont([0 90],[360 -90])
%
% SEE ALSO: MAPROTATE, SPHAREA, PHICURVE, RCENTER
%
% Last modified by fjsimons-at-alum.mit.edu, 06/01/2022

% Saved matrix as space-saving unsigned integer 
% - but that translates the NaN's into some  high number - take that out.
% Note how A==NaN does not return the NaN's!
% You'll have to maintain your own databases of continental outlines!

% Define default values
defval('res',0)
defval('ofs',0)
defval('pcol',grey)
defval('xyze',nan(1,3))
% Guyot Hall see guyotphysics(0)
defval('lolax',[-74.6548 40.3458])

switch res
 case 5
  % Overrride map for California
  defval('c11',[235.4320 43]);
  defval('cmn',[252.9990 22.8699]);
 case 6
  defval('c11',[90 10])
  defval('cmn',[180 -60])
  res=1;
 otherwise
  defval('c11',[0 90])
  defval('cmn',[360 -90])
end

% Where are the data kept? Create a directory $IFILES/COASTS
% with $IFILES set as an environment variable...
defval('ddir',fullfile(getenv('IFILES'),'COASTS'))

% Load the data sets
switch res
 case {0,7,9}
  fid=fopen(fullfile(ddir,'cont.mtl'),'r','b');
  cont=fread(fid,[5217 2],'uint16');
  fclose(fid);
 case {1,3,4,11}
  fid=fopen(fullfile(ddir,'cost.mtl'),'r','b');  
  cont=fread(fid,[9598 2],'uint16');
  fclose(fid);
 case 2
  load(fullfile(ddir,'conm'))
  hold on
  handl{1}=plot(conxm,conym,'k');
  axlim=[];
  handl{2}=plot(xbox,ybox,'k');
  XYZ=[conxm conym];
 case 10
  load(fullfile(ddir,'conmg'))
  hold on
  handl{1}=plot(conxm,conym,'k');
  axlim=[];
  handl{2}=plot(xbox,ybox,'k');
  XYZ=[conxm conym];
 case 5 
  cont=load(fullfile(ddir,'California'));
  cont=cont.cont;
 case 8
  cont=load(fullfile(ddir,'1to5M'));  
  cont=cont.cont;
  % Don't do this or else put in appropriate nans... see below for an
  % automatic way of doing this
  cont(:,1)=cont(:,1)+360*[cont(:,1)<0];
end

% Recast data in good form
switch res
  case {0,1,3,4,7,9,11}
   cont=cont/100-90;
   cont(cont==max(max(cont)))=NaN;
end

switch res
 case {0,1,5,8}
  % Restrict to the map requested
  lon=cont(:,1); lat=cont(:,2);
  lon(~(lon>=c11(1) & lon<=cmn(1)))=NaN;
  lat(~(lat<=c11(2) & lat>=cmn(2)))=NaN;
  % Plot these continental outlines
  hold on
  handl=plot(lon+ofs,lat,'k','Linewidth',1);
  axlim=[min(lon([~isnan(lon) & ~isnan(lat)]))...
	 max(lon([~isnan(lon) & ~isnan(lat)]))...
	 min(lat([~isnan(lon) & ~isnan(lat)]))...
	 max(lat([~isnan(lon) & ~isnan(lat)]))];
  XYZ=[lon+ofs lat];
 case 7
  % For the filled patches only
  % Collect the data
  lon=cont(:,1); lat=cont(:,2);
  % To use fill, eliminate all the NaNs
  fdem=find(isnan(lon));
  % Just need to know a lot about this data
  want=fdem(22)+1:fdem(23)-1;
  eant=fdem(95)+1:fdem(96)-1;
  % Redo the data set by fixing England and adding Antarctica last, etc
  lon=[lon([1:fdem(2)-1 fdem(2)+1:fdem(22) fdem(23)+1:fdem(95) ...
	    fdem(96)+1:fdem(104)-1 fdem(104)+1:end]) ; ...
       NaN; lon(eant) ; lon(want) ; ...
       360  ; 360 ; 0 ; 0];
  lat=[lat([1:fdem(2)-1 fdem(2)+1:fdem(22) fdem(23)+1:fdem(95) ...
	    fdem(96)+1:fdem(104)-1 fdem(104)+1:end]) ; ...
       NaN; lat(eant) ; lat(want) ; ...
       lat(want(end)) ; -90 ; -90 ; lat(eant(1))];
  % And partitioning again
  fdem=find(isnan(lon));
  % Take out last one
  tri=101;
  lon=lon([1:fdem(tri)-1 fdem(tri+1)+1:end]);
  lat=lat([1:fdem(tri)-1 fdem(tri+1)+1:end]);
  lon=lon+ofs;
  XYZ=[lon+ofs lat];
  % And partitioning again
  fdem=find(isnan(lon));
  % Now start with the patches
  beg=1; handl=nan(1,length(fdem));
  hold on
  for i=1:length(fdem)
    handl(i)=fill(lon(beg:fdem(i)-1),lat(beg:fdem(i)-1),pcol);
    beg=fdem(i)+1; 
  end
  fill(lon(beg:end),lat(beg:end),pcol)
  hold off
  axlim=[0 360 -90 90];
 case {3,4}
  % Convert to spherical coordinates
  lon=cont(:,1)/180*pi;
  lat=cont(:,2)/180*pi;
  rad=repmat(1.001,size(lat));
  % Convert to Cartesian coordinates
  [xx,yy,zz]=sph2cart(lon,lat,rad);
  XYZ=[xx yy zz];
  if res==4
    XYZ=XYZ(zz>0,:);
  end
  % Now need to take out the annoying connecting lines
  XYZ=penlift(XYZ);
  
  % And plot this, too
  handl=plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'k');
  axlim=[];
 case 9
  % Restrict to the map requested
  lon=cont(:,1); lat=cont(:,2);
  lon(~(lon>=c11(1) & lon<=cmn(1)))=NaN;
  lat(~(lat<=c11(2) & lat>=cmn(2)))=NaN;
  axlim=[min(lon([~isnan(lon) & ~isnan(lat)]))...
	 max(lon([~isnan(lon) & ~isnan(lat)]))...
	 min(lat([~isnan(lon) & ~isnan(lat)]))...
	 max(lat([~isnan(lon) & ~isnan(lat)]))];
  [xic,etac]=sphere2cube(lon,lat);
  hold on
  % Protect against unsightly junks
  for in=1:6
    [nx{in},ny{in}]=penlift(xic{in},etac{in});
  end
  [pc,pgc]=plotonchunk(nx,ny);
  hold off; set(pc,'Color','k','lines','-','marker','none')
  delete(cat(1,pgc{:})) 
  handl=pc;
  % XYZ=[lon+ofs lat];
  % Rather return the processed coordinates
  XYZ=[nx ; ny];
 case 11
  % This from LORIS1 and EOS1
  lon=lolax(1);
  lat=lolax(2);
  % Set view angles ahead of time as an explicit longitude and latitude
  [xv,yv,zv]=sph2cart(lon*pi/180,lat*pi/180,1);

  % Convert to spherical coordinates
  lonc=cont(:,1)/180*pi;
  latc=cont(:,2)/180*pi;
  rara=1.01;
  radc=repmat(rara,size(latc));
  % Convert to Cartesian coordinates
  [xx,yy,zz]=sph2cart(lonc,latc,radc);
  XYZ=[xx yy zz];

  % Inner product selectivity
  yes=[xv yv zv]*XYZ'>0; 
  % This protection from jumps is straight from PLOTCONT
  XYZ=penlift(XYZ(yes,1:3));
  % And then finally do it
  skl=0.99;
  handl(1)=plot3(XYZ(:,1)*skl,XYZ(:,2)*skl,XYZ(:,3)*skl,'k-');
  hold on
  % Plot a bit of the equator
  [xe,ye,ze]=sph2cart(linspace(0,2*pi,100),0,rara);
  xyze=[rotz(-lon*pi/180)*roty(-[90-lat]*pi/180)*[xe ; ye ; repmat(ze,1,length(ye))]]';
  handl(2)=plot3(xyze(:,1),xyze(:,2),xyze(:,3),'k');
  view([xv,yv,zv]); [AZ,EL]=view; % [xv yv zv];
  disp(sprintf('Azimuth: %i ; Elevation: %i',round(AZ),round(EL)))
  axlim=[-1 1 -1 1 -1 1]*rara*1.01;
  axis(axlim); axis off
  % Where is the North Pole?
  hold on
  pnp=plot3(0,0,1,'MarkerF','k','MarkerE','k','Marker','o');
  hold off
  % Where is the Viewing Axis?
  hold on
  pva=plot3(xv,yv,zv,'MarkerF','k','MarkerE','k','Marker','o');
  hold off
  delete(pnp)
  delete(pva)
end

defval('AZ',NaN)
defval('EL',NaN)

% Generate output
vars={axlim,handl,XYZ,xyze,AZ,EL};
varargout=vars(1:nargout);

% Last-minute cosmetic adjustment
axis equal
hold off
