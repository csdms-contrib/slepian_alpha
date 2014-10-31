function [Bl,dels,r,lon,lat,lmcosi]=geoboxcap(L,dom,N,degres,act,lonc,latc)
% [Bl,dels,r,lon,lat,lmcosi]=GEOBOXCAP(L,dom,N,degres,act,lonc,latc)
%
% Returns the spherical-harmonic POWER SPECTRUM in the UNIT-NORMALIZED
% spherical-harmonic basis of an all-or-nothing BOXCAR over a particular
% geographical region. It does that by calculating the power spectrum in
% the 4-pi-normalized basis and then adjusting the coefficients to pick up
% the power so that it matches that in the unit-normalized basis.
%
% INPUT:
% 
% L          Bandwidth of the output spectrum
% dom        'africa',  'eurasia',  'namerica', 'australia', 'greenland',
%            'samerica', 'amazon', 'orinoco', 'gpsnamerica', 'antarctica'
%            OR: [lon lat] an ordered list defining a closed curve [degrees]
% N          The splining smoothness of the geographical boundary [default: 10] 
% degres     The resolution of the underlying spatial grid [default: Nyquist]
% act        1 Actually perform the calculations [default]
%            0 Don't, but simply return the (rotated) mask function
% lonc,latc  Rotate coordinates by this amount 
%
% OUTPUT:
%
% Bl         The power spectrum 
% dels       The spherical harmonic degrees
% r          The "mask" with the "continent function"
% lon,lat    The grid on which  this is defined
% lmcosi     The complete set of spherical harmonic expansion coefficients
%
% EXAMPLE:
%
% [Bl1,dels1]=geoboxcap(18,'australia',[],[]);
% [Bl2,dels2]=geoboxcap(18,'australia',[],1);
%
% SEE ALSO:
%
% BPBOXCAP, KERNELC, GAMMAP
% 
% Last modified by fjsimons-at-alum.mit.edu, 05/13/2013

% Default inputs
defval('L',18)
defval('dom','australia')
defval('N',10)
defval('act',1)
defval('lonc',0)
defval('latc',0)
defval('xver',0)

if isstr(dom)
  % Run the named function to return the coordinates
  XY=eval(sprintf('%s(%i)',dom,N));
else
  % Get the coordinates as defined from the input in degrees
  XY=dom; 
end

% Make sure the coordinates make sense
XY(:,1)=XY(:,1)-360*[XY(:,1)>360];
% Not good for INPOLYGON is the next line
% XY(:,1)=XY(:,1)+360*[XY(:,1)<0];

% Make a grid of ones and zeroes depending on the desired resolution
degN=180/sqrt(L*(L+1));
defval('degres',degN);

% Default grid is all of the planet
defval('c11cmn',[0 90 360 -90]);

% The number of longitude and latitude grid points that will be computed
nlon=ceil([c11cmn(3)-c11cmn(1)]/degres+1);
nlat=ceil([c11cmn(2)-c11cmn(4)]/degres+1);

% Longitude grid vector in degrees
lon=linspace(c11cmn(1),c11cmn(3),nlon);
% Latitudex grid vector in degrees
lat=linspace(c11cmn(2),c11cmn(4),nlat);
% Make the input grid
r=repmat(0,nlat,nlon);
[LON,LAT]=meshgrid(lon,lat);

% Now decide if we're inside or outside of the region
% This isn't going to work for continents straddling the date line as in
% Africa, so, rotate such cases out of the way! E.g. for the Volta basin,
% stick on lonc=-10;
if lonc~=0 || latc~=0
  XYor=XY;
  [X,Y,Z]=sph2cart(XY(:,1)*pi/180,XY(:,2)*pi/180,1);
  [Xc,Yc,Zc]=sph2cart(lonc*pi/180,latc*pi/180,1);
  xyzp=[rotz(lonc*pi/180)*roty(-latc*pi/180)*[X(:) Y(:) Z(:)]']';
  [phi,piminth]=cart2sph(xyzp(:,1),xyzp(:,2),xyzp(:,3));
  lon=phi*180/pi; lat=piminth*180/pi;
  XY=[lon lat];
end
r(inpolygon(LON,LAT,XY(:,1),XY(:,2)))=1;

if act==1
  % And now do the spherical harmonic transform but only to L
  % This takes time!
  lmcosi=xyz2plm(r,L);
  
  if lonc~=0 || latc~=0
    % Now have the pleasure to rotate this back
    if latc==0
      % Only a longitudinal rotation, much faster!
      [C,S]=rotcof(lmcosi(:,3),lmcosi(:,4),-lonc*pi/180);
      lmcosi(:,3)=C;
      lmcosi(:,4)=S;
    else
      % A longitudinal and latitudinal rotation
      lmcosi=plm2rot(lmcosi,-lonc,latc,0);
    end
  end

  if xver==1
    % Make some plots
    figure(1)
    clf
    fridplot(LON,LAT)
    xlim(xpand(minmax(lon)))
    ylim(xpand(minmax(lat)))
    hold on
    plot(XY(:,1),XY(:,2),'r')
    hold off
    figure(2)
    clf
    rp=plm2xyz(lmcosi);
    imagefnan([0 90],[360 -90],rp)
    hold on
    % Good enough for plotting
    XYor(:,1)=XYor(:,1)+360*[XYor(:,1)<0];
    [X,Y]=penlift(XYor(:,1),XYor(:,2));
    plot(X,Y,'r')
    hold off
  end
  % Calculate the power spectral density in the 4pi basis
  [Bl,dels]=plm2spec(lmcosi,2);
else
  [Bl,dels,lmcosi]=deal(NaN);
end

% Check the B0 term which should equal the area^2 divided by (4pi)^2 in the
% 4pi-normalized basis, where the Y00 term equals 1
A1=4*pi*spharea(XY); A2=4*pi*areaint(XY(:,2),XY(:,1));
disp(sprintf('GEOBOXCAP A: %6.3f ; SPHAREA A: %6.3f ; AREAINT A: %6.3f',...
	4*pi*sqrt(Bl(1)),A1,A2))

% Make the adjustment so that this power spectrum is like the one from
% BPBOXCAP, i.e. so that it is also quoted in the unit-normalized
% harmonics.
Bl=Bl*4*pi;

% If these ever used to build loops etc they must be a row
dels=dels(:)';

% Provide output as desired
varns={Bl,dels,r,lon,lat,lmcosi};
varargout=varns(1:nargout);

