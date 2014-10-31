function varargout=polargrid(nx,ny,method)
% [PHI,THETA]=POLARGRID(nx,ny,method)
%
% Returns a ny-by-nx grid with the geographical North Pole a nonsingular
% point in its middle. Suitable for maps where the interest is
% polar. When you use this to interpolate geographic locations, you're
% looking down on the North Pole, the Greenwich meridian is down,
% colatitude increases radially outwards from the center and you
% increase in longitude anticlockwise from there. Run the topography
% example to orient yourself.
%
% INPUT:
%
% nx        Number of points in x; must be odd to contain South Pole
% ny        Number of points in y; must be odd to contain North Pole
%                                  must be quod to contain the Equator
%           Usually, you'll want nx about half the size of ny
% method    1 Much faster; rotates a lon-lat grid over two angles
%           2 Much slower; rotates a single meridian over successive angles
%             The results are equivalent though not identical in PHI:
%             differences are either modulo 2pi or else are at polar
%             latitudes, where they are immaterial.
%
% OUTPUT:
%
% PHI     Longitudes of the grid, in radians
% THETA   Colatitudes of the grid, in radians
%
% EXAMPLE:
%
% polargrid('demo1')
%
% Last modified by fjsimons-at-alum.mit.edu, 12/09/2010

defval('nx',2*20+1)

if ~isstr(nx)
  
  defval('ny',4*20+1)
  defval('method',1)

  % Check the input is odd, if not, the poles are not hit exactly 
  difer(mod(nx,2)-1,[],[],NaN) % Check "odd"
  difer(mod(ny,2)-1,[],[],NaN) % Check "odd"
  difer(mod(ny-1,4),[],[],NaN) % Check "quod"

  switch method
   case 1
    % Can I get the same result using just one rotate over two angles
    % NOTE THAT IF ALFA is 0 IN THE NEXT LINE, I PUT THE SOUTH POLE IN THE MIDDLE
    % NOTE THAT -GAMA CONTROLS THE LONGITUDES OF THE SINGULAR POINTS
    % Change these if you want the "pole" somewhere else
    [THETA,PHI]=rottp(linspace(0,pi,nx),linspace(0,2*pi,ny)',...
		      pi,pi/2,-pi/2);
    % Need to flip dimensions after rotation of course
    THETA=THETA';
    PHI=PHI';
   case 2
    % Construct a half a great circle hitting the South Pole first, and then
    % rotate it around the y axis, hitting the North Pole in the middle, and
    % finally, the South Pole, again. Start at [90,-90], over South Pole, back
    % up to [90,90], and rotate this
    th1=[linspace(pi/2,pi,(nx+1)/2) indeks(linspace(pi,pi/2,(nx+1)/2),'2:end')];
    % Claim the South Pole has longitude 0, which isn't technically necessary
    ph1=[repmat(-pi/2,1,(nx+1)/2-1) 0 repmat(pi/2,1,(nx-1)/2)];

    % Now find the coordinates that contain this thing through ny rotations
    rotem=linspace(0,2*pi,ny);
    [THETA,PHI]=deal(nan(ny,nx));
    PHI(1,:)=ph1; THETA(1,:)=th1;
    for ix=2:length(rotem)-1
      [THETA(ix,:),PHI(ix,:)]=rottp(th1,ph1,0,rotem(ix),0);
    end
    PHI(end,:)=ph1; THETA(end,:)=th1;
  end

  % Make sure that the PHI's are contained within zero and two pi
  PHI=PHI+2*pi*(PHI<0);

  % And order so as to preserve West to East as left to Right
  % And keep the front of the Earth in the front
  PHI=flipud(PHI);
  THETA=flipud(THETA);

  % Output
  vars={PHI,THETA};
  varargout=vars(1:nargout);
  
  % This is the physical grid on the sphere upon which you want to expand
  % or interpolate the function. You then just display the new function on
  % a rectangular grid with these PHI and TH rolled out. Should you try the
  % inverse transform, you'll find the North Pole smack in the middle of x-y,
  % and the South Pole is repeated on first and last lines, appearing at
  % somewhat random, since immaterial, x-distances. This suggested method 1
  % to me. However, I just need to display the result on a square grid:
  % where the North Pole will be in the center and the South Pole repeated
  % along the top and bottom lines, and the right and left edges are single
  % and singular points located on the former equator.

  % Inverse transform
  % [THETAX,PHIX]=rottp(THETA,PHI,pi/2,-pi/2,-pi);

  % As to the examples, so basically, you can use INTERP2 rather than
  % GRIDDATA, since the data is given on a rectangular grid in longitude
  % and latitude, and you want data at irregular locations in longitude and
  % latitude - knowing they transform to a regular grid in x and y, rather
  % than the other way round, as is the  solution I adopted in
  % PLOTPLM. This is much faster. [Our trick, of course, allows also to NOT
  % interpolate but rather EVALUATE the harmonics directly onto this
  % funky grid. Or should we? We may not have to.] However, what we lose
  % here is the same as what we lose  with the Mercator
  % projection... circular things are distorted progressively away from the
  % old North Pole, now a point on the new equator. But the distortion
  % isn't so bad if you watch the sizes of nx and ny. 

  % So the philosophy of old, in PLOTPLM, was to map the sphere to the
  % plane, which produces an irregular set of Cartesian points, and then
  % use GRIDDATA on a regular grid. Now it's to INTERPOLATE the regular
  % points in spherical coordinates on an irregular grid, which is then
  % mapped to produce a regular set. The latter approach is more efficient
  % yet, with the current POLARGRID, it produces a distorted
  % projection. Should perhaps fix that.

  % I suppose I should redo this whole thing and get a true polar grid for
  % evaluation of the harmonics.
elseif strcmp(nx,'demo1')
  %z=flipud(etopo(fullfile(getenv('IFILES'),'TOPOGRAPHY','EARTH','ETOPO1')));
  z=load('topo');
  z=z.topo;
  lon=linspace(0,2*pi,size(z,2));
  colat=linspace(0,pi,size(z,1));
  [PHI,THETA]=polargrid(2*120+1,4*120+1,1);
  zi=interp2(lon,colat,z,PHI,THETA);
  subplot(121); imagesc(zi); axis image ; demmap; title('Method 1')
  [PHI,THETA]=polargrid(2*120+1,4*120+1,2);
  zi=interp2(lon,colat,z,PHI,THETA);
  subplot(122); imagesc(zi); axis image ; demmap; title('Method 2')

  G=repmat(rindeks(galpha(10,18,1,linspace(0,pi,360),0),1),720,1)';
  Gi=interp2(linspace(0,2*pi,size(G,2)),...
	     linspace(0,pi,size(G,1)),G,PHI,THETA);
  subplot(122); imagesc(Gi); axis image ; kelicol
  fig2print(gcf,'portrait')
end

