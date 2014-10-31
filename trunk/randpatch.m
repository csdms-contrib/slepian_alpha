function varargout=randpatch(N,TH,phi0,theta0)
% [lon,lat,N]=randpatch(N,TH,phi0,theta0)
%
% Generates randomly distributed points on the sphere that are restricted
% to fall within a spherical patch
%
% INPUT:
%
% N        Number of points
% TH       Angular extent of the spherical cap (degrees)
% phi0     Longitude of the center (degrees)
% theta0   Colatitude of the center (degrees)
%
% OUTPUT:
%
% lon,lat  Vectors with random locations
% N        Actual number of points
%
% EXAMPLE:
%
% randpatch('demo1')
% randpatch('demo2')
%
% SEE ALSO:
%
% CAPLOC
% 
% Last modified by fjsimons-at-alum.mit.edu, 03/02/2010

% Supply default values
defval('N',50)

if ~isstr(N)
  defval('TH',15)
  defval('phi0',55)
  defval('theta0',70)

  % Figure out the fractional area of the patch
  A=spharea(TH,1);

  % Slight inflation factor
  infl=1.15;
  
  % Thus will need to generate roughly this many more points
  [lon,lat]=randsphere(ceil(infl*N/A));

  % Figure out the distance to the center
  [gcdkm,delta]=grcdist([phi0 90-theta0],[lon lat]);
  
  % And select only those within range
  lon=lon(delta<=TH);
  lat=lat(delta<=TH);
  
  % And now I get as close as I can to the requested number of points
  lon=lon(1:min(N,length(lon)));
  lat=lat(1:min(N,length(lat)));

  % Prepare output
  vars={lon,lat,length(lon)};
  varargout=vars(1:nargout);
  
elseif strcmp(N,'demo1')
  [lon1,lat1]=randpatch;
  [lon2,lat2]=caploc;

  disp(sprintf('Actual number of points is %i',length(lon1)))
  
  clf
  c=plot(lon2,lat2,'k'); hold on
  d=plot(lon1,lat1,'o'); axis image
  deggies(gca)
elseif strcmp(N,'demo2')
  Ntry=1000;
  Nact=nan(Ntry,1);
  for index=1:Ntry
    lon1=randpatch;
    Nact(index)=length(lon1);
  end
  
  clf
  hist(Nact)
  title('Number generated that fall within patch')
end
