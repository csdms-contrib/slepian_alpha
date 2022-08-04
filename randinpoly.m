function varargout=randinpoly(dom,N)
% [lon,lat]=randinpoly(dom,N) 
%
% Generates less than N equal-area random points within the region dom 
%
% INPUT:
%
% dom 		Named region such as eurasia, namerica, samerica, etc.
% N 		Number of points to generate. The resulting number will be
%            close but not exact
%
% OUTPUT:
%
% lon		longitudes of the random points
% lat 		latitudes of the random points
%
% EXAMPLE:
%
% randinpoly('demo1')
% 
% Last modified by charig@arizona.edu, 6/10/2022
% Last modified by plattner-at-alumni.ethz.ch, 6/27/2016

defval('dom','namerica')
defval('N',100)

if ~strncmp(dom,'demo',4)

    % Get the XY coordinates
    XY=eval(sprintf('%s()',dom));

    % Get the center of mass
    [lonc,latc] = rcenter(XY);

    % Calculate the spherical distance from the center of mass to each point on
    % the closed curve.
    [gcdkm,delta] = grcdist(XY,[lonc latc]);

    % Biggest distance in degrees?
    r = max(delta);

    % Figure out the fractional area of the patch
    A1 = spharea(dom);
    A2 = spharea(r,1);

    % Slight inflation factor
    %infl=1.15;
  
    % Thus will need to generate roughly this many more points
    [lon,lat]=randsphere(ceil(N/A1));

    % Sometimes the XY long goes beyond 360. Make lon go beyond by the same
    % amount
    if max(XY(:,1))>360
      shift=max(XY(:,1))-mod(max(XY(:,1)),360);
      lon2=lon+shift;
      lon=[lon;lon2];
      lat=[lat;lat];
    end

    % Now delete all points outside the region
    index=inpolygon(lon,lat,XY(:,1),XY(:,2));
    lon=lon(index);
    lat=lat(index);

    % And now I get as close as I can to the requested number of points
    lon=lon(1:min(N,length(lon)));
    lat=lat(1:min(N,length(lat)));

    % Prepare output
    vars={lon,lat};
    varargout=vars(1:nargout);

elseif strcmp(dom,'demo1')
    [lon,lat]=randinpoly('namerica',100);

    XY = namerica();
    figure
    plot(XY(:,1),XY(:,2),'k-'); axis equal; grid on
    hold on
    plot(lon,lat,'.');

end



