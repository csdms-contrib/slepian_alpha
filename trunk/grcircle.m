function [longclatgc,delta]=grcircle(lon1lat1,lon2lat2,longr,latgr)
% [longclatgc,delta]=GRCIRCLE([lon1 lat1],[lon2 lat2],longr,latgr)
%
% Calculates great circle points through a pair of points (maybe on a grid)
%
% INPUT:
%
% [lon1 lat1]  longitude and latitude of the start point [radians]
% [lon2 lat2]  longitude and latitude of the end point [radians]
% longr        regular-grid longitudes OR integer number of points to be calculated
% latgr        regular-grid latitudes [optional]
%
% OUTPUT:
%
% longclatgc   longitudes/latitudes of the great circle [radians]
% delta        incremental distance along the great circle
%
% If the longitudinal difference is 2pi inserts a NaN.
%
% EXAMPLE:
%
% LaJolla=[242.7340 32.8660]*pi/180; StLouis=[269.8160 38.6280]*pi/180; 
% Harvard=[288.4420 42.5060]*pi/180; LosAngeles=[241.6330 34.0830]*pi/180;
%% Calculate great circle on a grid
% longr=linspace(0,2*pi,100); latgr=linspace(-pi/2,pi/2,100);
% [lola,delta]=grcircle(Harvard,LosAngeles,longr,latgr);
% plot(lola(:,1)*180/pi,lola(:,2)*180/pi,'b+'); hold on
% [lo,la]=meshgrid(longr,latgr); fridplot(lo*180/pi,la*180/pi)
%% or calculate just a number of points
% [lola,delta]=grcircle(LaJolla,StLouis,100);
% plot(lola(:,1)*180/pi,lola(:,2)*180/pi,'ro');
% plotcont([220 60],[320 20]); hold off; axis image; axis([220 320 20 60])
%
% Last modified by fjsimons-at-alum.mit.edu, 03/10/2009

defval('longr',100);

[lon1,lat1]=deal(lon1lat1(1),lon1lat1(2));
[lon2,lat2]=deal(lon2lat2(1),lon2lat2(2));

lon1=lon1+(lon1<0)*2*pi;
lon2=lon2+(lon2<0)*2*pi;

[latmax,lonmax]=apex(lat1,lon1,lat2,lon2);

if nargin<=3 & round(longr)==longr 
  % Fixed number of points
  
  % Calculate total distance 
  totdis=grcdist([lon1 lat1]*180/pi,[lon2 lat2]*180/pi)/6371;
  %  'longr' points along total distance
  delta=linspace(0,totdis,longr)';
  [lon,lat]=distgrc([lon1 lat1],[lon2 lat2],delta);
  lon=lon+2*pi*(lon<0);  
  longclatgc=[lon lat];
  
else % Calculation on grid

  % Make sure 'latgr' and 'longr' are column vectors
  latgr=latgr(:); longr=longr(:);
  
  %-------------------------------------------------
  if abs(lon1-lon2)<pi
    l=find(longr>min(lon1,lon2) & longr < max(lon1,lon2));
  else
    l=find(longr<min(lon1,lon2) | longr > max(lon1,lon2));
  end
  
  lonmin=rem(latmax,2*pi);
  
  e=tan(latmax);
  
  if sin(lon1-lonmax)*sin(lon2-lonmax)<0
    
    if cos(lon1-lonmax)+cos(lon2-lonmax)>0
      w1=find(lat1<latgr & latmax>latgr);
      w2=find(lat2<latgr & latmax>latgr);
      w=[w1;w2];              
      if sin(lon1-lonmax)>0
	minus=[zeros(size(w1)); ones(size(w2))];
      else
	minus=[ones(size(w1)); zeros(size(w2))];
      end
    else
      w1=find(lat1>latgr & -latmax<latgr);
      w2=find(lat2>latgr & -latmax<latgr);
      w=[w1;w2];
      if sin(lon1-lonmax)>0
	minus=[zeros(size(w1)); ones(size(w2))];
      else
	minus=[ones(size(w1)); zeros(size(w2))];
      end
    end
    
  else
    w=find(min(lat1,lat2)<latgr & max(lat1,lat2)>latgr);
    
    if sin(lon1-lonmax)>0
      minus=zeros(size(w));
    else
      minus=ones(size(w));
    end
  end
  
  latca=atan(e.*cos(longr(l)-lonmax));
  
  lonca=(1-2*minus).*acos(tan(latgr(w))./e)+lonmax;
  
  latgc=[latgr(w);                        latca;     lat1; lat2];
  longc=[rem(lonca, 2*pi)+2*pi*(lonca<0); longr(l); lon1; lon2];
  
  delta=real(acos(sin(lat1).*sin(latgc)+cos(lat1).*cos(latgc).*cos(lon1-longc)));
  
  [delta,j]=sort(delta);
  
  longclatgc=[longc(j) latgc(j)];
end

% Should only be one number
lopro=find(abs(diff(longclatgc(:,1))-2*pi)<(2*pi/length(longclatgc)));
if ~isempty(lopro)
  longclatgc=[longclatgc(1:lopro,:) ; NaN NaN ; longclatgc(lopro+1:end,:)];
end


