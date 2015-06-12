function varargout=plotplm(data,lon,lat,meth,degres,th0,sres,cax)
% PLOTPLM(data,lon,lat,meth,degres,th0,sres,cax)
% PLOTPLM(lmcosi,[],[],meth,degres,th0,sres,cax)
% PLOTPLM(l,m,cosi,meth,degres,th0)
% [data,ch,ph,lon,lat]=PLOTPLM(...)
%
% Global plotting routine for fields coming out of PLM2XYZ,
% plotting REAL spherical harmonics normalized to 4pi.
%
% INPUT:
%
% data           2-D data/or LMCOSI matrix, OR
% lmcosi         Degree, Order, Cosine, Sine Coefficient matrix
%                (column size is either 4 or 6), OR
% l              A single angular degree AND
% m              A single angular order
% cosi           1 A unit expansion coefficient for the cosine [default] 
%                2 A unit expansion coefficient for the sine harmonic 
% lon,lat        2-D grid in radians [default (0,2pi) and (-pi/2,pi/2)]
% meth           1  Mollweide projection [default]
%                2  3-D sphere, no topography
%                3  Plots the spectrum (only if lmcosi specified)
%                4  Flat rectangular projection for entire globe
%                5  Plot looking down on the North Pole; longitudes off
%                6  Plot looking down on the South Pole; longitudes off
% degres         Resolution in degree (checked for Nyquist) [default: Nyquist]
% th0            Plot small circle at colatitude th0, in degrees
% sres           0 Take default for polar projection [default]
%                1 Preserve resolution (loosely) for polar projection
% cax            Put in color saturation for options 2, 5, and 6
%
% 'lon' and 'lat' are in radians.
%
% OUTPUT:
%
% data           Spatial data in case coefficients were given
%                Handles to spectral log-log plots for meth 3
% ch             Handle to the continents and the map edge:
%                ch{1}(1) is the handle to the continents
%                ch{2}(1) is the handle to the left border
%                ch{2}(2) is the handle to the right border
% ph             Handle to the plates or the small circles
% lon,lat        Useful if you want to plot contours on top of this
%
% EXAMPLE:
%
% load topo; topo=flipud(topo); plotplm(topo,[],[],2)
%
% See also PLOTONSPHERE, PLOTONEARTH, CPX2RSH, RSH2CPX, ADDCB.
%
% Last modified by fjsimons-at-alum.mit.edu, 02/20/2012
% Last modified by charig-at-princeton.edu, 05/14/2015

defval('meth',1)
defval('degres',[])
defval('th0',[])
defval('sres',0)

if (size(data,2)==4 | size(data,2)==6) & meth<7 & meth~=3
  [data,lon,lat]=plm2xyz(data,degres);
  lon=lon*pi/180;
  lat=lat*pi/180;
elseif length(data)==1 % Plot a single spherical harmonic order
  defval('lat',1) % This is the COSI component then
  [m,l,mz,lmcosi]=addmon(data);
  lmcosi(mz(data+1)+lon,2+lat)=1;
  [data,lon,lat]=plm2xyz(lmcosi,degres);
  lon=lon*pi/180;
  lat=lat*pi/180;
end

switch meth
 case 1
   defval('lon',linspace(0,2*pi,size(data,2)))
   defval('lat',linspace(pi/2,-pi/2,size(data,1)))
   % Need to make special provisions for PCOLOR compared to IMAGESC
   % The input values are PIXEL centered
   dlat=(lat(1)-lat(2))/2;
   dlon=(lon(2)-lon(1))/2;
   lat=[lat(1) lat(2:end)+dlat lat(end)];
   lon=[lon(1) lon(2:end)-dlon lon(end)];
   [lon,lat]=meshgrid(lon,lat);
   [xgr,ygr]=mollweide(lon,lat,pi);
   pc=pcolor(xgr,ygr,adrc(data)); shading flat
   [ax1,ch,XY1]=plotcont([],[],2);
   [ph,XY2]=plotplates([],[],2);
   axis image
   % OPENUP doesn't nearly do as good a job as this here
   ylim(ylim*1.05)
   axis off
 case 2
  % Make sphere for the data
  defval('lon',linspace(0,360,100)/180*pi);
  defval('lat',linspace(90,-90,100)/180*pi);
  [lon,lat]=meshgrid(lon,lat);
  rads=ones(size(lat));
  [x,y,z]=sph2cart(lon,lat,rads);
  % Color saturation
  defval('cax',minmax(data(:)));
  data(data>cax(end))=cax(end);
  data(data<cax(1))=cax(1);
  surface(x,y,z,'FaceColor','texture','Cdata',data);   
  hold on
  % Plot the continents
  [axlim,handl,XYZ]=plotcont([0 90],[360 -90],3);  
  delete(handl)
  hold on
  ch=plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'k-','LineWidth',1.5);
  axis image
  shading flat
  view(140,30)
  set(gca,'color',[.9 .9 .9]-.1)
  axis('off')
  colormap default
  ph=[];
 case 3
  [sdl,l,bta,lfit,logy,logpm]=plm2spec(data,2);
  % The following is now well taken care of by LIBBRECHT
  % l=l(l<256 & l>0); 
  % sdl=sdl(l<256 & l>0);
  l=l(l>0);
  sdl=sdl(l>0);
  disp(sprintf('PLOTPLM Spectrum plot from %i to %i',...
	       l(1),l(end)))
  clear data
  data(:,1)=loglog(l,sdl);
  set(data(:,1),'Marker','+')
  grid on
  tix=sort([l(1) 10.^[1 2 3] l(end)]);
  tix=tix(tix<=l(end)); 
  xlim([l(1) l(end)])
  set(gca,'xtick',tix)
  tiy=round(minmax(log10(sdl(~~sdl))));
  set(gca,'ytick',10.^[tiy(1):1:tiy(2)],...
	  'yminortick','off','yminorgrid','off')
  hold on
  data(:,2)=loglog(lfit,logy);
  % Do not use axis tight or openup on loglog
  ylim([0.9 1.1].*minmax([sdl(:) ; logy(:)]))  
  set(data(:,2),'Color','r')
  longticks(gca)
  data(:,3)=title(sprintf('%s= %8.3f','\beta',bta));
  data(:,4)=bta;
  [ch,ph]=deal([]);
 case 4
   imagef([0 90],[360 -90],data)
   [ax1,ch,XY1]=plotcont([],[],1);
   [ph,XY2]=plotplates([],[],1);
   axis image
 case 5
  % Project northern hemisphere only
  lon=linspace(0,2*pi,size(data,2));
  lat=linspace(pi/2,0,floor(size(data,1)/2));
  [LON,LAT]=meshgrid(lon,lat);
  % Radius from 0 to 1; longitude is azimuth
  r=cos(LAT);
  x=r.*cos(LON);
  y=r.*sin(LON);
  % Resolution of upper hemispheric projection
  if sres==0
    % Only project the upper hemisphere
    X=linspace(-1,1,500);
    Y=linspace(1,-1,500);
  else
    X=linspace(-1,1,length(lat));
    Y=linspace(1,-1,length(lat));
  end
  warning off MATLAB:griddata:DuplicateDataPoints
  % Watch out: the POLE is a new problem with the latest version of
  % Matlab; need to fake this entirely. See POLARGRID and LORIS1. 
  y(1,:)=y(2,:)/2;
  x(1,:)=x(2,:)/2;
  Z=griddata(x,y,data(1:length(lat),:),X,Y(:));
  disp('Data gridded by PLOTPLM')
  warning on MATLAB:griddata:DuplicateDataPoints
  % imagefnan([-1 1],[1 -1],Z,gray(10),minmax(Z(:)))  
  % colormap(gray(10)); hold on
  % IMAGEFNAN indexes the color map directly
  defval('cax',minmax(Z(:)));
  imagefnan([-1 1],[1 -1],Z,kelicol,cax)
  colormap(kelicol); hold on
  ph(1)=circ(1); 
  ph(2)=circ(sin(th0*pi/180));
  set(ph(1:2),'LineW',1)
  set(ph(2),'LineS','--')
  axis([-1.0100    1.0100   -1.0100    1.0100])
  axis off
  [axlim,handl,XYZ]=plotcont([0 90],[360 -90],4);
  delete(handl)
  hold on
  ch=plot(XYZ(:,1),XYZ(:,2),'k-','LineWidth',1);
  data=Z;
 case 6
  % Project southern hemisphere only
  data = flipud(data);
  data = fliplr(data);
  lon=linspace(0,2*pi,size(data,2));
  lat=linspace(pi/2,0,floor(size(data,1)/2));
  [LON,LAT]=meshgrid(lon,lat);
  % Radius from 0 to 1; longitude is azimuth
  r=cos(LAT);
  x=r.*cos(LON+pi/2);
  y=r.*sin(LON+pi/2);
  % Resolution of upper hemispheric projection
  if sres==0
    % Only project the upper hemisphere
    X=linspace(-1,1,500);
    Y=linspace(1,-1,500);
  else
    X=linspace(-1,1,length(lat));
    Y=linspace(1,-1,length(lat));
  end
  warning off MATLAB:griddata:DuplicateDataPoints
  % Watch out: the POLE is a new problem with the latest version of
  % Matlab; need to fake this entirely
  y(1,:)=y(2,:)/2;
  x(1,:)=x(2,:)/2;
  Z=griddata(x,y,data(1:length(lat),:),X,Y(:));
  disp('Data gridded by PLOTPLM')
  warning on MATLAB:griddata:DuplicateDataPoints
  % imagefnan([-1 1],[1 -1],Z,gray(10),minmax(Z(:)))  
  % colormap(gray(10)); hold on
  % IMAGEFNAN indexes the color map directly
  defval('cax',minmax(Z(:)));
  imagefnan([-1 1],[1 -1],Z,kelicol,cax)
  colormap(kelicol); hold on
  ph(1)=circ(1); 
  ph(2)=circ(sin(th0*pi/180));
  set(ph(1:2),'LineW',1)
  set(ph(2),'LineS','--')
  axis([-1.0100    1.0100   -1.0100    1.0100])
  axis off
  [axlim,handl,XYZ]=plotcont([0 90],[360 -90],11);
  delete(handl)
  hold on
  ch=plot(XYZ(:,1),XYZ(:,2),'k-','LineWidth',1);
  data=Z;

otherwise
  error('Not a valid method')
end

% Prepare desired output
varns={data,ch,ph,lon,lat};
varargout=varns(1:nargout);

