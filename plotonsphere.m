function varargout=plotonsphere(data,rang,mygrid,conts)
% pc=PLOTONSPHERE(data,rang,mygrid,conts)
%
% Plots data as raised topography and colors onto a sphere. Has the option
% to add grid lines and continent outlines.
%
% INPUT:
%
% data     Standard 2D geographic data; i.e. the first column is the
%          Greenwich meridian in the x-z plane.
% rang     Exaggeration of topography compared to 1 [default: 0]
% mygrid   Vector of grid spacing in the longitude and latitude 
%          dimensions (e.g. [30 15]) [default: empty, no grid]
% conts    Switch if you want continents drawn [default: 0]
%
% OUTPUT:
%
% pc       Handle to the continental outlines if wanted them
%
% See PLOTONEARTH, PLOTPLM
%
% EXAMPLE:
%
% plotonsphere('demo1')
% plotonsphere('demo2')
% plotonsphere('demo3') % Figure of filtered free-air anomaly
%
% Last modified by fjsimons-at-alum.mit.edu, 01/21/2021
% Last modified by charig-at-princeton.edu, 6/17/2016

defval('conts',0)
defval('mygrid',[])

if ~isstr(data)
  defval('rang',0);
  
  % Test periodicity over the longitudes
  data=reduntest(data);
  polestest(data)
  
  [ny,nx]=size(data);
  [x,y,z]=sphere(min([nx ny])-1);
  % Circumvent weird hold behavior
  p=plot(1,1); delete(p)

  % This to make it right - compare standard definitions DT. p 832.
  data=flipud(data);
  Cdata=interp2([1:nx],[1:ny]',data,linspace(1,nx,ny-1),[1:ny-1]');
  s=surface(-x,-y,z,'FaceColor','texture','Cdata',Cdata);
  axis equal

  % Interpolate the data on the sphere's grid
  if nx>ny 
    data=interp2([1:nx],[1:ny]',data,linspace(1,nx,ny),[1:ny]');
  elseif nx<ny
    data=interp2([1:nx],[1:ny]',data,[1:nx],linspace(1,ny,nx)');
  end

  % Do the topography scaling
  data=1+scale(data/max(abs(data(:))),[-1 1]*rang);

  % This is tricky business but it works
  x=-x.*data;      
  y=-y.*data;
  z=z.*data;
  set(s, 'xdata', x, 'ydata', y, 'zdata', z)
  % So x and y face at you
  view(60,20)
  set(gca, 'dataaspectratio', [1 1 1], 'cameraviewangle', 7)
  axis image
  shading flat
  axis off
  
  % The data are interpolated in 2D and mapped to a grid on the sphere. 
  % In order to draw the vector data on top of this grid we need to 
  % interpolate it to the same grid, and give the appropriate radius.
  % Otherwise we might have vectors that are drawn with depth below the
  % data surface.
  
  if conts==1
    % Plot the continents
    lx=xtraxis;
    [jk1,jk2,cont]=plotcont([0 90],[360 -90]);
    delete(lx)

    [lat,lon]=interpm(cont(:,2),cont(:,1),.1);
    cont = [lon lat];
    [ny,nx]=size(data);
    data2=flipud(data);
    % Make the elevations
    Vq=interp2([0:nx-1],[0:nx-1],data2,...
	       cont(:,1)*(nx-1)/360,(cont(:,2)-90)*-1*(ny-1)/180,'linear');
    [ind,colnr,rownr]=cor2ind(cont(:,1),cont(:,2),...
			      [-0.01 90.01],[360.01 -90.01],nx,ny);
    
    lonc=cont(:,1)/180*pi;
    latc=cont(:,2)/180*pi;
    rad=repmat(1,size(latc));
    
    indone=ind; 
    indone(isnan(ind))=1;
    data2=flipud(data);

    [xx,yy,zz]=sph2cart(lonc,latc,Vq);
    % Plot
    hold on
    pc=plot3(xx,yy,zz,'k-','LineWidth',1.5);
    hold off
  else
    pc=NaN;
  end

  % If you want a lon lat grid drawn
  if ~isempty(mygrid)
    if length(mygrid)==1
      degreslon = mygrid;
      degreslat = mygrid;
    elseif length(mygrid)==2
      degreslon = mygrid(1);
      degreslat = mygrid(2);
    end
    [ny,nx]=size(data);
    data2=flipud(data);
    degres=15;
    lon=linspace(0,360,721);
    lat=linspace(90,-90,361);
    % First one here: longitude coordinates of the longitude lines
    latlon=repmat(lon',1,180/degreslat+1);
    latlat=repmat([90:-degreslat:-90],length(lon),1);
    lonlon=repmat([0:degreslon:360],length(lat),1);
    lonlat=repmat(lat',1,360/degreslon+1);
    % Calculate radius for lines of lattitude
    Vq1=interp2([0:nx-1],[0:nx-1],data2,...
		latlon*(nx-1)/360,(latlat-90)*-1*(ny-1)/180,'linear');
    [xx1,yy1,zz1]=sph2cart(latlon/180*pi,latlat/180*pi,Vq1);
    % Calculate radius for lines of longitude
    Vq2 = interp2([0:nx-1],[0:nx-1],data2,...
		  lonlon*(nx-1)/360,(lonlat-90)*-1*(ny-1)/180,'linear');
    [xx2,yy2,zz2]=sph2cart(lonlon/180*pi,lonlat/180*pi,Vq2);
    % Plot
    hold on
    pc=plot3(xx1,yy1,zz1,'k-','LineWidth',1);
    pc=plot3(xx2,yy2,zz2,'k-','LineWidth',1);
    hold off
  end
  % Optional output
  varns={pc};
  varargout=varns(1:nargout);
elseif strcmp(data,'demo1')
  d70=plotplm(7,0,1,2,5); 
  d72=plotplm(7,2,1,2,5);
  d74=plotplm(7,4,1,2,5);
  d76=plotplm(7,6,1,2,5);
  d77=plotplm(7,7,1,2,5);
  clf
  % ah=krijetem(subnum(1,5));
  ah=krijetem(subnum(2,2));
  axes(ah(1))
  colormap jet
  plotonsphere(d70,0.1); shading faceted
  axes(ah(2))
  plotonsphere(d72,0.1); shading faceted
  axes(ah(3))
  plotonsphere(d74,0.1); shading faceted
  axes(ah(4))
  plotonsphere(d76,0.1); shading faceted
  % axes(ah(5))
  % plotonsphere(d77,0.1); shading faceted
  fig2print(gcf,'landscape')
  figdisp('Psevens')
  disp('User PAINTERS')
elseif strcmp(data,'demo2')
  clf
  % [E,V,Lmax,TH,C,K,V0,unc,sqz]=wieczorek(40,5,0,180/5,1);
  [E,V]=sdwcap(40,22,0,180/5);
  ah=krijetem(subnum(1,5));
  axes(ah(1))
  plotonsphere(repmat(E(:,1),1,length(E(:,1))*2),0.2); shading faceted
  axes(ah(2))
  plotonsphere(repmat(E(:,2),1,length(E(:,1))*2),0.2); shading faceted
  axes(ah(3))
  plotonsphere(repmat(E(:,3),1,length(E(:,1))*2),0.2); shading faceted
  axes(ah(4))
  plotonsphere(repmat(E(:,4),1,length(E(:,1))*2),0.2); shading faceted
  axes(ah(5))
  plotonsphere(repmat(E(:,5),1,length(E(:,1))*2),0.2); shading faceted
  fig2print(gcf,'landscape')
  figdisp('Psevens')
  disp('User PAINTERS')
elseif strcmp(data,'demo3')
  wat=2;
  [data,ah,cb,t,lmcosif]=gravifilt([2 40],[],0,wat,[]);
  clf
  plotonsphere(data,0.1,[15 15],1)
  caxis([-80 80])
  %view(100,30)
  colormap jet
  hold on
  axes('position',[0,0,1,1]);  % Define axes for the text.
  htext=text(0.52,0.92,...
	     ['Free-air anomaly filtered between L = 2 and 40'],...
	     'FontSize',18);
  set(htext,'HorizontalAlignment','center');
  set(gca,'Visible','off');
  [cb,xcb]=addcb('hor',[-80 80],[-80 80],'jet',20);
  shrink(cb); movev(cb,0.05);
  moveh(cb,.02)
  cb.XLabel.String = '[mgal]';
  cb.FontSize = 16;
  myf=gcf; myf.InvertHardcopy='off';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grd=reduntest(grd)
% Tests if last longitude repeats last (0,360)
% and removes last data column
if sum(abs(grd(:,1)-grd(:,end))) >= size(grd,2)*eps*10
  disp(sprintf('Data violate wrap-around by %8.4e',...
	       sum(abs(grd(:,1)-grd(:,end)))))
end
grd=grd(:,1:(end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function polestest(grd)
% Tests if poles (-90,90) are identical over longitudes 
var1=var(grd(1,:));
var2=var(grd(end,:));
if var1>eps*10 | var2>eps*10
  disp(sprintf('Poles violated by %8.4e and %8.4e',var1,var2))
end
