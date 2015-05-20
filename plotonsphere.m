function plotonsphere(data,rang)
% PLOTONSPHERE(data,rang)
%
% Plots data as topography and colors onto a sphere.
%
% INPUT:
%
% data     Standard 2D geographic data; i.e. the first column is the
%          Greenwich meridian in the x-z plane.
% rang     Exaggeration of topography compared to 1 [default: 0]
%
% See PLOTONEARTH, PLOTPLM
%
% EXAMPLE:
%
% plotonsphere('demo1')
% plotonsphere('demo2')
%
% Last modified by fjsimons-at-alum.mit.edu, 12/02/2008

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

  s=surface(-x,-y,z,'FaceColor','texture','Cdata',data);
  axis equal

  % Interpolate the data on the sphere's grid
  if nx>ny 
    data=interp2([1:nx],[1:ny]',data,linspace(1,nx,ny),[1:ny]');
  elseif nx<ny
    data=interp2([1:nx],[1:ny]',data,[1:nx],linspace(1,ny,nx)');
  end

  data=1+scale(data/max(abs(data(:))),[-1 1]*rang);

  % This is tricky business but it works
  x=-x.*data;      
  y=-y.*data;
  z=z.*data;
  set(s, 'xdata', x, 'ydata', y, 'zdata', z)
  view(140,30) % So x and y face at you
  set(gca, 'dataaspectratio', [1 1 1], 'cameraviewangle', 7)
  axis image

  shading flat

  axis off
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
  [E,V,Lmax,TH,C,K,V0,unc,sqz]=wieczorek(40,5,0,180/5,1);
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
