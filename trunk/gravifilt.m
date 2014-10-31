function varargout=gravifilt(L,wlen,pono,wat,lmcosi,degres)
% [data,ah,cb]=gravifilt(L,wlen,pono,wat,lmcosi,degres)
%
% Makes filtered map of the EGM2008 geoid expansion
%
% INPUT:
%
% L         Bandpass corner frequency [default: [20 50]]
% wlen      Filter length [default: 5]
% pono      0 Data will not be plotted
%           1 Data will be plotted [default]
% wat       1 gravitational potential [J/kg]
%           2 free-air gravity anomaly [m/s^2]
%           3 approximate geoid anomaly [m]
% lmcosi    Coefficients to be filtered [loaded if not supplied]
% degres    Grid cell size [default: 1/4 degree]
%
% OUTPUT:
%
% data      The data on a grid with requested resolution
% ah,cb,t   Axis handles to graph, color bar and title
% lmcosi    The filtered spherical harmonics coefficients
%
% SEE ALSO: TOPOFILT
%
% EXAMPLE:
%
% gravifilt('demoX'), where X=1,2,3,4,5
% gravifilt('demo6',[2 100])
%
% Last modified by fjsimons-at-alum.mit.edu, 02/20/2012

% Supply defaults
defval('L',[20 50])

if ~isstr(L)
  defval('wlen',5)
  defval('degres',1/4)
  defval('pono',1)
  defval('wat',3)
  
  % If the variable doesn't exist
  if exist('lmcosi','var')~=1
    % Read the file
    fnpl=fullfile(getenv('IFILES'),...
		  'EARTHMODELS','EGM2008','EGM2008_zerotide_720.mat');
    % If the file doesn't exist
    if exist(fnpl,'file')~=2
      % Make the file
      lmcosi=fralmanac('EGM2008_ZeroTide','SHM');
      % Restrict it to more manageable size
      lmcosi=lmcosi(1:addmup(720)-addmup(lmcosi(1)-1),:);
      save(fnpl,'lmcosi')
    else
      load(fnpl)
    end
  end
  defval('percs',[0 5 25 50 75 95 100]);

  % Change to [3] geoidal coefficients in m
  % Change to [2] free-air gravity coefficients in mgal
  lmcosi=plm2pot(lmcosi,[],[],[],wat);
  if wat==2
    lmcosi(:,3:4)=lmcosi(:,3:4)*1e5;
    legs='EGM2008 free-air gravity anomaly [mgal]';
  elseif wat==3
    legs='EGM2008 geoidal undulation [m]';
  end

  % Perform the filtering
  lmcosif=plmfilt(lmcosi,L,wlen);

  % Plot the spectrum to be sure
  clf
  [sdl,l]=plm2spec(lmcosif);
  semilogy(l,sdl)
  set(gca,'xtick',unique([0 1 2 L]),'xgrid','on')
  pause(1)
    
  % Perform the expansion and/or plot
  if pono==0
    data=plm2xyz(lmcosif,degres);
    [ah,cb,t]=deal(NaN);
  else
    clf
    ah=gca;
    
    % Plot the topography
    axes(ah(1))
    data=plotplm(lmcosif,[],[],4,degres);
    if length(L)==2
      t(1)=title(sprintf(...
	  '%s filtered between L = %i and %i',...
	  legs,L(1),L(2)));
    else
      t(1)=title(sprintf(...
	  '%s filtered to below L = %i',legs,L));
    end
    cb(1)=colorbar('hor');
    vmd=prctile(data(:),percs);
    tix=unique([0 vmd([1 find(percs==50) length(vmd)])]);
    set(cb(1),'xtick',tix,'xtickl',round(tix),...
	      'xlim',vmd([1 length(vmd)]))
    % Cosmetics
    fig2print(gcf,'landscape')
    movev(cb,-0.075)
    movev(t,0.075)
    set(t,'FontS',15)
    longticks([ah cb],2)
    if length(L)==1
      figdisp([],sprintf('%i_%i',wat,L(1)))
    else
      figdisp([],sprintf('%i_%i_%i',wat,L(1),L(2)))
    end
  end

  % Output only if desired
  vars={data,ah,cb,t,lmcosif};
  varargout=vars(1:nargout);
elseif strcmp(L,'demo1')
  c11cmn=[175 30 220 10];
  wat=3;
  [data,ah,cb,t]=gravifilt([2 180],[],[],wat);
  axes(ah); pause(3); axis(c11cmn([1 3 4 2]))
  shrink(ah); caxis([-20 20]); delete(cb)
  cb=colorbar('hor'); undeggies(ah)
  set(ah,'xtick',c11cmn([1 3]),'xtickl',c11cmn([1 3]),...
	 'ytick',c11cmn([4 2]),'ytickl',c11cmn([4 2]))
  deggies(ah); shrink(cb); movev(cb,0.1);
  axes(cb); t=xlabel(get(t,'string'));
elseif strcmp(L,'demo2')
  c11cmn=[175 30 220 10];
  wat=3;
  [data,ah,cb,t]=gravifilt([20 60],[],[],wat);
  axes(ah); pause(3); axis(c11cmn([1 3 4 2]))
  shrink(ah); caxis([-5 5]); delete(cb)
  cb=colorbar('hor'); undeggies(ah)
  set(ah,'xtick',c11cmn([1 3]),'xtickl',c11cmn([1 3]),...
	 'ytick',c11cmn([4 2]),'ytickl',c11cmn([4 2]))
  deggies(ah); shrink(cb); movev(cb,0.1);
  axes(cb); t=xlabel(get(t,'string'));
elseif strcmp(L,'demo3')
  c11cmn=[175 30 220 10];
  wat=2;
  [data,ah,cb,t]=gravifilt([2 360],[],[],wat);
  axes(ah); pause(3); axis(c11cmn([1 3 4 2]))
  shrink(ah); caxis([-60 60]); delete(cb)
  cb=colorbar('hor'); undeggies(ah)
  set(ah,'xtick',c11cmn([1 3]),'xtickl',c11cmn([1 3]),...
	 'ytick',c11cmn([4 2]),'ytickl',c11cmn([4 2]))
  deggies(ah); shrink(cb); movev(cb,0.1);
  axes(cb); t=xlabel(get(t,'string'));
elseif strcmp(L,'demo4')
  c11cmn=[175 30 220 10];
  wat=2;
  [data,ah,cb,t]=gravifilt([20 60],[],[],wat);
  axes(ah); pause(3); axis(c11cmn([1 3 4 2]))
  shrink(ah); caxis([-40 40]); delete(cb)
  cb=colorbar('hor'); undeggies(ah)
  set(ah,'xtick',c11cmn([1 3]),'xtickl',c11cmn([1 3]),...
	 'ytick',c11cmn([4 2]),'ytickl',c11cmn([4 2]))
  deggies(ah); shrink(cb); movev(cb,0.1);
  axes(cb); t=xlabel(get(t,'string'));
elseif strcmp(L,'demo5')
  % This inspired by Chase (1985), Annual Reviews
  c11cmn=[0 90 360 -90];
  wat=3;
  [data,ah,cb,t]=gravifilt(20,[],[],wat);
  movev(t,10)
  axes(ah)
  hold on
  lons=linspace(c11cmn(1),c11cmn(3),size(data,2));
  lats=linspace(c11cmn(4),c11cmn(2),size(data,1));
  contours=[-40 0 40];
  stile={'--',':','-'};
  cols={'b','k','r'};
  for index=1:length(contours)
    [c,hh]=contour(lons,lats,...
		   flipud(data),[contours(index) contours(index)]); 
    set(hh,'EdgeC',cols{index},'LineW',2,'LineS',stile{index})
  end
  hold off
  axes(cb)
  hold on
  for index=1:length(contours)
    plot([contours(index) contours(index)],ylim,...
	 'Color',cols{index},'LineW',2,'LineS','-')
  end
  newx=unique([contours get(cb,'xtick')]);
  set(cb,'xtick',newx,'xtickl',round(newx))
  hold off
elseif strcmp(L,'demo6')
  % Spherical harmonics range
  defval('wlen',[])
  ELS=wlen;
  defval('ELS',[2 100])
  % Which property for PLM2POT
  wat=2;
  % Calculate, don't plot
  [~,~,~,~,lmcosi]=gravifilt(ELS,[],0,wat);
  % Rendering resolution
  degres=1/4;
  % Color range
  cax=[-50 50];
  % Color scheme
  cosch='kelicol';
  % String for color bar
  cstr='free-air gravity anomaly (mgal)';
  
  clf
  [data,ch,ph,lon,lat]=plotplm(lmcosi,[],[],2,degres,[],[],cax);
  % The true center of the world
  lat=50.844850;
  lon=4.349754;
  % Set view angles ahead of time as an explicit longitude and latitude
  [xv,yv,zv]=sph2cart(lon*pi/180,lat*pi/180,1);
  % Set, then verify
  view([xv,yv,zv]); [AZ,EL]=view;
  % Use, but invert colorbar, keep the inversion for ADDCB
  colormap(flipud(eval(cosch))); invt=1;

  % Title
  t=title(sprintf('EGM2008 between l= %i and %i',ELS(1),ELS(2)));
  
  % Color scheme and bar
  [cb,xcb]=addcb('hor',cax,cax,cosch,25,invt);
  set(xcb,'string',cstr)
  shrink(cb,3,2)
  % Cosmetics
  fig2print(gcf,'portrait')
  % Weird, but let it go
  moveh(t,.5)
  movev(cb,.1)
  % Print command
  popts='-zbuffer'',''-r300';
  if length(ELS)==1
    figna=figdisp([],sprintf('%i_%i',wat,ELS(1)),popts,1);
  else
    figna=figdisp([],sprintf('%i_%i_%i',wat,ELS(1),ELS(2)),popts,1);
  end
  system(sprintf('epstopdf %s.eps',figna));
end
