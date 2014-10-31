function varargout=topofilt(L,wlen,pono,lmcosi)
% [data,data2]=topofilt(L,wlen,pono,lmcosi)
%
% Makes filtered map of the EGM2008 topography expansion
%
% INPUT:
%
% L         Lowpass corner frequency [default: 180]
% wlen      Filter length [default: 5]
% degres    Grid cell size [default: 1/4 degree]
% pono      0 Data will not be plotted
%           1 Data will be plotted [default]
% lmcosi    Coefficients to be filtered [loaded if not supplied]
%
% OUTPUT:
%
% data      The data on a grid with requested resolution
% data2     The norm of the gradient of the data
%
% SEE ALSO: GRAVIFILT
% 
% Last modified by fjsimons-at-alum.mit.edu, 12/06/2009

% Supply defaults
defval('L',180)
defval('wlen',5)
defval('degres',1/4)
defval('pono',1)
% If the variable doesn't exist
if exist('lmcosi','var')~=1
  % Read the file
  fnpl=fullfile(getenv('IFILES'),...
		'EARTHMODELS','EGM2008','EGM2008_topo_720.mat');
  % If the file doesn't exist
  if exist(fnpl,'file')~=2
    % Make the file
    lmcosi=fralmanac('EGM2008_Topography','SHM');
    % Restrict it to more manageable size
    lmcosi=lmcosi(1:addmup(720),:);
    save(fnpl,'lmcosi')
  else
    load(fnpl)
  end
end
defval('percs',[0 5 25 50 75 95 100]);

% Perform the filtering
lmcosif=plmfilt(lmcosi,L,wlen);

% Perform the expansion and/or plot
if pono==0
  data=plm2xyz(lmcosif,degres);
else
  clf
  [ah,ha,H]=krijetem(subnum(2,1));
  % Plot the topography
  axes(ah(1))
  data=plotplm(lmcosif,[],[],4,degres);
  t(1)=title(sprintf(...
      'Earth''s topography and bathymetry filtered to L = %i',L));
  cb(1)=colorbar('hor');
  vmd=prctile(data(:),percs);
  tix=unique([0 vmd([1 find(percs==50) length(vmd)])]);
  set(cb(1),'xtick',tix,'xtickl',round(tix),'xlim',vmd([1 length(vmd)]))

  % Calculate the norm of the surface gradient, DT (B.167)
  lmcosifdelnorm=lmcosif;
  [dems,dels]=addmon(lmcosif(end,1));
  % Divide by radius of Earth in meters for scale
  R=6378136.3;
  lmcosifdelnorm(:,3)=lmcosif(:,3).*sqrt(dels.*(dels+1))/R*100;
  lmcosifdelnorm(:,4)=lmcosif(:,4).*sqrt(dels.*(dels+1))/R*100;
  % Plot norm of surface gradient which is now in percent
  axes(ah(2))
  data2=plotplm(lmcosifdelnorm,[],[],4,degres);
  t(2)=title(sprintf('Signed norm of the surface gradient (%s)','%'));
  cb(2)=colorbar('hor');
  vmd2=prctile(data2(:),percs);
  tix=unique([vmd2([1 find(percs==50) length(vmd2)])]);
  set(cb(2),'xtick',tix,'xtickl',round(tix*10)/10,...
	    'xlim',vmd2([1 length(vmd2)]))
end

% Cosmetics
fig2print(gcf,'tall')
movev(cb,-0.075)
movev(t,0.075)
set(t,'FontS',15)
longticks([ah cb],2)
figdisp([],L)

% Output only if desired
vars={data,data2};
varargout=vars(1:nargout);

