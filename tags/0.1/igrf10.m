function lmcosi=igrf10(yr,yir)
% lmcosi=IGRF10(yr,yir)
%
% Interface to load the International Geomagnetic Reference Field
% and pass it on to other subroutines.
% 
% INPUT:
%
% yr        The year of interest (out of 1900:5:2005) [default: 2005]
%           OR: a string with the demo numbre
% yir       The year of interest if the first argument is a demo string
%
% OUTPUT:
%
% lmcosi    The tradition matrix with the ordered real coefficients
%           for the potential - however, in units of nT (nanoTesla)
%
% EXAMPLE:
%
% igrf10('demo1')
% igrf10('demo2') % The radial non-dipolar field, in nanoTesla
% igrf10('demo3') % The radial field, contoured, in nanoTesla
% igrf10('demo4') % The radial non-dipolar field, contoured, in nanoTesla
% igrf10('demo5') % The radial non-dipolar field
%
% SEE ALSO: 
%
% PLM2MAG
%
% Last modified by fjsimons-at-alum.mit.edu, 10/11/2010

% See plates at: http://pubs.usgs.gov/sim/2007/2964/

defval('yr',2005)

if ~isstr(yr)
  
 % Make it a single year - perhaps fix later
  if prod(size(yr))~=1
    error('Only a single year at the time can be requested for now')
  end
  
  % Open file
  fid=fopen(fullfile(getenv('IFILES'),...
		     'EARTHMODELS','IGRF-10','igrf10coeffs.txt'));
  
  % Define formats
  fmt1=['%s %s %s' repmat('%n',1,22) '%s'];
  fmt2=['%s' repmat('%n',1,25)];
  
  % The maximum expansion
  lmax=[repmat(10,1,20) 13 13];

  % Read the first line % TEXTSCAN better than TEXTREAD
  d=textscan(fid,fmt1,1);

  % Read the rest - supply zeroes where unavailable
  e=textscan(fid,fmt2,'emptyvalue',0);

  % Close file
  fclose(fid);

  % Available years
  years=[d{4:25}];
  % Spherical harmonic degrees
  EL=e{2};
  % Spherical harmonic orders
  EM=e{3};

  % Secular variation
  SV=e{26}(~isnan(e{26}));

  % The actual field expansion coefficients
  % Work from the back - you've got to know what's going on
  cosi=[e{4:25}];
  % Work-around the known knowns and the known unknowns
  ilmax=addmoff(lmax)-1;
  % This belongs with the last two
  cosi(ilmax(1)+1:end,21)=cosi(ilmax(1)+1:end,1);
  cosi(ilmax(1)+1:end,22)=cosi(ilmax(2)+1:end,2);
  cosi(ilmax(1)+1:end,1)=0;
  cosi(ilmax(2)+1:end,2)=0;

  % Now extract the data
  [C,iy]=intersect(years,yr);
  prepar=cosi(:,iy);
  % Stick in the non-existing zeros for degree and order zero
  prepar=[zeros(1,size(prepar,2)) ; prepar];

  % Reordering sequence
  [dems,dels,mz,lmcosi,mzi,mzo,bigm,bigl,rinm,ronm,demin]=...
      addmon(max(lmax));
  % What we have from the file is, effectively
  % [bigl(2:end) bigm(2:end)]
  % and what we want is lmcosi with the coefficients in the right position 
  % This is the output
  lmcosi(mzo+2*size(lmcosi,1))=prepar;
elseif strcmp(yr,'demo1')
  clf
  yir=2005;
  h=igrf10(yir);

  % Change Schmidt to full normalization for use in PLOTPLM
  % This converts the COEFFICIENTS to be multiplied with Schmidt to the 
  % COEFFICIENTS to be multiplied with 4pi-normalized harmonics, which
  % are Schmidt*sqrt(2l+1), see PLM2XYZ and note the TYPO in Blakely.
  h(:,3:4)=h(:,3:4)./repmat(sqrt(2*h(:,1)+1),1,2);
  % Make sure it is the RADIAL component of this at the surface
  h(:,3:4)=repmat(h(:,1)+1,1,2).*h(:,3:4);
  
  d=plotplm(h,[],[],4,1);
  longticks(gca,2)
  
  t(1)=title(sprintf('IGRF-10 magnetic field, year %i, degrees %i-%i',yir,...
		h(min(find(h(:,3))),1),h(end,1)));
  movev(t,5)
  
  cb=colorbar('hor');
  shrink(cb,2,2)
  axes(cb)
  longticks(cb,2)
  xlabel('radial component (nT)')
  
  movev(cb,-.1)
  
  fig2print(gcf,'portrait')
  figna=figdisp([],sprintf('%s-%i',yr,yir),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna));
  
elseif strcmp(yr,'demo2')
  defval('yir',2005);
  h=igrf10(yir);

  % Change Schmidt to full normalization for use in PLOTPLM
  h(:,3:4)=h(:,3:4)./repmat(sqrt(2*h(:,1)+1),1,2);

  % The nondipole field, as Blakely p170 and eq. (8.20)
  h(1:3,3:4)=0;
  
  % Make sure it is the RADIAL component of this at the surface
  h(:,3:4)=repmat(h(:,1)+1,1,2).*h(:,3:4);
  
  clf
  d=plotplm(h,[],[],4,1);
  axis image
  longticks(gca,2)
  t(1)=title(sprintf('IGRF-10 magnetic field, year %i, degrees %i-%i',yir,...
		h(min(find(h(:,3))),1),h(end,1)));
  movev(t,5)

  cb=colorbar('hor');
  shrink(cb,2,2)
  axes(cb)
  longticks(cb,2)
  xlabel('non-dipolar radial component (nT)');  
  movev(cb,-.1)
  
  fig2print(gcf,'portrait')
  figna=figdisp([],sprintf('%s-%i',yr,yir),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna));
elseif strcmp(yr,'demo3')
  clf
  yir=2005;
  h=igrf10(yir);
  
  % Change Schmidt to full normalization for use in PLOTPLM
  h(:,3:4)=h(:,3:4)./repmat(sqrt(2*h(:,1)+1),1,2);
  
  % Make sure it is the RADIAL component of this at the surface
  h(:,3:4)=repmat(h(:,1)+1,1,2).*h(:,3:4);

  d=plotplm(h,[],[],4,1); clf

  lons=linspace(0,360,size(d,2));
  lats=linspace(-90,90,size(d,1));
  
  % Don't forget to flip up down for contouring!
  [c,hh]=contour(lons,lats,...
		flipud(d),[-65000:5000:-5000]); 
  set(hh,'EdgeC','r')
  hold on
  [c,hh]=contour(lons,lats,...
		flipud(d),[5000:5000:65000]); 
  set(hh,'EdgeC','b')
  [c,hh]=contour(lons,lats,...
		flipud(d),[0 0]); 
  set(hh,'EdgeC','k','LineW',2)
  
  plotcont; axis image; ylim([-90 90])
  defval('dlat',45)
  set(gca,'ytick',[-90:dlat:90])
  set(gca,'xtick',[0:90:360])
  deggies(gca)

  longticks(gca,2)
  t(1)=title(sprintf('IGRF-10 magnetic field, year %i, degrees %i-%i',yir,...
		     h(min(find(h(:,3))),1),h(end,1)));
  movev(t,5)

  xl=xlabel(sprintf('minimum %i nT ; maximum %i nT',round(min(d(:))), ...
		  round(max(d(:)))));
  movev(xl,-10)

  fig2print(gcf,'portrait')
  figna=figdisp([],sprintf('%s-%i',yr,yir),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna));

elseif strcmp(yr,'demo4')
  clf
  yir=2005;
  h=igrf10(yir);
  
  % The nondipole field, as Blakely p170 and eq. (8.20)
  h(1:3,3:4)=0;

  % Change Schmidt to full normalization for use in PLOTPLM
  h(:,3:4)=h(:,3:4)./repmat(sqrt(2*h(:,1)+1),1,2);
  
  % Make sure it is the RADIAL component of this at the surface
  h(:,3:4)=repmat(h(:,1)+1,1,2).*h(:,3:4);

  d=plotplm(h,[],[],4,1); clf
  lons=linspace(0,360,size(d,2));
  lats=linspace(-90,90,size(d,1));

  clf
  % Don't forget to flip up down for contouring!
  [c,hh]=contour(lons,lats,...
		flipud(d),[-20000:1000:-1000]); 
  set(hh,'EdgeC','r')
  hold on
  [c,hh]=contour(lons,lats,...
		flipud(d),[1000:1000:20000]); 
  set(hh,'EdgeC','b')
  [c,hh]=contour(lons,lats,...
		flipud(d),[0 0]); 
  set(hh,'EdgeC','k','LineW',2)
  
  plotcont; axis image; ylim([-90 90])
  defval('dlat',45)
  set(gca,'ytick',[-90:dlat:90])
  set(gca,'xtick',[0:90:360])
  deggies(gca)

  longticks(gca,2)
  t(1)=title(sprintf('IGRF-10 magnetic field, year %i, degrees %i-%i',yir,...
		     h(min(find(h(:,3))),1),h(end,1)));
  movev(t,5)

  xl=xlabel(sprintf('minimum %i nT ; maximum %i nT',round(min(d(:))), ...
		  round(max(d(:)))));
  movev(xl,-10)

  fig2print(gcf,'portrait')
  figna=figdisp([],sprintf('%s-%i',yr,yir),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna));
elseif strcmp(yr,'demo5')
  clf
  defval('yir',2005);
  h=igrf10(yir);
  
  % The nondipole field, as Blakely p170 and eq. (8.20)
  h(1:3,3:4)=0;

  % Change Schmidt to full normalization for use in PLOTPLM
  h(:,3:4)=h(:,3:4)./repmat(sqrt(2*h(:,1)+1),1,2);
  
  % Make sure it is the RADIAL component of this at the surface
  h(:,3:4)=repmat(h(:,1)+1,1,2).*h(:,3:4);

  d=plotplm(h,[],[],4,1); 
  lons=linspace(0,360,size(d,2));
  lats=linspace(-90,90,size(d,1));

  % Don't forget to flip up down for contouring!
  [c,hh]=contour(lons,lats,...
		flipud(d),[-20000:2000:-2000]); 
  set(hh,'EdgeC','r')
  hold on
  [c,hh]=contour(lons,lats,...
		flipud(d),[2000:2000:20000]); 
  set(hh,'EdgeC','b')
  [c,hh]=contour(lons,lats,...
		flipud(d),[0 0]); 
  set(hh,'EdgeC','k','LineW',2)
  
  plotcont; axis image; ylim([-90 90])
  defval('dlat',45)
  set(gca,'ytick',[-90:dlat:90])
  set(gca,'xtick',[0:90:360])
  deggies(gca)

  longticks(gca,2)
  t(1)=title(sprintf('IGRF-10 magnetic field, year %i, degrees %i-%i',yir,...
		     h(min(find(h(:,3))),1),h(end,1)));
  movev(t,5)

  xl=xlabel(sprintf('minimum %i nT ; maximum %i nT',round(min(d(:))), ...
		  round(max(d(:)))));
  movev(xl,-10)

  fig2print(gcf,'portrait')
  figna=figdisp([],sprintf('%s-%i',yr,yir),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna));
end


