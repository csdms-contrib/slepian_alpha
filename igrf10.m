function lmcosi=igrf10(yr,yir)
% lmcosi=IGRF10(yr,yir)
%
% Interface to load the International Geomagnetic Reference Field
% and pass it on to other subroutines. Supplanted by IGRF which, to make
% it backwards compatible, calls this version when IGRF(10) is requested.
% 
% INPUT:
%
% yr        The year of interest (out of 1900:5:2005) [default: 2005]
%           OR: a string with the demo number
% yir       The year of interest if the first argument is a demo string
%
% OUTPUT:
%
% lmcosi    The tradition matrix with the ordered real coefficients
%           for the potential - however, in units of nT (nanoTesla)
%
% EXAMPLE:
%
% igrf10('demo1') % The radial field, in nanoTesla
% igrf10('demo2') % The radial non-dipolar field, in nanoTesla
% igrf10('demo3') % The radial field, only contoured, in nanoTesla
% igrf10('demo4') % The radial non-dipolar field, only contoured, in nanoTesla
% igrf10('demo5') % The radial non-dipolar field, also contoured
% igrf10('demo6') % The radial non-dipolar secular variation since 1900
%
% SEE ALSO: 
%
% PLM2MAG, IGRF
%
% Tested on 8.3.0.532 (R2014a)
%
% Last modified by fjsimons-at-alum.mit.edu, 03/17/2020

% See plates at: http://pubs.usgs.gov/sim/2007/2964/
% See /u/fjsimons/CLASSES/GEO371/2008/Images/igrf-10-Bv.gif

defval('yr',2005)

if ~isstr(yr)
  
 % Make it a single year - perhaps fix later
  if prod(size(yr))~=1
    error('Only a single year at the time can be requested for now')
  end
  
  % Open file
  fname=fullfile(getenv('IFILES'),'EARTHMODELS','IGRF-10','igrf10coeffs.txt');
  fid=fopen(fname);
  
  % Define formats
  fmt1=['%s %s %s' repmat('%n',1,22) '%s'];
  fmt2=['%s' repmat('%n',1,25)];

  % The maximum expansion for every year, you have to know what's going
  % on! Try determining it from the file itself, first, shunt otherwise
  %[~,NF]=system(sprintf('awk ''{print NF}'' %s',fname)); NF=str2num(NF);
  lmax=[repmat(10,1,20) repmat(13,1,2)];

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
  SV=e{26};

  %%% This applies when the unkowns are empties in the file %%%%%%%%%%
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
  if isempty(iy)
    error(sprintf('\n%s: Musty specify valid model year',upper(mfilename)));
  end
  prepar=cosi(:,iy);
  % Stick in the non-existing zeros for degree and order zero
  prepar=[zeros(1,size(prepar,2)) ; prepar];
  %%% This applied when the unkowns are empties in the file %%%%%%%%%%
  
  % Reordering sequence
  [dems,dels,mz,lmcosi,mzi,mzo,bigm,bigl,rinm,ronm,demin]=...
      addmon(max(lmax));
  % What we have from the file is, effectively
  % [bigl(2:end) bigm(2:end)]
  % and what we want is lmcosi with the coefficients in the right position 
  % This is the output
  lmcosi(mzo+2*size(lmcosi,1))=prepar;
  % Check that your inkling was right or you need to go back
  diferm(max(lmcosi(~~sum(lmcosi(:,3:4),2),1)),lmax(years==yr))
elseif strcmp(yr,'demo1')
  defval('yir',2005);
  h=igrf10(yir);

  % Plot and print
  plotandprint(h,yr,yir,0,0)
  
elseif strcmp(yr,'demo2')
  defval('yir',2005);
  h=igrf10(yir);

  % Plot and print
  plotandprint(h,yr,yir,1,0,[-20000:1000:-1000],[1000:1000:20000])

elseif strcmp(yr,'demo3')
  defval('yir',2005);
  h=igrf10(yir);
  
  % Plot and print
  plotandprint(h,yr,yir,0,1,[-65000:5000:-5000],[5000:5000:65000])

elseif strcmp(yr,'demo4')
  defval('yir',2005);
  h=igrf10(yir);

  % Plot and print
  plotandprint(h,yr,yir,1,1)
elseif strcmp(yr,'demo5')
  clf
  defval('yir',2005);
  h=igrf10(yir);
  
  % Plot and print
  plotandprint(h,yr,yir,1,2,[-20000:2000:-2000],[2000:2000:20000])
elseif strcmp(yr,'demo6')
  clf
  for yir=1900:5:2005
    h=igrf10(yir);
    % Plot and print
    plotandprint(h,yr,yir,1,2,[-20000:2000:-2000],[2000:2000:20000])
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction to plot and print 
function plotandprint(h,yr,yir,zro,cnt,negcont,poscont)
% Zero out the degree-1 component
defval('zro',0)
% Whether contours are plotted
defval('cnt',0)

% Change Schmidt to full normalization for use in PLOTPLM
% This converts the COEFFICIENTS to be multiplied with Schmidt to the 
% COEFFICIENTS to be multiplied with 4pi-normalized harmonics, which
% are Schmidt*sqrt(2l+1), see PLM2XYZ and note the TYPO in Blakely.
h(:,3:4)=h(:,3:4)./repmat(sqrt(2*h(:,1)+1),1,2);

if zro==1
  % The nondipole field, as Blakely p170 and p169 eq. (8.20)
  h(1:3,3:4)=0;
  xlab='non-dipolar radial component (nT)';
else
  xlab='radial component (nT)';
end

% Make sure it is the RADIAL component of this at the SURFACE
h(:,3:4)=repmat(h(:,1)+1,1,2).*h(:,3:4);

clf
% This resolution parameter will change the quoted maxima and minima
degres=1;
d=plotplm(h,[],[],4,degres);
kelicol

% The title string
ztit=sprintf('IGRF-10 magnetic field, year %i, degrees %i-%i',yir,...
		     h(min(find(h(:,3))),1),max(h(~~sum(h(:,3:4),2),1)));
switch cnt
  case 0
   % Just a color plot
  axis image
  longticks(gca,2)
  t(1)=title(ztit);
  movev(t,5)
  
  cb=colorbar('hor');
  shrink(cb,2,2)
  axes(cb)
  longticks(cb,2)
  xlabel(xlab)
  movev(cb,-.1)
 case {1,2}
  % A judicious contour plot
  if cnt==1
    % No overlay
    clf
  end
  
  % Negative and positive contour intervals
  defval('negcont',[-20000:1000:-1000])
  defval('poscont',[  1000:1000:20000])
  % Geographic grid
  lons=linspace(0,360,size(d,2));
  lats=linspace(-90,90,size(d,1));

  % Don't forget to flip up down for contouring!
  [c,hh]=contour(lons,lats,flipud(d),negcont); 
  set(hh,'EdgeC','r')
  hold on
  [c,hh]=contour(lons,lats,flipud(d),poscont); 
  set(hh,'EdgeC','b')
  [c,hh]=contour(lons,lats,flipud(d),[0 0]); 
  set(hh,'EdgeC','k','LineW',2)
  % Finalize
  plotcont; axis image; ylim([-90 90])
  defval('dlat',45)
  set(gca,'ytick',[-90:dlat:90])
  set(gca,'xtick',[0:90:360])
  if cnt==1
    % Otherwise you already had them
    deggies(gca)
  end
  longticks(gca,2)
  % Only quote the maximum degree where you actually have it
  t(1)=title(ztit);
  movev(t,5)
  xl=xlabel(sprintf('minimum %i nT ; maximum %i nT ; contour interval %i nT',...
                    round(min(d(:))),round(max(d(:))),...
                    unique([diff(negcont) diff(poscont)])));
  movev(xl,-10)
end

% Actual printing
fig2print(gcf,'portrait')
figna=figdisp('igrf10',sprintf('%s-%i',yr,yir),[],2);
% Maybe this...
% figna=figdisp([],sprintf('%s-%i',yr,yir),'-r300',1,'jpeg');
