function varargout=sdwtarget(co,TH,L,m,degres,SNm)
% [r,lon,lat,N,Nm,V,C]=SDWTARGET(co,TH,L,m,degres,SNm)
%
% Targets a region with a set of SINGLE-ORDER Slepian tapers
%
% INPUT
%
% co      [lon lat] coordinates of the center [degrees]
% TH      Radius of the window [degrees]
% L       Spherical-harmonic bandwidth of the window
% m       Angular order of the window, -L<=m<=L
% degres  The degree resolution of the output field [default 1]
% SNm     The number of tapers you require [defaulted]
%
% OUTPUT:
%
% r       A spatial field with the taper right where you want it
%         in a matrix where the third dimension is the taper number
% lon     The longitudes at which r is evaluated, in degrees
% lat     The latitudes at which r is evaluated, in degrees
% N       The full Shannon number for this problem
% Nm      The partial Shannon number for this problem
% V       The calculated eigenvalues belonging to this set of tapers
% C       The spherical harmonic coefficients of the polar windows
%
% EXAMPLE:
% 
% sdwtarget('demo1')
%
% SEE ALSO:
%
% PTOSLEP, PLM2ROT
%
% Last modified by fjsimons-at-alum.mit.edu, 11/19/2010

if ~isstr(co) 
  
  defval('co',[0 0])
  defval('TH',20)
  defval('L',18)
  defval('m',0)

  defval('degres',1)

  if m>L
    error('Order cannot exceed degree')
  end
  
  % Calculate the taper coefficients centered on the pole
  [E,Vg,th,C,T,V]=grunbaum(TH,L,m,0);
  
  % Calculate the full Shannon number
  N=(L+1)^2*(1-cos(TH/180*pi))/2;

  % Calculate the asymptotic partial Shannon number
  [Nm,Nsum]=nsubm(N,L,1,L);
  
  % Supply a default number of required tapers
  defval('SNm',round(Nm(abs(m)+1)))
  
  % This is for the output only
  Nm=Nm(abs(m)+1);

  % Make a blank array of indices and coefficients
  [em,el,mzero]=addmon(L);
  Cb=[el em repmat(0,length(em),2)];

  if SNm>0
    % Make a blank window array 
    r=repmat(NaN,[180/degres+1 360/degres+1 SNm]);
    
    % Find the appropriate rotation angles, see PLM2ROT
    alpa=0;       % Around z axis
    % sign flip 03/11/2010
    bita=co(2)-90; % To a desired colatitude
    % sign flip 03/11/2010
    gama=-co(1);   % To a desired longitude
    % which is identical to bita=90-co(2) and gama=180-co(1)
    
    % Do the following for every taper
    for index=1:SNm
      % Put the zonal coefficients in the right place
      Cb(mzero(abs(m)+1:end)+abs(m),3+(m<0))=C(:,index);
      
      % Rotate the tapers to a position on the sphere
      Cr=plm2rot(Cb,alpa,bita,gama);
      
      % Calculate the spatial functions
      [r(:,:,index),lon,lat]=plm2xyz(Cr,degres);
      
      % Adjust the sign
      if max(max(abs(r(:,:,index))))==abs(min(min(r(:,:,index))))
	r(:,:,index)=-r(:,:,index);
      end
    end
  else
    warning('You are getting nothing because your Shannon number is too low')
    [r,lon,lat]=deal([]);
  end
    
  % Output
  varn={r,lon,lat,N,Nm,V,C};
  varargout=varn(1:nargout);
elseif strcmp(co,'demo1')
  % Specify the centers
  col=round(rand*180);
  lon=round(rand*360); 
  % Specify the other parameters
  TH=20;
  L=round(36*rand(1));
  m=round(L*rand(1));
  r=sdwtarget([lon 90-col],20,L,m,1,4);
  % Calculate window location
  [lon2,lat2]=caploc([lon 90-col],TH);
  clf
  [ah,ha,H]=krijetem(subnum(2,2));
  for index=1:4
    axes(ah(index))
    imagef([],[],r(:,:,index)); axis image
    set(gca,'xtick',[0 lon 360],'xtickl',[0 lon 360])
    set(gca,'ytick',[-90 90-col 90],'ytickl',[-90 90-col 90])
    hold on
    plot(lon2,lat2,'k-')
  end
  movev(ah(1:2),-0.025)
  longticks(ah)
  set(ah,'xgrid','on','ygrid','on')
  serre(H',[],'down')
  supertit(ah(1:2),sprintf('L = %i ; m = %i ; lon= %i; lat= %i',...
			   L,m,lon,90-col));
  fig2print(gcf,'landscape')
  figdisp
end
