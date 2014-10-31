function varargout=boxcaptarget(co,TH,L,degres,ms)
% [r,lon,lat,N,Nm,V,C]=BOXCAPTARGET(co,TH,L,degres,ms)
%
% Targets a region with a set of axisymmetric of boxcar tapers
%
% INPUT
%
% co      [lon lat] coordinates of the center, in degrees, default [0 90]
% TH      Radius of the window [default: 180/sqrt(L(L+1))]
% L       Bandwidth of the window
% degres  The degree resolution of the output field [default 1]
% V        The "eigenvalue" - in this case the energy leakage parameter
% ms       0 Coefficients normalized to unit power [default]
%          1 Coefficients normalized according to Mark Simons
%
% OUTPUT:
%
% r       A spatial field with the taper right where you want it
%         in a matrix where the third dimension is the taper number
% lon     The longitudes at which r is evaluated, in degrees
% lat     The latitudes at which r is evaluated, in degrees
% N       The full Shannon number for this problem
% Nm      The partial Shannon number for this problem
% C       The spherical harmonic coefficients of the polar window
%
% EXAMPLE:
% 
% boxcaptarget('demo')
%
% Mark Simons, GJI 1996, Figure 6 (Note his RMS divides by 2l+1)
%
% plot(0:180,kindeks(boxcaptarget([],[],4,[],1),1)/4); hold on
% plot(0:180,kindeks(boxcaptarget([],[],8,[],1),1)/4);
% plot(0:180,kindeks(boxcaptarget([],[],16,[],1),1)/4); 
% axis([0 90 -10 150]); grid on; hold off
%
% Last modified by fjsimons-at-alum.mit.edu, 04.05.2005

if ~isstr(co) 
  defval('co',[0 90])
  defval('L',18)
  defval('TH',180/sqrt(L*(L+1)))
  defval('degres',1)
  defval('ms',0)
  
  disp(sprintf('TH= %8.3f',TH))
  
  % Calculate the tapers centered on the pole
  [jk,C,V]=boxcap(TH,L,ms);
  
  % Calculate the full Shannon number
  N=(L+1)^2*(1-cos(TH/180*pi))/2;

  % Calculate the asymptotic partial Shannon number
  [Nm,Nsum]=nsubm(N,L,1,L);
  
  % This is for the output only
  Nm=Nm(1);

  % Make a blank array of indices and coefficients
  [em,el,mzero]=addmon(L);
  Cb=[el em repmat(0,length(em),2)];

  % Make a blank array 
  r=repmat(NaN,[180/degres+1 360/degres+1]);

  % Find the appropriate rotation angles
  alpha=0; % Around z axis
  beta=90-co(2); % To a desired colatitude
  gamma=co(1); % To a desired longitude

  % Put the zonal coefficients in the right place
  Cb(mzero,3)=C;
  % Rotate the tapers to a position on the sphere
  Cr=plm2rot(Cb,alpha,beta,gamma);
  
  % Calculate the spatial functions
  [r,lon,lat]=plm2xyz(Cr,degres);

  % Output
  varn={'r','lon','lat','N','Nm','V','C'};
  for index=1:nargout
    varargout{index}=eval(varn{index});
  end
else
  col=round(rand*180); 
  lon=round(rand*360); 
  TH=20; L=12;
  r=boxcaptarget([lon 90-col],TH,L);
  clf
  imagef([],[],r)
  set(gca,'xtick',[0 lon 360],'xtickl',[0 lon 360])
  set(gca,'ytick',[-90 90-col 90],'ytickl',[-90 90-col 90])
  title(sprintf('Lon= %i; Lat= %i',lon,90-col));
  longticks(gca)
  set(gca,'xgrid','on','ygrid','on')
end

