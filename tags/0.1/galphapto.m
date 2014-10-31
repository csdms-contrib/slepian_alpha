function varargout=...
    galphapto(TH,L,phi0,theta0,omega,theta,phi,J,irr,Glma,V,N,EL,EM)
% [Gar,V,N,J,phi0,theta0,omega,theta,phi,TH,L,Glma,EL,EM]=...
%          GALPHAPTO(TH,L,phi0,theta0,omega,theta,phi,J,irr,Glma,V,N,EL,EM)
%
% This function returns an (alpha)X(r) matrix with the spatially expanded
% BANDLIMITED or PASSBAND Slepian functions of the SINGLE CAP as rotated
% to a desired location and azimuthally by a certain amount. The
% evaluation points can be a regular grid or an irregular collection of
% points. In the first case, the column dimension is
% length(theta)*length(phi) and in the second case, it is length(theta). 
% The normalization is as in Simons, Wieczorek and Dahlen, eq. (3.15)/(5.7)
%
% INPUT:
%
% TH       Angular extent of the spherical cap (degrees)
% L        Bandwidth (maximum angular degree) or passband (two degrees)
% phi0     Longitude of the center (degrees)
% theta0   Colatitude of the center (degrees)
% omega    Anticlockwise azimuthal rotation (degrees) [default: 0]
% theta    Colatitude vector (0<=theta<=pi) [default: 2TH around cap]
% phi      Longitude vector (0<=phi<=2*pi) [default: 2TH around cap]
%          Unless irr=1, we assume you mean a 2-D (theta,phi) grid.
%          But if irr=1, length(theta(:)) must equal length(phi(:)).
% J        How many of the best eigenfunctions do you want? [default: N]
% irr      0 Regular grid, no matter how you input theta and phi [default]
%          1 Irregular grid, input interpreted as distinct pairs of theta, phi
% Glma     The spectral eigenfunctions in case you already have them
% V        The spectral eigenvalues in case you already have them
% N        The Shannon number in case you already have it
% EL       The degrees in question if you already have them
% EM       The orders in question if you already have them
%
% OUTPUT:
% 
% Gar      The spatial eigenfunctions
% V        All concentration eigenvalues, sorted globally and descending
% N        The Shannon number
% J        The truncation level
% phi0     Longitude of the center (degrees)
% theta0   Colatitude of the center (degrees)
% omega    Anticlockwise azimuthal rotation (degrees)
% theta    The colatitudes that you put in or received (radians)
% phi      The longitudes that you put in or received (radians)
% TH       Angular extent of the spherical cap (degrees)
% L        Bandwidth (maximum angular degree) or passband (two degrees)
% Glma     The spectral eigenfunctions, for use in, e.g. PLOTSLEP
% EL       The degrees in question
% EM       The orders in question
%
% EXAMPLES:
%
% galphapto('demo1') % Plots spatial functions on a subgrid
% galphapto('demo2') % Extracts profiles from the spatial plots
% galphapto('demo3') % Plots bandpass spatial functions on a subgrid
% galphapto('demo4') % Plots an example from Harig et al. 
% galphapto('demo5') % Illustrates local sensitivity of pre-Shannon coefficients
%
% SEE ALSO:
%
% GALPHA, GLMALPHA, GLMALPHAPTO, PTOSLEP, SDWCAP
%
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012

% Supply default values
defval('TH',15)

if ~isstr(TH)

  defval('L',36)
  defval('phi0',15)
  defval('theta0',70)
  defval('omega',10)
  defval('theta',[theta0-2*TH:1:theta0+2*TH]*pi/180)
  defval('phi',[phi0-2*TH:1:phi0+2*TH]*pi/180)
  defval('irr',0)


  % Basic error check
  if irr==1 & ~all(size(theta)==size(phi))
    error('Input arrays must have the same dimensions for irregular grids')
  end
  
  % Figure out if it's lowpass or bandpass
  lp=length(L)==1;
  bp=length(L)==2;

  maxL=max(L);
  
  if exist('Glma')~=1 ||  exist('V')~=1 ||  exist('N')~=1 ...
	|| exist('EL')~=1 ||  exist('EM')~=1
    % Construct the rotated basis in the spectral domain
    [Glma,V,EL,EM,N]=glmalphapto(TH,L,phi0,theta0,omega);
  end
  
  % Default truncation is at the Shannon number
  defval('J',round(N))
  
  % Note that the transpose of the orthogonality is lost in case you no
  % longer have a full matrix...

  % Sort the tapers globally and potentially restrict their number
  % Note that GRUNBAUM individually will be inversely sorted, and that,
  % depending on numerical precision, the new sorting maybe slightly
  % counterintuitive. 
  [V,i]=sort(V,'descend');
  Glma=Glma(:,i); Glma=Glma(:,1:J); 
  % Don't cut these off but return them all, see SPIE2009_1 and SPIE2009_9
  % V=V(1:J);

  % Get the real spherical harmonics in the right places
  [Y,t,p]=ylm([0 maxL],[],theta,phi,[],[],[],irr);

  % Could, in here, also allow for the expansion using PLM2XYZ or
  % course. Should write routine to go back and forth Glma -> lmcosi
  % as sometimes, due to memory restrictions, this will be the only thing 
  
  % Perform the expansion, watching out for the phase factor
  % See the demos in GLMALPHAPTO for an alternative
  % This should be done in GALPHA also though it can wait
  Glmap=Glma.*repmat((-1).^EM,1,size(Glma,2));
  Gar=Glmap'*Y;
  
  % I'm now going to call the alternative GLM2LMCOSI
  
  % Provide output - with the original Glma!
  varns={Gar,V,N,J,phi0,theta0,omega,theta,phi,TH,L,Glma,EL,EM};
  varargout=varns(1:nargout);
elseif strcmp(TH,'demo1')
  TH=15;
  phi0=15;
  theta0=70;
  L=36; L=72;

%   TH=15;
%   phi0=35;
%   theta0=35;
%   L=36;
  
  [Gar,V,N,J,phi0,theta0,omega,theta,phi,TH,L,Glma]=...
      galphapto(TH,L,phi0,theta0);   

  % Make a decent plot
  clf
  [ah,ha,H]=krijetem(subnum(3,3));
  for index=1:length(ah)
    axes(ah(index))
    % Extract the data directly from the expansion
    data=reshape(Gar(index,:),length(theta),length(phi));
    c11cmn=[phi(1) pi/2-theta(1) phi(end) pi/2-theta(end)]*180/pi;
    imagefnan(c11cmn(1:2),c11cmn(3:4),setnans(data))
    ploco(c11cmn,theta0,phi0,TH)
    title(sprintf('%s = %8.3f','\lambda',V(index)))
  end
  % Cosmetic adjustments
  cosmo(ah,ha,H)
  % One more annotation
  axes(ha(6))
  xlabel(sprintf('%s = %i%s ; L = %i','\omega',...
		 omega,str2mat(176),L))
  figdisp([],'demo1')
elseif strcmp(TH,'demo2')
  TH=15;
  phi0=15;
  theta0=70;
  L=36;

%   TH=15;
%   phi0=35;
%   theta0=35;
%   L=36;

  % Try to extract the spatial eigenfunctions along a great circle going
  % through the center
  incr=linspace(-TH,TH,100)*fralmanac('DegDis')/1000;
  [lon2,lat2]=grcazim([phi0 90-theta0],incr,10);

  % Calculate the spatial eigenfunctions on this profile
  [Gar,V,N,J,phi0,theta0,omega,theta,phi,TH,L,Glma]=...
      galphapto(TH,L,phi0,theta0,[],[90-lat2]*pi/180,lon2*pi/180,[],1);
  
  % Make a decent plot
  clf
  [ah,ha,H]=krijetem(subnum(3,3));
  for index=1:length(ah)
    axes(ah(index))
    % Extract the data directly from the expansion
    plot(Gar(index,:));
    
    title(sprintf('%s = %8.3f','\lambda',V(index)))
    drawnow
  end
  figdisp([],'demo2')
elseif strcmp(TH,'demo3')
  TH=15;
  phi0=15;
  theta0=70;
  L=[17 72];

  [Gar,V,N,J,phi0,theta0,omega,theta,phi,TH,L,Glma]=...
      galphapto(TH,L,phi0,theta0);   

  % Make a decent plot
  clf
  [ah,ha,H]=krijetem(subnum(3,3));
  for index=1:length(ah)
    axes(ah(index))
    % Extract the data directly from the expansion
    data=reshape(Gar(index,:),length(theta),length(phi));
    c11cmn=[phi(1) pi/2-theta(1) phi(end) pi/2-theta(end)]*180/pi;
    imagefnan(c11cmn(1:2),c11cmn(3:4),setnans(data))
    ploco(c11cmn,theta0,phi0,TH)
    title(sprintf('%s = %8.3f','\lambda',V(index)))
  end
  % Cosmetic adjustments
  cosmo(ah,ha,H)
    % One more annotation
  axes(ha(6))
  xlabel(sprintf('%s = %i%s ; L = %i-%i','\omega',...
		 omega,str2mat(176),L(1),L(2)))
  figdisp([],'demo3')
elseif strcmp(TH,'demo4')
  TH=40;
  phi0=134;
  theta0=90--24;
  L=5;

  [Gar,V,N,J,phi0,theta0,omega,theta,phi,TH,L,Glma]=...
      galphapto(TH,L,phi0,theta0,0);   

  % Make a decent plot
  clf
  [ah,ha,H]=krijetem(subnum(2,2));
  for index=1:length(ah)
    axes(ah(index))
    % Extract the data directly from the expansion
    data=reshape(Gar(index,:),length(theta),length(phi));
    c11cmn=[phi(1) pi/2-theta(1) phi(end) pi/2-theta(end)]*180/pi;
    imagefnan(c11cmn(1:2),c11cmn(3:4),setnans(data))
    ploco(c11cmn,theta0,phi0,TH)
    title(sprintf('%s = %8.3f','\lambda',V(index)))
  end
  % Cosmetic adjustments
  cosmo(ah,ha,H)
    % One more annotation
  axes(ha(6))
  xlabel(sprintf('%s = %i%s ; L = %i-%i','\omega',...
		 omega,str2mat(176),L(1),L(2)))
  figdisp([],'demo4')
elseif strcmp(TH,'demo5')
  TH=20;
  % This must be commensurate with the number of samples, if not you're
  % going to have trouble inverting the matrix and making sense of it all
  L=18;
  phi0=120;
  theta0=60;
  omega0=0;
  % Some scattered points close to the region
  [lon2,lat2]=randpatch(250,TH,phi0,theta0);
  % Some more ones everywhere
  [lon1,lat1]=randsphere(750-length(lon2));
  % The complete data set
  lon=[lon1 ; lon2];
  lat=[lat1 ; lat2];
  % Split into the inside and outside parts
  [gcdkm,delta]=grcdist([phi0 90-theta0],[lon(:,1) lat(:,1)]);
  keepit=delta<=TH;
  init=find(keepit(:));
  outof=find(~keepit(:));
  % Generat the sampled Slepian matrix
  [Gar,V,N,J,phi0,theta0,omega,theta,phi,TH,L,Glma]=...
      galphapto(TH,L,phi0,theta0,0,(90-lat)*pi/180,lon*pi/180,(L+1)^2,1);
  % Reorder the columns to reflect the geometry
  Gar=Gar(:,[init(:)' outof(:)']);
  % Make the generalized inverse as if you were going to invert it
  Garinv=pinv(Gar');
  clf
  ah=krijetem(subnum(1,2));
  axes(ah(1))
  % Maybe truncate for color display?
  imagefnan([1 1],[size(Gar')],abs(Gar),...
	    [],halverange(abs(Gar),10)); grid on; axis ij normal
  ylabel(sprintf('rank %s','\alpha'))
  xlabel('spatial position')
  title(sprintf('G_{%s%s}','\alpha','r'))
  axes(ah(2))
  % Maybe truncate for color display?
  imagefnan([1 1],[size(Garinv')],abs(Garinv),...
	    [],halverange(abs(Garinv),10)); grid on; axis ij normal
  set(ah,'xtick',[1 sum(keepit) size(Gar,2)],...
	  'ytick',[1 round(N) (L+1)^2])
  ylabel(sprintf('rank %s','\alpha'))
  xlabel('spatial position')
  title(sprintf('G^{-1}_{%s%s}','\alpha','r'))
  longticks(ah,2)
end

% Some subroutines useful for the demos above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploco(c11cmn,theta0,phi0,TH)
% Plot the continents - note the Greenwich trick
yesorno=360*[c11cmn(1)<0 0];
[ax,f,lola1]=plotcont(c11cmn(1:2)+yesorno,...
		      c11cmn(3:4)+yesorno,[],-360);
[ax,f,lola2]=plotcont([0 c11cmn(2)],c11cmn(3:4));
axis(c11cmn([1 3 4 2]))
% Plot the circle of concentration
hold on
[lon2,lat2]=caploc([phi0 90-theta0],TH);
f=plot(lon2-360*[lon2>180],lat2,'k-');
% Plot the center
d=plot(phi0,90-theta0,'o','MarkerF','w','MarkerE','k');
set(gca,'xlim',phi0+2*TH*[-1 1],'ylim',90-theta0+2*TH*[-1 1],...
	'xtick',phi0+[-2*TH -TH 0 TH 2*TH],...
	'ytick',90-theta0+[-2*TH -TH 0 TH 2*TH])
grid on
hold off
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=setnans(data)
% Cut out the smallest values for ease of visualization
data(abs(data)<max(abs(data(:)))/1000)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cosmo(ah,ha,H)
% Cosmetic adjustments
longticks(ah); deggies(ah)
fig2print(gcf,'landscape')
nolabels(ah(1:6),1)
nolabels(ha(4:9),2)
serre(H,1,'across')
serre(H',1/2,'down')
