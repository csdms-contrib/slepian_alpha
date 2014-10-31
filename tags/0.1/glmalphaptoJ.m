function varargout=glmalphaptoJ(TH,L,phi,theta,omega,J)
% [G,V,EL,EM,N,MTAP,IMTAP]=GLMALPHAPTOJ(TH,L,phi,theta,omega,J)
%
% Returns an (lm)X(alpha) matrix with spherical harmonic coefficients of
% the bandlimited or bandpass Slepian functions of the SINGLE CAP rotated
% to a desired location and azimuthally by a certain amount. Only J
% functions are being calculated, and these will correspond to the best
% concentrated functions in the rough sense that the best of every order
% will be calculated first, then the second best of every order,
% etc. There is NO guarantee that you are getting the same result as when
% you ask for the J best of every order, and in general this function
% could behave badly so use at your own risk.
%
% INPUT:
%
% TH        Angular extent of the spherical cap (degrees)
% L         Bandwidth (maximum angular degree) or passband (two degrees)
% phi       Longitude of the center (degrees)
% theta     Colatitude of the center (degrees)
% omega     Anticlockwise azimuthal rotation (degrees) [default: 0]
% J         The number of eigenfunctions that are being asked (and saved)
%
% OUTPUT:
%
% G        The unitary matrix of localization coefficients
% V        The eigenvalues in this ordering (not automatically sorted)
% EL       The vector with spherical harmonic degrees as first index of G
% EM       The vector with spherical harmonic orders as first index of G
% N        The Shannon number
% MTAP     The original uncorrupted order of the eigentapers
% IMTAP    The rank within that particular order of the eigentapers
%
% EXAMPLE:
%
% glmalphaptoJ('demo')
%
% Last modified by fjsimons-at-alum.mit.edu, 01/11/2011

% Supply defaults
defval('TH',30)

if ~isstr(TH)
  defval('L',18)
  defval('phi',78)
  defval('theta',78)
  defval('omega',10)
  defval('J',3)

  % Figure out if it's lowpass or bandpass
  lp=length(L)==1;
  bp=length(L)==2;
  maxL=max(L);

  % The spherical harmonic dimension
  ldim=(L(2-lp)+1)^2-bp*L(1)^2;

  % If angles are integers, save the results
  if ~mod(phi,round(phi)) && ~mod(theta,round(theta)) ...
	&& ~mod(omega,round(omega)) 
    if lp
      fname=fullfile(getenv('IFILES'),'GLMALPHAPTOJ',...
		     sprintf('glmalphaptoJ-%i-%i-%i-%i-%i-%i.mat',...
			     TH,L,phi,theta,omega,J));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GLMALPHAPTOJ',...
		     sprintf('glmalphaptoJbl-%i-%i-%i-%i-%i-%i-%i.mat',...
			     TH,L,phi,theta,omega,J));
    end
  else
    fname='neveravailable';
  end

  if exist(fname,'file')==2 
    load(fname)
    disp(sprintf('Loading %s',fname))
  else
    % Initialize matrices
    G=repmat(0,(maxL+1)^2,J);
    V=repmat(0,1,J);
    % Initialize ordering matrices
    MTAP=repmat(0,1,J);
    IMTAP=repmat(0,1,J);

    % Find indices into G belonging to the orders
    [EM,EL,mz,blkm]=addmout(maxL);
    
    % Find increasing column index; that's how many belong to this order
    % alpha=cumsum([1 L+1 gamini(L:-1:1,2)]);
    alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
  		gamini(L(2-lp)-bp*(L(1)-1),bp*2*L(1)) ...
  		gamini(L(2-lp)-bp*L(1):-1:1,2)]);
  
    % For AXISYMMETRIC REGIONS
    % For the SINGLE CAP ONLY
    % This is rearranged from the regular GLMALPHAPTO
    oldm=-L:L;
    [ii,jj]=sort(abs(oldm));
    oldm=oldm(jj);
    % Don't do too much work - maximum number per order calculated,
    % noting that this works worst when J is close to (J+1)^2, in which
    % case you really need them all. There is no telling how the
    % eigenvalues will interleave until you've done them all
    numm=ceil(J/length(oldm));
    mmin=1;
    for m=oldm(1:min(length(oldm),J))
      % Get the coefficients of the rotated bases from rotating GRUNBAUM
      % if lowpass or SDWCAP if bandpass
      [lm,cosi,C,EpL,EM,Vp]=ptoslep(phi,theta,omega,TH,L,m,numm);

      mmax=mmin+size(C,2)-1;

      % Add this over the newly small matrix
      G(:,mmin:mmax)=C;
      V(mmin:mmax)=Vp(1:min(numm,length(Vp)));
      MTAP(mmin:mmax)=m;
      % It's all neatly ordered here, downgoing within every order
      IMTAP(mmin:mmax)=1:min(numm,length(Vp));
	
      % Update the index
      mmin=mmax+1;
      
      if mmin>J
	break
      end
    end

    % At the end you very well may end up needing to sort again
    
    % Calculate the Shannon number
    N=ldim*(1-cos(TH/180*pi))/2;
    
    % Save the results
    if ~strcmp(fname,'neveravailable') 
      % If the variable is HUGE you must use the -v7.3 flag, if not, you
      % can safely omit it and get more backwards compatibility
      save(fname,'-v7.3','G','V','EL','EM','N','MTAP','IMTAP')
    end
  end

  % Provide output
  varns={G,V,EL,EM,N,MTAP,IMTAP};
  varargout=varns(1:nargout);
elseif strcmp(TH,'demo')
  clf
  L=72;
  J=3;
  lon1=180+round(rand*60);
  lat1=round(rand*30);
  omg=round(rand*360);
  TH=20;
  [G,V,EL,EM,N,MTAP,IMTAP]=glmalphaptoJ(TH,L,lon1,90-lat1,omg,J);
  % Collect the output into a format that PLM2XYZ knows how to interpret
  [a1,a2,a3,lmcosi,a5,mzo,a7,a8,rinm,ronm]=addmon(L);
  % Create the blanks
  cosi=lmcosi(:,3:4);
  % Stick in the coefficients of the 1st eigentaper
  wot=ceil(rand*J);
  cosi(ronm)=G(:,wot);
  % Construct the full matrix
  lmcosi(:,3:4)=cosi;
  plotplm(lmcosi,[],[],4,1)
  undeggies(gca)
  set(gca,'xtick',lon1,'xtickl',lon1)
  set(gca,'ytick',lat1,'ytickl',lat1); grid on
  deggies(gca)
  longticks(gca)
  hold on
  [lon2,lat2]=caploc([lon1 lat1],TH);
  plot(lon2,lat2) 
  title(sprintf(...
      'L = %i ; %s = %i ; %s = %i ; N = %i; i = %i ; m = %i ; V = %5.2f',...
      L,'\Theta',TH,'\omega',omg,round(N),wot,MTAP(wot),V(wot)))
end
