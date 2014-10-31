function varargout=glmalphapto(TH,L,phi,theta,omega)
% [G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=GLMALPHAPTO(TH,L,phi,theta,omega)
%
% Returns an (lm)X(alpha) matrix with spherical harmonic coefficients of
% the bandlimited or bandpass Slepian functions of the SINGLE CAP rotated
% to a desired location and azimuthally by a certain amount.
%
% INPUT:
%
% TH        Angular extent of the spherical cap (degrees)
% L         Bandwidth (maximum angular degree) or passband (two degrees)
% phi       Longitude of the center (degrees)
% theta     Colatitude of the center (degrees)
% omega     Anticlockwise azimuthal rotation (degrees) [default: 0]
%
% OUTPUT:
%
% G        The unitary matrix of localization coefficients
% V        The eigenvalues in this ordering (not automatically sorted)
% EL       The vector with spherical harmonic degrees as first index of G
% EM       The vector with spherical harmonic orders as first index of G
% N        The Shannon number
% GM2AL    The sum over all orders of the squared coefficients, i.e. the
%          TOTAL power, NOT the power spectral density
% MTAP     The original uncorrupted order of the eigentapers
% IMTAP    The rank within that particular order of the eigentapers
%
% EXAMPLE:
%
% glmalphapto('demo1') % For some basic consistency checking
% glmalphapto('demo2') % For a few quick plots
% glmalphapto('demo3') % One way of doing the spatial expansion
% glmalphapto('demo4') % Another way of doing the spatial expansion
% 
% SEE ALSO: GLMALPHA, PTOSLEP, GALPHA
%
% Last modified by fjsimons-at-alum.mit.edu, 01/17/2010

% This becomes troublesome at large (>100) spherical-harmonic
% bandwidths. Obviously, there we enter the domain of the Cartesian
% modeling. But we should still be able to provide computational and
% memory savings by rewriting this routine in a better way.

% Supply defaults
defval('TH',30)

if ~isstr(TH)
  defval('L',18)
  defval('phi',78)
  defval('theta',78)
  defval('omega',10)

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
      fname=fullfile(getenv('IFILES'),'GLMALPHAPTO',...
		     sprintf('glmalphapto-%i-%i-%i-%i-%i.mat',...
			     TH,L,phi,theta,omega));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GLMALPHAPTO',...
		     sprintf('glmalphaptobl-%i-%i-%i-%i-%i-%i.mat',...
			     TH,L,phi,theta,omega));
    end
  else
    fname='neveravailable';
  end

  % Initialize ordering matrices
  MTAP=repmat(0,1,ldim);
  IMTAP=repmat(0,1,ldim);

 
  if exist(fname,'file')==2 
    load(fname)
    disp(sprintf('Loading %s',fname))
  else
    % Initialize matrices
    G=repmat(0,(maxL+1)^2,ldim);
    V=repmat(0,1,ldim);
    
    % Find indices into G belonging to the orders
    [EM,EL,mz,blkm]=addmout(maxL);
    
    % Find increasing column index; that's how many belong to this order
    % alpha=cumsum([1 L+1 gamini(L:-1:1,2)]);
    alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
  		gamini(L(2-lp)-bp*(L(1)-1),bp*2*L(1)) ...
  		gamini(L(2-lp)-bp*L(1):-1:1,2)]);
  
    % For AXISYMMETRIC REGIONS
    % For the SINGLE CAP ONLY
    for m=-L:L
      disp(sprintf(...
	  'Figure out how to work with abs(m) instead - now %i',m))

      % Get the coefficients of the rotated bases from rotating GRUNBAUM
      % if lowpass or SDWCAP if bandpass
      [lm,cosi,C,EL,EM,Vp]=ptoslep(phi,theta,omega,TH,L,m);

      % Distribute this at the right point in the huge matrix
      if m<0
	% Here you do the originally negative orders, now all over
	G(:,alpha(2*abs(m)):alpha(2*abs(m)+1)-1)=C;
	V(alpha(2*abs(m)):alpha(2*abs(m)+1)-1)=Vp;
	MTAP(alpha(2*abs(m)):alpha(2*abs(m)+1)-1)=m;
	% It's all neatly ordered here, downgoing within every order
	IMTAP(alpha(2*abs(m)):alpha(2*abs(m)+1)-1)=1:length(Vp);
      else
	% And here you do the originally positive orders, now all over
	G(:,alpha(2*m+1):alpha(2*m+2)-1)=C;
	V(alpha(2*m+1):alpha(2*m+2)-1)=Vp;
	MTAP(alpha(2*m+1):alpha(2*m+2)-1)=m;
	% It's all neatly ordered here, downgoing within every order
	IMTAP(alpha(2*m+1):alpha(2*m+2)-1)=1:length(Vp);
      end
    end
    % Calculate the Shannon number and compare it to the theory
    N=sum(V);
    difer(N-ldim*(1-cos(TH/180*pi))/2,[],[],NaN);
    
    % Compute the sum over all orders of the squared coefficients
    % Thus works when they have not been blocksorted yet. 
    GM2AL=repmat(0,ldim,maxL+1);
    for l=0:maxL
      b=(l-1+1)^2+1;
      e=(l+1)^2;
      GM2AL(:,l+1)=sum(G(b:e,:).^2,1)';
    end
    % Make sure that the sum over all degrees is 1
    difer(sum(GM2AL,2)-1,[],[],NaN)

    % Save the results
    if ~strcmp(fname,'neveravailable') 
      % If the variable is HUGE you must use the -v7.3 flag, if not, you
      % can safely omit it and get more backwards compatibility
      save(fname,'-v7.3','G','V','EL','EM','N','GM2AL','MTAP','IMTAP')
    end
  end
  
  % Provide output
  varns={G,V,EL,EM,N,GM2AL,MTAP,IMTAP};
  varargout=varns(1:nargout);
elseif strcmp(TH,'demo1')
  % Matrices not orthogonal; total power etc identical
  TH=30; L=18;
  [G1,V1,EL1,EM1,N1,GM2AL1,MTAP1,IMTAP1]=glmalphapto(TH,L,78,78,10);
  [G2,V2,EL2,EM2,N2,GM2AL2,MTAP2,IMTAP2]=glmalphapto(TH,L,0,0,0);
  [G3,V3,EL3,EM3,N3,GM2AL3,MTAP3,IMTAP3]=glmalpha(TH,L);
  difer(GM2AL1-GM2AL2); difer(V1-V2); difer(MTAP1-MTAP2)
  difer(IMTAP1-IMTAP2); difer(N1-N2); difer(EL1-EL2); difer(EM1-EM2)
  difer(G2-G3); difer(GM2AL2-GM2AL3); difer(V2-V3); difer(MTAP2-MTAP3)
  difer(IMTAP2-IMTAP3); difer(N3-N3); difer(EL3-EL3); difer(EM3-EM3)
  difer(G1'*G1-eye(size(G1)))
  difer(G2'*G2-eye(size(G2)))
  difer(G3'*G3-eye(size(G3)))
  % We lost the block ordering in the rotation, it is still instructive
  % to sort as if it were - remember the columns are always order by order
  [EM,EL,mz,blkm]=addmout(L);
  imagesc(G3(blkm,:))
elseif strcmp(TH,'demo2')
  clf
  L=20;
  p=200;
  t=90;
  o=0;
  [G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=glmalphapto(20,L,p,t,o);
  % Collect the output into a format that PLM2XYZ knows how to interpret
  [~,~,~,lmcosi,~,mzo,~,~,rinm,ronm]=addmon(L);
  % Create the blanks
  cosi=lmcosi(:,3:4);
  % Stick in the coefficients of the 1st eigentaper
  wot=max(abs(guess(1)),1);
  cosi(ronm)=G(:,wot);
  % Construct the full matrix
  lmcosi(:,3:4)=cosi;
  plotplm(lmcosi,[],[],4,1)
  undeggies(gca)
  set(gca,'xtick',p,'xtickl',p)
  set(gca,'ytick',90-t,'ytickl',90-t); grid on
  deggies(gca)
  title(sprintf('L = %i ; i = %i',L,wot))
elseif strcmp(TH,'demo3')
  TH=15;
  phi0=15;
  theta0=70;
  L=36;
  
  defval('meth',1);
  defval('degres',1);

  % Construct the rotated basis
  [G,V,EL,EM,N]=glmalphapto(TH,L,phi0,theta0);
  
  % Sort the taper and take the first N good ones
  [V,i]=sort(V,'descend');
  G=G(:,i); G=G(:,1:round(N));
  
  % Get the spherical harmonics on the desired grid
  % Plot the map
  clf
  [ah,ha]=krijetem(subnum(3,3));
  
  % Define the region of interest for plotting
  c11cmn=[phi0-2*TH (90-theta0)+2*TH phi0+2*TH (90-theta0)-2*TH];
  
  % The next thing is now part of GALPHAPTO
  switch meth
   case 1 % FIRST METHOD
      % Generate indices that PLOTPLM/PLM2XYZ know how to interpret
      [a1,a2,a3,lmcosi,a5,mzo,a7,a8,rinm,ronm]=addmon(L);
      % Create the blanks
      cosi=lmcosi(:,3:4);
      disp(sprintf('\nMethod I\n'))
   case 2 % SECOND METHOD
    theta=[90-c11cmn(2):degres:90-c11cmn(4)]*pi/180;
    phi=[c11cmn(1):degres:c11cmn(3)]*pi/180;
    [XYlmr,t,p,dem,del]=ylm([0 L],[],theta,phi);
    % There is a need for an additional phase factor here - because my
    % PLM2ROT/PLM2XYZ/KERNELC don't include it whereas YLM/XLM do.
    G=G.*repmat((-1).^EM,1,size(G,2));
    Gspace=[G'*XYlmr]';
  end
  
  for index=1:length(ah)
    axes(ah(index))
    switch meth
     case 1 % FIRST METHOD
        % Stick in the coefficients of the indexth eigentaper
	cosi(ronm)=G(:,index);
	% Construct the full matrix
	lmcosi(:,3:4)=cosi/sqrt(4*pi);
	% Note that we need to remove the 4 pi factor if we are going to
        % use PLM2XYZ which will be putting it in!
	% Expand to space around the area of interest
	[data,lond,latd]=plm2xyz(lmcosi,degres,c11cmn);
     case 2 % SECOND METHOD
	% Extract the data directly from the expansion
	data=reshape(Gspace(:,index),length(theta),length(phi));
    end
    
    % Is this identical? Compare the two methods offline. Yes, is the answer.
        
    % Cut out the smallest values for ease of visualization
    data(abs(data)<max(data(:))/1000)=NaN;
    imagefnan(c11cmn(1:2),c11cmn(3:4),data)
    % Plot the center
    hold on
    d=plot(phi0,90-theta0,'o','MarkerF','w','MarkerE','k');
    % Plot the circle of concentration
    [lon,lat]=caploc([phi0 90-theta0],TH);
    e=plot(lon-360*[lon>180],lat,'k-');
    % Plot the continents - note the Greenwich trick
    yesorno=360*[c11cmn(1)<0 0];
    [ax,f,lola1]=plotcont(c11cmn(1:2)+yesorno,...
			  c11cmn(3:4)+yesorno,[],-360);
    [ax,f,lola2]=plotcont([0 c11cmn(2)],c11cmn(3:4));
    axis(c11cmn([1 3 4 2]))
    title(sprintf('%s = %8.3f','\lambda',V(index)))
    drawnow
  end
  
  % Cosmetic adjustments
  longticks(ah); deggies(ah)
  fig2print(gcf,'landscape')
  nolabels(ah(1:6),1)
  nolabels(ha(4:9),2)
  figdisp
elseif strcmp(TH,'demo4')
  TH=15;
  phi0=35;
  theta0=35;
  L=36;

  defval('meth',2);
  defval('degres',1);

  % Construct the rotated basis
  [G,V,EL,EM,N]=glmalphapto(TH,L,phi0,theta0);
  
  % Sort the taper and take the first N good ones
  [V,i]=sort(V,'descend');
  G=G(:,i); G=G(:,1:round(N));
  
  % Get the spherical harmonics on the desired grid
  % Plot the map
  clf
  [ah,ha]=krijetem(subnum(3,3));
  
  % Define the region of interest for plotting
  c11cmn=[phi0-2*TH (90-theta0)+2*TH phi0+2*TH (90-theta0)-2*TH];
  
  % The next thing is now part of GALPHAPTO
  switch meth
   case 1 % FIRST METHOD
     % Generate indices that PLOTPLM/PLM2XYZ know how to interpret
     [a1,a2,a3,lmcosi,a5,mzo,a7,a8,rinm,ronm]=addmon(L);
     % Create the blanks
     cosi=lmcosi(:,3:4);
   case 2 % SECOND METHOD
     theta=[90-c11cmn(2):degres:90-c11cmn(4)]*pi/180;
     phi=[c11cmn(1):degres:c11cmn(3)]*pi/180;
     [XYlmr,t,p,dem,del]=ylm([0 L],[],theta,phi);
     % There is a need for an additional phase factor here - because my
     % PLM2ROT/PLM2XYZ/KERNELC don't include it whereas YLM/XLM do.
     G=G.*repmat((-1).^EM,1,size(G,2));
     Gspace=[G'*XYlmr]';
     disp(sprintf('\nMethod II\n'))
  end
  
  for index=1:length(ah)
    axes(ah(index))
    switch meth
     case 1 % FIRST METHOD
        % Stick in the coefficients of the indexth eigentaper
	cosi(ronm)=G(:,index);
	% Construct the full matrix
	lmcosi(:,3:4)=cosi/sqrt(4*pi);
	% Note that we need to remove the 4 pi factor if we are going to
        % use PLM2XYZ which will be putting it in!
	% Expand to space around the area of interest
	[data,lond,latd]=plm2xyz(lmcosi,degres,c11cmn);
     case 2 % SECOND METHOD
        % Extract the data directly from the expansion
	data=reshape(Gspace(:,index),length(theta),length(phi));
    end
    
    % Cut out the smallest values for ease of visualization
    data(abs(data)<max(data(:))/1000)=NaN;
    imagefnan(c11cmn(1:2),c11cmn(3:4),data)
    % Plot the center
    hold on
    d=plot(phi0,90-theta0,'o','MarkerF','w','MarkerE','k');
    % Plot the circle of concentration
    [lon,lat]=caploc([phi0 90-theta0],TH);
    e=plot(lon-360*[lon>180],lat,'k-');
    % Plot the continents - note the Greenwich trick
    yesorno=360*[c11cmn(1)<0 0];
    [ax,f,lola1]=plotcont(c11cmn(1:2)+yesorno,...
			  c11cmn(3:4)+yesorno,[],-360);
    [ax,f,lola2]=plotcont([0 c11cmn(2)],c11cmn(3:4));
    axis(c11cmn([1 3 4 2]))
    title(sprintf('%s = %8.3f','\lambda',V(index)))
    drawnow
  end
  
  % Cosmetic adjustments
  longticks(ah); deggies(ah)
  fig2print(gcf,'landscape')
  nolabels(ah(1:6),1)
  nolabels(ha(4:9),2)
  figdisp
end
