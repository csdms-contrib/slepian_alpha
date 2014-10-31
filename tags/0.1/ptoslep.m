function varargout=ptoslep(phi,theta,omega,TH,L,m,amax)
% [lm,cosi,glma,EL,EM,V]=PTOSLEP(phi,theta,omega,TH,L,m,amax)
%
% Computes SINGLE-CAP FIXED-ORDER azimuthally rotated bandlimited or
% passband Slepian tapers centered at a desired geographical location
% from an initial North-Polar location, that is.
%
% INPUT:
%
% phi        Longitude of the new center (degrees)
% theta      Colatitude of the new center (degrees)
% omega      Clockwise azimuthal rotation (degrees)
%            This is consistent with PLM2ROT(lmcosi,omega,-theta,-phi)
% TH         Radius of the concentration region (degrees)
% L          Bandwidth of the Slepian basis, or passband (two degrees)
% m          Single angular order of the Slepian basis function, -l=<m<=l
% amax       The number of tapers you want [default: all available]
%
% OUTPUT:
%
% lm         The degrees and orders of the expansion, an (L+1)^2 X 2 array
% cosi       Real spherical harmonic expansion coefficients suited for
%            plotting by PLM2XYZ; this is a (L+1)^2 X 2 X amax array
% glma       The same as if they came out of GLMALPHA, (L+1)^2 X amax
%            After rotation of course these aren't single-order anymore,
%            so block-ordering is not an option here
% EL,EM      An array with degrees and orders belonging to the latter
% V          The eigenvalues or quadratic concentration fractions
%
% EXAMPLE:
%
% ptoslep('demo1') % Some random rotations
% ptoslep('demo2') % Rotations by pi/2/m
% ptoslep('demo3') % Various rotations over omega
% ptoslep('demo4') % Check that unrotated polar caps are same as GLMALPHA
% ptoslep('demo5') % Check that rotations over pi/2/m are as GLMALPHA
% ptoslep('demo6') % Check the results for changed order signs
%
% See also: SDWTARGET, GLMALPHAPTO
%
% Last modified by fjsimons-at-alum.mit.edu, 06/20/2011

% Should I be able to predict the relation between -m (cos) and +m (sin)
% expansion coefficients - AFTER MAKING THESE COMPLEX? Remember, GRUNBAUM
% and SDWCAP are the same for +/- m. But then they go in at the right
% orders into lmcosi. Then rotated. If then made complex, the relation
% should be easy. But should look at UMMP again to see exactly in what
% order we need to present the coefficients.

defval('phi',78)

if ~isstr(phi)
  defval('theta',78)
  defval('omega',10)
  defval('TH',20)
  defval('L',18)
  defval('m',0)
  
  % Figure out if it's lowpass or bandpass
  lp=length(L)==1;
  bp=length(L)==2;
  maxL=max(L);

  % The spherical harmonic dimension
  ldim=(L(2-lp)+1)^2-bp*L(1)^2;
    
  defval('amax',maxL-max(abs(m),bp*min(L))+1);

  % Amax can't be more than however many are available
  amax=min(amax,maxL-max(abs(m),bp*min(L))+1);
  
  % Calculate the taper coefficients centered on the pole
  if lp
    [E,Vg,th,C,T,V]=grunbaum(TH,L,m,0);
  elseif bp
    % Note that the small-eigenvalue eigenfunctions might be
    % numerically degenerate and thus not as using Grunbaum - if
    % you need to compare, compare where the eigenvalues are "big"
    [E,V,Np,th,C]=sdwcap(TH,L,m,0,-1);
  else
    error('The degree range should be either one or two numbers')
  end

  % Make a blank array of indices and coefficients
  [em,el,mzero,jk1,jk2,jk3,jk4,jk5,jk6,ronm]=addmon(maxL);
  lm=[el em];
  % The final arrays
  cosi=repmat(0,[length(em) 2 amax]);
  glma=repmat(0,[(maxL+1)^2 amax]);
  
  % The intermediate array
  lmcosi=[el em repmat(0,length(em),2)];
  
  % Do it!
  for index=1:amax
    % Put the zonal coefficients in the right place
    % Negative orders should be the cosine tapers
    lmcosi(mzero(abs(m)+1:end)+abs(m),3+(m>0))=C(:,index);

    % Rotate the tapers to a position on the sphere
    tempry=kindeks(plm2rot(lmcosi,omega,-theta,-phi),3:4);

    % Put them into the lmcosi type matrix
    cosi(:,:,index)=tempry;

    % Put them into the not-block-ordered GLMALPHA type matrix
    glma(:,index)=tempry(ronm);
  end

  % Supply the indices for the GLMALPHA type matrix
  [EM,EL]=addmout(maxL);

  % Assign output
  varns={lm,cosi,glma,EL,EM,V};
  varargout=varns(1:nargout);
elseif strcmp(phi,'demo1')
  phi=round(rand*180)+90;
  theta=round(rand*90)+45;
  omega=round(rand*360);
  TH=30; L=18; m=1;
  [lm,cosi]=ptoslep(phi,theta,omega,TH,L,m,4);
  clf
  % Calculate cap coordinates
  [lon1,lat1]=caploc([phi 90-theta],TH);
  % Calculate azimuth coordinates, define EAST from SOUTH
  [lon2,lat2]=grcazim([phi 90-theta],...
		      linspace(0,TH*fralmanac('DegDis')/1000,100),180+omega);
  for index=1:4
    ah(index)=subplot(2,2,index);
    r=plm2xyz([lm cosi(:,:,index)],1);
    imagef([],[],r)
    set(gca,'xtick',[0 phi 360],'xtickl',[0 phi 360])
    set(gca,'ytick',[-90 90-theta 90],'ytickl',[-90 90-theta 90])
    hold on
    p1(index)=plot(lon1,lat1,'k','LineW',1);
    p2(index)=plot(lon2,lat2,'k','LineW',1);
  end
  supertit(ah(1:2),sprintf('lon= %i ; lat= %i ; rot = %i',phi,90-theta,omega));
  movev(ah(1:2),-0.025)
  longticks(ah)
  set(ah,'xgrid','on','ygrid','on')
elseif strcmp(phi,'demo2')
  phi=round(rand*180)+90; phi=180;
  theta=round(rand*90)+45; theta=90;
  % Pick parameters; positive order m is
  TH=30; L=12; m=-2;
  % Make a "sine taper" from a "cosine taper"
  omega=180/2/abs(m);
  [lm,cosi1]=ptoslep(phi,theta,omega,TH,L,m,2);
  % Make the "sine taper" by flipping the sign of the order
  [lm,cosi2]=ptoslep(phi,theta,0,TH,L,-m,2);
  % Calculate cap coordinates
  [lon1,lat1]=caploc([phi 90-theta],TH);
  % Calculate azimuth coordinates, define EAST from SOUTH
  [lon2,lat2]=grcazim([phi 90-theta],...
		      linspace(0,TH*fralmanac('DegDis')/1000,100),180-omega);
  % Start plot subdivision
  clf
  [ah,ha]=krijetem(subnum(2,2));
  % Plot the rotated cosine taper
  for index=1:2
    axes(ah(index))
    r=plm2xyz([lm cosi1(:,:,index)],1);
    imagef([],[],r)
    set(gca,'xtick',[0 phi 360],'xtickl',[0 phi 360])
    set(gca,'ytick',[-90 90-theta 90],'ytickl',[-90 90-theta 90])
    hold on
    p1(index)=plot(lon1,lat1,'k','LineW',1);
    p2(index)=plot(lon2,lat2,'k','LineW',1);
  end
  hold off
  % Plot the sine taper
  for index=1:2
    axes(ah(index+2))
    r=plm2xyz([lm cosi2(:,:,index)],1);
    imagef([],[],r)
    set(gca,'xtick',[0 phi 360],'xtickl',[0 phi 360])
    set(gca,'ytick',[-90 90-theta 90],'ytickl',[-90 90-theta 90])
    hold on
    p1(index)=plot(lon1,lat1,'k','LineW',1);
  end
  hold off
  supertit(ah(1:2),sprintf('lon= %i ; lat= %i ; rot = %i',phi,90-theta,omega));
  movev(ah(1:2),-0.025)
  longticks(ah)
  set(ah,'xgrid','on','ygrid','on')
elseif strcmp(phi,'demo3')
  phi=round(rand*180)+90; phi=180;
  theta=round(rand*90)+45; theta=90;
  % Pick parameters; positive order m is the cosine taper
  TH=30; L=12; m=-1;
  % Make a "sine taper" from a "cosine taper"
  omega=linspace(0,360,8);
  clf
  % Calculate cap coordinates
  [lon1,lat1]=caploc([phi 90-theta],TH);
  % Start plot subdivision
  [ah,ha]=krijetem(subnum(2,4));
  % Plot the rotated cosine taper
  for index=1:length(omega)
    axes(ah(index))
    % Make the taper
    [lm,cosi]=ptoslep(phi,theta,omega(index),TH,L,m,1);
    % Calculate azimuth coordinates, define EAST from SOUTH
    [lon2,lat2]=grcazim([phi 90-theta],...
			linspace(0,TH*fralmanac('DegDis')/1000,100),...
			180-omega(index));
    % Perform the expansion
    r=plm2xyz([lm cosi],1);
    imagef([],[],r)
    set(gca,'xtick',[0 phi 360],'xtickl',[0 phi 360])
    set(gca,'ytick',[-90 90-theta 90],'ytickl',[-90 90-theta 90])
    hold on
    p1(index)=plot(lon1,lat1,'k','LineW',1);
    p2(index)=plot(lon2,lat2,'k','LineW',1);
    x(index)=xlabel(sprintf('rot = %i',round(omega(index))));
  end
  supertit(ah(1:length(omega)/2),sprintf('lon= %i ; lat= %i',phi,90-theta));
  movev(ah(1:length(omega)/2),-0.025)
  longticks(ah)
  set(ah,'xgrid','on','ygrid','on')
  fig2print(gcf,'landscape')
elseif strcmp(phi,'demo4')
  % Check that output from GLMALPHA equals PTOSLEP at North Pole and omega=0 
  TH=30; L=18; m=-L:L;
  % Figure out where the equivalent sits in the output of GLMALPHA
  alpha=cumsum([1 L+1 gamini(L:-1:1,2)]);
  % See GLMALPHA for a bandpass update
  for index=1:length(m)
    [lm,cosi,glma,EL,EM]=ptoslep(0,0,0,TH,L,m(index),1);
    [G,V,EL2,EM2]=glmalpha(TH,L,1);
    difer(glma-G(:,alpha(2*abs(m(index))+(m(index)>=0))));
  end
elseif strcmp(phi,'demo5')
  % The first cosine coefficient, rotated by omega+180/(2|m|)...
  m=-1; omega=30; L=18; phi=78; theta=78; TH=30;
  [a,b,G1]=ptoslep(phi,theta,omega+180/2/abs(m),TH,L,m,1);
  % ... should be the same as the sine coefficient, rotated by just omega
  m=1;
  [a,b,G2]=ptoslep(phi,theta,omega,TH,L,m,1);
  difer(G1-G2)
  % Next check this with GLMALPHAPTO
elseif strcmp(phi,'demo6')
  % What should be the relation between +/-m?
  % Well, the result is NOT simply a rotation around the z-axis, so much
  % is clear, but rather an azimuthal rotation around the center
  m=3;
  [lm1,cosi1,glma1,EL1,EM1,V1]=ptoslep(40,70,10,30,12,m,1);
  [lm2,cosi2,glma2,EL2,EM2,V2]=ptoslep(40,70,10+180/2/abs(m),30,12,-m,1);
  difer(lm1-lm2); difer(EL1-EL2); difer(EM1-EM2); difer(V1-V2)
  difer(glma1-glma2); difer(cosi1-cosi2)
  clf
  subplot(211); plotplm([lm1 cosi1(:,:,1)],[],[],4,1)
  subplot(212); plotplm([lm2 cosi2(:,:,1)],[],[],4,1)
end

