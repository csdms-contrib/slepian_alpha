function varargout=rcenter(lola,method)
% [lonc,latc,A,lonm,latm]=RCENTER(lola,method)
%
% Determines the center of mass of a surface on the sphere bounded by a
% closed curve given as longitude/latitude in degrees, i.e. the
% geographic centroid. This only works as long as the patch is contained
% within one hemisphere.
%
% INPUT:
% 
% lola    [lon(:) lat(:)] coordinates of a closed curve [degrees]
% method  'GL' or 'Green'
%         'Matlab' (area only)
%
% OUTPUT:
% 
% lonc,latc  The coordinates of the center of the region [degrees]
% A          The area of the region on the unit sphere [steradians]
% lonm,latm  The geographic mean of the coordinates of the curve [degrees]
%
% EXAMPLE:
%
% rcenter('demo1') % Random spherical caps
% rcenter('demo2') % Australia
% rcenter('demo3') % North American geological domains
%
% SEE ALSO:
%
% SPHAREA, AREAINT, CURVECHECK, CAPLOC
%
% Last modified by fjsimons-at-alum.mit.edu, 12/24/2009

defval('plots',0)
defval('method','GL')

% Initialize
[lonc,latc]=deal(NaN);

if ~isstr(lola)
  % You should get a first idea of the center and work relative to that,
  % if not you may end up getting strange results for very large and for
  % near-polar regions - assuming the curve is "uniformly densely" defined
  [latm,lonm]=meanm(lola(:,2),lola(:,1));
  
  switch method
    case 'GL'
     % Rotate to the best-guess origin, don't just subtract
     [tp,pp]=rottp((90-lola(:,2))*pi/180,lola(:,1)*pi/180,...
		   lonm*pi/180,-latm*pi/180,0);
     if plots
       [x1,y1,z1]=sph2cart(lola(:,1)*pi/180,pi/180*lola(:,2),...
			   ones(size(lola(:,1))));
       [x2,y2,z2]=sph2cart(pp,pi/2-tp,...
			   ones(size(lola(:,1))));
       plot3(x1,y1,z1,'bo'); hold on
       plot3(x2,y2,z2,'rv'); hold off; axis image
       view(0,0); pause
     end
     % Note that this remapping appears to change the area ever so slightly
     lola(:,1)=pp*180/pi;
     lola(:,2)=(pi/2-tp)*180/pi;
     
     % Find the bounding points of the region, convert to radians
     thN=90-max(lola(:,2)); thN=thN*pi/180;
     thS=90-min(lola(:,2)); thS=thS*pi/180;
     
     if plots
       nidit=abs(lola(:,1))>90;
       lola(nidit,1)=lola(nidit,1)-sign(lola(nidit,1))*180;
       plot(lola(:,1),lola(:,2),'o'); hold on
       plot(0,0,'bo','MarkerF','b')
       plot(0,90-thN*180/pi,'kv','MarkerF','k')
       plot(0,90-thS*180/pi,'kv','MarkerF','k')
     end
     
     % Set up the Gaussian integration
     intv=cos([thS thN]);
     % The number of Gauss-Legendre points can actually be fairly low in
     % order to integrate a sine curve
     nGL=32;
     % These are going to be the required colatitudes
     [w,x,N]=gausslegendrecof(nGL,[],intv);
     % Get the integration info for the domain - treat as "flat"
     [phint,thh,phh]=phicurve([90-lola(:,2) lola(:,1)],acos(x)*180/pi);
     % Riemann-sum equivalent?
     % xR=linspace(cos(thS-eps),cos(thN+eps),50)';
     % [phiR,thhR,phhR]=phicurve([90-lola(:,2) lola(:,1)],...
     %        acos(xR)*180/pi);
     if plots
       plot(phh,90-thh,'k-'); hold off ; pause
     end
     % Convert back to radians
     phint=phint*pi/180;
     % phiR=phiR*pi/180;
     % Do the longitudinal integrals - avoid COSCOS yet it is like it
     % See SPHAREA, in other words
     % The unweighted one is part of the th- and unweighted th-integrand
     IPH1=sum(phint(:,2:2:end)-phint(:,1:2:end),2);
     % Riemann-sum equivalent?
     % IPH1R=sum(phiR(:,2:2:end)-phiR(:,1:2:end),2);
     % The phi-weighted one is part of the unweighted th-integrand
     IPH2=sum(phint(:,2:2:end).^2-phint(:,1:2:end).^2,2)/2;
     
     % The area of the region - GL better than a straight Riemann sum
     A=w'*IPH1;
     % AR=sum(IPH1R*indeks(diff(xR),1));
     % The first moment of the colatitude
     TH0=w'*(acos(x).*IPH1);
     % The first moment of the longitude
     PH0=w'*IPH2;
     
     % Return the properly normalized and unitized output
     [thc,phc]=rottp(TH0/A,PH0/A,0,latm*pi/180,0-lonm*pi/180);
   case 'Green'
    % This inspired by AREAINT and see RB VIII p 103--104
    % As a line integral, see GEO371 notes also
    % Take the midpoints, put away the sign
    % Any constant stuck in front vanishes!!
    th=pi/2-lola(:,2)*pi/180;
    newth=th(1:end-1)+diff(th)/2;
    phi=wrapTo2Pi(lola(:,1)*pi/180);
    newphi=phi(1:end-1)+diff(phi)/2;
    % Intermediate step for the integrand and the step
    intm=-cos(newth).*diff(phi);
    A=abs(sum(intm));
    % Not quite right yet, though I wonder why, check sign of area and
    % see how the poles need to be flipped into submission
    phc=sum(newphi.*intm);
    thc=sum((sin(newth)-newth.*cos(newth)).*intm);
    % So it's going to depend on the resolution of the graph, the
    % evenness of the spacing of the graph along the curve, and so on 
    % Work with the splined boundaries, thus
    % and now do the center of mass, the same way
   case 'Matlab'
    % Compare with AREAINT, which uses a different origin
    A=areaint(lola(:,2),lola(:,1))*4*pi;
  end

  % Output
  lonc=phc*180/pi;
  latc=90-thc*180/pi;
  vars={lonc,latc,A,lonm,latm};
  varargout=vars(1:nargout);
elseif strcmp(lola,'demo1')
  % Calculate positions of a spherical patch
  TH=ceil(rand*90);
  [lon1,lat1]=randsphere(1);
  [lon2,lat2]=caploc([lon1 lat1],TH,256);
  [lonc,latc,A,lonm,latm]=rcenter([lon2 lat2]);
  % Make plotting adjustments
  plot(lon2,lat2,'+'); hold on
  plot(lon1,lat1,'o')
  plot(lonc,latc,'rv');
  plot(lonc+360,latc,'rv');
  plot(lonm,latm,'g^');
  axis([lon1-2*TH lon1+2*TH lat1-2*TH lat1+2*TH]); hold off
  title(sprintf('(%8.3f,%8.3f) vs (%8.3f,%8.3f) for %s=%i and A/A=%12.5e',...
		lonc,latc,lon1,lat1,'\Theta',TH,A/spharea(TH,1)/4/pi))
  fig2print(gcf,'portrait')
elseif strcmp(lola,'demo2')
  lola=australia;
  plot(lola(:,1),lola(:,2))
  hold on
  [lonc,latc,A,lonm,latm]=rcenter(lola);
  plot(lonc,latc,'rv'); 
  plot(lonc+360,latc,'rv');
  plot(lonm,latm,'g^');
  axis([90+360 180+360 -60 10]); hold off
elseif strcmp(lola,'demo3')
  dirname=fullfile(getenv('IFILES'),'GEOLOGY','NORTHAMERICA');
  load(fullfile(dirname,'tapestry.mat'))
  names={'ColoradoPlateaus','ColumbiaPlateau'};
  par=2-round(rand); 
  % Find the curve defining the regions in GEOCENTRIC LONGITUDE and LATITUDE
  lola=[lon.(names{par})' lat.(names{par})'];
  % Figure out the center of mass for this projection
  [lonc,latc,A,lonm,latm]=rcenter(lola);
  % Compare with the 'means' of the longitude and the latitude
  clf
  plot(lola(:,1),lola(:,2))
  hold on
  plot(mean(lola(:,1)),mean(lola(:,2)),'o')
  [mla,mlo]=meanm(lola(:,2),lola(:,1));
  plot(mlo,mla,'+')
  plot(lonc,latc,'rv')
end

% In Maple:
% with(VectorCalculus):
% SetCoordinates( 'spherical'[r,t,p] );
% F := VectorField( <0,-p*sin(t)+p*(cos(t)+sin(t)^2-cos(t)^2),(1-cos(t))>);
% F := VectorField( <0,p*sin(t)*(2*cos(t)-1),sin(t)>);
% F := VectorField( <0,p*(1-2*sin(t)^2-sin(t)),cos(t)>);
% simplify(eval(Curl(F),r=1));
% All of these fields have a radial curl that is equal to one, and thus
% the integral of the radial curl component is the surface area under the
% sphere 


