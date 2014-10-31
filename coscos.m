function [I,phint,th,pars]=coscos(th,m1,m2,dom,pars)
% [I,phint,th,pars]=coscos(th,m1,m2,dom,pars)
%
% Calculates the longitudinal integral 
% \int_{A}{B}\cos(m1\phi)\cos(m2\phi)\,d\phi
% over a specified domain at some (Gauss-Legendre) integration point(s).
%
% INPUT:
%
% th         Colatitude(s) that you may need to define the patch
% m1,m2      Angular order of the cosine factors
% dom        A matrix with the 'phint' hatchings from PHICURVE, or
%            'patch'   spherical patch [default] with the following specs
%            'england' England, Scotland and Wales
% pars       [th0,ph0,thR] for 'patch', with 
%                 th0  Colatitude of the cap center, in radians
%                 ph0  Longitude of the cap center, in radians
%                 thR  Radius of the cap, in radians
%            N for 'england' with N the smoothness of the splining
%
% EXAMPLE:
%
% coscos('demo1')   % Illustrates the area of a random patch
% coscos('demo2')   % Compares result with GL integration
% [a,b]=coscos('demo2',0); % Compares result with GL integration
%
% Last modified by fjsimons-at-alum.mit.edu, 04/20/2009

defval('th',[])
defval('m1',10) % Not zero (see demo2)
defval('m2',4)
defval('dom','patch')    

if ~isstr(th)
  if isstr(dom)
    switch dom
     case 'patch'
      % Get the parameters of the dom
      defval('pars',[pi/4 pi/2 pi/9]);
      th0=pars(1); ph0=pars(2); thR=pars(3);
      defval('th',linspace(th0-thR,th0+thR,100));
      phint=dphpatch(th,thR,th0,ph0);
     case 'england' 
      defval('N',10)
      % Now we may have multiple pairs
      phint=dphengland(th*180/pi,N);
      phint=phint*pi/180;
     otherwise
      error('Specify valid domain')
    end
  else
    phint=dom;
  end       

  I=repmat(0,size(phint,1),1);
  for index=1:size(phint,2)/2
    A=phint(:,2*index-1); B=phint(:,2*index);
    if m1~=m2 & m1~=-m2
      I=I+1/2*(sin((m1-m2)*B)*m1+...
	     sin((m1-m2)*B)*m2+...
	     sin((m1+m2)*B)*m1-...
	     sin((m1+m2)*B)*m2-...
	     sin((m1-m2)*A)*m1-...
	     sin((m1-m2)*A)*m2-...
	     sin((m1+m2)*A)*m1+...
	     sin((m1+m2)*A)*m2)/...
	((m1-m2)*(m1+m2));
    elseif m1~=0
      m=m1;
      I=I+1/2*(+cos(m*B).*sin(m*B)+...
	     m*(B-A)-cos(m*A).*sin(m*A))/m;
    else
      I=I+(B-A);
    end
  end
  
elseif strcmp(th,'demo1')
  % Center of the circle, latitude, in degrees
  defval('cent',[10+rand*340 90-(rand*180)]);
  % Colatitude/Longitude of center, in radians
  th0=pi/2-cent(2)/180*pi;
  ph0=cent(1)/180*pi;
  % Radius of the circle, in degrees
  defval('rad',rand*min([90-cent(2) cent(1) abs(360-cent(1))]));
  % Radius of the circle, in radians
  thR=rad/180*pi;
  % Longitude/latitude of the circle, in degrees
  [lon,lat]=caploc(cent,rad,101,1);
  % North and South latitude, in degrees
  lan=cent(2)+rad; 
  las=cent(2)-rad; 
  % North and South colatitude, in radians
  thN=max(th0-thR,0);
  thS=min(th0+thR,pi);

  % Some generic sampling
  nth=abs(150/pi/2*(thS-thN));
  % A set of latitudes from South to North
  la=linspace(las,lan,nth)';
  % A set of colatitudes from North to South 
  th=linspace(thN,thS,nth)';
  % The corresponding longitudes; all others are wrong
  [Id,phintd]=coscos(th,1,1,'patch',[th0 ph0 thR]);

  % Make plot
  plot(lon,90-lat,'k','LineW',2); hold on
  plot(cent(1),90-cent(2),'o','MarkerE','k','MarkerF','g')
  a=phintd(:,1)*180/pi;
  b=phintd(:,2)*180/pi;
  % Don't try to fix the longitude mess-up, it will show up as an exterior
  % domain 
  plot([a b]',[th th]'*180/pi,'k'); 
  hold off
  axis ij image
  ylim([0 180])
  xlim([min(0,min(a)) max(max(b),360)])
  grid on
  xlabel('Longitude \phi');
  ylabel('Colatitude \theta');
  longticks(gca)
  title(sprintf('(%s,%s,%s) = (%i,%i,%i)',...
		'\theta_0','\phi_0','\theta_R',...
		round(90-cent(2)),round(cent(1)),round(rad)))

  figdisp
elseif strcmp(th,'demo2')
  if m1==0; yesplot=0; else yesplot=1; end
  defval('cent',[rand*360 90-(rand*180)])
  th0=pi/2-cent(2)/180*pi;
  ph0=cent(1)/180*pi;
  defval('rad',rand*90);
  [lon,lat]=caploc(cent,rad,101,1);
  thR=rad/180*pi;
  mmax=60;
  m1=round(rand*mmax-mmax/2);
  m2=round(rand*mmax-mmax/2);
  
  thN=max(th0-thR,0);
  thS=min(th0+thR,pi);
  
  thd=thN+rand*[thS-thN];  
  [Id(1),phintd]=coscos(thd,m1,m2,'patch',[th0 ph0 thR]);
  Id(2)=sinsin(thd,m1,m2,'patch',[th0 ph0 thR]);
  Id(3)=sincos(thd,m1,m2,'patch',[th0 ph0 thR]); 
  
  % Remember the dot product!
  cc=inline(sprintf('cos(%i*x).*cos(%i*x)',m1,m2));
  ss=inline(sprintf('sin(%i*x).*sin(%i*x)',m1,m2));
  sc=inline(sprintf('sin(%i*x).*cos(%i*x)',m1,m2));
  
  % This is quite plainly not enough
  GL1=abs(m1)+abs(m2);
  % This needs to be really, really high to be correct
  GL=200;
  % It's a good thing we are not doing this
  J(1)=gausslegendre(phintd,cc,GL);
  J(2)=gausslegendre(phintd,ss,GL);
  J(3)=gausslegendre(phintd,sc,GL);
  
  phin=linspace(phintd(1),phintd(2),10000);
  K(1)=trapeze(phin,feval(cc,phin));
  K(2)=trapeze(phin,feval(ss,phin));
  K(3)=trapeze(phin,feval(sc,phin));
  
  if yesplot==1
    ph=linspace(min(0,phintd(1)),max(2*pi,phintd(2)),10000);
    
    clf
    ah=krijetem(subnum(3,1));
    
    axes(ah(1))
    p(1)=plot(ph,feval(cc,ph)); grid on; hold on
    yl=ylim; f(1)=fillbox([phintd yl(2) yl(1)]);
    yb(1)=ylabel(sprintf('cos(m_1%s)cos(m_2%s)','\phi','\phi'));
    
    axes(ah(2))
    p(2)=plot(ph,feval(ss,ph)); grid on; hold on
    yl=ylim; f(2)=fillbox([phintd yl(2) yl(1)]);
    yb(2)=ylabel(sprintf('sin(m_1%s)sin(m_2%s)','\phi','\phi'));
    
    axes(ah(3))
    p(3)=plot(ph,feval(sc,ph)); grid on; hold on
    yl=ylim; f(3)=fillbox([phintd yl(2) yl(1)]);
    yb(3)=ylabel(sprintf('sin(m_1%s)cos(m_2%s)','\phi','\phi'));
    
    set(p,'LineW',2); longticks(ah,2)
    for ind=1:3
      axes(ah(ind))
      set(ah(ind),'Children',[p(ind) ; f(ind)]); 
      t(ind)=title(sprintf(...
	  'I_%i= %8.3f ; J_%i= %8.3f ; K_%i= %8.3f ; |%s_{JI}|= %8.3e',...
	  ind,Id(ind),ind,J(ind),ind,K(ind),'\Delta',abs(J(ind)-Id(ind))));
    end 
    set(ah,'xlim',[min(0,phintd(1)) max(2*pi,phintd(1))])
    nolabels(ah(1:2),1)
    axes(ah(3))
    xlabel(sprintf('m_1= %i ; m_2= %i ; GL= %i',m1,m2,GL))
  end
  if any(abs(J(:)-Id(:))>1e-12); error('Something wrong') ; end
  % Make decent output
  if nargout==1
    varargout{1}=abs(J(:)-Id(:));
    varargout{2}=GL1;
  end
end
