function [I,dph,th,thR,th0,ph0]=sinsin(th,m1,m2,dom,pars)
% [I,dph,th,nthR,th0,ph0]=sinsin(th,m1,m2,dom,pars)
%
% Calculates the longitudinal integral 
% \int_{A}^{B}\sin(m1\phi)\sin(m2\phi)\,d\phi
% over a specified domain at some (Gauss-Legendre) integration point(s).
%
% INPUT:
%
% th         Colatitude(s) that you may need to define the patch
% m1,m2      Angular order of the cosine factors
% dom        A matrix with the hatchings from PHICURVE, or
%            'patch'   Spherical patch [default] with the following specs
%            'england' England, Scotland and Wales
% pars       [th0,ph0,thR] for 'patch', with 
%                 th0  Colatitude of the cap center, in radians
%                 ph0  Longitude of the cap center, in radians
%                 thR  Radius of the cap, in radians
%            N for 'england' with N the smoothness of the splining
%
% EXAMPLE:
%
% sinsin('demo1');
%
% Last modified by fjsimons-at-alum.mit.edu, 04/20/2009

defval('th',[])
defval('m1',10);
defval('m2',4);
defval('dom','patch')

if isstr(dom)
  switch dom
   case 'patch'
    defval('pars',[pi/4 pi/2 pi/9]);
    th0=pars(1); ph0=pars(2); thR=pars(3);
    defval('th',linspace(th0-thR,th0+thR,100));
    phint=dphpatch(th,thR,th0,ph0);
   case 'england'     
    defval('N',10)
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
	   sin((m1-m2)*B)*m2-...
	   sin((m1+m2)*B)*m1+...
	   sin((m1+m2)*B)*m2-...
	   sin((m1-m2)*A)*m1-...
	   sin((m1-m2)*A)*m2+...
	   sin((m1+m2)*A)*m1-...
	   sin((m1+m2)*A)*m2)/...
      ((m1-m2)*(m1+m2));
  elseif m1==-m2 & m1~=0
    m=m1;
    I=I-1/2*(-cos(m*B).*sin(m*B)+...
	    m*(B-A)+cos(m*A).*sin(m*A))/m;
  elseif m1~=0
    m=m1;
    I=I+1/2*(-cos(m*B).*sin(m*B)+...
	   m*(B-A)+cos(m*A).*sin(m*A))/m;
  else
    I=I+0;
  end
end

