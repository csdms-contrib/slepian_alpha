function [PpL,PLm1,PL,PLp1]=legendrediff(L,x,norma)
% [PpL,PLm1,PL,PLp1]=legendrediff(L,x,norma)
%
% Computes the derivative of the Legendre polynomial with respect to its
% argument x=cos(theta). Compare LIBBRECHT, where it is wrt acos(x).
%
% INPUT: 
%
% L              Degree of spherical harmonic < 256
% x              Evaluation point(s)
% norma          'sch' Schmidt-normalized [default]
%                'fnr' Fully normalized real
%
% OUTPUT:
%
% PpL            Schmidt-normalized derivative of P_L(x), m=0
% PLm1           Schmidt-normalized polynomial P_{L-1}(x), m=0
% PL             Schmidt-normalized polynomial P_L(x), m=0
% PLp1           Schmidt-normalized polynomial P_{L+1}(x), m=0
%
% See Wolfram under Legendre-Gauss Quadrature.
%
% EXAMPLE:
%
% legendrediff('demo1') % Comparison with LIBBRECHT and DIFF
%
% SEE ALSO: LIBBRECHT, YLM
%
% Last modified by fjsimons-at-alum.mit.edu, 05/17/2011

defval('L','demo1')
defval('norma','sch')

if ~isstr(L)
  switch norma
   case 'sch'
    fac1=1;
    fac2=1;
    fac3=1;
   case 'fnr'
    fac1=sqrt(2*(L-1)+1);
    fac2=sqrt(2*L+1);
    fac3=sqrt(2*(L+1)+1);
   otherwise
    error('Specify valid normalization')
  end
  % Must build in what it means to be -1
  PLm1=rindeks(legendre(L-1,x,'sch'),1);  
  PL=rindeks(legendre(L,x,'sch'),1);
  warning off
  PpL=L*(x(:).*PL(:)-PLm1(:))./(x(:).^2-1);
  warning on

  % From Boyd (2001)
  PpL(x==1)=L*(L+1)/2;
  PpL(x==-1)=(-1)^(L-1)*L*(L+1)/2;
  PpL=fac2*PpL(:)';

  PLm1=PLm1*fac1;
  PL=PL*fac2;

  if nargout==4
    PLp1=fac3*rindeks(legendre(L+1,x,'sch'),1);  
  end
elseif strcmp(L,'demo1')
  clf
  theta=linspace(0,pi,500);
  x=cos(theta); more off
  deg=round(rand*20);deg=9
  m=0;
  p=rindeks(legendre(deg,x,'sch'),1);
  [pp,jk,p2]=legendrediff(deg,x,'sch');
  [pl,dpl]=libbrecht(deg,x,'sch',[]);
  pl=rindeks(pl,1);
  dpl=rindeks(dpl,1);
  subplot(211)
  plot(x,p,'b-','LineW',2); hold on
  plot(x,pl,'y-','LineW',1); 
  if max(abs(p2(:)-p(:)))>0; error('Something wrong?'); end
  title(sprintf(...
      'Legendre functions (Schmidt); l=%i (m= %i)',...
      deg,m),'FontS',12)  
  grid on; axis tight ; openup(gca,6); 
  nolabels(gca,1)
  yl=ylabel('N_l^m\timesP_l^m(cos\theta)');
  movev(gca,-.1)
  l1=legend('LEGENDRE','LIBBRECHT');
  longticks(gca,2)
  subplot(212)
  plot(x,pp,'r-','LineW',2); 
  hold on
  % Watch out since x is not equally spaced.
  plot(x(2:end)-indeks(diff(x),1)/2,...
       diff(p)./diff(x),'y')
  warning off
  plot(x,-dpl./sin(theta),'k--','LineW',1);   
  warning on
  yl=ylabel('dP_l^m(cos(\theta))/dcos\theta');
  xl=xlabel('cos(\theta)');
  axis tight ; grid on; openup(gca,6)
  l=legend('LEGENDREDIFF','DIFF','LIBBRECHT');
  longticks(gca,2)
  fig2print(gcf,'portrait'); id; axes(l)
  figdisp
end
