function [K11,K21,K22]=cdlkernel(x1,x2,L,m,spd)
% [K11,K21,K22]=CDLKERNEL(x1,x2,L,m,spd)
%
% Calculates a Christoffel-Darboux kernel with Legendre functions
% whose harmonics are normalized to 4pi (option 'fnr' elsewhere).
% Kernels are usually calculated at Gauss-Legendre integration points, at
% a mixture of GL points and a complete set, and at the complete set.
%
% INPUT:
%
% x1        Vector with first set of points
% x2        Vector with second set of points
% L         Maximum degree of the expansion
% m         Single angular order (set to zero)
% spd       1 Fast method using the Christoffel-Darboux identity
%           2 Slow method involving all L+1 summations
%
% OUTPUT:
%
% K11       The first kernel at point set x1 (GL points)
%           \suml_{l=0}^L P_{l}(x1)P_{l}(x1)
% K21       The interpolating kernel mixing x1 and x2
%           \suml_{l=0}^L P_{l}(x2)P_{l}(x1)
% K22       The full kernel at point set x2 (complete set)
%           \suml_{l=0}^L P_{l}(x2)P_{l}(x2)
%
% EXAMPLE:
%
% cdlkernel('demo1')     % Compares the two methods
%
% Last modified by fjsimons-at-alum.mit.edu, June 3rd, 2004

if ~isstr(x1)
  defval('m',0')
  defval('spd',1)
  if m~=0
    spd=2;
  end

  switch spd
   case 1
    % Using the Christoffel-Darboux identity
    % First, get P'(L+1), P(L) and P(L+1)
    % Since we're dividing polynomials don't worry about sqrt(2l+1)
    [Pdp1,P1,Pp1]=legendrediff(L+1,x1,'sch');
    [Pdp2,P2,Pp2]=legendrediff(L+1,x2,'sch');
    % Then, also get P'(L)
    Pd1=legendrediff(L,x1,'sch');
    Pd2=legendrediff(L,x2,'sch');

    X1=x1(:)*ones(1,length(x1)); X2=x2(:)*ones(1,length(x2));
    X21=x2(:)*ones(1,length(x1)); X12=x1*ones(1,length(x2));

    % For the non-diagonal elements
    warning off
    K11=(L+1)*(Pp1(:)*P1(:)'-P1(:)*Pp1(:)')./(X1-X1');
    K22=(L+1)*(Pp2(:)*P2(:)'-P2(:)*Pp2(:)')./(X2-X2');
    K21=(L+1)*(Pp2(:)*P1(:)'-P2(:)*Pp1(:)')./(X21-X12');
    warning on
    % For the diagonal elements
    K11(ondiag(K11))=(L+1)*(Pdp1.*P1-Pd1.*Pp1);
    K22(ondiag(K22))=(L+1)*(Pdp2.*P2-Pd2.*Pp2);

    if ~~sum(sum((X21-X12')==0))
      error('Diagonal elements not properly accounted for');
    end
   
   case 2
    % Straightforward brute-force calculation
    % Initialize kernels
    K11=zeros(length(x1),length(x1));
    K21=zeros(length(x2),length(x1));
    K22=zeros(length(x2),length(x2));
    for l=m:L
      % Make sure Plm(1)=sqrt(2l+1) so Ylm normalized to 4\pi
      Plm=(rindeks(legendre(l,x1(:)','sch')*sqrt(2*l+1),m+1));
      Plmint=(rindeks(legendre(l,x2(:)','sch')*sqrt(2*l+1),m+1));
      % Kernel at Gauss-Legendre points
      K11=K11+Plm(:)*Plm(:)';
      % Kernel to go from Gauss-Legendre to full resolution
      K21=K21+Plmint(:)*Plm(:)';
      % Kernel at full resolution
      K22=K22+Plmint(:)*Plmint(:)';
    end
  end
elseif strcmp(x1,'demo1')
  defval('th0',40)
  defval('m',0)
  defval('SN',5)
  defval('nth',720)
  nth0=ceil(th0/180*nth);
  TH=linspace(0,th0/180*pi,nth0);
  xint=cos(TH);
  th0=th0*pi/180;
  index=0;
  els=[0:50];
  for L=els
    index=index+1;
    [w,xGL]=gausslegendrecof(2*L,[],[cos(th0) 1]);
    tic
    [KGL1,Kint1,K1]=cdlkernel(xGL,xint,L,m,1);
    tic1(index)=toc;
    tic
    [KGL2,Kint2,K2]=cdlkernel(xGL,xint,L,m,2);
    tic2(index)=toc;
    erro(index)=mean(abs(K1(:)-K2(:)));
  end
  ah(1)=subplot(121);
  p{1}=plot(2*els,erro,'+'); grid on; title('Accuracy')
  xlabel('Product degree'); ylabel('Difference')
  ah(2)=subplot(122); 
  p{2}=plot(2*els,tic1,'bo',2*els,tic2,'kv');
  title('Time Cost'); xlabel('Product degree') ; ylabel('Seconds')
  yll=ylim; ylim([0 yll(2)])
  set(p{1},'MarkerF','r','MarkerE','r')
  set(p{2}(1),'MarkerF','b','MarkerE','b')
  set(p{2}(2),'MarkerF','r','MarkerE','r')
  set(ah(2),'YScale','Log')
  l=legend('Christoffel-Darboux formula','Complete summation',4);
  grid on
  longticks(ah)
  figdisp
else
  error('Specify valid option')
end



