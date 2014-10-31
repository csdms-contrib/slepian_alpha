function [HNP,HSP,epsi,m2l2]=hilbxlm(l,m,X,norma)
% [HNP,HSP,epsi]=HILBXLM(l,m,X,norma)
%
% Hilb's uniformly valid asymptotic expressions for the fully normalized
% associated Legendre functions, by Dahlen (1980).
% We get two; one near the North pole and one near the South pole.
%
% INPUT:
%
% l       Angular degree
% m       Angular order
% X       Argument X=cos(theta), so -1 < X < 1
% norma   'sch' Schmidt normalization
%         'fnr' Fully normalized real
%         'fnc' Fully normalized complex
%
% OUPUT:
%
% HNP     North Pole approximation, valid for [0->pi-epsi] if m==0
% SNP     South Pole approximation, valid for [epsi->pi] if m==0
% epsi    Validity parameter; the error in this range is of O(m^2/l^2)
% m2l2    The value of m^2/l^2
%
% EXAMPLE:
%
% hilbxlm('demo1')
%
% See also BACKUS and LIBBRECHT
%
% Last modified by fjsimons-at-alum.mit.edu, Feb 12th, 2004

if ~isstr(l)
  if any(X>=1 | X<=-1)
    error('X must contain real values between -1 < X < 1')
  end
  if m>l
    error('The order m must be smaller than the degree l')
  end

  lbd=l+1/2; epsi=1/lbd; m2l2=m^2/l^2;
  k12=sqrt(lbd/2/pi);
  
  X=acos(X(:));

  % Extra normalization factors
  switch norma
   case 'sch'
    fac=(-1)^m/sqrt((2*l+1)/4/pi)*(sqrt(2-(m==0)));
   case 'fnc'
    fac=1;
   case 'fnr'
    fac=(-1)^m*sqrt(4*pi)*(sqrt(2-(m==0)));
   otherwise
    error('Specify valid normalization')
  end

  HNP=fac*(-1)^m*k12*sqrt(X./sin(X)).*besselj(m,lbd*X);
  HSP=fac*(-1)^l*k12*sqrt((pi-X)./sin(pi-X)).*besselj(m,lbd*(pi-X));  

elseif strcmp(l,'demo1')
  deg=round(rand*30);deg=23
  ord=round(rand*deg/2); % Works best for (m/l)<<1
  ord=4
  % For the cosine, just doing eps is not enough
  theta=linspace(0+1e-7,pi-1e-7,500);
  X=cos(theta);
  norms='sch';
  P=rindeks(libbrecht(deg,X,norms,[],ord),1);
  [NP,SP,epsi,m2l2]=hilbxlm(deg,ord,X,norms);
  yli=minmax([P(:)]); 
  if yli(1)==yli(2)
    yli=yli+[-1 1]*1e-0;
  end
  clf
  ah(1)=subplot(221);
  r(1)=plot(acos(X),P); hold on
  a(1)=plot(acos(X),NP,'y-');
  axis tight; grid on
  ylim(yli*1.1)
  xlim([0 pi])
  pp(1)=plot([pi-epsi pi-epsi],yli*1.1,'k');
  yl(1)=ylabel('X_l^m(\theta)=N_l^m\timesP_l^m(cos\theta)');
  pilabels(gca)
  nolabels(gca,1)
  ah(2)=subplot(222);
  r(2)=plot(acos(X),P); hold on
  a(2)=plot(acos(X),SP,'y-');
  axis tight; grid on
  ylim(yli*1.1)
  xlim([0 pi])
  pp(2)=plot([epsi epsi],yli*1.1,'k');
  yl(2)=ylabel('X_l^m(\theta)=N_l^m\timesP_l^m(cos\theta)');
  movev(ah(1:2),-.1)
  l=legend('LIBBRECHT','HILB');
  pilabels(gca)
  nolabels(gca,1)
  ah(3)=subplot(223);
  ylo=[1e-5 5];
  e(1)=plot(acos(X),abs(P(:)-NP(:)));
  xl(1)=xlabel('\theta');
  yl(3)=ylabel('Absolute Error');
  xlim([0 pi]); grid on; hold on
  pp(3)=plot([pi-epsi pi-epsi],ylo,'k');
  pp(4)=plot([0 pi],[m2l2 m2l2],'k');
  pilabels(gca)
  pilabels(gca)
  ah(4)=subplot(224);
  e(2)=plot(acos(X),abs(P(:)-SP(:)));
  xlim([0 pi]); grid on; hold on
  pp(5)=plot([epsi epsi],ylo,'k');
  pp(6)=plot([0 pi],[m2l2 m2l2],'k');
  xl(2)=xlabel('\theta');
  yl(4)=ylabel('Absolute Error');
  pilabels(gca)
  st=supertit(ah(1:2),sprintf(...
	'Associated Legendre functions (Schmidt); l=%i, m= %i',...
	deg,ord),12)
  movev(st,.15)
  set([r e],'LineW',2)
  set(ah(3:4),'yscale','log','ylim',ylo)
  set(pp,'Color',grey)
  longticks(ah,2)
  fig2print(gcf,'portrait'); id; axes(l)
  figdisp
end
