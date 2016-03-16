function [ap1,th0]=dahlen(L,m,X,norma)
% [ap1,th0]=DAHLEN(L,m,X,norma)
%
% Dahlens JWKB asymptotic approximation to the associated Legendre
% polynomials of large degree L.
%
% INPUT:
%
% L          Degree of the associated Legendre polynomial
% m          Order of the associated Legendre polynomial
% X          Argument X=cos(theta), so -1 < X < 1
% norma      'sch' Schmidt-normalized polynomials 
%            'fnr' Fully-normalized real
%            'fnc' Fully-normalized complex
%
% OUTPUT:
%
% ap1        The approximation due to DT (B.84)
% th0        The validity range is [th0 pi-th0]
%
% See Dahlen and Tromp (1998), Theoretical Global Seismology,
% DT (X.NN) refer to their numbered equations.
%
% EXAMPLE:
%
% dahlen('demo1')
%
% See also BACKUS, HILBXLM.
%
% Last modified by fjsimons-at-alum.mit.edu, 03/16/2016

if ~isstr(L)
  defval('norma','fnc')

  % DT B.80)
  th0=asin(abs(m)/(sqrt(L*(L+1))));

  if any(X>=1 | X<=-1)
    error('X must contain real values between -1 < X < 1')
  end
  if m>L
    error('The order m must be smaller than the degree l')
  end

  th=acos(X(:));

  % Extra normalization factors
  switch norma
   case 'sch'
    fac=(-1)^m/sqrt((2*L+1)/4/pi)*(sqrt(2-(m==0)));
   case 'fnc'
    fac=1;
   case 'fnr'
    fac=(-1)^m*sqrt(4*pi)*(sqrt(2-(m==0)));
   otherwise
    error('Specify valid normalization')
  end

  % The approximation DT (B. 84)
  warning off
  ap1=fac/pi*(sin(th).^2-sin(th0)^2).^(-1/4).*...
      cos((L+1/2).*acos(cos(th)/cos(th0))-...
	   m*acos(cot(th)/cot(th0))+m*pi-1/4*pi);
  warning on
  % This is only real withing th0 and pi-th0
  ap1(abs(imag(ap1))>0)=NaN;
  
elseif strcmp(L,'demo1')
  deg=round(rand*30);
  ord=round(rand*deg);
  % For the cosine, just doing eps is not enough
  theta=linspace(0+1e-7,pi-1e-7,500);
  X=cos(theta);
  norms='sch';
  P=rindeks(libbrecht(deg,X,norms,[],ord),1);
  [apd,th0]=dahlen(deg,ord,X,norms);
  yli=minmax([P(:)]); 
  if yli(1)==yli(2)
    yli=yli+[-1 1]*1e-0;
  end
  clf
  ah(1)=subplot(211);
  r(1)=plot(acos(X),P); hold on
  a(1)=plot(acos(X),apd,'y-');
  movev(ah(1),-.1)
  l=legend('LIBBRECHT','DAHLEN');
  axis tight; grid on
  ylim(yli*1.1)
  xlim([0 pi])
  pp(1)=plot([th0 th0],yli*1.1,'k');
  pp(2)=plot([pi-th0 pi-th0],yli*1.1,'k');
  yl(1)=ylabel('X_l^m(\theta)=N_l^m\timesP_l^m(cos\theta)');
  st=title(sprintf(...
      'Associated Legendre functions (Schmidt); l=%i, m= %i',...
      deg,ord),'FontS',12)
  pilabels(ah(1))
  nolabels(gca,1)
  ah(2)=subplot(212);
  ylo=[1e-5 5];
  e(1)=plot(acos(X),abs(P(:)-apd(:)));
  xl(1)=xlabel('Colatitude (\theta)');
  yl(2)=ylabel('Absolute Error');
  xlim([0 pi]); grid on; hold on
  pp(3)=plot([th0 th0],ylo,'k');
  pp(4)=plot([pi-th0 pi-th0],ylo,'k');
  set([r e],'LineW',2)
  set(ah(2),'yscale','log','ylim',ylo)
  set(pp,'Color',grey)
  pilabels(ah(2))
  longticks(ah,2)
  fig2print(gcf,'portrait'); id
  figdisp
end
