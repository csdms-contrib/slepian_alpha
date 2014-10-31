function [knum,av]=isav(mat,kx,ky)
% [knum,av]=isav(mat,kx,ky)
%
% Makes isotropic average in wave number domain.
% ZERO is in the center of the matrix 'mat'.
%
% Like RADAVG but for actual wave numbers, i.e. bounded by the lowest
% wavenumber (the wavelength associated with the physical box size).
%
% [knum,av]=isav(data,linspace(-1/2,1/2,size(data,2)),linspace(-1/2,1/2,size(data,1)));
%
% Used last in LOADING
%
% SEE ALSO:
%
% RADAVG
%
% EXAMPLE:
%
% isav('demo');
%
% Last modified by fjsimons-at-alum.mit.edu, 08/09/2013

if ~isstr(mat)
  mlk=min(length(kx),length(ky));
  % This is where you want it, will be used for interpolation.
  % knum=linspace(0,min(kx(end),ky(end)),mlk/2);
  % Actually, make sure the Nyquist is IN it - always from the left
  knum=linspace(0,min(max(abs(kx)),max(abs(ky))),mlk/2);
  
  kth=indeks(linspace(0,2*pi,mlk*2+1),1:mlk*2);

  [KNUM,KTH]=ndgrid(knum,kth);
  kxi=KNUM.*cos(KTH);
  kyi=KNUM.*sin(KTH);

  % The original kx and ky are correct and include 0
  zi=interp2(kx,ky,mat,kxi,kyi);

  % Take mean over all angles at constant radius
  av=nanmean(zi');
elseif strmatch(mat,'demo')
  clf
  m=50; nlev=0.4;
  AA=(1-nlev*rand(1,m)).*sin(linspace(pi/2,pi,m));
  subplot(121) ; plot(AA)
  revolved=revolve(AA);
  subplot(122) ; surf(revolved) ; shading flat
  for index=1:100    
    revolv2(:,index)=interp(revolved(:,index),2); 
  end
  [knum,av]=isav(revolv2,linspace(-50,50,size(revolv2,2)),...
		 linspace(-50,50,size(revolv2,1)));
  subplot(121); hold on
  plot(knum,av,'r-')
end


