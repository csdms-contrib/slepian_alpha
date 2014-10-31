function [E,V,th,C,ngl1,ngl2,K]=sdwcap2(TH,L,m,nth,vcut,grd,eo)
% [E,V,th,C,ngl1,ngl2,K]=SDWCAP2(TH,L,m,nth,vcut,grd,eo)
%
% Spherical harmonic localization to a DOUBLE spherical polar cap: 
% bandlimited and optimally spatially concentrated solutions.
%
% INPUT:
%
% TH          Angular extent of both spherical caps, in degrees
% L           Bandwidth (maximum angular degree)
% m           Angular order of the required data window, -l<m<l
% nth         Number of colatitudes to evaluate
%             If nth=0 then only does the eigenvalues. 
% vcut        Cut-off eigenvalue [default= eps*10]
%             One must realize that calculating eigenfunctions with
%             near-zero eigenvalues is next to impossible to do
%             correctly; what we get is an arbitrary orthogonal system
% grd         1 Colatitudes only; returns matrix E [default]
%             2 Colatitude/Longitude; returns cell E
% eo          1 even/odd eigenvalue problems solved separately [default]
%             2 regular algorithm treating all degrees together
%
% OUTPUT:
%
% E           Optimally concentrated tapers, expanded to space domain 
% V           Eigenvalues, sorted
% th          Colatitudes at which the functions are evaluated, in degrees
% C           Optimally concentrated tapers, spherical harmonics
%             coefficients, normalized to unity on the unit sphere
% ngl1, ngl2  Number of GL points on unit sphere and domain, respectively
% K           The matrix that is diagonalized to produce the eigenfunctions
%
% See also DOUBLECAP, GRUNBAUM2
%
% Last modified by fjsimons-at-alum.mit.edu, 03/27/2009

defval('TH',30)
defval('L',18)
defval('m',0)
defval('nth',720)
defval('vcut',eps*10)
defval('grd',1)
defval('eo',1)

% Work with the absolute value of m
mor=m;
m=abs(m);

if(m>L)
  error('Order cannot exceed degree')
end

% Filename of saved data
fnpl=sprintf('%s/SDW-%i-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDWCAP2'),TH,L,nth,m);

if exist(fnpl,'file')==2  && ~(L==180 & m==0) && vcut>0
  eval(sprintf('load %s',fnpl))
%  disp(sprintf('%s loaded by SDWCAP2',fnpl))
else
  % Work with the results applicable to the SINGLE cap
  [E,V,N,th,C,ngl1,ngl2,unc,com,sdl,K]=sdwcap(TH,L,m,nth,vcut,grd);
  % And construct the double cap from the symmetry relation
  Ks=meshgrid(1:L-m+1);
  K(~~mod(Ks'+Ks,2))=0;
  % But it's either 0 or 2
  K=K*2;

  if eo==1
    % Split in even and odd degrees; this guarantees
    % exactly even/odd results, which is preferable
    lodd=even(L-m+1);
    Ke=K(lodd,lodd);
    Ko=K(~lodd,~lodd);
    [Ce,Ve]=eig(Ke);
    [Co,Vo]=eig(Ko);
    C=repmat(0,L-m+1,L-m+1);
    C(lodd,lodd)=Ce;
    C(~lodd,~lodd)=Co;
    V=repmat(0,L-m+1,1);
    V(lodd)=diag(Ve);
    V(~lodd)=diag(Vo);
    V=diag(V);
  else 
    [C,V]=eig(K);
  end
  
  % Convert to radians
  TH=TH*pi/180;

  [ngl1,ngl2]=orthocheck(C,V,TH,m,2);

  % Order eigenvalues and eigenfunctions downgoing
  [V,isrt]=sort(sum(V,1),'descend');
  % V=fliplr(V); C=C(:,fliplr(isrt));
  C=C(:,isrt);

  % Only return nonzero "useful" eigenvalues
  C=C(:,V>vcut); 
  V=V(V>vcut);

  % Compute spatial functions, colatitudinal part only
  if nth~=0
    % Zonal functions only 
    if m==0
      % Make spatial functions
      % This is SDW (2005) equation (5.7) combined 
      % with the sqrt(2-dom) of (5.9) already included!
      [E,th]=pl2th(C,nth,1);
      th=th*180/pi;
      nlon=2*nth-1;
    else
      % This is SDW (2005) equation (5.7) combined 
      % with the sqrt(2-dom) of (5.9) already included!
      [E,nlon,lat]=plm2th(C,nth,m,1);
      th=linspace(0,180,size(E,1));
    end    
    % Make E start with a positive lobe and ajust C too
    % Don't take first sample as for m~=0 it is 0
    for index=1:size(E,2)
      C(:,index)=C(:,index)*sign(E(2,index));
      E(:,index)=E(:,index)*sign(E(2,index));
    end
  else
    E=0;
    th=0;
    nlon=0;
  end
  % Output in degrees
  save(fnpl,'E','V','L','N','TH','C','th','nlon','ngl1','ngl2','K')
end

if nth~=0 & grd==2
  % Output on full grid; watch the sign of m
  if mor<=0
    EE=E; clear E
    for index=1:size(EE,2)
      E{index}=EE(:,index)*cos(m*linspace(0,2*pi,nlon));
    end
  end
  if mor>0
    EE=E; clear E
    for index=1:size(EE,2)
      E{index}=EE(:,index)*sin(m*linspace(0,2*pi,nlon));
    end
  end
end
