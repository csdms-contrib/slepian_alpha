function [E,V,th,C]=sdwcapt2(TH,L,m,nth,vcut,grd)
% [E,V,th,C]=SDWCAPT2(TH,L,m,nth,vcut,grd)
%
% Spherical harmonic localization to DOUBLE spherical polar cap:
% spatially concentrated and optimally bandlimited solutions.
%
% INPUT:
%
% TH          Angular extent of the spherical cap, in degrees
% L           Bandwidth
% m           Angular order of the required data window, -l<m<l
% nth         Number of colatitudes to evaluate
% vcut        Cut-off eigenvalue [default= eps*10]
% grd         1 Colatitudes only; returns matrix E
%             2 Colatitude/Longitude; returns cell E
%
% OUTPUT:
%
% E           Optimally concentrated tapers, expanded to space domain 
% V           Eigenvalues, sorted
% N           Sum of all the eigenvalues
% th          Colatitudes at which the functions are evaluated, in degrees
% C           Optimally concentrated tapers; coefficients
%
% Last modified by fjsimons-at-alum.mit.edu, 10/03/2008

defval('TH',40)
defval('L',18)
defval('m',0)
defval('nth',720)
defval('vcut',eps*10)
defval('grd',1)

% Work with the absolute value of m
mor=m;
m=abs(m);

if(m>L)
  error('Order cannot exceed degree')
end

% Filename of saved data
fnpl=sprintf('%s/SDW-%i-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDWCAPT2'),TH,L,nth,m);

if exist(fnpl,'file')==2
  eval(sprintf('load %s',fnpl))
  disp(sprintf('%s loaded by SDWCAPT2',fnpl))
else 
  [E,V,th,C,ngl1,ngl2,K]=sdwcap2(TH,L,m,nth,vcut,grd);
  % Use the "trick": expand truncated basis functions
  % Put on the black lines
  belt=[th>TH]&[th<180-TH];
  E(belt,:)=0;
  if m==0
    % Expand to Nyquist; watch for the factors!
    C2=th2pl(E,size(E,1)-1)*sqrt(4*pi);
  else
    C2=th2plm(E,size(E,1)-1,m);
  end
  % Don't need to do this, but it should return the zeroed input E!
  % Reexpand
  % E2=pl2th(C2,[],1); % Watch for the factors!
  % Now E with the zeroed portions and E2 should be identical
  % and the C2's are still orthogonal
  C=C2;
  % Output in degrees
  eval(sprintf('save %s E V L TH C th',fnpl))
end

if grd==2
  % Output on full grid
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
