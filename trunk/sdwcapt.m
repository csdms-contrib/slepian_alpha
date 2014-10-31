function [E,V,N,th,C,E0]=sdwcapt(TH,L,m,nth,vcut,grd)
% [E,V,N,th,C,E0]=SDWCAPT(TH,L,m,nth,vcut,grd)
%
% Spherical harmonic localization to a spherical polar cap:
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
% E0          Tapers with zeros added to make them of length nth
% 
% SEE ALSO:
%
% CDLKERNEL
%
% Last modified by fjsimons-at-alum.mit.edu, 05/21/2009

defval('TH',40)
defval('L',18)
defval('m',0)
defval('nth',720)
defval('vcut',eps*10)
defval('grd',1)

% Work with the absolute value of m
mor=m;
m=abs(m);

fnpl=sprintf('%s/SDW-%i-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDWCAPT'),TH,L,nth,m);

if exist(fnpl,'file')==2 % & 1==3
  eval(sprintf('load %s',fnpl))
  disp(sprintf('%s loaded by SDWCAPT',fnpl))
else 
  N=(L+1)^2*(1-cos(TH/180*pi))/2;
  
  % Approximate Nyquist degree is $\pi/d\theta$
  if nth < (L+1) & nth~=0
    error('Sample finer to avoid aliasing')
  end

  % Regular sampling at number of requested points
  nTH=ceil(TH/180*nth);
  % But you sample only until TH
  th=linspace(0,TH/180*pi,nTH);

  % Convert to radians
  thor=th*180/pi;
  TH=TH*pi/180;

  % The degree of the integrand will be 2*Lmax
  % Gauss-Legendre integration points and weights scaled on the interval
  % [cos(TH) 1] after transformation of variables.
  ngl=max(2*L,200*(m~=0));
  % If you want to be fancy do this
  % [jk1,V2,jk2,jk3,jk4,ngl1,ngl2]=sdwcap(TH*180/pi,L,m);
  % ngl=ngl2
  [w,xGL,nsel]=gausslegendrecof(ngl,[],[cos(TH) 1]);

  disp(sprintf('Gauss-Legendre with %i points',nsel))

  % Later on, we will need points over the entire interval from cos([0,pi])
  % Interpolation points at the requested number of TH's
  xint=cos(th);

  % Construct the Christoffel-Darboux kernel
  [KGL,Kint]=cdlkernel(xGL,xint,L,m);
  
  % The interval rolls out for Gauss-Legendre quadrature
  % but you need the scaling; contains longitudinal integral
  normz=(1+(m==0))/4;
  KGL=normz*KGL;
  Kint=normz*Kint;
  
  % Calculate eigenvectors and values of symmetrized kernel
  [EGL,V]=eig(diag(sqrt(w))*KGL*diag(sqrt(w)));
  
  % Need to unscale E
  EGL=EGL./repmat(sqrt(w(:)),1,size(EGL,2));
     
  % With this, inner product over domain becomes V; no dot!
  EGL=EGL*sqrt(V)/sqrt((1+(m==0))*pi);

  % Should you ever want to plot E at the GL values.
  th=acos(xGL);

  % Order eigenvalues and eigenfunctions downgoing
  [V,isrt]=sort(sum(V,1),'descend');
  % V=fliplr(V); EGL=EGL(:,fliplr(isrt));
  EGL=EGL(:,isrt);
  
  % We know how many there should be
  EGL=EGL(:,1:L-m+1);
  V=V(:,1:L-m+1);

  % Only return nonzero and real eigenvalues (order matters)
  EGL=EGL(:,V>vcut & imag(V)==0);
  V=V(V>vcut & imag(V)==0);

  % But usually we "interpolate" them using the same kernel at the desired
  % values.
  E=(Kint*diag(w))*EGL*diag(1./V);

  % Make E start with a positive lobe
  for index=1:size(E,2)
    E(:,index)=E(:,index)*sign(E(2,index));
  end

  % Verify automatic orthonormality over restricted interval
  % orv1=w(:)'*[EGL.^2]*(1+(m==0))*pi;
  % Calculate the latitudinal and longitudinal integral
  orv2=[EGL'*diag(w)*EGL]*(1+(m==0))*pi; 
  disp(sprintf('Orthonormality: Mean absolute error %8.3e',...
	       mean(mean(abs(orv2-diag(V)))))) 

  % Output coefficients C
  E0=[E ; zeros(nth-size(E,1),size(E,2))];  
 
  if m==0    
    % These are normalized to the AREA of the unit sphere
    [C,rnk]=th2pl(E0,[],'im');
    % Now they're normalized to unit area for the unit sphere
    C=C*sqrt(4*pi);
  else
    % Make it -2 today, August 11th, 2004
    [C,rnk]=th2plm(E0,nth-2,m);
    C=[C ; repmat(NaN,1,size(C,2))];
  end
  nlon=ceil(2*nth-1);

  % Output in degrees
  TH=TH*180/pi;
  th=thor;
  eval(sprintf('save %s E V L N TH C th nlon E0',fnpl)) 
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
