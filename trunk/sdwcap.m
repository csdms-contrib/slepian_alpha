function [E,V,N,th,C,ngl1,ngl2,unc,com,sdl,K]=sdwcap(TH,L,m,nth,vcut,grd,method)
% [E,V,N,th,C,ngl1,ngl2,unc,com,sdl,K]=SDWCAP(TH,L,m,nth,vcut,grd,method)
%
% Spherical harmonic localization to a SINGLE fixed-order spherical polar
% cap: bandlimited (lowpass or bandpass) and optimally spatially
% concentrated solutions.
%
% INPUT:
%
% TH          Angular extent of the spherical cap, in degrees
% L           Bandwidth (maximum angular degree), or passband (two degrees)
% m           Angular order of the required data window, -l<m<l
% nth         Number of colatitudes to evaluate (none if 0) [default: 720]
% vcut        Cut-off eigenvalue [default= eps*10]
%             If vcut==-1 automatically recalculates everything
% grd         1 Colatitudes only; returns matrix E
%             2 Colatitude/Longitude; returns cell E
% method      'gl' for Gauss-Legendre via LEGENDREPRODINT [default]
%             'paul' for Paul via PAUL 
%             'siam' should write one here exactly as SDW (2006) eq. (5.5)
%
% OUTPUT:
%
% E           Optimally concentrated tapers, expanded to space domain 
% V           Eigenvalues, sorted
% N           Sum of ALL the eigenvalues (Shannon number for ALL orders)
%             compare with the sum that you got using NSUBM(N,m)...
% th          Colatitudes at which the functions are evaluated, in degrees
% C           Optimally concentrated tapers, spherical harmonics
%             coefficients, normalized to unity on the unit sphere
% ngl1, ngl2  Number of GL points on unit sphere and domain, respectively
% unc         Uncertainty product
% com         Center of mass
% sdl         Spherical harmonics standard deviation
% K           The matrix that is diagonalized to produce the eigenfunctions
%
% See also DOUBLECAP, GRUNBAUM, KERNELC
%
% Last modified by plattner-at-princeton.edu, 05/25/2011
% Last modified by fjsimons-at-alum.mit.edu, 05/29/2012

% For bandpass and possibly Grunbaum, see Morrison63 and SenGupta+2012

defval('TH',40)
defval('L',18)
defval('m',0)
defval('nth',720)
defval('vcut',eps*10)
defval('grd',1)
defval('method','gl')

% Work with the absolute value of m
mor=m;
m=abs(m);

if(m>max(L))
  error('Order cannot exceed degree')
end

% Figure out if it's lowpass or bandpass
lp=length(L)==1;
bp=length(L)==2;
maxL=max(L);

% The spherical harmonic dimension
ldim=(L(2-lp)+1)^2-bp*L(1)^2;

% Filename of saved data
dirname=fullfile(getenv('IFILES'),'SDWCAP');
if lp
  fnpl=fullfile(dirname,sprintf(...
      'SDW-%i-%i-%i-%i.mat',TH,L,nth,m));
elseif bp
  fnpl=fullfile(dirname,sprintf(...
      'SDWBL-%i-%i-%i-%i-%i.mat',TH,L(1),L(2),nth,m));
else
  error('The degree range should be either one or two numbers')
end

% In 7.0.0.19901 (R14), in very rare cases (*-180-0-0) this
% save was buggy, and could not be subsequently loaded.
% Thus used  ~(L==180 & m==0) as an extra condition, now gone.
% Used matzerofix to remediate this.
if exist(fnpl,'file')==2 && (vcut>0) & 1==3
  load(fnpl)
  disp(sprintf('%s loaded by SDWCAP',fnpl))
else
  % Get the full Shannon number
  N=ldim*(1-cos(TH/180*pi))/2;
  
  % Approximate Nyquist degree is $\pi/d\theta$
  if nth < (maxL+1) & nth~=0
    error('Sample finer to avoid aliasing')
  end
  
  % Convert to radians
  TH=TH*pi/180;

  % Initialize kernel
  K=zeros(maxL+1-max(m,bp*min(L)),maxL+1-max(m,bp*min(L)));

  % Construct kernel: integral of Ylm over the cap
  lmin=max(m,bp*min(L));

  switch method
   case 'gl'
    tic
    for lr=lmin:maxL
      for lc=lr:maxL
        % Orthonormalization of Ylm is to unity over unit sphere
        % When TH=180, K should be the identity matrix
        % The pi*(1+(m==0)) is the longitudinal integral
        % Note that this line ALSO would work with 'paul' but then it
        % would be very inefficient in not reusing the database
        K(lr+1-lmin,lc+1-lmin)=...
            legendreprodint(lr,m,lc,m,cos(TH),method)...
            *sqrt(2*lr+1)*sqrt(2*lc+1)/(4*pi)*pi*(1+(m==0));
        % Symmetrize
        K(lc+1-lmin,lr+1-lmin)=...
            K(lr+1-lmin,lc+1-lmin);
      end
    end
    toc
   case 'paul'
    tic
    % Load all the wigner0j symbols up to level 2*maxL
    [~,C0,S0,L0]=zeroj(0,0,2*maxL);
    % Load all the wigner3j symbols up to level 2*maxL
    % If S3 is an empty you just got the symmetrical variety back
    [~,C3,S3,L3]=threej(0,0,2*maxL);
    
    % Get all the integrals of single functions
    Itab=paul(2*maxL,cos(TH));

    % We won't load the wigner3j data base since we've symmetrized the
    % storage 
    for lr=lmin:maxL
      for lc=lr:maxL
        % The required degree range
        ELL=(max(abs(lr-lc),2*m):(lr+lc));
        % Note that the zero-bottom symbol needs to have the top row sum
        % even. This is never empty, I think, but if it is, see LEGENDREPRODINT
        ELL=ELL(~mod(ELL+lr+lc,2));
        lELL=length(ELL); ELL=ELL(:);
        
        % The actual wigner0j symbols needed
        w0=zeroj(lr,lc,ELL,L0,2,C0,S0); w0=w0(:);
        % The actual wigner3j symbols needed
        if isempty(S3)
          % You have the symmetrized version, proceed as in THREEJ
          % Find where they're being kept
          [CC,oddperm,phasefix]=wignersort(ELL,lr,lc,-2*m,m,m);         
          % Do the initial evaluation from the loaded variable
          wm=full(C3(CC));
          % Fix the phase
          wm(oddperm)=wm(oddperm).*phasefix;
          % Now fix the triangle condition violations
          wm=wm.*triangle(repmat(lr,lELL,1),repmat(lc,lELL,1),ELL);
          % Now fix the order violations
          wm=wm.*~[lr<abs(m) | lc<abs(m) | ELL(:)<abs(-2*m)];
        else
          % You have the full sparse database ready to pass on to THREEJ
          wm=threej(ELL,lr,lc,-2*m,m,m,L3,[],C3,S3); wm=wm(:);
        end

        % And this here is as in LEGENDREPRODINT
        Q=(-1)^(m+m)*(2*ELL+1).*wm.*w0*sqrt(2-(m==0))*sqrt(2-(m==0))...
          /sqrt(2-[(m+m)==0]);
        % Figure out the right indices according to the way the Itab is set up
        indices=ELL.*(ELL+1)/2+m+m+1;       
	% CHECK WHAT ORDERS COME OUT!
        
        % The Paul-Gaunt integral with extra adaptation to the YLM
        K(lr+1-lmin,lc+1-lmin)=Q(:)'*Itab(indices,:)...
            *sqrt(2*lr+1)*sqrt(2*lc+1)/(4*pi)*pi*(1+(m==0));
        % Symmetrize
        K(lc+1-lmin,lr+1-lmin)=...
            K(lr+1-lmin,lc+1-lmin);
      end
    end
  end
  toc
  
  % Calculate eigenvalues and eigenvectors; C'*C=I
  [C,V]=eig(K);
  
  if lp && vcut<=0
    % Check the partial Shannon number
    difer(sum(diag(V))-indeks(nsubm(N,m,3,L),'end'))
  elseif bp
    disp('Still need to fix NSUBM for bandpass functions')
  end

  % Fill in the missing degrees in case it was bandpass
  C=[zeros(lmin-m,size(C,2)) ; C];
  
  % Check normalization and get number of GL points
  [ngl1,ngl2,com,Vc,nofa,zmean]=orthocheck(C,V,TH,m,1,0,bp);
  
  % If it's bandpass, check the mean
  if bp ; difer(zmean,[],[],'Zero-mean eigenfunctions'); end
  
  if length(C)>1 & m==0 & lp
    % Calculate uncertainty relation
    sdl=sqrt(sum(repmat([0:L]'.*[1:L+1]',1,size(C,2)).*C.^2,1));
    NN=1;
    sdv=sqrt(NN*(1+(1-2/NN)*com.^2));
    unc=sdv./com.*sdl;
  else
    com=0;
    unc=0;
    sdl=0;
    sdv=0;
  end

  % Order eigenvalues and eigenfunctions downgoing
  [V,isrt]=sort(sum(V,1),'descend');
  % V=fliplr(V); C=C(:,fliplr(isrt));
  C=C(:,isrt);
  unc=fliplr(unc); sdl=fliplr(sdl); com=fliplr(com); sdv=fliplr(sdv); 

  % Only return nonzero "useful" eigenvalues
  C=C(:,V>vcut); 
  V=V(V>vcut);

  % Compute spatial functions, colatitudinal part only
  if nth~=0
    % Zonal functions only 
    if m==0
      % Make spatial functions
      % This is SDW (2005) equation (5.10) combined with the sqrt(2-dom) of
      % (5.12) already included!
      [E,th]=pl2th(C,nth,1);
      th=th*180/pi;
      nlon=2*nth-1;
    else
      % This is SDW (2005) equation (5.10) combined with the sqrt(2-dom) of
      % (5.12) already included!
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
  % In 7.0.0.19901 (R14), in very rare cases (3-180-0-0) this
  % save is buggy, and cannot be subsequently loaded
  % Output in degrees
  save(fnpl,...
       'E','V','L','N','TH','C','th','nlon','ngl1','ngl2','unc','com','sdl','K')
end

if nth~=0 & grd==2
  % Output on full grid; watch the sign of m
  % Note that sqrt(2-dom) is already part of the harmonic
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
