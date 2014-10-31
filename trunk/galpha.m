function varargout=...
    galpha(TH,L,sord,theta,phi,srt,upco,resc,blox,J,irr,Glma,V,N,EL,EM)
% [G,V,EM,GK,VK,NA,N,theta,phi,Glma,EL]=...
%           GALPHA(TH,L,sord,theta,phi,srt,upco,resc,blox,J,irr)
%
% Constructs a matrix of axisymmetric-polar-cap Slepian eigenfunctions
% evaluated at a set of spatial locations. Normalization is to UNITY over
% the surface of the entire sphere.
%
% INPUT:
%
% TH          Angular extent of the spherical cap, in degrees
% L           Bandwidth (maximum angular degree) or passband (two degrees)
% sord        1 Single polar cap of radius TH [default]
%             2 Double polar cap, each of radius TH
%             3 Equatorial belt of width 2TH
% theta       colatitude (0 <= theta <= pi) [default: 720 linearly spaced]
% phi         longitude (0 <= theta <= 2*pi)  [default: 0], if NaN, get Xlm
% srt         'global' Global sorting of eigenvalues across all orders
%             'local'  Sorting of cap eigenvalues within orders [default]
%             'belt'   Global sorting of eigenvalues on the belt
% upco        +ve fraction of unit radius for upward continuation [default: 0]
%             -ve fraction of unit radius for downward continuation
% resc        0 Not rescaled [default, unit-normalization is in effect]
%             1 Rescaled to retain unit integral over the unit sphere
%               (this is only relevant for the down/upward continued functions)
% blox        1 or 0 but should not affect any output (see 'demo1')
% J           How many of the best eigenfunctions do you want? [default: all]
% irr         0 Regular grid, no matter how you input theta and phi [default]
%             1 Irregular grid, input interpreted as distinct pairs of theta, phi
% Glma        The spectral eigenfunctions in case you already have them
% V           The spectral eigenvalues in case you already have them
% N           The Shannon number in case you already have it
% EL          The degrees in question if you already have them
% EM          The orders in question if you already have them
%
% OUTPUT:
%
% G           A matrix with dimensions
%               (L+1)^2 x (length(theta)*length(phi))
% V           A vector with the eigenvalues 
% EM          If axisymmetric, the order on which the taper is based
%             (otherwise, this makes no sense and all orders are involved)
% NA          N/A, which must be equal to diag(G'*G)
% N           N, the Shannon number, rounded to the nearest integer
% GK          Matrix G truncated to Shannon-number proportions
% VK          Vector V truncated to Shannon-number proportions
% theta       The colatitudes that you put in or received
% phi         The longitudes that you put in or received
% Glma        The spectral eigenfunctions, for use in, e.g. PLOTSLEP
% EL          The degrees in question
%
% EXAMPLE:
%
% galpha('demo1') % should return no error messages
% galpha('demo2') % makes some quick plots at the Greenwich meridian
% galpha('demo3') % makes some quick 2D plots, checks normalization
% galpha('demo4') % makes some quick 2D plots, with a rotation
%
% G=galpha(40,18,1,linspace(0,pi,180),linspace(0,2*pi,360),'global');
% 
% See also DOUBLECAP
%
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012

% If the output was on a Driscoll-Healey or HEALPIX grid (and maybe I
% should think about doing that), this would be an orthogonal matrix, but
% now it isn't, of course.... Should adapt for irregular grids as GALPHAPTO.

% Note to self:
%
% To find how many caps at different orders m are required for a total
% number of well concentrated function N, you could either run the
% fixed-order problem, and take the first round(Nm) of each of those,
% knowing that you have to take m as well as -m. But, although
% N0+2*sum(N1...NL) equals N exactly, it is not true that round(N) equals
% round(N0+2*sum(N1...NL))... So you might be making a small mistake if
% you go via route 1 above. Instead, you could run this code 
% [G,V,EM,GK,VK,NA,N]=galpha(TH,L,1,[],NaN,[],'global');
% and figure out which orders to take by looking at EM(1:N)... the only
% trouble there, then, is, that you might be splitting a multiplet +/- m,
% and it will be up to you to decide whether or not to take another one
% or one fewer... Seems like the round(Nm) approach after all would be OK.

defval('TH',40)

if ~isstr(TH)
  defval('L',18)
  defval('sord',1)
  defval('theta',linspace(0,pi,720))
  defval('phi',0)
  defval('srt','local')
  defval('upco',0)
  defval('resc',0)
  % Note: the block sorting only matters in here and not for the output
  defval('blox',0)
  defval('irr',0)

  % Basic error checks
  if irr==1 & ~all(size(theta)==size(phi))
    error('Input arrays must have the same dimensions for irregular grids')
  end
  if isempty(strcmp(lower(srt),'global')) & ...
	isempty(strcmp(lower(srt),'local'))& ...
	isempty(strcmp(lower(srt),'belt'))
    error('Specify valid option ''srt''')
  end

  % Figure out if it's bandlimited or bandpass
  lp=length(L)==1;
  bp=length(L)==2;
  maxL=max(L);
  
  % The spherical-harmonic dimension
  ldim=(L(2-lp)+1)^2-bp*L(1)^2;
  
  % Adjust the calling sequence, then revert
  if sord==3
    sordo=2; THo=90-TH;
  else
    sordo=sord; THo=TH;
  end
  % Get the coefficients and the orders etc
  % The EM's are for the tapers! the EMrows disappear later
  if exist('Glma')~=1 ||  exist('V')~=1 ||  exist('N')~=1 ...
	|| exist('EL')~=1 ||  exist('EM')~=1
    [Glma,V,EL,EMrow,N,GAL,EM]=glmalpha(THo,L,sordo,blox,upco,resc);
  end
  
  % Calculate N/A and rounded Shannon number
  NA=ldim/(4*pi);
  
  % Slight quirk in the spharea function
  if sord==1 || sord==3
    N=round(ldim*spharea(TH,sord));
  elseif sord==2
    N=round(ldim*spharea(90-TH,sord));
  end

  % I suppose here later we can put the geographical regions, 
  % but, the lon/lat parameterization here in GALPHA is not ideal.

  % Get the unit-normalized real spherical harmonics
  if length(phi)==1 && isnan(phi)
    XYlmr=xlm([0 maxL],[],theta,0,0,blox);
  else
    XYlmr=ylm([0 maxL],[],theta,phi,0,0,blox,irr);
  end
  % I believe this should get the (-1)^m in there also because if not it
  % won't conform to how PLM2XYZ would do it should you use GLM2LMCOSI...
  % But this can wait - not a critical function. GALPHAPTO already does
  % it. 

  % Default truncation is none at all, so the 'srt' option has power
  defval('J',ldim)
  % Eigenvalue sorting? If 'local' doesn't do anything
  if strcmp(lower(srt),'belt') || sord==3
    V=1-V;
  end
  if strcmp(lower(srt),'global') || strcmp(lower(srt),'belt')
    [V,i]=sort(V,'descend');
    Glma=Glma(:,i); EM=EM(i);
  end
  % Truncation? If past Shannon number, do it from here on, if not, wait
  % But doing the spectral truncation before the spatial expansion saves time
  if J~=ldim && J>=N
    Glma=Glma(:,1:J);
  elseif J~=ldim && J<N
    Glma=Glma(:,1:N);
  end
  
  % Expand - it's just here that the block ordering matters
  % i.e. for the EMrow which we won't need anymore. The EM are always block
  % sorted to begin with, since the GLMALPHA matrix is filled order by order.
  G=Glma'*XYlmr;

  [GK,VK]=deal([]);
  if nargout>=4
    if strcmp(lower(srt),'global') || strcmp(lower(srt),'belt')
      GK=G(1:N,:);
      VK=V(1:N);
    end
  end
  
  % And now do the final truncation that you wanted in the first place
  G=G(1:J,:);
  V=V(1:J);
  EM=EM(1:J);

  % What do we predict the sum of the eigenfunctions squared is?
  % See RB VI p89
  % And check that it's all nicely done
  if resc==0
    if upco==0 && J==ldim
      % This doesn't work if it's no longer full-dimensionsal
      difer(abs(diag(G'*G)-NA),[],[],NaN);
    elseif upco>0
      a=upco;
      p=(-((2*a+a^2+4*(L+1)*a+2*(L+1)*a^2+2)/...
	   ((1+a)^2)^(L+1))+(2*a+a^2+2))/(a^2*(2+a)^2)/4/pi;
      if J==ldim
	difer(abs(diag(G'*G)-p),10,0);
      end
    elseif upco<0
      a=abs(upco);
      p=(((-2*a-a^2+4*(L+1)*a+2*(L+1)*a^2-2)*(1+a)^(2*L+4))-...
	 (-2*a-a^2-2)*(1+a)^2)/(a^2*(2+a)^2)/4/pi;
      % Watch out, these numbers grow very large and accumulate errors!
      % But the formula is, however, correct
      if J==ldim
	difer(abs(diag(G'*G)-p),10,0);
      end
    end
  end
  
  % Output
  varns={G,V,EM,GK,VK,NA,N,theta,phi,Glma,EL};
  varargout=varns(1:nargout);
elseif strcmp(TH,'demo1')
  disp('Checking that internal block sorting has no effect whatsoever')
  for srt={'global' 'local' 'belt'}
    disp(srt)
    [G0,V0,EM0,GK0,VK0,NA0,N0,theta0]=galpha([],[],1,[],[],srt,[],[],0);
    [G1,V1,EM1,GK1,VK1,NA1,N1,theta1]=galpha([],[],1,[],[],srt,[],[],1);
    difer([G0(:);V0(:);EM0(:);GK0(:);VK0(:);NA0(:);N0(:);theta0(:)]-...
	  [G1(:);V1(:);EM1(:);GK1(:);VK1(:);NA1(:);N1(:);theta1(:)])
    TH=44; L=23; 
    [G,V,EM,GK,VK,NA,N,th]=galpha(TH,L,1,[],[],srt);
    % This should be the Shannon number everywhere
    difer(diag(G'*G)-NA)
    % This should be close to the Shannon number in the region
    plot(th*180/pi,diag(G'*diag(V)*G))
    hold on
    plot([TH TH],[0 NA],'k--')
    axis tight
    disp('Hit return to go on...')
    pause
  end
  hold off
elseif strcmp(TH,'demo2')
  sord=3;
  TH=40; L=18; [G,V,EM,GK,VK,NA,N,th]=galpha(TH,L,sord,[],[],'global');
  % To give a nice visual
  clf
  for index=1:(L+1)^2; 
    plot(th,G(index,:),'linew',2); 
    axis tight ; grid on ;
    title(sprintf('%s= %i ; order m= %i ; %s = %8.6f',...
		  '\alpha',index,EM(index),'\lambda',V(index))); 
    ylim([-3.3 3.3]) ; disp('Hit return to go on...') ; pause
  end
elseif strcmp(TH,'demo3')
  sord=3;
  TH=40; L=12; [G,V,EM,GK,VK,NA,N,~,~,Glma]=...
     galpha(TH,L,sord,linspace(0,pi,50),linspace(0,2*pi,100),'global');
  % Approximate normalization
  [lon,lat]=Fibonacci_grid;
  Gf=galpha(TH,L,sord,pi/2-lat*pi/180,lon*pi/180,'global',[],[],[],[],1);
  
  % To give a nice visual
  clf
  for index=1:(L+1)^2; 
    imagesc(reshape(G(index,:),50,100)); 
    % Check the normalization
    norma=sum(Gf(index,:).*Gf(index,:))*[4*pi/length(Gf(index,:))];
    disp(' '); disp(sprintf('Normalization over the sphere is %f',norma)); disp(' ')
    % Check again using GLM2LMCOSI which now contains the renormalization
    Gff=plm2xyz(glm2lmcosi(Glma,index),lat,lon);
    norma=sum(sum(Gff.*Gff))*[4*pi/length(Gff)];
    disp(' '); disp(sprintf('Normalization over the sphere is %f',norma)); disp(' ')

    axis tight ; grid on ;
    title(sprintf('%s= %i ; order m= %i ; %s = %8.6f',...
		  '\alpha',index,EM(index),'\lambda',V(index))); 
    disp('Hit return to go on...') ; pause
  end
elseif strcmp(TH,'demo4')
  sord=3;
  TH=40; L=12; [G,V,EM,GK,VK,NA,N,th]=...
     galpha(TH,L,sord,linspace(0,pi,50),linspace(0,2*pi,100),'global');
  % Apply rotation
  R=rots(L,V,EM,pi/6);
  G=R*G;
  clf
  % To give a nice visual
  for index=1:(L+1)^2; 
    imagef([],[],reshape(G(index,:),50,100)); 
    axis tight ; grid on ;
    title(sprintf('%s= %i ; order m= %i ; %s = %8.6f',...
		  '\alpha',index,EM(index),'\lambda',V(index))); 
    disp('Hit return to go on...') ; pause ; end
else
  error('Specify valid demo string, see HELP GALPHA')
end
