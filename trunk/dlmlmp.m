function D=dlmlmp(TH,L,Lmax,sord,xver)
% D=DLMLMP(TH,L,Lmax,sord,xver)
%
% Constructs the (overcomplete) localization kernel for axisymmetric
% concentration SINGLE-CAP domains. 
%
% INPUT:
%
% TH       Angular extent of the spherical cap, in degrees
% L        Bandwidth (maximum angular degree) of the BANDLIMITED functions
% Lmax     Requested bandwidth of the SPACELIMITED functions (note that
%          they won't be truly spacelimited unless Lmax is infinity
% sord     1 Single cap of diameter 2TH [default]
%          2 Double cap left by subtracting belt of width 2TH
% xver     1 Verifies GRUNBAUM per order and GLMALPHA for all orders
%          0 No excessive verification
% 
% OUTPUT:
%
% D        The (Lmax+1)^2 x (L+1)^2 axisymmetric localization kernel,
%          Simons, Dahlen, Wieczorek 2006, eq. (5.1). Automatically block
%          sorted per order as in GLMALPHA. Note that this matrix
%          contains zeros for all those orders exceeding L.
%
% EXAMPLE:
%
% D=dlmlmp(30,18,18,1,0); subplot(211); imagefnan(D); % Real harmonics
% U=flipud(fliplr(ummp(18))); subplot(212); imagefnan(U'*D*U); % Complex
%
% EXAMPLE:
%
% D=dlmlmp(40,18,18,1); [EM,EL,a,b,dblk]=addmout(18); 
% DD=D(dblk,dblk); imagesc(DD); mean(mean(DD-DD*DD))
%
% If the kernel is (L+1)^2 x (L+1)^2, its eigenfunctions are the
% optimally spatially concentrated bandlimited functions g; the
% coefficients out to Lmax of the optimally spectrally concentrated
% spatially limited eigenfunctions h can be obtained by multiplying the
% overcomplete (Lmax+1)^2 x (L+1)^2 kernel D with the (L+1)^2
% coefficients of g. Lmax will be twice the maximum bandwidth of the
% spectral estimate sought.
%
% See also: GLMALPHA, HLMALPHA, KERNELC
%
% Last modified by fjsimons-at-alum.mit.edu, 04/11/2007

defval('TH',30)
defval('L',30)
defval('Lmax',200)
defval('sord',1)
defval('xver',0)

if Lmax<L
  error('Lmax must be bigger than L for this to make sense')
end

% Load precomputed functions if they exist
fname=fullfile(getenv('IFILES'),'DLMLMP',...
	       sprintf('DLMLMP-%i-%i-%i-%i.mat',TH,L,Lmax,sord));

if exist(fname)==2
  load(fname)
else
  % Initialize huge matrix
  D=repmat(0,(Lmax+1)^2,(L+1)^2);
  
  % Construct vectors with the degree and the order
  [ML,LL,mz,b]=addmout(L); ML=ML(b); LL=LL(b);
  [MLmax,LLmax,mz,b]=addmout(Lmax); MLmax=MLmax(b); LLmax=LLmax(b);
  
  % Find increasing column index; that's how many belong to this order
  % See GLMALPHA for a bandpass update
  alphaL=cumsum([1 L+1 gamini(L:-1:1,2)]);
  alphaLmax=cumsum([1 Lmax+1 gamini(Lmax:-1:1,2)]);
  
  h=waitbar(0);
  % Construct kernel: integral of Ylm over the cap; only same orders matter
  for m=0:L
    waitbar(m/L,h,sprintf('DLMLMP: building order %i of %i',m,L))
    % Initialize kernel
    Dm=repmat(NaN,Lmax+1-m,L+1-m);
    for lc=m:L
      for lr=lc:Lmax
	% See also SDWCAP and SDWCAP2... granted, could do a better job
        % here and simply construct the sord=2 TH cap from the sord=1
        % 90-TH result - so only the sord=1 would have to be stored
	% Constructs one block matrix; lower "triangular" part first
	% But it'd would have to have a big array with lml'm' in here
	Dm(lr+1-m,lc+1-m)=[1+(sord==2)*(-1)^(lc+lr)]*...
	    legendreprodint(lr,m,lc,m,...
	    cos(pi/2*(sord==2)+(-1)^(sord+1)*TH*pi/180))...
	    *sqrt(2*lr+1)*sqrt(2*lc+1)/(4*pi)*pi*(1+(m==0));	
	% Symmetrize if needed
	if lr<=L
	  Dm(lc+1-m,lr+1-m)=Dm(lr+1-m,lc+1-m);
	end
      end
    end
    % Now stick this thing into the huge matrix at the right spot, just
    % like in GLMALPHA
    if m>0
      D(alphaLmax(2*m):alphaLmax(2*m+1)-1,alphaL(2*m):alphaL(2*m+1)-1)=Dm;
    end
    % Duplicate for the other order in case the region is axisymmetric
    D(alphaLmax(2*m+1):alphaLmax(2*m+2)-1,alphaL(2*m+1):alphaL(2*m+2)-1)=Dm;
    if xver==1 % Excessive verification, of commutation relation!
      [E,Vg,th,C,T,V]=grunbaum(TH,L,m,0);
      % Be slightly less critical in evaluating the commutativity
      difer(T*Dm(1:L+1-m,1:L+1-m)-Dm(1:L+1-m,1:L+1-m)*T,9)
%      imagefnan([],[],D,[],[],[],[],1)
%      pause
    end
  end

  close(h)
  
  % Now check whether or not this obeys the relations we think they should
  if xver==1
    disp('Using excessive verification')
    [G,V,EL,EM]=glmalpha(TH,L,1,1);
    if L==Lmax
      difer(D*G-G*diag(V),9);
    else
      % The first L elements of the new eigenfunctions must agree with the
      % rescaled versions of the old eigenfunctions; beyond that, we have
      % no control; we check this for the zeroth order only, so we don't
      % have to go hunt for the dumb indices... it works, anyway
      for index=1:L+1
	difer(indeks(D*G(:,index),1:L+1)-G(1:L+1,index)*V(index),9)
      end
    end
  end
  eval(sprintf('save %s D TH L Lmax',fname))
end
