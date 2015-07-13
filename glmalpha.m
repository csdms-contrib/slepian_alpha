function varargout=glmalpha(TH,L,sord,blox,upco,resc,J,anti,rotb)
% [G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=GLMALPHA(TH,L,sord,blox,upco,resc,J,anti,rotb)
%
% Returns an (lm)X(alpha) matrix with unit-normalized spherical harmonic
% coefficients of the BANDLIMITED or PASSBAND Slepian functions of the
% SINGLE or DOUBLE polar cap, or of a geographical region of
% interest. Only in the geographical case are the eigenvalues automatically
% sorted; if not, the column dimension is always block-ordered by virtue of the
% construction. The matrix G is orthogonal, G'*G is the identity. In column
% truncated form, G(:,1:J)*G(:,1:J)' is not the identity also, but rather
% a projection with eigenvalues 0 and 1. 
%
% Should put an option to save only the essentials up to a certain truncation
%
% INPUT:
%
% TH       Angular extent of the spherical cap, in degrees OR
%          'england', 'eurasia',  'namerica', 'australia', 'greenland', 
%          'africa', 'samerica', 'amazon', 'orinoco', 'antarctica', 
%          'contshelves', 'alloceans',
%          OR: [lon lat] an ordered list defining a closed curve [degrees]
%          OR: {'region' buf} where buf is the distance in degrees that 
%          the region outline will be enlarged by BUFFERM
% L        Bandwidth (maximum angular degree), or passband (two degrees)
% sord     1 Single polar cap of radius TH [default]
%          2 Double polar cap, each of radius TH
%          N Splining smoothness for geographical regions [default: 10]
%
% The following options are only for axisymmetric polar caps:
%
% blox     0 Standard (lm) row ordering, l=0:L, m=-l:l as ADDMOUT [default]
%          1 Block-diagonal row ordering, m=[0 -1 1 -2 2 ... -L L], l=abs(m):L
% upco     +ve fraction of unit radius for upward continuation [default: 0]
%          -ve fraction of unit radius for downward continuation
% resc     0 Not rescaled [default]
%          1 Rescaled to have a unit integral over the unit sphere
%            (this is only relevant for the down/upward continued functions)
%
% Though there are three more options that hold only for the geographic cases
%
% J        The number of eigenfunctions that are being asked (and saved);
% anti     1 get the opposite of the region you specify 
%          0 get exactly the region that you specify [default]
% rotb     0 nothing special happens
%          1 eigenfunctions of rotated kernels are rotated back to the
%          correct positions, as recommended for 'antarctica',
%          'contshelves'. This does not need to be recorded by the
%          filename as this is what any reasonable person should want
%
% OUTPUT:
%
% G        The unitary matrix of localization coefficients; note how
%          LOCALIZATION delivers these as LMCOSI arrays into PLM2XYZ
% V        The eigenvalues in this ordering (not automatically sorted)
% EL       The vector with spherical harmonic degrees as first index of G
% EM       The vector with spherical harmonic orders as first index of G
% N        The Shannon number
% GM2AL    The sum over all orders of the squared coefficients, i.e. the
%          TOTAL power, NOT the power spectral density
% MTAP     The order of the eigentapers, if the region is axisymmetric
% IMTAP    The rank within that particular order of the eigentapers
%
% EXAMPLE:
%
% glmalpha('demo1') % Illustrates the block sorting and checks unitarity
% glmalpha('demo2') % Makes a coupling kernel a la BCOUPLING
%
% SEE ALSO:
%
% GLMALPHAPTO, ADDMOUT, ADDMON, KERNELC, LOCALIZATION, GALPHA, DLMLMP, GLM2LMCOSI
%
% Last modified charig-at-princeton.edu, 06/16/2015
% Last modified by fjsimons-at-alum.mit.edu, 06/05/2013

% Should be able to update this to retain the rank order per m as well as
% the global ordering. Does this work for the whole-sphere? In that case,
% should really want G to be the identity - all though of course,
% anything works, too. You don't get necessarily the spherical harmonics
% back...

defval('TH',30)

if ~(ischar(TH) && ~isempty(strfind(TH(:)','demo')))

  defval('L',18)
  defval('dom',[]);
  % This is only relevant for the axisymmetric cap
  defval('blox',0);
  defval('upco',0);
  defval('resc',0);
  defval('anti',0);
  defval('rotb',1);

  defval('mesg','GLMALPHA Check passed')
  % Hold all messages
  mesg=NaN;

  % Figure out if it's lowpass or bandpass
  lp=length(L)==1;
  bp=length(L)==2;
  maxL=max(L);

  % The spherical harmonic dimension
  ldim=(L(2-lp)+1)^2-bp*L(1)^2;

  % Just get the file names here
  if upco==0 && resc==0
    if ~isstr(TH) && ~iscell(TH) && length(TH)==1 % POLAR CAPS
      defval('sord',1) % SINGLE OR DOUBLE CAP
      if lp
        fname=fullfile(getenv('IFILES'),'GLMALPHA',...
	      sprintf('glmalpha-%i-%i-%i-%i.mat',TH,L,sord,blox));
      elseif bp
        fname=fullfile(getenv('IFILES'),'GLMALPHA',...
          sprintf('glmalphabl-%i-%i-%i-%i-%i.mat',...
          TH,L(1),L(2),sord,blox));
      else
        error('The degree range is either one or two numbers')       
      end

      % Initialize ordering matrices
      MTAP=repmat(0,1,ldim);
      IMTAP=repmat(0,1,ldim);
    else % GEOGRAPHICAL REGIONS and XY REGIONS
      defval('sord',10) % SPLINING SMOOTHNESS
      defval('buf',0) % BUFFER REGION
      % We'll put in a Shannon number based on the area only, not based on
      % an actual sum of the eigenvalues
      defval('J',ldim)
      % Not the next line, though we can change our minds
      % defval('J',ldim*spharea(TH)) % beware, this currently breaks for buffers
      if isstr(TH) % Geographic (keep the string)
        h=TH; dom=TH;
      elseif iscell(TH) % Geographic + buffer
        if TH{2}==0; h=TH{1}; else h=[TH{1} num2str(TH{2})]; end
        %h=[TH{1} num2str(TH{2})];
        dom=TH{1}; buf=TH{2};
      else % Coordinates (make a hash)
        h=hash(TH,'sha1');
      end
      if lp
        fname=fullfile(getenv('IFILES'),'GLMALPHA',...
		     sprintf('glmalpha-%s-%i-%i.mat',h,L,J));
      elseif bp
        fname=fullfile(getenv('IFILES'),'GLMALPHA',...
		     sprintf('glmalphabl-%s-%i-%i-%i.mat',h,L(1),L(2),J));
      else
        error('The degree range is either one or two numbers')       
      end
      defval('GM2AL',NaN) % If not, calculate order per taper
      defval('MTAP',NaN) % If not, calculate order per taper
      defval('IMTAP',NaN) % And rank ordering within that taper
      defval('xver',0) % For excessive verification of the geographical case
    end
  else
    fname='neveravailable';
    defval('xver',1) % For excessive verification of the upco'd case
  end

  if anti==1
    % Update the file name to reflect the complimentarity of the region
    fname=sprintf('%s-1.mat',pref(fname)); 
  end

  if exist(fname,'file')==2
    load(fname)
    disp(sprintf('Loading %s',fname))
  else
    % Initialize matrices
    G=repmat(0,(maxL+1)^2,ldim);
    V=repmat(0,1,ldim);
    
    % Find row indices into G belonging to the orders
    [EM,EL,mz,blkm]=addmout(maxL);
    
    % Find increasing column index; that's how many belong to this order
    % alpha=cumsum([1 L+1 gamini(L:-1:1,2)]);
    % The middle bit is twice for every nonzero order missing
    % alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
    %   		gamini(L(2-lp)-bp*(L(1)-1),bp*2*(L(1)-1)) ...
    %   		gamini(L(2-lp)-bp*(L(1)-1):-1:1,2)]);
    % This should be the same for L and [0 L]
    alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
      gamini(L(2-lp)-bp*(L(1)-1),bp*2*L(1)) ...
      gamini(L(2-lp)-bp*L(1):-1:1,2)]);
    
    % For GEOGRAPHICAL REGIONS or XY REGIONS
    if isstr(TH) || length(TH)>1
      % Calculates the localization kernel for this domain
      % See if we can run this calculation in parallel
      tl = license('test','distrib_computing_toolbox'); % license?
      if tl
        if verLessThan('matlab', '8.2')
            % For MATLAB older than MATLAB 8.2, we need to check if the pool is open
            s = matlabpool('size');
            if s
              disp('Running KERNELCP (parallel)');
              [Klmlmp,XY]=kernelcp(maxL,TH,sord);
            else
              disp('No open matlabpool.  Running KERNELC (non-parallel).');
              [Klmlmp,XY]=kernelc(maxL,TH,sord);
            end    
        else
            % For MATLAB 8.2 and newer, a parpool should start automatically
            disp('Running KERNELCP (parallel)');
            [Klmlmp,XY]=kernelcp(maxL,TH,sord);
        end
      else
        disp('No Parallel Computing License.  Running KERNELC (non-parallel).');
        [Klmlmp,XY]=kernelc(maxL,TH,sord);  
      end
      
      if anti==1
        % Get the complimentary region
        Klmlmp=eye(size(Klmlmp))-Klmlmp;
      end

      if bp
        % So far we only have wanted to remove small portions of the
        % kernel.  Therefore at the moment here, we make the whole thing,
        % and the apply the bp after the fact.  However, in the future if
        % you want a bp at large L, we should modify kernelcp to only make
        % a partial kernel.
        
        % Remove the beginning section of the kernel 
        rem=bp*L(1)^2;
        Klmlmp(:,1:rem)=[];
        Klmlmp(1:rem,:)=[];
      end
      % Calculates the eigenfunctions/values for this localization problem
      if lp 
        [G,V]=eig(Klmlmp);
      elseif bp
        [Gbp,V]=eig(Klmlmp);
        G(L(1)^2+1:end,:) = Gbp;
      end
      [V,isrt]=sort(sum(real(V),1));
      V=fliplr(V);
      G=G(:,fliplr(isrt));
    
      [~,~,~,~,~,~,ems,els,R1,ronm]=addmon(maxL);
      % This indexes the orders of G back as 0 -101 -2-1012 etc
      G=G(R1,:);
      % Check indexing
      difer(els(R1)-EL,[],[],mesg)
      difer(ems(R1)-EM,[],[],mesg)
    
      % Calculate Shannon number and compare this with the theory
      N=sum(V);

      if lp
        % Is the Shannon number right? Need the area of the region
        difer(ldim*Klmlmp(1)-N,[],[],mesg)
      elseif bp
        difer(ldim*spharea(TH)-N,[],[],mesg)
      end

      % Check if the expansion of a basis function is indeed either 1 or 0
      if xver==1
        disp('Excessive verification')
        % Is the area right?
        difer(Klmlmp(1)-spharea(TH),4,[],mesg)

        % This is a bit double up... but it's only for excessive verification
        [V1,C]=localization(L,TH,sord);
        difer(V-V1',[],[],mesg)
        for index=1:length(C)
	  salpha=G'*C{index}(ronm);
	  % Only one of these functions should get "hit"
	  difer(sum(abs(salpha)>1e-9)-1,[],[],mesg)
        end
      end
    
      % If it's various parts of Antarctica or ContShelves, need to rotate 
      % it back in shape
      % Note: currently the "serial" version of ROTATEGP is itself since it
      % defaults to a regular for loop
      if [strcmp(dom,'antarctica') || strcmp(dom,'eantarctica')...
            || strcmp(dom,'eantarcticaCoasts1') || strcmp(dom,'eantarcticaCoasts2')...
            || strcmp(dom,'eantarcticaInt') || strcmp(dom,'eantarcticaIntG')...
            ] && rotb==1
        defval('pars',10)
        % Return the rotation parameters also, to rotate G
        % Each of these coordinate files return lonc and latc which are not
        % necessarilly the center of the region, but are appropriate centers
        % for rotation (e.g. [0 -90])
        [~,lonc,latc]=eval(sprintf('%s(%i,%f)',dom,pars,buf));
        [Grot] = rotateGp(G,lonc,latc);
        G = Grot;
      elseif [strcmp(dom,'antarcticaGP') || strcmp(dom,'eantarcticaCoasts1OceanBuf')...
            || strcmp(dom,'eantarcticaCoasts2OceanBuf')...
            || strcmp(dom,'eantarcticaIntGOceanBuf')] && rotb==1
        defval('pars',0)
        [~,lonc,latc]=eval(sprintf('%s(%i,%f)',dom,pars,buf));
        [Grot] = rotateGp(G,lonc,latc);
        G = Grot;
      elseif strcmp(dom,'contshelves') && rotb==1
        defval('pars',10)
        [~,lonc,latc]=eval(sprintf('%s()',dom));
        [Grot] = rotateGp(G,lonc,latc);
        G = Grot;
      end
      % You can plot this here, if you want, by doing, e.g.
      % cosi = lmcosi(:,3:4);
      % cosi(ronm)=Grot(:,1);
      % plotplm([lmcosi(:,1:2) cosi],[],[],2,0.5); view(145,-35)
      % This now does show up in the right spot

      % Here's it's going to be almost trivial to get the integral over the
      % region of the Slepian eigenfunctions, given that it is a linear
      % combination of a scaled row of Klmlmp
      % sarea=G*Klmlmp(:,1);
      % Which you could verify in the space domain using galpha

      % Here I should save the actual eigenfunctions
      % defval('J',round(N))
      % Truncate to the smaller amount of eigenfunctions and -values
      G=G(:,1:J);
      V=V(1:J);
      save(fname,'-v7.3','G','V','EL','EM','N')
    else
      % For AXISYMMETRIC REGIONS
      if blox~=0 && blox~=1
        error('Specify valid block-sorting option ''blox''')
      end
      % For the SINGLE or DOUBLE POLAR CAPS
      for m=0:maxL
        % Same results for +/- m; set nth=0 thus no sign correction!
        if sord==1
	  if lp
	    [E,Vg,th,C,T,Vp]=grunbaum(TH,L,m,0);
	  elseif bp
	    % Note that the small-eigenvalue eigenfunctions might be
        % numerically degenerate and thus not as using Grunbaum - if
        % you need to compare, compare where the eigenvalues are "big"
	    [E,Vp,Np,th,C]=sdwcap(TH,L,m,0,-1);
      end
        elseif sord==2
	  if lp
	  [E,Vg,th,C,T,Vp]=grunbaum2(TH,L,m,0);
	  elseif bp
	  error('Bandpass double-cap tapers not ready yet')
      end
        else
	  error('Specify single or double polar cap')
        end
      
      if upco~=0
	if upco>0
	  % The upward continuation matrix
	  A=diag((1+upco).^[-(m:L)-1]);
	elseif upco<0
	  % The downward continuation matrix
	  A=diag((1+abs(upco)).^[(m:L)+1]);
	end
	
	% Comparisons with Grunbaum only make sense for lowpass
	if xver==1 & lp
	  % This should give the same result, more or less, less accurate 
	  if sord==1
	    [a,Vs,c,d,Cs,e,f,g,h,j,D]=sdwcap(TH,L,m,0,-1);
	  else
	    [a,Vs,c,Cs,e,f,D]=sdwcap2(TH,L,m,0,-1);
	  end
	  % This should give the eigenvalues again, which we'd had from
      % orthocheck 
	  warning off
	  % Check difference integration and kernel eigenvalues
	  difer(Vp(:)-diag((C'*D*C)./(C'*C)),[],[],mesg)
	  % Check difference integration and diagonalization eigenvalues
	  difer(Vp(:)-Vs(:),[],[],mesg)
	  % Check difference between the eigenfunctions barring the sign
	  % and only wherever the eigenvalues amount to anything
	  difer(abs(Cs(:,Vp>1e-6))-abs(C(:,Vp>1e-6)),[],[],mesg)
	  warning on
	  Vc=diag((C'*A*C*diag(Vp)*C'*A*C));
	  Vp0=Vp;
	end
	
	% Upward continuation from 1 to 1+a or from 1+a to 1:
	% New eigenfunctions, same name
	C=A*C;
	% Calculate new eigenvalues, same name
	[jk1,jk2,jk3,Vp,nofa]=orthocheck(C,[],TH/180*pi,m,sord,1);

	% Make sure these are sorted, since that's not automatically the case
	% [Vp,ind]=sort(Vp,'descend');
	% C=C(:,ind);
	% Current thinking is: do NOT resort, as you'll want to compare the
	% best at a=0 with whatever it becomes later!
	
	if xver==1
	  warning off
	  % Check difference integration eigenvalues and those from kernel
	  difer(Vp(:)-diag((C'*D*C)./(C'*C)),[],[],mesg)
	  warning on
	  % Check how many Vp>Vp0 
	  disp(sprintf('%i/%i eigenvalues greater',sum(Vp(:)>Vp0(:)), ...
		       length(Vp0)))
	  disp(sprintf('Shannon number greater: %i',sum(Vp)>sum(Vp0)))
	end
	if resc==1
	  % Rescale these functions to have an integral to unity over the
	  % sphere; note: this doesn't make the set orthonormal of course
	  C=C*diag(1./nofa);
	end
      end

        % Distribute this at the right point in the huge matrix
        if m>0
          % Here you supply the negative orders
          G(EM==-m,alpha(2*m):alpha(2*m+1)-1)=C;
          V(alpha(2*m):alpha(2*m+1)-1)=Vp;
          MTAP(alpha(2*m):alpha(2*m+1)-1)=-m;
          % It's all neatly ordered here, downgoing within every order
          IMTAP(alpha(2*m):alpha(2*m+1)-1)=1:length(Vp);
        end
        % Duplicate for the positive order in case the region is axisymmetric
        G(EM==m,alpha(2*m+1):alpha(2*m+2)-1)=C;
        V(alpha(2*m+1):alpha(2*m+2)-1)=Vp;
        MTAP(alpha(2*m+1):alpha(2*m+2)-1)=m;
        % It's all neatly ordered here, downgoing within every order
        IMTAP(alpha(2*m+1):alpha(2*m+2)-1)=1:length(Vp);
      end

      % Calculate the Shannon number and compare it to the theory
      N=sum(V);
      if upco==0
        difer(N-ldim*sord*(1-cos(TH/180*pi))/2,[],[],mesg);
      end
    
      % Compute the sum over all orders of the squared coefficients
      % Thus works when they have not been blocksorted yet. 
      GM2AL=repmat(0,ldim,maxL+1);
      for l=0:maxL
        b=(l-1+1)^2+1;
        e=(l+1)^2;
        GM2AL(:,l+1)=sum(G(b:e,:).^2,1)';
      end

      % Make sure that the sum over all degrees is 1 - but I forgot why
      difer(sum(GM2AL,2)-1,[],[],mesg)

      % This is not blockdiagonal, unless you make it thus
      if blox==1
        G=G(blkm,:);
        EM=EM(blkm);
        EL=EL(blkm);
      end
      if ~strcmp(fname,'neveravailable') 
        % Save the results if it isn't a geographical region
        % If the variable is HUGE you must use the -v7.3 flag, if not, you
        % can safely omit it and get more backwards compatibility
        save(fname,'-v7.3','G','V','EL','EM','N','GM2AL','MTAP','IMTAP')
      end
    end
  end

  % Provide output
  varns={G,V,EL,EM,N,GM2AL,MTAP,IMTAP};
  varargout=varns(1:nargout);

elseif strcmp(TH,'demo1')
  
  % Note that using ADDMOUT you can get this back to block-diagonal form
  G=glmalpha; [~,~,~,bl]=addmout(18); imagefnan(G(bl,:)); axis tight
  difer(G'*G-eye(size(G)))

elseif strcmp(TH,'demo2')
  % Make a demo for Saarimaki et al.
  L=36;
  [G0,V0,EL0,EM0]=glmalpha(30,L,1,0);
  [G1,V1,EL1,EM1]=glmalpha(30,L,1,1);
  [i1,j1]=sort(V1,'descend');
  [i0,j0]=sort(V0,'descend');
  N0=round(sum(V0));
  N1=round(sum(V1));
  imagefnan(G0(:,1:N0)*G0(:,1:N0)')
  imagefnan(G1(:,1:N1)*G1(:,1:N1)')

  % Not block sorted
  [G0,V0,EL0,EM0]=glmalpha('africa',L,[],0);
  % Block sorted
  [G1,V1,EL1,EM1]=glmalpha('africa',L,[],1);

  % Fake a coupling kernel a la BCOUPLING
  M=G0(:,1:N0)*G0(:,1:N0)';
  M=M.^2;
  for l=0:L
    b=l^2+1;
    e=(l+1)^2;
    % Construct the symmetric part
    for lp=l:L
      bp=lp^2+1;
      ep=(lp+1)^2;
      % The symmetric expression
      ML(l+1,lp+1)=sum(sum(M(b:e,bp:ep)));
      ML(lp+1,l+1)=ML(l+1,lp+1);
    end
  end
  % Forget about the asymmetric scaling and look at a certain cut
  for i=1:L+1
    plot(ML(i,:)/ML(i,i))
    j=indeks(find(diff(ML(i,:)/ML(i,i)>0.05)),1:2);
    set(gca,'xtick',sort([i j]),'xgrid','on',...
	    'ytick',[0 0.05 1],'ygrid','on')
    pause
  end
end
