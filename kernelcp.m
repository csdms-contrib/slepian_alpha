function varargout=kernelcp(Lmax,dom,pars,ngl,rotb)
% [Klmlmp,XY,K1,K]=KERNELCP(Lmax,dom,pars,ngl,rotb)
%
% Parallel version of KERNELC.
% Calculation of the localization matrix for some domain on the sphere.
% NOT FOR POLAR PATCHES! AND NOT GOOD FOR NEAR-POLAR PATCHES! (See GRUNBAUM)
% NOT WITHOUT MODIFICATIONS FOR REGIONS CONTAINING THE NORTH POLE OR THE
% SOUTH POLE! (For that, see GLMALPHA). Unit normalization as in YLM.
%
% INPUT:
%
% Lmax       Maximum angular degree (bandwidth)
% dom        'patch'   spherical patch [default], with specs in 'pars'
%                      NOTE: better use GRUNBAUM / PLM2ROT in this case
%            'sqpatch' square patch with [thN thS phW phE] in 'pars'
%            'england', 'eurasia',  'namerica', 'australia', 'greenland'
%            'africa', 'samerica', 'amazon', 'orinoco', 'gpsnamerica',
%            'antarctica', 'contshelves', 'alloceans', with specs in 'pars'
%            OR: [lon lat] an ordered list defining a closed curve [degrees]
%            OR: {'region' buf} where buf is the distance in degrees that 
%            the region outline will be enlarged by BUFFERM
% pars       [th0,ph0,thR] for 'patch'
%                 th0  Colatitude of the cap center, in radians
%                 ph0  Longitude of the cap center, in radians
%                 thR  Radius of the cap, in radians
%            N  splining smoothness for geographical regions [default: 10]
%            OR: the string with the name you want the result saved as
% ngl        The degree of the Gauss-Legendre integration [default: 200] OR
%            'alternative' for an alternative calculation method [not parallel]
% rotb       0 That's it, you got it [default: 0]
%            1 For, e.g., 'antarctica', 'contshelves', if you rotated coordinates
%            to make the integration procedure work, this option makes
%            sure that the kernel matrix reflects this. If not, you have
%            to apply counterrotation after diagonalizing in LOCALIZATION. 
%
% OUTPUT:
%
% Klmlmp     The localization kernel whose eigenfunctions we want,
%            indexed as: degree  [0  1  1  1  2  2  2  2  2]
%                        order   [0  0 -1  1  0 -1  1 -2  2]
%            The function LOCALIZATION later reindexes this in LMCOSI
%            fashion. Note: you can use ADDMOUT and ADDMON to modify, and
%            see, e.g. PLOTSLEP and KLMLMP2ROT for some implementations
% XY         The outlines of the region into which you are localizing
% K1         An intermediate result useful when rotb=1, see KLMLMP2ROT
% K          An verification result useful when rotb=1, see KLMLMP2ROT
% 
% EXAMPLE:
%
% L=19;
% [Klmlmp,XY]=kernelc(L,'australia'); 
% Klmlmp2=kernelc(L,'australia',[],'alternative'); 
% and then DIFER, EIG, PLOTSLEP, etc, to evaluate the difference
%
% kernelc('demo1') % For an illustration of the Antarctica matrix
% kernelc('demo2') % For an illustration of the Antarctica functions
% kernelc('demo3') % For a show of Australia
% kernelc('demo4') % For a demonstration of the rotation of the kernel
% kernelc('demo5') % For a demonstration of kernelcp
%
% See also LOCALIZATION, SDWREGIONS, GLMALPHA, DLMLMP, KERNELC2D, 
%          PLOTSLEP, PLM2AVG, KERNELC, LEGENDREPRODINT, DLMLMP
%
% Last modified by charig-at-princeton.edu, 05/14/2015
% Last modified by fjsimons-at-alum.mit.edu, 06/26/2012

t0=clock;
defval('Lmax',12); 
defval('dom','greenland')
defval('ngl',200)
defval('rotb',0)
defval('K1',NaN)

if ~isstr(Lmax)
  % Generic path name that I like
  filoc=fullfile(getenv('IFILES'),'KERNELCP');
  if isstr(dom)
    switch dom
      % If the domain is a square patch
     case 'sqpatch'
      fnpl=sprintf('%s/%s-%i-%i-%i-%i-%i.mat',filoc,dom,Lmax,...
		   round(pars(1)*180/pi),round(pars(2)*180/pi),...
		   round(pars(3)*180/pi),round(pars(4)*180/pi));
      % If the domain is a spherical patch
     case 'patch'
      fnpl=sprintf('%s/%s-%i-%i-%i-%i.mat',filoc,dom,Lmax,...
		   round(pars(1)*180/pi),round(pars(2)*180/pi),...
		   round(pars(3)*180/pi));
      % If the domain is a named region or a closed contour
     otherwise
      fnpl=sprintf('%s/WREG-%s-%i.mat',filoc,dom,Lmax);
      % For some of the special regions it makes sense to distinguish
      % It it gets rotb=1 here, it doesn't in LOCALIZATION
       if [strcmp(dom,'antarctica') || strcmp(dom,'contshelves')] && rotb==1 
	 fnpl=sprintf('%s/WREG-%s-%i-%i.mat',filoc,dom,Lmax,rotb);
       end
    end
  elseif iscell(dom)
    % We have a cell, and we expect the format to be {'greenland' 1.0}
    % where we have a region and want to add a buffer region around it.
    % However, if dom{2} turns out to be zero, we should ignore it.
    if dom{2}==0; h=dom{1}; else h=[dom{1} num2str(dom{2})]; end
    buf=dom{2};
    % However, if dom{2} turns out to be zero, we should ignore it.
    fnpl=sprintf('%s/WREG-%s-%i.mat',filoc,h,Lmax);
  else
    % If, instead of a string or cell, we have closed form coordinates, 
    % then make a hash from the coordinates and use it as the filename.
    h=hash(dom,'sha1');
    fnpl=sprintf('%s/%s-%i.mat',filoc,h,Lmax);  
  end
  
  if exist(fnpl,'file')==2 && ~isstr(ngl)
    load(fnpl)
    disp(sprintf('%s loaded by KERNELCP',fnpl))
 else
    if strcmp(dom,'patch')
      defval('pars',[90 75 30]*pi/180);
      % For future reference 
      th0=pars(1); ph0=pars(2); thR=pars(3);
      if th0==0
	disp('Really, should be putting in the GRUNBAUM call here')
	error('Not for polar caps! Use GRUNBAUM or SDWCAP instead')
	% BUT IN COMPARING, NOTE THAT THE SIGN MAY BE OFF
      end
      if thR>th0
	disp('Really, should be putting in the GRUNBAUM call here')
	error('Not for near-polar caps! Use GRUNBAUM, SDWCAP, then rotate')
      end
      % Convert all angles to degrees for CAPLOC only
      [lon,lat]=caploc([ph0 pi/2-th0]/pi*180,thR/pi*180,100,1);
      % Northern and Southern points, in radians
      thN=(th0-thR);
      thS=(th0+thR);
      XY=[lon lat];
    elseif strcmp(dom,'sqpatch')
      defval('pars',[30 90 10 90]*pi/180);
      thN=pars(1); thS=pars(2); phW=pars(3); phE=pars(4);
      XY=[phW pi/2-thN ; phW pi/2-thS ; phE pi/2-thS ; ...
	  phE pi/2-thN ; phW pi/2-thN]*180/pi;
    else
      defval('buf',0);
      if isstr(dom)
	% If it's a named geographical region (+buffer?) or a coordinate boundary
	defval('pars',10);
	% Run the named function to return the coordinates
	if [strcmp(dom,'antarctica') || strcmp(dom,'contshelves')] && rotb==1
	  % Return the rotation parameters also, to undo later
	  [XY,lonc,latc]=eval(sprintf('%s(%i)',dom,pars));
	else
	  % Don't, the result will be the kernel for the rotated dom
	  XY=eval(sprintf('%s(%i)',dom,pars));
	end
      elseif isstr(pars)
	defval('pars','supplyyourownfilename');
	XY=dom;
	% Use the input to define the file name that will be created
	fnpl=sprintf('%s/WREG-%s-%i.mat',filoc,pars,Lmax);
      elseif iscell(dom)
    XY=eval(sprintf('%s(%i,%f)',dom{1},pars,buf));
      else
	XY=dom;
      end
      thN=90-max(XY(:,2)); thN=thN*pi/180;
      thS=90-min(XY(:,2)); thS=thS*pi/180;
    end

    % No more set-up after this point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Introduce and dimensionalize variables and arrays
    [dems,dels,mz,lmc,mzi,mzo,bigm,bigl]=addmon(Lmax);
    dimK=(Lmax+1)^2; lenm=length(dems);
    Klmlmp=repmat(NaN,dimK,dimK);

    if isstr(ngl) && strcmp(ngl,'alternative')
      % Hold on and see if we can go more quickly here using GEOBOXCAP
      fax=2^5; 
      % Try oversampling the Nyquist degree by a certain factor
      degres=180/sqrt(Lmax*(Lmax+1))/fax;
      % Calculate the mask function
      [Bl,dels,r,lon,lat]=geoboxcap(Lmax,dom,[],degres);
      % Prepare the reindexing arrays
      % We're not fully using the recursion here, so there is wastage  
      % Perform the masked spherical harmonic transform
      h=waitbar(0,'KERNELCP: Loop over all degrees and orders');
      % With the recursions as they are, we are not yet taking full
      % advantage of this method. See Mark Wieczorek's Fortran code which
      % presumably works better for this case.
      for l=0:Lmax
	% Remember the normalization conventions
	theYplus=2*sqrt(pi)*ylm(l,0:l,(90-lat)/180*pi,lon/180*pi);
	theYmins=2*sqrt(pi)*ylm(l,-1:-1:-l,(90-lat)/180*pi,lon/180*pi);
	for m=0:l
	  waitbar((addmup(l)+m)/addmup(Lmax),h);
	  % Return the expansion coefficients in "standard" real
          % harmonics order 
	  lmcosiplus=xyz2plm((-1)^m*theYplus(:,:,m+1).*r,Lmax);
	  % Reindex the coefficients to "standard" localization kernel order
	  % and put them into where the positive m sits
	  posm=addmoff(l-1)+2*m+1;
	  % Redundant check that we are at the right order and degree
	  % difer(bigl(posm)-l)
	  % difer(bigm(posm)-m)
	  Klmlmp(posm,:)=lmcosiplus(2*size(lmcosiplus,1)+mzo)';
	  if m>0
	    % Return the expansion coefficients in "standard" real
            % harmonics order 
	    lmcosimins=xyz2plm((-1)^m*theYmins(:,:,m).*r,Lmax);
	    % Also do the negatives which come right before the positives
	    Klmlmp(posm-1,:)=lmcosimins(2*size(lmcosimins,1)+mzo)';
	  end
	end
      end
      delete(h)

      % NOTE : THIS PIECE OF THE CODE IS REPEATED VERBATIM BELOW
      % By whichever way you obtained the kernel, now check if you might want
      % to rotate it back so its eigenvectors are "right", right away,
      % e.g. for Antarctica or ContShelves without needing to rotate as
      % part of LOCALIZATION 
      if rotb==1
	disp(' ')
	disp('The input coordinates were rotated. Kernel will be unrotated,')
	disp('so its eigenfunctions will show up in the right place')
	disp(' ')
	% Get the rotation parameters for this particular region
	[XY,lonc,latc]=eval(sprintf('%s(%i)',dom,pars));
	if nargout<4
	  % Rotate the kernel, properly
	  [Klmlmp,K1]=klmlmp2rot(Klmlmp,lonc,latc);
	else
	  % Some extra verification in here
	  [Klmlmp,K1,K]=klmlmp2rot(Klmlmp,lonc,latc);
	end
      else
	[lonc,latc,K1,K]=deal(0);
      end
      
      % Do not save this way of calculating the kernels
      fnpl='neveravailable';
    else
      % Regular Gauss-Legendre method here (not alternative)  
        
      % Calculating different Gauss-Legendre points for all possible product
      % degrees is not a good idea since they get multiplied by more
      % functions of theta 
      intv=cos([thS thN]);
      nGL=max(ngl,2*Lmax);
      % These are going to be the required colatitudes - forget XY
      [w,x,N]=gausslegendrecof(nGL,[],intv);
      disp(sprintf('%i Gauss-Legendre points and weights calculated',N))
      
      % First calculate the Legendre functions themselves
      % Note that dimK==sum(dubs)
      dubs=repmat(2,lenm,1); dubs(mz)=1; comdi=[];
      % First, calculate all the Legendre functions themselves
      Xlm=repmat(NaN,length(x),lenm);
      
      % Calculate the Legendre polynomials
      ind=0;
      for l=0:Lmax  
	Xlm(:,ind+1:ind+l+1)=[legendre(l,x(:)','sch')*sqrt(2*l+1)]';
	ind=ind+l+1;
      end
      
      % Note: Xlmlmp is length ((lenm^2)+lenm)/2 because the Legendre
      % products have a redundant half (almost), and can be ordered 
      % as [0 11 222].  In order to use this with the kernel, which is
      % length (dimK^2+dimK)/2, we need an indexing array of the same
      % length which is filled with indices to Xlmlmp.  This array (bigo)
      % fills Xlmlmp back out to the ordering used in the kernel, [0 111
      % 22222].  Each Legendre polynomial has the shortened ordering, so
      % Xlmlmp essentially has redundancy in two dimensions.  When "coss"
      % is made, this will expand Xlmlmp in basically one dimension.  The
      % second pass is made when "ins" is inserted at certain points into
      % "coss" to form "bigo."  Later, the Legendre products will be multiplied by
      % different "I" matrices, representing the sine and cosine products 
      % from the longitudinal integrals.
      % 
      % For more information on this or other functions, see the 
      % Simons' group wiki page.
      
      % for lm1dex=1:lenm
      %  comdi=[comdi ; dubs(lm1dex:lenm)];
      % end
      
      % In our ordering, the -1 precedes 1 and stands for the cosine term 
      % comdex=[1:((lenm^2)+lenm)/2]';
      % coss=gamini(comdex,comdi);
      % Need a vector of length "index" that points to the right
      % combination in XlmXlmp for the next array we are
      % designing. First, find the positions we've been missing 
      h=[dimK:-1:1']'; k=find(dems); kk=k+[1:length(k)]';
      % Where to insert other elements
      inpo=[indeks(cumsum(skip(h,kk)),k)+1]';
      % How many elements to insert
      inel=h(kk);
      % Which elements to insert
      beg=inpo-h(k)+[1:length(inel)]';
      ent=inpo-h(k)+inel+[0:length(inel)-1]';
      ins=[];
      % for ind=1:length(beg)
      %  ins=[ins coss(beg(ind):ent(ind))];
      % end
      % And how to do it
      %bigo=insert(coss,ins,gamini(inpo,inel));
      % Get the longitudinal integration info for the domain
      if isstr(dom)
	switch dom
	 case 'patch'
	  % Get the parameters of the dom
	  phint=dphpatch(acos(x),thR,th0,ph0);
	 case 'sqpatch'
	  % Always the same longitudinal integration interval
	  phint=repmat([phW phE],length(x),1);
	 otherwise
	  defval('Nk',10);
	  % Now we may have multiple pairs
	  phint=dphregion(acos(x)*180/pi,Nk,dom);
	  phint=phint*pi/180;
	end
      else
	% Now we may have multiple pairs
	phint=dphregion(acos(x)*180/pi,[],XY); % Changed "dom" to "XY" here CTH
	phint=phint*pi/180;
      end
      
      % The number of elements that will be calculated is
      nel=(dimK^2+dimK)/2;
      parfor lm1dex=1:dimK
        l1=bigl(lm1dex);
        m1=bigm(lm1dex);
        % Can only use the loop variable once per index.  So for Klmlmp,
        % also use index=lm1dex
        index=lm1dex;
        ondex=0;
        I=repmat(NaN,length(x),dimK-index+1);
        % Instead of counting up andex and undex, write expressions for
        % them analytically for the row of Klmlmp we are calculating
        countdown = [dimK:-1:1];
        andex = 1 + sum(countdown(1:(index-1)));
        undex = sum(countdown(1:index));
        % We know bigo, andex, and undex, so just calculated exactly which
        % parts of XlmXlmp you need for this specific iteration
        smalll1 = abs(l1);
        smallm1 = abs(m1);
        pos1=1/2*(smalll1)^2+1/2*smalll1+smallm1+1;
        XlmXlmp = 0;
        for lm2dex=lm1dex:dimK
          l2=bigl(lm2dex);
          m2=bigm(lm2dex);
          smalll2 = abs(l2);
          smallm2 = abs(m2);
          pos2=1/2*(smalll2)^2+1/2*smalll2+smallm2+1;
          ondex=ondex+1;
          XlmXlmp(1:length(x),ondex)=Xlm(:,pos1).*Xlm(:,pos2);
          % Now evaluate the longitudinal integrals at the GL points
          if m1>0 & m2>0
            I(:,ondex)=sinsin(acos(x),m1,m2,phint);
          elseif m1<=0 & m2<=0
            I(:,ondex)=coscos(acos(x),m1,m2,phint);
          elseif m1>0 & m2<=0 % Got rid of redundant ,pars below here
            I(:,ondex)=sincos(acos(x),m1,m2,phint);
          elseif m1<=0 & m2>0
            I(:,ondex)=sincos(acos(x),m2,m1,phint);
          end
        end
        
        % Do the calculation and set as a temp variable
        temprow = (w(:)'*(XlmXlmp.*I));
        % Pad the temp variable with the appropriate zeros out front
        temprow = [zeros(1,(index-1)) temprow];
        % Now we can distribute over the kernel.  We need to do it this way
        % because if you slice Klmlmp with the loop variable (lm1dex) then
        % all other indicies need to be constant, or ':', or 'end.'
        Klmlmp(lm1dex,:)=temprow;
      end %parfor
    
      % Symmetrize the Kernel
      Klmlmp = Klmlmp + Klmlmp' - eye(size(Klmlmp)).*Klmlmp;
	
      % Verify that the first value correctly gives the area of
      % the patch 
      if strcmp(dom,'patch')
        parea=2*pi*(1-cos(thR));
        apo=abs(parea-Klmlmp(1))/parea;
        disp(sprintf(...
            'Area of the patch approximated to within %5.2f %s',...
            apo*100,'%'))
        if apo*100>1;
          error('Something wrong with the area element: radians/degrees ?')
        end
      elseif strcmp(dom,'sqpatch')
        parea=[cos(thN)-cos(thS)]*[phE-phW];
        apo=abs(parea-Klmlmp(1))/parea;
        disp(sprintf(...
            'Area of the patch approximated to within %5.2f %s',...
            apo*100,'%'))
        if apo*100>1;
          error('Something wrong with the area element: radians/degrees ?')
        end
      else
        disp(sprintf('Area of the domain approximated as %8.3e',...
                     Klmlmp(1)))
      end
      
      % To make this exactly equivalent to Tony's \ylm, i.e. undo what we
      % did above here, taking the output of YLM and multiplying
      Klmlmp=Klmlmp/4/pi;
      
      % This then makes Klmlmp(1) the fractional area on the sphere
      
      % Save this now
      if isstr(pars)
	dom=pars;
      end
      % This is where the save statement used to be
    end
    disp(sprintf('KERNELCP (Matrix)  took %8.4f s',etime(clock,t0)))

    % NOTE : THIS PIECE OF THE CODE IS REPEATED VERBATIM ABOVE
    % By whichever way you obtained the kernel, now check if you might want
    % to rotate it back so its eigenvectors are "right", right away, 
    % e.g. for Antarctica or ContShelves without needing to rotate as
    % part of LOCALIZATION 
    if rotb==1
      disp('The input coordinates were rotated. Kernel will be unrotated,')
      disp('so its eigenfunctions will show up in the right place')
      disp(' ')
      % Get the rotation parameters for this particular region
      [XY,lonc,latc]=eval(sprintf('%s(%i)',dom,pars));
      if nargout<4
	% Rotate the kernel, properly
	[Klmlmp,K1]=klmlmp2rot(Klmlmp,lonc,latc);
	K=0;
      else
	% Some extra verification in here
	[Klmlmp,K1,K]=klmlmp2rot(Klmlmp,lonc,latc);
      end
    else
      [lonc,latc,K1,K]=deal(0);
    end
    
     % This is only saved when it's not the alternative calculation method
     if ~strcmp(fnpl,'neveravailable')
       save(fnpl,'Lmax','Klmlmp','dom','ngl','XY',...
 	   'lonc','latc','K1','K')
     end
  end
  
  defval('K1',NaN)
  defval('K',NaN)
  
  % Provide output where requested
  varns={Klmlmp,XY,K1,K};
  varargout=varns(1:nargout);
elseif strcmp(Lmax,'demo1')
  L=12;
  % Construct the kernel for Antarctica rotated to the equator and not
  % rotated back
  K0=kernelc(L,'antarctica',[],[],0);
  % Excessive verification
  xver=0;
  % Construct the kernel for Antarctica in physical space which means
  % the kernel has been rotated back
  if xver==1
    % With extra steps
    [K2,XY,K1,K]=kernelc(L,'antarctica',[],[],1);
    % This intermediate result better be very nearly the same:
    difer(K-K0)
  else
     % Without the extra steps
    [K2,XY,K1]=kernelc(L,'antarctica',[],[],1);
  end

  % Compare the kernels visually
  clf
  strx='rank | degree and order';
  stry1='degree | degree and order';
  stry2='degree | degree and order';
  stry3='order | degree and order';
  [ah,ha,H]=krijetem(subnum(1,3));
  
  c11=[1 1];
  cmn=[(L+1)^2 (L+1)^2];

  axes(ah(1))
  % Regular view, degrees indicated
  imagefnan(c11,cmn,K0); axis ij
  xl(1)=xlabel(strx);
  yl(1)=ylabel(stry1); 
  t(1)=title('The original kernel');
    
  axes(ah(2))
  % Regular view, degrees indicated
  % Should look "zonal", almost all on m=0
  imagefnan(c11,cmn,K1); axis ij
  xl(2)=xlabel(strx);
  yl(2)=ylabel(stry2); 
  t(2)=title('Once rotated');

  axes(ah(3))
  % Sort orders and degrees from KERNELC to ADDMOUT
  [dems,dels,mz,lmc,mzin,mzo,Km,Kl,rinm]=addmon(L);
  % Now [Kl(rinm) Km(rinm)] is like ADDMOUT
  [EM,EL,mz,blkm,dblk]=addmout(L);
  % Now [Kl(rinm(blkm)),Km(rinm(blkm))] is block sorted
  imagefnan(c11,cmn,K2(rinm(blkm),rinm(blkm))); axis ij
  blox=Km(rinm(blkm));
  ordshew=[1; find(diff(Km(rinm(blkm))))+1];
  xl(3)=xlabel(strx);
  yl(3)=ylabel(stry3); 
  t(3)=title('Fully rotated');

  % Cosmetics
  longticks(ah)
  set(ah(1:2),'ytick',addmoff([0:L]-1)+1,'ytickl',[0:L])
  set(ah(3),'ytick',ordshew,'ytickl',num2cell(Km(ordshew)))
  set(ah,'xtick',[1 (L+1)^2],'xtickl',[1 (L+1)^2])
  fig2print(gcf,'landscape')
  figdisp([],sprintf('antarctica_%i',L))
elseif strcmp(Lmax,'demo2')
  L=18;
  % Construct the kernel for Antarctica in rotated space
  K0=kernelc(L,'antarctica',[],[],0);

  % Construct the kernel for Antarctica in physical space
  [K2,XY,K1]=kernelc(L,'antarctica',[],[],1);

  % Now check the diagonalization of the first kernel 
  [C0,V0]=eig(K0);
  [V0,isrt]=sort(sum(real(V0),1),'descend');
  C0=C0(:,isrt(1:addmoff(L)));
 
  % Now check the diagonalization of the second kernel
  [C2,V2]=eig(K2);
  [V2,isrt]=sort(sum(real(V2),1),'descend');
  C2=C2(:,isrt(1:addmoff(L)));

  % Compare the eigenvalues
  difer(V0-V2)
  
  % Compare the eigenfunctions!
  % Clean this up!
  J=5;
  [dems,dels,mz,lmc,mzin]=addmon(L);
  [XY,lonc,latc]=antarctica;

  CC=cellnan(J,length(dels),2);
  for index=1:J
    % This statement reappears exactly in LOCALIZATION
    % And the same arguments appear exactly in KLMLMP2ROT
    CC{index}=plm2rot(...
	[dels dems reshape(insert(C0(:,index),0,mzin),2,length(dems))'],...
	-lonc,latc,0);
    % Note that the rotation may lead to some roundoff errors and
    % imaginary parts in the eigenfunctions and eigenvalues!
  end
    
  % Do the plotting
  clf
  [ah,ha,H]=krijetem(subnum(J,3));

  for index=1:J
    axes(ha(index))
    % The non-polar eigenfunctions for Antarctica
    plotslep(C0,index,2);
    axes(ha(index+J))
    % The near-polar eigenfunctions for Antarctica
    p=plotslep(C2,index,2);
    axes(ha(index+2*J))
    % The rotated eigenfunctions of the non-polar kernel
    d=plotplm(CC{index},[],[],4,1);
    set(ah,'clim',halverange(d,100))
    % Let us appraise the comparison also - knowing that the sign remains
    % arbitrary
    a=minmax(abs(d-p)./d);
    b=minmax(abs(d+p)./d);
    disp(sprintf('The min/max relative difference is %8.3e | %8.3e',...
	 [max([a(1) b(1)]) min([a(2) b(2)])]))
  end  

  % Cosmetics
  nolabels(ha(J+1:end),2)
  nolabels(ah(1:end-3),1)
  longticks(ah)
  fig2print(gcf,'tall')
elseif strcmp(Lmax,'demo3')
  % Reconcile the indexing of KERNELC | KLMLMP and GLMALPHA
  Klmlmp=kernelc(18,'australia');
  % Actually, just see inside GLMALPHA or LOCALIZATION, that will do
elseif strcmp(Lmax,'demo4')
  L=18;
  [Klmlmp,XY,K1,K]=kernelc(L,'antarctica',[],[],1);
  % Diagonalize
  [C,V]=eig(Klmlmp);
  % What's up with the eigenvalues? They have small imaginary parts when
  % their magnitude is small. 
  difer(eye(size(C))-(C'*C))
  % If it fails the test, it will be due to this numerical degeneracy!

  % However, the large-eigenvalue ones should still look good!
  % These rotated kernels (whose eigenfunctions should be right on
  % Antarctica without any further ado!) should be orthonormal!
  [dems,dels,mz,lmc,mzin]=addmon(L);
  % Preallocate/initialize cell array
  CC=cellnan(J,length(dels),2);
  % See in LOCALIZATION and GLMALPHA
  for index=1:size(C)
    % Note that you can undo this using the rinm variable in ADDMON
    CC{index}=reshape(insert(C(:,index),0,mzin),2,length(dems))';
    plotplm([dels dems CC{index}],[],[],2,0.5); view(145,-65)
    pause
  end
elseif strcmp(Lmax,'demo5')
  Lmax=22;
  dom='greenland';
  K2=kernelc(Lmax,dom);
  K1=kernelcp(Lmax,dom);
  difer(K2-K1)
end

