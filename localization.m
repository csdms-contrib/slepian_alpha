function varargout=localization(L,dom,N,J,rotb,anti)
% [V,C,dels,dems,XY,Klmlmp,G]=LOCALIZATION(L,dom,N,J,rotb,anti)
%
% Returns bandlimited spectral eigenfunctions and their associated
% eigenvalues localized to a closed domain on the unit sphere.
%
% INPUT:
%
% L         Bandwidth, maximum angular spherical harmonic degree
% dom       'patch', 'sqpatch', 'africa', 'eurasia', 'namerica', 
%           'australia', 'greenland', 'samerica', 'amazon', 'orinoco',
%           'gpsnamerica', 'antarctica', 'alloceans', 'contshelves'
%           OR: [lon lat] an ordered list defining a closed curve [degrees]
% [N |pars] The splining smoothness [default: 10] or
%           The 3 or 4 parameters for 'patch' and 'sqpatch' passed to KERNELC,
%           go look there for the definitions
% J         Loads only so many eigenfuctions as required [default: all]
% rotb      0 That's it, nothing else needs to be said here
%           1 The obtained eigenfunctions are subsequently rotated by an
%             amount that is specified by the region name. Currently this
%             is only sometimes necessary and even then only possible for
%             'antarctica' and 'contshelves'. The alternative approach
%             might start from KERNELC and KERNELCP with the rotation in
%             it but then this uses KLMLMP2ROT which renders the resulting
%             kernel slightly non-Hermitian so the option where
%             KERNELC(P) is called without rotation (as it is now inside
%             this function) and LOCALIZATION is called with, is preferred
%             for these regions here. You might notice the difference.
% anti      1 Get the region complementary to the requested region
%           0 Do not, proceed as usual
%
% OUTPUT:
%
% V        A list of all the eigenvalues 
% C        A cell array with cosine/sine coefficients eigenfunctions;
%          note how GLMALPHA will unwrap these into a neat matrix
% dels     Spherical harmonic degrees 
% dems     Spherical harmonic orders 
% XY       Coordinates of the outlines of the region
% Klmlmp   The localization matrix whose eigenfunctions are G and thus C
% G        The entire matrix of spectral eigenfunctions as from GLMALPHA
%
% EXAMPLES:
%
% localization('demoX') % where X=1,2,3,4,5,6, or 7
%
% SEE ALSO:
%
% KERNELC, KERNELCP, SDWCAP, SDWCAP2, CDLKERNEL, SDWCAPT, SDWCAPT2,
% LOCALIZATION2D, PLOTPLM, PLM2XYZ, PLOTSLEP, KLMLMP2ROT, GLMALPHA,
% ROTATEGP 
% 
% Last modified by fjsimons-at-alum.mit.edu, 1/22/2016 

% Study covariance at some point?

% Default inputs
defval('L',18)

if ~isstr(L)
  defval('dom','australia')
  defval('N',10)
  defval('rotb',0)
  defval('anti',0)
  
  % If it's been pre-done this isn't strictly necessary anymore - is it?
  % You could thus avoid calling KERNELC(P) if you don't want this
  % output and if for some reason you threw the kernel away. 
  
  % Calculates the localization kernel for this domain
  if nargout>4
    % Note that the "rotb" parameter is always set to 0 by default and that
    % you then may specify whether or not to counter-rotate the
    % coefficients. If you want to rotate the kernel, don't change this
    % option here but rather call KERNELC or KERNELCP independently.
    try
      [Klmlmp,XY]=kernelcp(L,dom,N);
    catch
      [Klmlmp,XY]=kernelc(L,dom,N);
    end
    
    % We're going to here interpret the call 'alloceans' as needing to
    % further subtract some continental kernels
    if strcmp(dom,'alloceans')
      regs={'samerica','africa','australia'};
      for ind=1:length(regs)
        disp(sprintf('Removing also %s',regs{ind}))
        try 
          Klmlmp=Klmlmp-kernelcp(L,regs{ind},N);
        catch
          Klmlmp=Klmlmp-kernelc(L,regs{ind},N);
        end
      end
    end
    defval('J',length(Klmlmp))
  else
    defval('XY',NaN)
    defval('Klmlmp',NaN)
    defval('G',NaN)
  end
  defval('J',(L+1)^2)
  
  % If XY is just coordinates we come up with a hash to store the eigenfunctions.
  if ~isstr(dom)
    doms=hash(dom,'sha1');
  else
    doms=dom;
  end

  % See if the diagonalization has been done before
  switch doms
   case 'sqpatch'
    fnpl=sprintf('%s/%s-%i-%i-%i-%i-%i.mat',...
		 fullfile(getenv('IFILES'),'LOCALIZE'),doms,L,...
		 round(N(1)*180/pi),round(N(2)*180/pi),...
		 round(N(3)*180/pi),round(N(4)*180/pi),J);
   case 'patch'
    fnpl=sprintf('%s/%s-%i-%i-%i-%i-%i.mat',...
		 fullfile(getenv('IFILES'),'LOCALIZE'),doms,L,...
		 round(N(1)*180/pi),round(N(2)*180/pi),...
		 round(N(3)*180/pi),J);
   otherwise
    fnpl=sprintf('%s/%s-%i-%i-%i.mat',...
		 fullfile(getenv('IFILES'),'LOCALIZE'),doms,L,N,J);
  end
  if anti==1
    Klmlmp=eye(size(Klmlmp))-Klmlmp;
    fnpl=sprintf('%s/%s-%i-%i-%i-%i.mat',...
		 fullfile(getenv('IFILES'),'LOCALIZE'),doms,L,N,J,anti);
    disp(sprintf('Remember that you requested the anti-%s region',doms))
  end
  
  if exist(fnpl,'file')==2
    disp(sprintf('Loading %s',fnpl))
    load(fnpl)
  else
    try
      % Calculates the eigenfunctions/values for this localization problem
      if J<length(Klmlmp)-1 & 1==3
	% Note that this often does not work very well at all, so don't try
	OPTS.disp=0;
	[C,V]=eigs(Klmlmp,J,'LA',OPTS);
	[V,isrt]=sort(sum(real(V),1),'descend');
	C=C(:,isrt(1:J));
      else
	% This is slow but better
	[C,V]=eig(Klmlmp);
	[V,isrt]=sort(sum(real(V),1),'descend');
	% Global sorting
	C=C(:,isrt(1:J));
	save(fnpl,'C','V','dom','L','N','J')
      end
    catch
      error('Better call KERNELC or KERNELCP so that you have a saved kernel')
    end
  end
  
  % Prepare for differently-ordered degree and order output
  [dems,dels,mz,lmc,mzin]=addmon(L);
  % Preallocate/initialize cell array
  CC=cellnan(J,length(dels),2);
  % Sticks the cosine/sine coefficients back
  % into the right place of LMCOSI
  % See also PLOTSLEP or, ultimately, KLM2LMCOSI
  parfor index=1:J
    % Note that you can undo this using the rinm variable in ADDMON
    CC{index}=reshape(insert(C(:,index),0,mzin),2,length(dems))';
  end

  % You can plot this here, if you want, by doing, e.g.
  % plotplm([dels dems CC{10}],[],[],4,1); hold on
  % [X,Y,Z]=sph2cart(XY(:,1)*pi/180,XY(:,2)*pi/180,1);
  % plot3(X,Y,Z)

  % Also check out ROTATEGP
  
  % It it's Antarctica or ContShelves, need to rotate it back in shape
  if [strcmp(dom,'antarctica') || strcmp(dom,'contshelves')] && rotb==1
    % But see now also KERNELC and KLMLMP2ROT, and ANTARCTICA
    % Return the rotation parameters also, to undo later
    [XY,lonc,latc]=eval(sprintf('%s(%i)',dom,N));
    % h=waitbar(0);
    try
      matlabpool open 
    end
    parfor index=1:J
      % This here was changed 10/18/2010 to reflect the changed
      % conventions in PLM2ROT
      CC{index}=kindeks(plm2rot([dels dems CC{index}],...
				-lonc,latc,0),3:4);
      %waitbar(index/J,h,...
      %        sprintf('Rotated eigenfunction %i/%i',index,J)) 
    end
  end
  % You can plot this here, if you want, by doing, e.g.
  % plotplm([dels dems CC{10}],[],[],2,0.5); view(145,-35)
  % This now does show up in the right spot

  % Prepare for standard output
  V=V(:);

  % And assign it - keep the cell structure and the big matrix
  varns={V,CC,dels,dems,XY,Klmlmp,C};
  varargout=varns(1:nargout);
elseif strcmp(L,'demo1')
  L=36/2;
  % Must get nargout to force running KERNELCP
  [V,C,dels,dems,XY,Klmlmp]=localization(L,'africa',[],[],[],1);
  subplot(121)
  plotplm([dels dems C{1}],[],[],[],1)
  % They look weird - you are in a degenerate eigenspace
  % While we had them (just didn't return them) let us reconstruct the
  % original eigenfunctions of the localization matrix by undoing the
  % degree-order ordering scheme that made us get this output
  subplot(122)
  [~,~,~,~,~,~,~,~,~,rinm]=addmon(L);
  bigC=nan(length(C),length(C));
  for i=1:length(C)
    bigc(:,i)=C{i}(rinm);
  end
  imagefnan(bigc'*bigc)
elseif strcmp(L,'demo2')
  [V,C,dels,dems,XY,Klmlmp]=localization(18,'amazon');
  plotplm([dels dems C{1}])
elseif strcmp(L,'demo3')
  th0=114; ph0=134; thR=40;
  [V,C,dels,dems,XY,Klmlmp]=localization(6,'patch',[th0 ph0 thR]*pi/180);
  plotplm([dels dems C{1}],[],[],[],2); hold on; plot(ph0,90-th0,'w')
elseif strcmp(L,'demo4')
  L=60/3; J=100; rotb=1;
  [V,C,dels,dems,XY,Klmlmp]=localization(L,'antarctica',[],J,rotb);
  subplot(121)
  plotplm([dels dems C{12}],[],[],2,0.5)
  view(145,-35)
  subplot(122)
  [~,~,~,~,~,~,~,~,~,rinm]=addmon(L);
  bigC=nan(length(C),length(C));
  for i=1:length(C)
    bigc(:,i)=C{i}(rinm);
  end
  imagefnan(bigc'*bigc)
elseif strcmp(L,'demo5')
  % Verify that LOCALIZATION and GLMALPHA give essentially the same
  % information 
  L=18;
  [V,C,dels,dems,XY,Klmlmp,GG]=localization(L,'samerica');
  [G,VG,EL,EM,N,GM2AL,MTAP,IMTAP]=glmalpha('samerica',L);
  difer(V(:)-VG(:))
  [~,~,~,~,~,~,~,~,R1,R2]=addmon(L);
  % It's only a matter of indexing and ordering
  difer(GG(R1,:)-G)
elseif strcmp(L,'demo6')
  L=72/3;
  % What will be the Shannon number in this construction?
  N=[spharea('alloceans')-spharea('africa')-spharea('samerica')-...
    spharea('australia')]*(L+1)^2;
  % The difference lies in the second polar gap
    
  % Computes all Slepian functions for the oceans
  [V,C,dels,dems,XY,Klmlmp,G]=localization(L,'alloceans');

  % Compare this to what the procedure tells us
  sum(V)

  % Make some plots
  for index=1:10:size(C,2)
    clf
    [d,c,p]=plotplm([dels dems C{index}],[],[],[],1);
    delete(p)
    title(sprintf('Slepian function %i',index))
    disp('Press enter to continue...')
    pause
  end

  % Compute the eigenvalue-weighted sum of squares
  [F,lo,la,Plms]=plm2xyz([dels dems C{1}],1); F=F.^2*V(1);
  for index=2:size(C,2)
    F=F+plm2xyz([dels dems C{index}],1,[],[],[],Plms).^2*V(index);
  end
  % Then plot it all
  clf
  [d,c,p]=plotplm(F,[],[],4,1);
  delete(p)
  title(sprintf('Eigenvalue-weighted sum of %i Slepian functions',index))

  % Look at the localization matrix, after order-ordering as GLMALPHA
  % Find row indices into G belonging to the orders
  [EM,EL,mz,blkm]=addmout(L);
  imagefnan(Klmlmp)
  
  % Look at the orthogonality of the matrices when split
  imagefnan(G'*G)
  N=round(N);
  imagefnan(G(1:N,1:N)'*G(1:N,1:N))
  imagefnan(G(N+1:end,N+1:end)'*G(N+1:end,N+1:end))
elseif strcmp(L,'demo7')
  L=18;
  % Computes all Slepian functions for the continental shelves
  [V,C,dels,dems,XY,Klmlmp,G]=localization(L,'contshelves',[],[],1);

  % Compute the eigenvalue-weighted sum of squares
  [F,lo,la,Plms]=plm2xyz([dels dems C{1}],1); F=F.^2*V(1);
  parfor index=2:size(C,2)
    F=F+plm2xyz([dels dems C{index}],1,[],[],[],Plms).^2*V(index);
  end
  
  % Then plot it all
  clf
  [d,c,p]=plotplm(F,[],[],4,1);
  delete(p)
  title(sprintf('Eigenvalue-weighted sum of %i Slepian functions',size(C,2)))
  hold on
  % Just use CONTSHELVES([],1)
  whereitsat=fullfile(getenv('IFILES'),'COASTS');
  XY=load(fullfile(whereitsat,'ContShelves.txt'));
  XY(:,1)=XY(:,1)+[XY(:,1)<0]*360;
  [X,Y]=penlift(XY(:,1),XY(:,2));
  plot(X,Y,'Linew',2,'Color','k');
  hold off
else
  error('Specify valid demo number')
end


