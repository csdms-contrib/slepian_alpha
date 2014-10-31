function [K2,K1,K]=klmlmp2rot(Klmlmp,lonc,latc)
% [K2,K1,K]=KLMLMP2ROT(Klmlmp,lonc,latc)
% 
% INPUT:
%
% Klmlmp     A localization kernel coming out of, e.g. KERNELC, or SDWCAP
% lonc       A "longitudinal" rotation parameter 
% latc       A "latitudinal" rotation parameter
%
% OUTPUT:
%
% K2         The same kernel rotated such that its eigenfunctions are
%            rotated versions of the eigenfunctions of the unrotated kernel
% K1         The same kernel in an intermediate step that should reveal
%            the geometry of the concentration region when plotted
% K          This should be the same as the original kernel down to eps
%
% When a region (e.g. Antarctica or ContShelves) is rotated to a nonpolar
% position in order for standard KERNELC to be able to compute a
% localization matrix, its solutions will need to be rotated back, as in
% LOCALIZATION. However, here in KLMLMP2ROT we will rotate the kernel
% itself so its solutions give the required near-polar eigenfunctions
% right away. The  advantage here is that we can combine this with
% localization kernels for any other region in an unrotated coordinate
% frame, e.g. to add the localization kernels for Africa and
% Antarctica. The other advantage is that we can use the decomposition
% directly by matrix multiplication the spherical harmonics
% coefficients. The disadvantage is that finite-precision rotation will
% lead to a slight non-Hermiticity of the matrix, which might show.
%
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012

% Get the maximum degree for this square, symmetric kernel
Lmax=sqrt(size(Klmlmp,1))-1;
% Get the right indices ready
[dems,dels,mz,lmc,mzin,mzo]=addmon(Lmax);

% Apply the rotation on the ROWS
% Signs flipped on 10/19/2010 to ensure consistency with
% ANTARCTICA, LOCALIZATION and KERNELC 
alp=-lonc; bta=latc; gam=0;

% Do this once
disp('First rotation')
K1=rotateonce(Klmlmp,alp,bta,gam,dems,dels,mzin,mzo);

% For a near-zonal function the structure should now be apparent
% Note that norm3d is just
% f=sqrt(sum(sum(sum(x.^2))));

% Check the norm, it should not have changed
difer(norm3d(K1)-norm3d(Klmlmp),[],[],NaN)
% Nor should the norm per column have changed
difer(sum(K1.^2,1)-sum(Klmlmp.^2,1),[],[],NaN)

% Apply the rotation on the COLUMNS thus the transpose
disp('Second rotation')
K2=rotateonce(K1',alp,bta,gam,dems,dels,mzin,mzo);

% Check the norm, it should not have changed
difer(norm3d(K2)-norm3d(K1),[],[],NaN)
% Nor should the norm per column have changed
difer(sum(K2.^2,1)-sum(K1'.^2,1),[],[],NaN)

% Check the symmetry of the output
difer(K2-K2',[],[],NaN)

% Only if you want to test
if nargout>2
  % This still needs to be fixed as of 10/20/10, though it's likely to be
  % unimportant just how as we know it does work. Simple sign flip?
  % Apply a DIFFERENT rotation on the ROWS of K1
  alp=0;
  bta=latc;
  gam=-lonc;
  disp('Third rotation (verification round)')
  K=rotateonce(K1,alp,bta,gam,dems,dels,mzin,mzo);
  % And check that you get the input back
  difer(K-Klmlmp,[],[],NaN)
  % but it ain't going to be identical, for sure
  % minmax(K-Klmlmp)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K=rotateonce(K,alp,bta,gam,dems,dels,mzin,mzo)
% Initialize temporary structure
Ktmp=cellnan(size(K,1),length(dels),2);

tic
% First, rearrange the ROWS into LMCOSI format for each COLUMN
parfor index=1:size(K,1)
  Ktmp{index}=reshape(insert(K(:,index),0,mzin),2,length(dems))';
end
% Now perform the counterrotation as per alp, bta and gam
% Should really write something that avoids having to loop here
%h=waitbar(0,sprintf('Performing spherical harmonic rotations'));

parfor index=1:size(K,1)
  % It is this PLM2ROT call that we need for the functions also, see
  % e.g. LOCALIZATION, where it is applied to the eigenfunctions instead 
  %  Ktmp{index}=kindeks(plm2rot([dels dems Ktmp{index}],...
  %			      alp,bta,gam),3:4);
  % And then stick these guys back into the same format as they came in
  % See PLOTSLEP where I've done this before
  %  K(:,index)=Ktmp{index}(mzo);
  % Or all at once, really, which is slightly faster (but not much)
  K(:,index)=indeks(kindeks(plm2rot([dels dems Ktmp{index}],...
			     alp,bta,gam),3:4),mzo);
  %waitbar(index/size(K,1),h)
end
%close(h)
toc
