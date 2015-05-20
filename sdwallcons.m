function [F,V,K,Vc,Kc]=sdwallcons(dom,L,mkmov)
% [F,V,K,Vc,Kc]=SDWALLCONS(dom,L,mkmov)
%
% Simons, Dahlen and Wieczorek (2006): Localization to the continents
%
% INPUT:
%
% dom     Cell with domains (strings) you want included [default: all]
% L       Maximal spherical harmonic bandwidth [default: 18]
% mkmov   1 Make a movie of this
%         0 Don't [default]
%
% OUTPUT:
%
% F       The eigenvalue-weighted sum of the energies
% V       The eigenvalues
% K       The localization kernel 
% Vc      The complement of V
% Kc      The complement of K
%
% Last modified by fjsimons-at-alum.mit.edu, 08/27/2009
% Last modified by charig-at-princeton.edu, 05/14/2015

defval('dom',...
       {'africa','eurasia','namerica','australia','greenland', ...
	'samerica'});
defval('L',18);
defval('mkmov',0);

if ~iscell(dom)
  dom={dom};
end

% Initialize the kernel matrix
K=repmat(0,(L+1)^2,(L+1)^2);
for index=1:length(dom)
  disp(sprintf('Working on %s',dom{index}))
  [Klmlmp,XY]=kernelcp(L,dom{index});
  K=K+Klmlmp;
  clear Klmlmp
end

% Stolen from LOCALIZATION
[C,V]=eig(K);
[V,isrt]=sort(sum(real(V),1));
V=fliplr(V);
C=C(:,fliplr(isrt));

% Sticks the cosine/sine coefficients back
% into the right place of LMCOSI
[dems,dels,mz,lmc,mzin]=addmon(L);
for index=1:size(C,2)
  CC{index}=reshape(insert(C(:,index),0,mzin),2,length(dems))';
end

% Eigenvalues and eigenfunctions
V=V(:);
C=CC;

if nargout>3
  Kc=eye(size(K))-K;
  [Cc,Vc]=eig(Kc);
  [Vc,isrt]=sort(sum(real(Vc),1));
  Vc=fliplr(Vc);
  Cc=Cc(:,fliplr(isrt));
  for index=1:size(Cc,2)
    CCc{index}=reshape(insert(Cc(:,index),0,mzin),2,length(dems))';
  end
  
  % Eigenvalues and eigenfunctions
  Vc=Vc(:);
  Cc=CCc;
end

% Provide the progressive spatial sum
F=0;

if mkmov==1
  % MAKE SURE THE FIGURE FILE IS UNINTERRUPTED VIEW
  mov=avifile('sdwmovie.avi');
end
for index=1:(L+1)^2
  fnpl=sprintf('%s/SDWALLCONS-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'KERNELCP'),L,index);
  if exist(fnpl,'file')==2
    FF=load(fnpl);
    F=F+FF.F;
    % This doesn't work here
    % F(F<max(F(:))/100)=NaN;
    imagefnan([0 90],[360 -90],F,[],minmax(F))
    plotcont
    title(sprintf('%3.3i',index)) 
    axis tight
    % pause
    if mkmov==1
      % Make and save movie
      % mov(index)=getframe;
      % movie(mov)
      mov=addframe(mov,getframe(gcf));
    end
  else
    % I could be using the precomputation - see LOCALIZATION demo5
    F=F+plm2xyz([dels dems C{index}],1).^2*V(index);
    % eval(sprintf('save %s F',fnpl))
    save(fnpl,'F')
  end
end
if mkmov==1
  mov=close(mov);
end

