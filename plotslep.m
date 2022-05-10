function varargout=plotslep(G,i,fmt,degres)
% data=PLOTSLEP(G,i,fmt,degres)
%
% Plots Slepian functions coming out of GLMALPHA, GLMALPHAPTO, or those
% that are the direct result of diagonalizing the kernel from KERNELC
%
% INPUT:
%
% G          Output from GLMALPHA or GLMALPHAPTO [default], OR
%            Output from diagonalizing the result of KERNELC
% i          The rank of the function you want plotted [default: 1]
% fmt        1 The input is the output of GLMALPHA|PTO [default]
%            2 The input is from diagonalizing KERNELC
% degres     Pixelation of the plot as input to PLOTPLM
%
% OUTPUT:
%
% data       Spatial data in case coefficients were given
%
% SEE ALSO: GLM2LMCOSI, KLM2LMCOSI
%
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012

% Set defaults
defval('i',1)
defval('fmt',1)

% Produce and empty array
[dems,dels,~,lmcosi,mzin,~,~,~,~,ronm]=addmon(sqrt(length(G))-1);

switch fmt
 case 1
  % Construct from the GLMALPHA standard-ordered coefficients
  % See GLM2LMCOSI
  lmcosi(2*length(lmcosi)+ronm)=G(:,i);
 case 2
  % Construct from the KERNELC standard-ordered coefficients
  % See LM2LMCOSI
  lmcosi(:,3:4)=reshape(insert(G(:,i),0,mzin),2,length(dems))';
end

% Do the standard plotting routine
defval('meth',4)
defval('degres',1)
% This should be renormalized!!
data=plotplm(lmcosi,[],[],meth,degres);

% Prepare output if needed
varns={data};
varargout=varns(1:nargout);
