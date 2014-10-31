function varargout=rotateGp(G,lonc,latc)
% [G2]=rotateGp(G,lonc,latc)
%
% Rotates a G matrix normally returned by GLMALPHA, in parallel, so that 
% it can be used to project into a Slepian basis.
%
% INPUT:
%
% G       A G matrix of eigenvectors made from some kernel.  The format and
%          ordering of this is assumed to be the same as what is returned
%          from GLMALPHA
% lonc    The amount by which you need to rotate it back over z
% latc    The amount by which you need to rotate it back over y.  Both lonc
%          and latc are as produced from ANTARCTICA.m
%
% OUTPUT:
%
% G2        The unitary matrix of localization coefficients
%
% SEE ALSO:
%
% KLMLMP2ROT, PLM2ROT, ADDMOUT, ADDMON, KERNELC, LOCALIZATION, GALPHA, DLMLMP
%
% Last modified by charig-at-princeton.edu, 06/26/2012
% Last modified by fjsimons-at-alum.mit.edu, 07/13/2012

defval('G','glmalpha(''antarctica'',20)');
defval('xver',0)  
defval('lonc','antarctica(10)')

if isstr(G)
  % Evaluate the specified expression
  [G] = eval(G);
end
if isstr(lonc)
  % Evaluate the specified expression
  [XY,lonc,latc] = eval(lonc);
end

defval('L',(sqrt(size(G,2))-1));
defval('J',size(G,2))

% See if we can run this calculation in parallel
try
  matlabpool open
  disp('Running rotateGp in parallel')
catch
  error('Run ROTATEG instead or close your open pool')
end

% Prepare for differently-ordered degree and order output
[dems,dels,mz,lmcosi,mzin,mzo,bigm,bigl,rinm,ronm,demin]=addmon(L);
% Preallocate/initialize cell array
CC=cellnan(J,length(dels),2);  
  
% Collect the eigenvector output into a format that PLM2XYZ knows how to interpret
for j=1:size(G,2)
  % Create the blanks
  cosi=lmcosi(:,3:4);
  % Stick in the coefficients of the 1st eigentaper
  cosi(ronm)=G(:,j);
  % Construct the full matrix
  CC{j} = cosi; 
end

% You can plot this here, if you want, by doing, e.g.
% plotplm([dels dems CC{1}],[],[],4,1);

% Rotate the eigenfunctions to the new location
% NOTE: this will error if you have not run plm2rot before for this degree.
%  This is because the first time it is run, it creates a file like
%  $IFILES/DLMB/dlmb-005.mat, and if you run a parfor loop before this file
%  is made, then some function calls try to load it while it is in the
%  process of being made.  Hence it errors, and we fix this by calling
%  index=1 beforehand to guarantee the file is made properly.

% See under LOCALIZATION also
CC{1}=kindeks(plm2rot([dels dems CC{1}],-lonc,latc,0),3:4);
parfor index=2:J
  % This here was changed 10/18/2010 to reflect the changed
  % conventions in PLM2ROT
  CC{index}=kindeks(plm2rot([dels dems CC{index}],-lonc,latc,0),3:4);
end
  
% You can plot this here, if you want, by doing, e.g.
% plotplm([dels dems CC{1}],[],[],2,1); view(145,-35)

% Now we go back to the G matrix format
for j=1:size(G,2)
  % Get the coefficients
  cosi = CC{j};
  % Remove the m=0 sine coefficients
  cosinozero = cosi(mzo);
  % Reorder from [0 0-11 0-11-22] to [0 -101 -2-1012]
  Grot(:,j) = cosinozero(rinm);
end

% Close your matlab pool
try
  matlabpool close
end

% Collect output
varns={Grot};
varargout=varns(1:nargout); 
