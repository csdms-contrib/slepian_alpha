function varargout=imagefdir(c11,cmn,matrix)
% IMAGEFDIR(c11,cmn,matrix)
% H=IMAGEFDIR(...)
%
% Like IMAGEF but NOT SCALED - indexes colormap directly
% Also works with 3-D matrix rgb input.
%
% Last modified by fjsimons-at-alum.mit.edu, Feb 7th, 2001

cm1=[c11(1) cmn(2)];
c1n=[cmn(1) c11(2)];

h=image([cm1(1) c1n(1)],[cm1(2) c1n(2)],flipdim(matrix,1));  axis xy

if nargout
  varargout{1}=h;
end



