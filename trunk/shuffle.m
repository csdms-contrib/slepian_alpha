function [outmat,unshuf]=shuffle(inmat)
% [outmat,unshuf]=SHUFFLE(inmat)
% 
% Shuffles the rows of a matrix or vector
%
% SEE ALSO: UNSHUFFLE
%
% EXAMPLE:
%
% inmat=outmat(unshuf,:)
%
% Last modified by fjsimons-at-alum.mit.edu, 11/17/2008

[m,n]=size(inmat);

% Make sure the dimensions are right
if m==1
  inmat=inmat(:);
  [m,n]=size(inmat);
end

inmat=[randn(m,1) inmat];
[outmat,unshuf]=sortrows(inmat);
outmat=outmat(:,2:1+n);

unshuf=[unshuf [1:m]'];
[junk,unshuf]=sortrows(unshuf);
  
if ~all(outmat(unshuf,:)==inmat(:,2:1+n)) ; error ; end
  
  
