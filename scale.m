function scalmat=scale(inmat,ranscale)
% scalmat=SCALE(inmat,ranscale)
% 
% Rescales a REAL-VALUED matrix or vector to a new range
%
% INPUT:
%
% inmat      Input data
% ranscale   Scaling (default: [-1 1])
%
% EXAMPLE:
%
% plot(scale(scale(sin(0:0.01:3*pi),[1898 9100]),[-1 1]),'r-')
%
% Last modified by fjsimons-at-alum.mit.edu, 08/15/2012

defval('ranscale',[-1 1])

if diff(ranscale)<=0
  ranscale=fliplr(ranscale);
end

reensj=max(inmat(:))-min(inmat(:));
minim=min(inmat(:));

scalmat=ranscale(1)+(inmat-minim)*range(ranscale)/reensj;

