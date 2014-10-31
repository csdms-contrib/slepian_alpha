function [w,wl,wr]=fhanning(n)
% [w,wl,wr]=fhanning(n)
%
% Calculates Hanning windows of a certain length
%
% INPUT:
%
% n       The required length of the window
%
% OUTPUT:
%
% w       The Hanning window for bandpass
% wl      The left half of the window for lowpass
% wr      The right half of the window for lowpass
%
% Last modified by fjsimons-at-alum.mit.edu, 08/05/2014

if ~rem(n,2)
   % Even length window
   half = n/2;
   wl = channing(half,n);
   w = [wl; wl(end:-1:1)];
else
   % Odd length window
   half = (n+1)/2;
   wl = channing(half,n);
   w = [wl; wl(end-1:-1:1)];
end

if nargout>=3
  wr=flipud(wl);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = channing(m,n)
w = .5*(1 - cos(2*pi*(1:m)'/(n+1)));


