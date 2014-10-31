function varargout=fftaxis(fsize,ftsize,spun)
% [xfaks,yfaks,fnx,fny,xsint,ysint]=FFTAXIS(fsize,ftsize,spun)
%
% Makes linear frequency axis going through zero [at floor((dim+1)/2)]
% -f_N < f <= f_N with f_N=1/(2Dt)
%
% INPUT:
%
% fsize     is the SIZE of the field (Y X)
% ftsize    is the SIZE of its transform (kY kX)
% spun      is the physical dimension of the field (lY lX)
%
% OUTPUT:
%
% fnx,fny   is the Nyquist rate in both dimensions
%
% SEE ALSO:
%
% KNUM2, FFTAXIS1D, which are actually more like Matlab does it. 
% So this function should start becoming OBSOLETE.
%
% For use after FFTSHIFT, for non-geographical data, where the "-x,-y"
% or "first" is in the upper left, as opposed to geographical data, where
% the "-x,+y" or "first" is in the upper left. Thus, for geographical
% data, this needs to be FLIPUD. For 1-D and 2-D but see also FFTAXIS1D
% and KNUM2. The advantage is that the frequencies are symmetric around
% the zero-center, as they should. The magnitude of the Fourier transform
% of a discrete-time signal is always an even function. Analogously, for
% 2D, the center is a symmetry point. The lowest frequency resolvable
% that is not the DC-component is the Rayleigh frequency, given by
% 1/T=1/NDt with T the data length. 
%
% Last modified by fjsimons-at-alum.mit.edu, 04/16/2010

% See Percival and Walden, 1993

% Figure sampling interval in both directions
if nargout>3
  xsint=spun(2)/(fsize(2)-1);
  ysint=spun(1)/(fsize(1)-1);
  % Calculate centered frequency axis (PW p 112)
  M=ftsize(1);
  N=ftsize(2);
  intvectX=linspace(-floor((N+1)/2)+1,N-floor((N+1)/2),N);
  intvectY=linspace(-floor((M+1)/2)+1,M-floor((M+1)/2),M);
  xfaks=intvectX/N/xsint;
  yfaks=intvectY/M/ysint;
  % Calculate Nyquist frequencies
  fnx=1/2/xsint;
  fny=1/2/ysint;
  varns={xfaks,yfaks,fnx,fny,xsint,ysint};
else
  xsint=spun/(max(fsize)-1);

  % Calculate centered frequency axis (PW p 112)
  N=max(ftsize);
  intvectX=linspace(-floor((N+1)/2)+1,N-floor((N+1)/2),N);
  % In other words the above is equal to
  % -floor((N+1)/2)+1+[0:N-1]
  xfaks=intvectX/N/xsint;
  
  % Calculate Nyquist frequencies
  fnx=1/2/xsint;
  varns={xfaks,fnx,xsint};
end

% So if we put in
% [a,b,c,d,e,f]=fftaxis([101 51],[100 111],[100 50])
% we see how the frequency goes from -1/2 to 1/2;
% alternatively the angular frequency goes from -pi to pi;
% we can also plot the frequency normalized by the Nyquist:
% in that case, this scale should go from -1 to 1.

% Prepare output
varargout=varns(1:nargout);
