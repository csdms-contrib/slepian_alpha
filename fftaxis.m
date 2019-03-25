function varargout=fftaxis(fsize,ftsize,pl)
% [xfaks,yfaks,fnx,fny,dx,dy]=FFTAXIS(fsize,ftsize,pl)
% [xfaks,fnx,dx]=FFTAXIS(fsize,ftsize,pl)
%
% Makes a one-or-two-dimensional linear frequency (not angular wavenumber) axis
% going through zero at floor((dim+1)/2), with -fnx < xfaks <= fnx and fnx=1/2/dx
% If fewer than three outputs are requested, makes a one-dimensional axis
% using the larges size and dimensions and physical lengths input
%
% INPUT:
%
% fsize     is the SIZE of the two-dimensional field (Y X)
% ftsize    is the SIZE of its Fourier transform (fY fX)
% pl        is the physical dimension of the field (lY lX)
%
% OUTPUT:
%
% xfaks, yfaks   the frequency axis
% fnx, fny       the Nyquist frequencies
% dx, dy         the sampling intervals
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
% EXAMPLE:
% 
%% So... if we put in
% [xfaks,yfaks,fnx,fny,dx,dy]=fftaxis([101 51],[100 111],[100 50])
%% we see how the frequency ranges between from -1/2 to 1/2;
%% alternatively the angular frequency goes from -pi to pi;
%% we can also plot the frequency normalized by the Nyquist:
%% in that case, this scale should go from -1 to 1.
%
% Last modified by fjsimons-at-alum.mit.edu, 03/19/2019

% See Percival and Walden, 1993

% Figure sampling interval in both directions
if nargout>3
  dy=pl(1)/(fsize(1)-1);
  dx=pl(2)/(fsize(2)-1);
  % Calculate centered frequency axis (PW p 112)
  M=ftsize(1);
  N=ftsize(2);
  % Makes the axis
  intvectX=linspace(-floor((N+1)/2)+1,N-floor((N+1)/2),N);
  intvectY=linspace(-floor((M+1)/2)+1,M-floor((M+1)/2),M);
  % Scales the axis
  xfaks=intvectX/N/dx;
  yfaks=intvectY/M/dy;
  % Calculate Nyquist frequencies
  fnx=1/2/dx;
  fny=1/2/dy;
  % Prepare output
  varns={xfaks,yfaks,fnx,fny,dx,dy};
else
  % One-dimensional case
  dx=pl/(max(fsize)-1);

  % Calculate centered frequency axis (PW p 112)
  N=max(ftsize);
  intvectX=linspace(-floor((N+1)/2)+1,N-floor((N+1)/2),N);
  % In other words the above is equal to
  % -floor((N+1)/2)+1+[0:N-1]
  xfaks=intvectX/N/dx;
  
  % Calculate Nyquist frequencies
  fnx=1/2/dx;
  % Prepare output
  varns={xfaks,fnx,dx};
end

% Prepare output
varargout=varns(1:nargout);
