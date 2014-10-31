function [fax,selekt,lax]=fftaxis1D(sig,nfft,physl)
% [fax,selekt,lax]=fftaxis1D(sig,nfft,physl)
%
% Returns a vector of physical frequencies suitable to interpret the
% Fourier transform of a real or complex one-dimensional time-series.
% If the signal is real, the axis goes from 0 to <f_N, the Nyquist
% frequency, whereas if it is complex, it goes from 0 to <twice that (as
% an alternative to going from -f_N to f_N). The wavelengths go, not from
% the impractical infinity, but from the inverse of the Rayleigh
% frequency, i.e. the physical length of the signal, to <1/f_N.
%
% INPUT:
% 
% sig       A vector (time series) with signal values, real or complex
% nfft      The number of frequencies requested
% physl     The physical length scale of the entire signal, in UNIT
%
% OUTPUT:
%
% fax       Physical frequency axis, in UNIT^{-1}
% selekt    Index vector to pick out the non-redundant portions
% lax       Physical wavelength axis, in UNIT
%
% EXAMPLE:
%
% Use without FFTSHIFT, as in:
%
% FX=fft(S,nfft); SX=abs(FX).^2; plot(fax,SX(selekt)
%
% The "complete" axis is again, for a power-of-two signal
% [kax ; -flipud(kax(2:end-1))]
%
% The equivalence of FFTAXIS1D and KNUM2 is explicit:
% N=256;
% [fax,selekt]=fftaxis1d(rand(N,1),N,N-1);
% [K,kx]=knum2([2 N],[2 N]-1);
% fx=-fliplr(indeks(kx,selekt)/2/pi);
% difer(fx(:)-fax(:))
%
% SEE ALSO: KNUM2
%
% Last modified by fjsimons-at-alum.mit.edu, 10/22/2012

xsint=physl/(length(sig)-1);

if isreal(sig)
  selekt = [1:floor(nfft/2)+1];
else  
  selekt = 1:nfft; 
end
fax = (selekt - 1)'/xsint/nfft;
warning off MATLAB:divideByZero
lax=1./fax; lax(1)=physl;
warning on MATLAB:divideByZero

