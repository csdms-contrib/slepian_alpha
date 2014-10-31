function [K,kx,ky,dci,dcn,dx,dy]=knum2(mn,pl)
% [K,kx,ky,dci,dcn,dx,dy]=knum2([M N],[lY lX])
%
% Wavenumber axis of a 2D Fourier transform, suitable after FFTSHIFT and
% for wavenumber manipulation to ensure Hermitian symmetry with IFFT. The
% wavenumbers increase symmetrically away from zero, but the Nyquist
% frequency is only reached on the top and/or left parts of the
% dimensions that have an even number of samples.
%
% INPUT:
% 
% M,N       Sample size of FFT matrix, equal to size of original matrix 
%           [M   N] [rows columns]. The need to be larger than one, obviously.
% lY,lX     Physical size of original matrix
%           [lY lX] [down across]
%
% OUTPUT:
%
% K         The wavenumber matrix (the norm of the wave vectors)
% kx,ky     The components of the wave vectors
% dci       The (m,n) indices to the DC component [at floor(dim/2)+1]
% dcn       The (m,n) indices to the components at exactly the Nyquist
% dx,dy     The sampling interval in the space domain
%
% The "corner Nyquist" rate is (1,1) and the DC rate is the "center
% point", dci, which is in the LR quadrant for even lengths; the way
% IFFT/FFTSHIFT expect it. Note that this is not for geographical data!
% The negative x and y frequencies are in the TOP LEFT corner -> FLIPUD.
%
% EXAMPLE:
%
% N=60+round(rand); h=peaks(N); H=fftshift(fft2(h)); F=randn*10
% hf=ifft2(ifftshift(H.*exp(-F*knum2([N N],[100 100])))); imagef(hf)
% 
% EXAMPLE:
%
% H=peaks(200); Hk=fftshift(fft2(H));
% [K,kx,ky,dci,dcn,dx,dy]=knum2(size(Hk),[1000 1000]);
%% Where and what is the "DC" component?
% difer(sum(H(:))-Hk(dci(1),dci(2)))
%% How is Parseval's theorem satisfied?
% difer(H(:)'*H(:)-Hk(:)'*Hk(:)/prod(size(Hk)),8)
%
% SEE ALSO: FFTAXIS, FFTAXIS1D, KNUM, RANDGPN, BRACEWELL
%
% Last modified by fjsimons-at-alum.mit.edu, 04/08/2014

% Supply defaults for testing
defval('mn',[32 32]+round(rand))
defval('pl',mn-1)

M=mn(1);
N=mn(2);

ly=pl(1);
lx=pl(2);

dx=lx/(N-1);
dy=ly/(M-1);

kx=2*pi*linspace(-floor(N/2),floor((N-1)/2),N)/N/dx;
ky=2*pi*linspace(-floor(M/2),floor((M-1)/2),M)/M/dy;

[KX,KY]=meshgrid(kx(:),ky(:));
K=sqrt(KX.^2+KY.^2);

% Check that the result is good for Hermitian - see the example perhaps 

% The index of the zero frequency, guaranteed
dci=[floor(M/2)+1 floor(N/2)+1];
difer(K(dci(1),dci(2)),[],[],NaN)

% The indices where the Nyquist is exactly reached, which is only the
% case for one side of dimensions that are even.
dcn=[1 dci(2) ; dci(1) 1];
isev=find(~mod(mn,2));
for ind=isev
  difer(K(dcn(ind,1),dcn(ind,2))/2/pi*pl(ind)/(mn(ind)-1)-1/2,[],[],NaN)
end
dcn=dcn(isev,:);

% There may be others that are close... by chance, let's ignore them
%[i,j]=find(abs(abs(k/2/pi/what*what)-1/2)<eps)

% This is only subtly different from FFTAXIS and can be obtained from
% it. Just fliplr the kx and ky axes from fftaxis. See notes. (Not yet in
% Notebook.)
