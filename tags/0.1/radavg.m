function [knum,rav,thet,thav,nrad,naz]=radavg(A)
% [knum,rav,thet,thav]=RADAVG(A)
%
% Returns radial and azimuthal averages of a "wavenumber"-matrix
%
% INPUT: 
%
% A           A (possibly rectangular) matrix with Hermitian symmetry,
%             with assumed IDENTICAL Nyquist frequency in both dimensions
%
% OUTPUT:
%
% knum        Wavenumbers (0 to 1/2) for averages over all angles... not needed
% rav         The average of the matrix values at those dimensionless wavenumbers
% thet        Angles (between 0 and 2*pi) for the averages over all wavenumbers
% thav        The average of the matrix values at those angles
%
% EXAMPLE:
%
% [k,kx,ky,dci,dcn]=knum2([412 412],[8 8]*1e6); 
% th=[0.025 1.5 23000];
% Sk=maternos(k,th);
% Sk(k>max(abs(kx)))=NaN;
%% Average of the matrix values
% [~,rav,thet,thav]=radavg(v2s(Sk));
%% Average of the wavenumber values
% [knumii,knum]=radavg(k);
%% Alternative computation
% [~,avi]=isav(v2s(Sk),kx,ky);
% [knumi,knums]=isav(k,kx,ky);
%% Plot average value at the average wavenumber using RADAVG
% semilogx(knum,rav,'v'); hold on; 
% Plot the exact value at the average wavenumber using MATERNOS
% semilogx(knum,maternos(knum,th),'o')
% Plot the average value at the exact wavenumber 
% semilogx(knumi,avi,'+');
% xlabel('average wavenumbers') ; ylabel('average power spectral density')
% % plot(thet,thav,'o')
% % plot(knumi./knumii)
% % plot(knums./knum')
% % plot(knumi./knum')
% % plot(rav./avi')
% % plot(rav./maternos(knum,th))
%
% SEE ALSO:
%
% RADAVG_DW, ISAV
%
% Contributions by Oded Aharonson, 01/01/2000
% Last modified by fjsimons-at-mit.edu, 04/11/2014

% Make the proper connection with a two-dimensional wavevector scale
[ny,nx]=size(A);
[K,kx,ky,dci,dcn]=knum2([ny nx],[ny nx]-1);

% Both are between -1/2 and 1/2 as we have assumed a unit sampling interval
kx=kx/2/pi;
ky=ky/2/pi;

% Number of radial points and azimuths
nr=ceil(min(nx,ny)/2);
nth=min(nx,ny)*2;

% Prepare to where you want to interpolate 
knum=linspace(0,0.5,nr);
thet=linspace(0,2*pi,nth+1);
thet=thet(1:end-1);

% One column of radii per angle you want to consider
kxi=knum(:)*cos(thet(:)');
% One row of angles per radius you want to consider
kyi=knum(:)*sin(thet(:)');

% Interpolate
zi=interp2(kx,ky,A,kxi,kyi);

% Average
rav=nanmean(zi,2);
thav=nanmean(zi,1);

% Check that at the centerpoint, which is unique, you should get the same
% thing coming in as going out (if the center value was zero)
try
  difer([A(dci(1),dci(2))-rav(1)]/rav(1),[],[],NaN)
catch
  difer([A(dci(1),dci(2))-rav(1)],[],[],NaN)
end
