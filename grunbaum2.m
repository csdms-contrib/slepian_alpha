function [E,Vg,th,C,T,V]=grunbaum2(TH,L,m,nth,grd)
% [E,Vg,th,C,T,V]=GRUNBAUM2(TH,L,m,nth,grd)
%
% Eigenfunctions of the DOUBLE POLAR CAP concentration problem.
%
% Calculates the matrix the way Grunbaum et al. (1982) propose.
% Orders the eigenfunctions in decreasing order.
%
% INPUT:
%
% TH          Angular extent of the spherical cap, in degrees
% L           Bandwidth (maximum angular degree)
% m           Angular order of the required data window, -l<m<l
% nth         Number of points sampling between 0 and pi
%             If nth=0, E will be zero and no sign-correction performed
% grd         1 Colatitudes only; returns matrix E [default]
%             2 Colatitude/Longitude; returns cell E
%
% OUTPUT:
%
% E           Optimally concentrated tapers, expanded to space domain 
% Vg          The Grunbaum eigenvalues
% th          The colatitudes at which the functions are evaluated
% C           The spherical harmonics coefficients
% T           The "tridiagonal" matrix in its full form
% V           The eigenvalues given by integration over the patch
%
% EXAMPLE:
%
% Compare Grunbaum's approach with our numerical one
% 
% TH=ceil(rand(1)*20); L=ceil(rand(1)*20); m=round(rand*(L-1));
% [E,Vg,th,C,T,V]=grunbaum2(TH,L,m);
% [E2,V2,th,C2,jk1,jk2,K]=sdwcap2(TH,L,m);
% difer(abs(T*K-K*T))
% [E,Vg,th,C,T,V]=grunbaum2(TH,L+1,m);
% [E2,V2,th,C2,jk1,jk2,K]=sdwcap2(TH,L+1,m);
% difer(abs(T*K-K*T))
%
% See also SDWCAP, BOXCAP
%
% Last modified by fjsimons-at-alum.mit.edu, 08/20/2008

defval('TH',30)
defval('L',18)
defval('m',0)
defval('nth',720)
defval('grd',1)

if m>L | m<-L
  error('Order cannot exceed the bandwidth')
end

mor=m;
m=abs(m);

b=cos(TH/180*pi);

% Even/odd degrees
ele=[m:2:L];   % Always an even function, for m even or odd
elo=[m+1:2:L]; % Always an odd function,  for m even or odd
% Compare to other, ridiculous formalism in RIDICOMP

% On/off-diagonal terms
Glle=ondiago(ele,max(ele),m,b);
Gl2le=offdiago(ele,max(ele),m);
if ~isempty(elo)
  Gllo=ondiago(elo,max(elo),m,b);
  Gl2lo=offdiago(elo,max(elo),m);
end

% Construct the even-degree tridiagonal matrix
Te=tridiag(Gl2le(1:end-1),Glle,Gl2le(1:end-1));
if ~isempty(elo)
  % Construct the odd-degree tridiagonal matrix
  To=tridiag(Gl2lo(1:end-1),Gllo,Gl2lo(1:end-1));
end

% Even/odd degree vector
lodd=even(L-m+1);

if nargout>=5
  % Construct the full matrix (not strictly required)
  T=repmat(0,L-m+1,L-m+1);
  if ~isempty(elo)
    T(lodd,lodd)=To;
  end
  T(~lodd,~lodd)=Te;
end

% Perform diagonalization separately for the even/odd matrices
[Ce,Ve]=eig(Te);
if ~isempty(elo)
  [Co,Vo]=eig(To);
end

% Construct the full eigensolution matrix
C=repmat(0,L-m+1,L-m+1);
if ~isempty(elo)
  % Collect the even-degree solutions
  C(lodd,lodd)=Co;
end
% Collect the odd-degree solutions
C(~lodd,~lodd)=Ce;
% Collect all the Grunbaum eigenvalues
Vg=repmat(0,L-m+1,1);
if ~isempty(elo)
  Vg(lodd)=diag(Vo);
end
Vg(~lodd)=diag(Ve);

% Output needs to be resorted according to the eigenvalues,
% but this is hard, so rely on the Sturm-Liouville eigenvalues
if ~mod(length(Vg),2)
  % This adjustment, however, is necessary for the even cases
  % I kind of forget why this is so, but hey, it works.
  res=reshape(flipud(reshape(1:length(Vg),2,length(Vg)/2)),length(Vg),1);
  Vg=Vg(res);
  C=C(:,res);
end

% Check normalization and calculate the eigenvalues
% from a straightforward GL integration
[ngl1,ngl2,com,V]=orthocheck(C,Vg,TH/180*pi,m,2);

% Compute spatial functions, colatitudinal part only
if nth~=0
  % Zonal functions only 
  if m==0
    % Make spatial functions
    % This is SDW (2006) equation (5.7) combined 
    % with the sqrt(2-dom) of (5.9) already included!
    [E,th]=pl2th(C,nth,1);
    th=th*180/pi;
    nlon=2*nth-1;
  else
    % This is SDW (2006) equation (5.7) combined 
    % with the sqrt(2-dom) of (5.9) already included!
    [E,nlon,lat]=plm2th(C,nth,m,1);
    th=linspace(0,180,size(E,1));
  end
  % Make E start with a positive lobe and ajust C too
  % Take the first NONZERO sample! Not just a numbered sample!
  % This was the source of a very nasty bug
  for index=1:size(E,2)
    C(:,index)=C(:,index)*sign(indeks(E(~~E(:,index),index),1));
    E(:,index)=E(:,index)*sign(indeks(E(~~E(:,index),index),1));
  end
else
  E=0;
  th=0;
  nlon=0;
end

if nth~=0 & grd==2
  % Output on full grid; watch the sign of m
  if mor<=0
    EE=E; clear E
    for index=1:size(EE,2)
      E{index}=EE(:,index)*cos(m*linspace(0,2*pi,nlon));
    end
  end
  if mor>0
    EE=E; clear E
    for index=1:size(EE,2)
      E{index}=EE(:,index)*sin(m*linspace(0,2*pi,nlon));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the on-diagonal terms
function Gll=ondiago(el,L,m,b)
Gll=-el.*(el+1)*b^2+2./(2*el+3).*((el+1).^2-m^2)+...
    ((el-2).*(el+1)-L*(L+3)).*...
    (1/3-2/3*(3*m^2-el.*(el+1))./(2*el+3)./(2*el-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Gl2l=offdiago(el,L,m)
% Computes the off-diagonal terms
Gl2l=(el.*(el+3)-L*(L+3))./(2*el+3).*...
     sqrt(((el+2).^2-m^2).*((el+1).^2-m^2)./(2*el+5)./(2*el+1));

