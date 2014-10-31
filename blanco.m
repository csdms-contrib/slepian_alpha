function varargout=blanco(L,bta)
% [DR,DC]=BLANCO(L,bta)
%
% Constructs spherical harmonics rotation matrices for use with spherical
% harmonics coefficients, using the algorithm of Blanco, Florez & Bermejo (1997).
%
% INPUT:
% 
% L        Maximum degree of the harmonics
% alp      Euler angle in degrees
% bta      Euler angle in degrees [default: 90]
% gam      Euler angle in degrees
%
% OUTPUT:
%
% DR       Rotation matrix for m>=0 only (lower right quarter)
% DC       Rotation matrix for m=-l:l (full matrix)
%
% EXAMPLE:
%
% Compare the difference; it works for the decomposition where it's just
% ninety 
%
% L=round(rand*180)
% d1=blanco(L,90); d2=dlmb(L); 
% for index=1:L+1; difer(d1{index}-d2{index}); end
%
% See also DLMB, PLM2ROT
%
% Last modified by fjsimons-at-alum.mit.edu, January 12th, 2003

t0=clock;
defval('bta',90)
bta=bta*pi/180;
if bta<0
  warning('Only for positive rotations')
end

% Eq. (48-52)
d{1}(1,1)=1;
d{2}(2,2)=cos(bta);
d{2}(3,1)=sin(bta/2)^2;
d{2}(3,2)=-1/sqrt(2)*sin(bta);
d{2}(3,3)=cos(bta/2)^2;
% Avoid round-off error
if abs(bta-pi/2)<eps
  d{2}(2,2)=0;
  d{2}(3,1)=1/2;
  d{2}(3,3)=1/2;
end
if abs(bta-pi)<eps
  d{2}(3,2)=0;
  d{2}(3,3)=0;
end

% Use l, m, and mp as degrees and
% li, mi, mpi as indices of d{l}
% li1, mi1, mpi1 as indices of d{l-1}
% li2, mi2, mpi2 as indices of d{l-2}
% Loop over degrees l
for l=2:L
  % Start with zero to initialize
  li=l+1;
  d{li}=repmat(0,2*l+1);
  li1=li-1;
  li2=li-2;
  % Use Eq. 64
  for m=0:l-2
    % Index of m for l
    mi =m+l+1;
    % Index of m for l-1
    mi1=m+l;
    % Index of m for l-2
    mi2=m+l-1;
    for mp=-m:m
      mpi =mp+l+1;
      mpi1=mp+l;
      mpi2=mp+l-1;
      fac1=l*(2*l-1)/sqrt((l^2-m^2)*(l^2-mp^2));
      fac2=d{2}(2,2)-m*mp/l/(l-1);
      fac3=sqrt(((l-1)^2-m^2)*((l-1)^2-mp^2))/(l-1)/(2*l-1);
      d{li}(mi,mpi)=fac1*(fac2*d{li1}(mi1,mpi1)...
			  -fac3*d{li2}(mi2,mpi2));
    end
  end
  % Eq. 65
  % Last index of current l, double
  m=l;  mi=m+l+1;
  % Last index of previous l, double
  mi11=m+l-1; 
  d{li}(mi,mi)=d{2}(3,3)*d{li1}(mi11,mi11);
  % Eq. 66
  % One but last index of current l
  m=l-1; mi=m+l+1; mi1=m+l;
  d{li}(mi,mi)=(l*d{2}(2,2)-l+1)*d{li1}(mi11,mi11); 
  % Eq. 67
  % Last row of the matrix except its last column
  m=l;
  mi=m+l+1;
  for mp=l-1:-1:-l
    mpi=mp+l+1;
    warning off
    fac1=-sqrt((l+mp+1)/(l-mp-1+1))...
	 *sqrt(d{2}(3,1)/d{2}(3,3));
    warning on
    if isinf(sqrt(d{2}(3,1)/d{2}(3,3)))
      if mp==-l
	d{li}(mi,mpi)=1;
      else
	d{li}(mi,mpi)=0;
      end
    else
      d{li}(mi,mpi)=fac1*d{li}(mi,mpi+1);
    end
    if isinf(fac1)
      d{li}(mi,mpi)=0;     
    end
    if isnan(fac1)
      d{li}(mi,mpi)=0;
    end
  end
  % Eq. 68
  % One but last row
  m=l-1;
  mi=m+l+1;
  for mp=(l-2):-1:(1-l)
    mpi=mp+l+1;
    warning off
    fac1=-(l*d{2}(2,2)-mp-1+1)/(l*d{2}(2,2)-mp-1)*...
	 sqrt((l+mp+1)/(l-mp-1+1))...
	 *sqrt(d{2}(3,1)/d{2}(3,3));
    warning on
    if isinf(sqrt(d{2}(3,1)/d{2}(3,3)))
      if mp==1-l
	d{li}(mi,mpi)=l-1;
      else
	d{li}(mi,mpi)=0;
      end
    else
      d{li}(mi,mpi)=fac1*d{li}(mi,mpi+1);
    end
    if isinf(fac1)
      d{li}(mi,mpi)=0;     
    end
    if isnan(fac1)
      d{li}(mi,mpi)=0;
    end
  end
end
% Something is wrong - either the l-1 elements are not good
% or we're not using enough of the symmetry properties.
% Symmetrize using mirror properties of Eq. 38
for l=1:L
  li=l+1;
  % First symmetry equation: flip across LLUR diagonal
  d{li}=(d{li}+flipud(fliplr(d{li}')))./(flipud(eye(2*l+1))+1);
  % Second symmetry relation
  mmp=repmat(-l:l,2*l+1,1)+repmat([-l:l]',1,2*l+1);
  d{li}=(d{li}+(-1).^mmp.*d{li}')./(eye(2*l+1)+1);
end

% Collect only lower right quarter 
for l=0:L
  li=l+1;
  dr{li}=d{li}(li:end,li:end);
end

varn={'dr','d'};
for index=1:nargout
  varargout{index}=eval(varn{index});
end
disp(sprintf('BLANCO took %8.4f s',etime(clock,t0)))
