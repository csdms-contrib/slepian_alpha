function [D,d]=dlmb(L)
% [D,d]=DLMB(L)
%
% INPUT: 
%
% L          Maximum angular degree
%
% OUTPUT:
%
% D          Lower right quarter of the Wigner D-matrix, for m>=0
% d          Masters' concatenated output 
%
% Computes matrix elements for spherical harmonic polar rotation around
% the y-axis over 90 degrees only. This splits the rotation into two
% parts, both of which only contain a constant -pi/2 polar rotation; the
% others are azimuthal rotations.
% 
% We need the rotation matrix 
%   D_{mm'}(a,b,g)=exp(-ima)d_{mm'}(b)exp(-im'g)
% but we factor the rotation itself into:
%    R(a,b,g)=R(a-pi/2,-pi/2,b)R(0,pi/2,g+pi/2)
% thus we only need to compute d_{mm'} for b=90.
%
% See also BLANCO, PLM2ROT
%
% After a code by T. Guy Masters.
% See also McEwen, 2006.
% Last modified by fjsimons-at-alum.mit.edu, 08/05/2008

fname=fullfile(getenv('IFILES'),'DLMB',sprintf('dlmb-%3.3i.mat',L));

% Need to figure out what the size of d will be
d=zeros(1,sum(([0:L]+1).^2));

if exist(fname,'file')==2
  load(fname)
else
  t0=clock;
  % Initialize using D&T C.115.
  % For l=0
  d(1)=1; 
  % For l=1 % Hold on, not necessary for lower ones
  if L>=1
    d(2)=0;  
    d(3)=1/sqrt(2);
    d(4)=-1/sqrt(2);
    d(5)=1/2;
  end
  ind=5;
  f1=1/2;
  Lwait=100;
  if L>Lwait
    h=waitbar(0,'Looping over the degrees');
  end
  for l=2:L
    lp1=l+1;
    knd=ind+lp1;
    fl2p1=l+lp1;
%      for i=1:l
%        f(i)=sqrt(i*(fl2p1-i));
%      end
    f=sqrt([1:l].*(fl2p1-[1:l]));

    f1=f1*(2*l-1)/(2*l);
    % For N=0
    d(knd)=-sqrt(f1);
    d(knd-1)=0;
    for i=2:l
      j=knd-i;
      d(j)=-f(i-1)*d(j+2)/f(i);
    end
    % Positive N (bottom triangle)
    f2=f1;
    g1=l;
    g2=lp1;
    for N=1:l
      knd=knd+lp1;
      en2=N+N;
      g1=g1+1.;
      g2=g2-1.;
      f2=f2*g2/g1;
      d(knd)=-sqrt(f2);
      d(knd-1)=d(knd)*en2/f(1);
      for i=2:l-N
	j=knd-i;
	d(j)=(en2*d(j+1)-f(i-1)*d(j+2))/f(i);
      end
    end
    % Fill upper triangle and fix signs
    for j=1:l
      for m=j:l
	d(ind+m*lp1+j)=d(ind+j*lp1+m-l);
      end
    end
    isn=1+mod(l,2);
    for n=0:l
      knd=ind+n*lp1; 
      for i=isn:2:lp1
	d(knd+i)=-d(knd+i);
      end
    end
    ind=ind+lp1*lp1;
    if L>Lwait
      waitbar(l/L,h)
    end
  end
  if L>Lwait
    close(h)
  end
  d=d(:);
  
  % Now let's rearrange the coefficients as 1x1, 2x2, 3x3 etc rotation
  % matrices.
  cst=1;
  for l=1:(L+1)
    % Start of coefficient sequence; need transpose!
    D{l}=reshape(d(cst:cst+l^2-1),l,l)';
    cst=cst+l^2;
  end
  disp(sprintf('DLMB took %8.4f s',etime(clock,t0)))
  eval(sprintf('save %s D d L',fname))
end

  
