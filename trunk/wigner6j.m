function [w6j,L,norma]=wigner6j(l2,l3,l4,l5,l6)
% [w6j,L,norma]=WIGNER6J(l2,l3,l4,l5,l6)
%
% Calculates Wigner 6j symbols by recursion, for all values up to L
% that are allowed in the expression (L  l2 l3)
%                                    (l4 l5 l6)
%  There is no truncation at any bandwidth - they are all returned.
%
% INPUT:
%
% l2, l3, l4, l5, l6   The quantities defined above, all scalar
%
% OUTPUT:
%
% w6j                  The Wigner 6j symbols (with prepended zeroes)
% L                    The degrees at which these are evaluated
% norma                Normalization check which passed muster
%
% EXAMPLES:
%
% wigner6j('demo1') % Just tests some values against Stone's code
% wigner6j('demo2') % Reproduces Table IV and Figure 3 of Schulten and Gordon
% wigner6j('demo3') % Reproduces Table III of Schulten and Gordon
% wigner6j('demo4') % Reproduces Figure 4 of Schulten and Gordon
% 
% COMPARE WITH:
%
% Fortran code WIGNER6J.F based on original by Schulten and Gordon
% g77 -o $PFILES/wigner6j $FFILES/wigner6j.f $FFILES/drc6j.f ; wigner6j
% sprintf('%23.16e\n',wigner6j(0,8,7,6.5,7.5,7.5))
%
% Last modified by fjsimons-at-alum.mit.edu, 1/11/2007

if ~isstr(l2)
  % Overflow parameters
  %huge=sqrt(1.79D+308/20);
  huge=sqrt(realmax/20);
  sqhuge=sqrt(huge);
  tiny=1/huge;
  srtiny=1/sqhuge;

  % Limits for L1
  l1min=max(abs(l2-l3),abs(l5-l6));
  l1max=min(l2+l3,l5+l6);  
  
  %  Check error condition 4.
  if mod(l1max-l1min,1)
    error('l1max-l1min not integer')
  end

  % disp(sprintf('Computation from %i to %i',l1min,l1max))
  % Check error conditions 1, 2, and 3.
  if mod(l2+l3+l5+l6,1) || mod(l4+l2+l6,1)
    % warning('l2+l3+l5+l6 or l4+l2+l6 not integer')
    w6j=repmat(0,l1max+1,1)';
    L=0:l1max; return
  elseif  l4+l2-l6<0 || l4-l2+l6<0 || -l4+l2+l6<0
    % warning('l2, l4, l6 triangular condition not satisfied')
    w6j=repmat(0,l1max+1,1)';
    L=0:l1max; return
  elseif  l4-l5+l3<0 || l4+l5-l3<0 || -l4+l5+l3<0
    % warning('l3, l4, l5 triangular condition not satisfied')
    w6j=repmat(0,l1max+1,1)';
    L=0:l1max; return
  end

  % Compute single coefficient directly by analytic formula
  if l1min<l1max
    nfin=floor(l1max-l1min+1);
  elseif l1min==l1max
    w6j=repmat(0,l1max+1,1)';
    %  This is reached in case that L1 can take only one value
    w6j(l1max+1)=(-1)^floor(l2+l3+l5+l6)/sqrt((2*l1min+1)*(2*l4+1));
    % Check normalization
    norma=(2*l1min+1)*(2*l4+1)*w6j(l1max+1)^2;
    difer(norma-1)
    L=0:l1max;
    return
  else
    % Check error condition 5
    error('l1min greater than l1max')
  end

  %  Start of forward recursion
  l1=l1min;
  newfac=0;
  c1=0;
  w6j(1)=srtiny;
  sum1=(2*l1+1)*tiny;

  lstep=1;
  while lstep~=nfin % WHILE LOOP 30
    lstep=lstep+1;
    l1=l1+1;
    oldfac=newfac;
    a1=(l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)*(-l1+l2+l3+1);
    a2=(l1+l5+l6+1)*(l1-l5+l6)*(l1+l5-l6)*(-l1+l5+l6+1);
    newfac=sqrt(a1*a2);

    if l1>1
      dv=2*(l2*(l2+1)*l5*(l5+1)+l3*(l3+1)*l6*(l6+1)-l1*(l1-1)*l4*(l4+1))...
	 -(l2*(l2+1)+l3*(l3+1)-l1*(l1-1))*(l5*(l5+1)+l6*(l6+1)-l1*(l1-1));
      denom=(l1-1)*newfac;
      if ~(lstep-2<=0)
	c1old=abs(c1);
      end
      c1=-(2*l1-1)*dv/denom;
    else
      %  If L1=1, (L1-1) has to be factored out of DV, hence
      c1=-2*(l2*(l2+1)+l5*(l5+1)-l4*(l4+1))/newfac;
    end

    if ~(lstep>2)
      % If L1=L1MIN+1, the third term in recursion equation vanishes
      x=srtiny*c1;
      w6j(2)=x;
      sum1=sum1+tiny*(2*l1+1)*c1^2;
      if lstep==nfin
	sumuni=sum1;
	%  Normalize 6j coefficients
	cnorm=1/sqrt((2*l4+1)*sumuni);
	%  Sign convention for last 6j coefficient determines overall phase
	sign1=sign(w6j(nfin));
	sign2=(-1)^floor(l2+l3+l5+l6);
	if sign1*sign2<=0
	  cnorm=-cnorm;
	end
	if abs(cnorm)<1
	  thresh=tiny/abs(cnorm);
	  for n=1:nfin
	    if abs(w6j(n))<thresh
	      w6j(n)=0;
	    end
	    w6j(n)=cnorm*w6j(n);
	  end
	else
	  for n=1:nfin
	    w6j(n)=cnorm*w6j(n);
	  end
	end

	% Check normalization
	norma=sum((2*l4+1)*(2*(l1min:l1max)+1).*w6j.^2);
	difer(norma-1)
	% Provide output
	w6j=[repmat(0,1,l1min) w6j];
	L=0:l1max;
	return
      end
    else
      c2=-l1*oldfac/denom;
      %  Recursion to the next 6j coefficient X
      x=c1*w6j(lstep-1)+c2*w6j(lstep-2);
      w6j(lstep)=x;
      
      sumfor=sum1;
      sum1=sum1+(2*l1+1)*x*x;
      if lstep==nfin
	break
      end
      if abs(x)<sqhuge
	%  As long as the coefficient ABS(C1) is decreasing, the recursion
	%  proceeds towards increasing 6j values and, hence, is numerically
	%  stable.  Once an increase of ABS(C1) is detected, the recursion
	%  direction is reversed.
	if c1old-abs(c1)<=0
	  break
	end
      else
	% %  This is reached if last 6j coefficient larger than sqhuge,
	% %  so that the recursion series W6J(1), ... ,W6J(LSTEP)
	% %  has to be rescaled to prevent overflow
	for i=1:lstep
	  if abs(w6j(i))<srtiny
	    w6j(i)=0;
	  end
	  w6j(i)=w6j(i)/sqhuge;
	end
	sum1=sum1/huge;
	sumfor=sumfor/huge;
	x=x/sqhuge;
      end
    end
  end % WHILE LOOP 30
  
  % WHILE BREAK 100
  % %  Keep three 6j coefficients around LMATCH for comparison later
  % %  with backward recursion.
  x1=x;
  x2=w6j(lstep-1);
  x3=w6j(lstep-2);

  % %  Starting backward recursion from L1MAX taking NSTEP2 steps, so
  % %  that forward and backward recursion overlap at the three points
  nfinp1=nfin+1;
  nfinp2=nfin+2;
  nfinp3=nfin+3;
  nstep2=nfin-lstep+3;
  l1=l1max;

  w6j(nfin)=srtiny;
  sum2=(2*l1+1)*tiny;
  l1=l1+2;

  % WHILE LOOP LABEL 110
  lstep=1;
  while lstep~=nstep2
    lstep=lstep+1;
    l1=l1-1;
    oldfac=newfac;
    a1s=(l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)*(-l1+l2+l3+2);
    a2s=(l1+l5+l6)*(l1-l5+l6-1)*(l1+l5-l6-1)*(-l1+l5+l6+2);
    newfac=sqrt(a1s*a2s);
    dv=2*(l2*(l2+1)*l5*(l5+1)+l3*(l3+1)*l6*(l6+1)-l1*(l1-1)*l4*(l4+1))...
       -(l2*(l2+1)+l3*(l3+1)-l1*(l1-1))*(l5*(l5+1)+l6*(l6+1)-l1*(l1-1));
    denom=l1*newfac;
    
    c1=-(2*l1-1)*dv/denom;
    if ~(lstep>2)
      y=srtiny*c1;
      w6j(nfin-1)=y;
      if lstep==nstep2
	break
      end
      sumbac=sum2;
      sum2=sum2+(2*l1-3)*c1^2*tiny;
    else
      c2=-(l1-1)*oldfac/denom;
      y=c1*w6j(nfinp2-lstep)+c2*w6j(nfinp3-lstep);
      if lstep==nstep2
	break
      else
	w6j(nfinp1-lstep)=y;
	sumbac=sum2;
	sum2=sum2+(2*l1-3)*y*y;
	if ~(abs(y)<sqhuge)   
	  for i=1:lstep
	    index=nfin-i+1
	    if abs(w6j(index))<srtiny
	      w6j(index)=0;
	    end
	    w6j(index)=w6j(index)/sqhuge;
	  end
	  sumbac=sumbac/huge;
	  sum2=sum2/huge;
	end
      end
    end
  end

  % WHILE BREAK 200
  y3=y;
  y2=w6j(nfinp2-lstep);
  y1=w6j(nfinp3-lstep);

  % Determine now RATIO such that YI=RATIO*XI  (I=1,2,3) holds
  % with minimal error.
  warning off
  ratio=(x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3);
  warning on
  nlim=nfin-nstep2+1;

  if abs(ratio)<1
    nlim=nlim+1;
    ratio=1/ratio;
    for n=nlim:nfin
      w6j(n)=ratio*w6j(n);
    end
    sumuni=sumfor+ratio^2*sumbac;
  else
    for n=1:nlim
      w6j(n)=ratio*w6j(n);
    end
    sumuni=ratio^2*sumfor+sumbac;
  end

  %  Normalize 6j coefficients
  cnorm=1/sqrt((2*l4+1)*sumuni);
  %  Sign convention for last 6j coefficient determines overall phase
  sign1=sign(w6j(nfin));
  sign2=(-1)^floor(l2+l3+l5+l6);
  if sign1*sign2<=0
    cnorm=-cnorm;
  end
  if abs(cnorm)<1
    thresh=tiny/abs(cnorm);
    for n=1:nfin
      if abs(w6j(n))<thresh
	w6j(n)=0;
      end
      w6j(n)=cnorm*w6j(n);
    end
  else
    for n=1:nfin
      w6j(n)=cnorm*w6j(n);
    end
  end
  % Check normalization
  norma=sum((2*l4+1)*(2*(l1min:l1max)+1).*w6j.^2);
  difer(norma-1)
  
  % Provide output:
  % Provide the right number of zeroes
  w6j=[repmat(0,1,l1min) w6j];
  L=0:l1max;
elseif strcmp('demo1',l2)
  % Useing Anthony Stones' RFF calculator for most of these
  difer(indeks(wigner6j(1,10,0,10,1),10:12)-...
	[sqrt(1/7)/3 -sqrt(1/7)/3 sqrt(1/7)/3])
  difer(wigner6j(2,2,3,2,3)-[0 sqrt(3/7)/5 -sqrt(3/2)/35 -1/10 -sqrt(1/3)/14])
  difer(wigner6j(1,1,2,1,2)-[0 -sqrt(1/5)/2 -1/10])
  difer(wigner6j(1,2,1,2,2)-[0 -1/10 sqrt(7/3)/10 -sqrt(2/3)/5])
  difer(wigner6j(1,1,1,0,0)-sqrt(1/3))
  difer(indeks(wigner6j(0,1,1,0,1),2)-1/3)
  difer(indeks(wigner6j(2,1,1,0,1),2)-1/3)
  difer(indeks(wigner6j(0,1,1,0,1),2)-1/3)
  difer(wigner6j(2,2,2,0,0)-sqrt(1/5))
  difer(indeks(wigner6j(3,2,2,0,1),2)-sqrt(1/15))
  difer(indeks(wigner6j(4,2,2,0,2),3)-1/5)
  difer(indeks(wigner6j(1,3,3,0,2),3)-sqrt(1/35))
  difer(indeks(wigner6j(3,3,3,0,2),3)-sqrt(1/35))
  difer(indeks(wigner6j(5,3,3,0,2),3)-sqrt(1/35))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),2)-(7/680)*sqrt(23/2))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),3)--(13/680)*sqrt(23/6))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),4)-(1147/167960)*sqrt(23/3))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),5)-(23/33592)*sqrt(115))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),6)--(181/33592)*sqrt(115/6))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),7)-(4343/167960)*sqrt(23/42))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),8)-(1941/83980)*sqrt(1/322))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),9)--(10993/83980)*sqrt(1/46))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),10)-(1709/33592)*sqrt(5/46))  
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),11)-(13/2584)*sqrt(55/46))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),12)--(691/12920)*sqrt(11/69))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),13)-(103/12920)*sqrt(13/69)) 
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),14)-(1621/12920)*sqrt(13/322)) 
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),15)-(533/2584)*sqrt(5/966))
  difer(indeks(wigner6j(8,7,13/2,15/2,15/2),16)-(13/1292)*sqrt(5/69))
  % An example from Rotenberg
  difer(indeks(wigner6j(8,8,8,8,8),9)-...
	rotenberg([-2 2 -2 0 0 -2 -2 -2 -2],-1)*4073)
elseif strcmp('demo2',l2)
  [w,L]=wigner6j(48,80,112,120,72);
  disp(sprintf('------------------------------'))
  disp(sprintf('L1  Values of 6j coefficients'))
  disp(sprintf('------------------------------'))
  disp(sprintf('%2.2i %23.16i\n',[L(1:8:end)' w(1:8:end)']'))
  disp(sprintf('------------------------------'))
  clf
  w(1:5)=w(1:5)*1e3;
  w=[w(1:5) NaN w(6:end)];
  L=[L(1:5) NaN L(6:end)];
  plot(L,w*1e4,'x-')
  ylim([-12 16.5])
  set(gca,'ytick',-12:2:16)
  xlim([45 130])
  shrink(gca,2,2)
  grid on
  figdisp('Schulten+75-fig3',[],[],1)
elseif strcmp('demo3',l2)
  [w,L]=wigner6j(8,7,13/2,15/2,15/2);
  disp(sprintf('------------------------------'))
  disp(sprintf('L1  Values of 6j coefficients'))
  disp(sprintf('------------------------------'))
  disp(sprintf('%2.2i %23.16i\n',[L(1:9)' w(1:9)']'))
  disp(sprintf('------------------------------'))
  plot(L,w*1e4,'+-')
elseif strcmp('demo4',l2)
  clf
  subplot(311)
  [w,L]=wigner6j(48,80,120,120,72);
  plot(L,w*1e4,'x-')
  ylim([0 16.5])
  axis tight
  set(gca,'ytick',0:2:16)
  shrink(gca,2,1)
  grid on
  subplot(312)
  [w,L]=wigner6j(48,80,115,120,72);
  plot(L,w*1e4,'x-')
  ylim([-12 16.5])
  axis tight
  set(gca,'ytick',-12:2:16)
  shrink(gca,2,1)
  grid on
  subplot(313)
  [w,L]=wigner6j(48,80,110,120,72);
  plot(L,w*1e4,'x-')
  ylim([-12 14])
  axis tight
  set(gca,'ytick',-12:2:14)
  shrink(gca,2,1)
  grid on
  figdisp('Schulten+75-fig4',[],[],1)
end

