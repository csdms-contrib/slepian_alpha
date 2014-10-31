function G=guseinov(l1,l2,l3,m1,m2,m3,wat)
% G=guseinov(l1,l2,l3,m1,m2,m3,wat)
%
% Direct calculation of Clebsch-Gordan and Gaunt coefficients.
% With funny phase convention; use for absolute values only.
%
% INPUT:
%
% wat       'gaunt' Calculates the Gaunt coefficient [default]
%           'clebsch' Calculates the Clebsch-Gordan coefficients
%
% SEE ALSO: GAUNT, WIGNER3JM
%
% EXAMPLES:
%
% guseinov('demo1') % Table 1 (left) in Guseinov, compared to WIGNER3JM
% guseinov('demo2') % Table 1 (right) in Guseinov, compared to WIGNER3JM 
% guseinov('demo3') % Table 2 in Guseinov, compared to WIGNER3JM 
% guseinov('demo4') % Table 1 in Sebilleau, compared to WIGNER3JM 
%
% Last modified by fjsimons-at-alum.mit.edu, 31.07.2006

% Using the method proposed by Guseinov et al. (1995a,b, 2005).
% See also Sebilleau (1998) and Xu (1996), among others.
% Switching all the bottom signs should leave the value intact...
% These computations are not very accurate; we know that some work, some
% really don't. Changing the sign, really, is a good check; if the values are
% very different, you're in trouble. Sometimes the one with the flipped sign
% is the one you need to reproduce the table more accurately.

% warning('Only QUADRUPLE precision will do this justice.')

% Gain speed by using recurrence relations for binomials, not NCHOOSEK

if ~isstr(l1)

  defval('wat','gaunt')

  warning off
  switch wat
   case 'clebsch'
    defval('m3',m1+m2)
    % Calculate the Clebsch-Gordan coefficient
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if m3~=[m1+m2] | ~triangle(l1,l2,l3) | abs(m3)>l3
      G=0; return ; end
      f2=(2*l3+1)^2/(2*l1+1)/(2*l2+1)*...
	 [nchoosek(l1+l2+l3+1,l1+l2-l3)/nchoosek(l1+l2+l3+1,l1-l2+l3)]*...
	 [nchoosek(2*l3,l3+m3)/nchoosek(l1+l2+l3+1,l2-l1+l3)]/...
	 nchoosek(2*l1,l1+m1)/nchoosek(2*l2,l2+m2);
      tmin=max([0 l2+m2-(l3+m3) l1-m1-(l3-m3)]);
      tmax=min([l1+l2-l3 l2+m2 l1-m1]);
      if tmax>=tmin; tt=tmin:tmax; else tt=tmin; end
      disp(sprintf('%i terms; m1=%i  m2=%i  m3=%i',length(tt),m1,m2,m3))
      f3=0;
      % Avoid alternating the signs doesn't help
      % for t=[tt(even(tt)) tt(~even(tt))]
      for t=tt
	f3=f3+(-1)^t*nchoosek(l1+l2-l3,t)*nchoosek(l3+m3,l2+m2-t)*...
	   nchoosek(l3-m3,l1-m1-t);
      end
      % Don't bother trying to get the right sign
      G=abs(sqrt(f2)*f3);
   case 'gaunt'
    defval('m3',m1-m2)
    % Calculate the Gaunt Coefficient
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if m3~=[m1-m2] | ~triangle(l1,l2,l3) | abs(m3)>l3 ...
	   sum(m1+m2+m3)==0 & mod(sum(l1+l2+l3),2); 
      G=0; return; end
      
      g=(l1+l2+l3)/2;
      p=(-1)^(g-(l2-m2)+(abs(m1)+abs(m2)+abs(m3))/2);
      f1=[nchoosek(2*g-l1-l2,g-l1)/nchoosek(2*g,2*l3)]*...
	 [nchoosek(g,l3)/(2*g+1)];
      f2=(2*l1+1)*(2*l2+1)*...
	 [nchoosek(l1+l2+m3,l1+m1)/nchoosek(l1+l2-m3,l1-m1)]*...
	 [nchoosek(2*l3+l1+l2+m3,l1+l2+m3)/...
	  nchoosek(l1+l2+2*l3+m3,l1+l2-m3)]/... 
	 nchoosek(2*l3,l3-m3)/nchoosek(2*l3+2*m3,l3+m3);
      tmin=max(0,l3-m1-l2);
      tmax=min([l1-abs(m1) l3-m3 l3-m1+l2]);
      if tmax>=tmin; tt=tmin:tmax; else tt=tmin; end
      disp(sprintf('%i terms; m1=%i  m2=%i  m3=%i',length(tt),m1,m2,m3))
      f3=0;
      % Avoid alternating the signs doesn't help
      % for t=[tt(even(tt)) tt(~even(tt))]
      for t=tt
	f3=f3+(-1)^t*nchoosek(l1+m1+t,t)*...
	   [nchoosek(l2-m2+l3-m3-t,l3-m3-t)/nchoosek(l1+l2+m3,l1+m1)]*...
	   nchoosek(l1+l2-l3,l1-m1-t)*nchoosek(l1+l2+m3,l1+l2-l3);
      end
      G=p*f1*sqrt(f2)*f3;
      
      % Note that this next factor was omitted from Guseinov's papers... 
      G=G*sqrt((2*l3+1)/4/pi);

      warning on
  end
elseif strcmp(l1,'demo1')
  % RELATION OF CLEBSCH-GORDAN TO WIGNER 3J-SYMBOLS:
  % Table I of Guseinov (1995), all in Guseinov phase
  
  L1=[20 15 20 20 25 35 25 25 40 37 40 40 60 58 60 60 80 77 80 80]; 
  L2=[15 20 9 15 35 25 40 35 37 40 49 37 58 60 7 58 77 80 6 77]; 
  L3=[9 9 15 9 40 40 35 40 49 49 37 49 7 7 58 7 6 6 77 6];
  M1=[-3 2 -3 3 12 -17 12 -12 -2 1 -2 2 3 2 3 -3 1 -3 1 -1];
  M2=[2 -3 1 -2 -17 12 5 17 1 -2 1 -1 2 3 -5 -2 -3 1 2 3];
  M3=M1+M2; clear comp
  for ind=1:length(L1)
    l1=L1(ind); l2=L2(ind); l3=L3(ind);
    m1=M1(ind); m2=M2(ind); m3=M3(ind);
    % COMPARE BOTH APPROACHES
    % In Guseinov's stupid phase convention
    G1=guseinov(l1,l2,l3,m1,m2,m3,'clebsch');
    % In Condon-Shortley phase convention
    G2=indeks(wigner3jm(l1,l2,l3,m1,m2,-m3),'end')...
       *sqrt(2*l3+1)*(-1)^(m3+l1-l2);
    % Every third in this table switches l3 and l2
    if ~mod(ind-3,4)
      G1=G1*sqrt(2*l2+1)/sqrt(2*l3+1);
      G2=G2*sqrt(2*l2+1)/sqrt(2*l3+1);
    end
    % And finally, compare
    comp(ind,:)=[abs(G1) abs(G2) abs(abs(G1)-abs(G2))];
  end
  comp
elseif strcmp(l1,'demo2')
  % RELATION OF CLEBSCH-GORDAN TO WIGNER 3J-SYMBOLS:
  % Table I of Guseinov (1995), all in Condon-Shortley phase

  L1=[20 15 20 20 25 35 25 25 40 37 40 40 60 58 60 60 80 77 80 80]; 
  L2=[15 20 34 15 35 25 60 35 37 40 77 37 58 60 116 58 77 80 78 77]; 
  L3=[34 34 15 34 60 60 35 60 77 77 37 77 116 116 58 116 78 78 77 78];
  M1=[-3 2 -3 3 12 -17 12 -12 -2 1 -2 2 3 2 3 -3 1 -3 1 -1];
  M2=[2 -3 1 -2 -17 12 5 17 1 -2 1 -1 2 3 -5 -2 -3 1 2 3];
  M3=M1+M2; clear comp
  for ind=1:length(L1)
    l1=L1(ind); l2=L2(ind); l3=L3(ind);
    m1=M1(ind); m2=M2(ind); m3=M3(ind);
    % In Guseinov's stupid phase convention
    G1=guseinov(l1,l2,l3,m1,m2,m3,'clebsch');
    % In Condon-Shortley phase convention
    G2=indeks(wigner3jm(l1,l2,l3,m1,m2,-m3),'end')...
       *sqrt(2*l3+1)*(-1)^(m3+l1-l2);
    % Every third in this table switches l3 and l2
    if ~mod(ind-3,4)
      G1=G1*sqrt(2*l2+1)/sqrt(2*l3+1);
      G2=G2*sqrt(2*l2+1)/sqrt(2*l3+1);
    end
    comp(ind,:)=[abs(G1) abs(G2) abs(abs(G1)-abs(G2))];
  end
  comp
elseif strcmp(l1,'demo3')
  L1=[20 20 20 20  25  25  25 40 40 40 40  60 60 38  2  80  80 80 80];
  L2=[15 15 17  9  35  35  35 37 37 21  5  58 58 58 58  77  77 77 77];
  L3=[35 31 15 15  60  48  38 75 59 37 37 118 58 60 60 155 131 83  5];
  M1=[-3  3 -3 -3  12  12  12  2  2 -2 -2   3  3  1  1   1   1  1  1];
  M2=[2  -2 -5 -5 -17 -17 -17 -1 -1 -3 -3   2  1 -2 -2  -3  -3 -3 -3];
  M3=M1-M2;
  clear comp
  for ind=1:length(L1)
    l1=L1(ind); l2=L2(ind); l3=L3(ind);
    m1=M1(ind); m2=M2(ind); m3=M3(ind);
    % In Guseinov's stupid phase convention
    G1=guseinov(l1,l2,l3,m1,m2,m3,'gaunt');
    % In Condon-Shortley convention
    G2=...
    indeks(wigner3jm(l1,l2,l3,m1,-m2,-m3),'end')*...
	sqrt(2*l1+1)*sqrt(2*l2+1)*sqrt(2*l3+1)/sqrt(4*pi)...
	*indeks(wigner3jm(l1,l2,l3,0,0,0),'end')*(-1)^m3;
    if abs(G1)>1e3; G1=NaN; end
    comp(ind,:)=[abs(G1) abs(G2) abs(abs(G1)-abs(G2))];
  end
  comp
elseif strcmp(l1,'demo4')
  L1=[10 12 20 29];
  L2=[10 15 20 29];
  L3=[12 5 40 34];
  M1=[ -9 -2 -1 -10];
  M2=[ 3   3 -1  -5];
  M3=[-12 -5  0  -5];
  clear comp
  for ind=1:length(L1)
    l1=L1(ind); l2=L2(ind); l3=L3(ind);
    m1=M1(ind); m2=M2(ind); m3=M3(ind);
    % In Guseinov's stupid phase convention; flip signs for stability 
    G1=guseinov(l1,l2,l3,-m1,-m2,-m3,'gaunt');
    % In Condon-Shortley convention; flip signs for the fun of it
    G2=...
    indeks(wigner3jm(l1,l2,l3,-m1,m2,m3),'end')*...
	sqrt(2*l1+1)*sqrt(2*l2+1)*sqrt(2*l3+1)/sqrt(4*pi)...
	*indeks(wigner3jm(l1,l2,l3,0,0,0),'end')*(-1)^m1;
    if abs(G1)>1e3; G1=NaN; end
    comp(ind,:)=[abs(G1) abs(G2) abs(abs(G1)-abs(G2))];
  end
  comp
end

