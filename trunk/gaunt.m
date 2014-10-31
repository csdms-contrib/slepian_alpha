function [G,l]=gaunt(l1,l2,l3,m1,m2,m3,meth,L,C0,S0,C3,S3)
% [G,l]=GAUNT(l1,l2,l3,m1,m2,m3,meth,L,C0,S0,C3,S3)
%
% Gaunt coefficients: the l1m1th expansion coefficient of the product
% Yl2m2 x Yl3m3 where Ylm is a Condon-Shortley complex spherical
% harmonic.
%
% INPUT:
%
% l1,l2,l3   Top row of the Wigner 3j symbol [may be vector for 'table']
% m1,m2,m3   Bottom row of the Wigner 3j symbol [may be vector for 'table']
% meth       'recursive' Recursive method and Wigner symbols (WIGNER3JM)
%            'gl'        Gauss-Legendre integration
%            'table'     Table look-up prestored from recursive (THREEJ, ZEROJ)
%            'guseinov'  Binomial method of Guseinov (no sign)
% L          Bandwidth of the table vectors (only if they are supplied next)
% C0,S0      The column/element vectors of the table for ZEROJ if available
% C3,S3      The column/element vectors of the table for THREEJ if available
% 
% OUTPUT:
%
% G       The scalar Gaunt coefficient
%                (for scalar input: for 'guseinov', 'gl' and 'table')
%         A vector of Gaunt coefficients 
%                (for scalar input: for 'recursive'; for vector inputs)
% l       The first degrees from 0 to the maximum allowed value
%
% EXAMPLE: 
%
% gaunt('demo1') % Generate Table 1 [Sebilleau, 1996]; compare approaches
% gaunt('demo2') % Generate Table 2 [Guseinov, 1995]; compare approaches
% gaunt('demo3') % Compare the recursive method with the prestored tables
% gaunt('demo4') % Compare the Gauss-Legendre method with the prestored tables
%
% Last modified by fjsimons-at-alum.mit.edu, 01/18/2007
 
% Evaluates the integral of the product of three complex fully
% orthonormalized spherical harmonics, with the Condon-Shortley phase... 
% l1,m1 belong to the harmonic whose complex conjugate appears in the
% integral... Thus this is the l1m1th expansion coefficient of the
% product Yl2m2 x Yl3m3 in the complex spherical harmonic basis Yl1m1.
%
% Note that an alternative formulation contains Wigner 3j-symbols with a
% bottom row, which must sum to one, of -m1 m2 m3 and that thus,
% m3=m1-m2. Flipping all the order signs leaves the overall sign intact.
% m1 regulates the sign by the factor (-1)^m1.

% SHOULD BUILD IN RULES BEFORE YOU EVEN START
% MUST CHECK! THREEJ and GAUNT DO NOT BEHAVE THE WAY I WANT THEM TO
% REGARDING REPEATED VALUES OF A DEGREE BUT NOT THE CORRESPONDING ORDERS!
% SHOULD REWRITE THIS TO BE NOTATIONALLY MORE UNIFORM

if ~isstr(l1)

  defval('meth','table')
  defval('L',[])
  defval('C0',[])
  defval('S0',[])
  defval('C3',[])
  defval('S3',[])
  
  %disp(sprintf('Using method %s',meth))
  
  switch meth
   case 'guseinov'
    if length(l1)>1||length(l2)>1||length(l3)>1; error('Need scalar input'); end
    G=abs(guseinov(l1,l2,l3,m1,m2,m3));
    Gp=abs(guseinov(l1,l2,l3,-m1,-m2,-m3));
    if abs(G-Gp)>1e-10
      warning('No stable result. Choose (n)either value.'); 
      G=[G Gp];
    end
    % Returns single value
    l=l1;
   case 'recursive'
    if length(l1)>1||length(l2)>1||length(l3)>1; error('Need scalar input'); end
    [wm,l]=wigner3jm(l1,l2,l3,-m1,m2,m3);
    % Note that WIGNER0J uses a different algorithm
    [w0,ll]=wigner3jm(l1,l2,l3,0,0,0);
    % Returns an array of all allowable values
    difer(l-ll)
    G=(-1)^m1*sqrt((2*l+1)*(2*l2+1)*(2*l3+1)/4/pi).*wm.*w0;
   case 'gl'
    if length(l1)>1||length(l2)>1||length(l3)>1; error('Need scalar input'); end
    % Construct function names first... remember the complex conjugate
    integrand=inline(sprintf(...
	['(-1)^%i*xlm(%i,%i,acos(x)).*'...
	 'xlm(%i,%i,acos(x)).*'...
	 'xlm(%i,%i,acos(x))'],m1,...
	l1,-m1,l2,m2,l3,m3));
    % And the complex exponential with zero sum of orders evaluates to
    % 2*pi, but to zero if the selection rules aren't applied... like so
    G=gausslegendre([-1 1],integrand,l1+l2+l3)*2*pi*~[-m1+m2+m3];
    l=l1;
   case 'table'
    % Used to be THREEJ, now have ZEROJ
    w0=zeroj(l1,l2,l3,L,[],C0,S0);
    % disp('Done zero-bottom row')
    % The result may be a whole vector
    % Could probably load them all at the same time to save time
    % Put in condition that if they're all zero you don't have to anymore
    if all(m1==0) && all(m2==0) && all(m3==0)
      wm=w0;
    else
      wm=threej(l1,l2,l3,-m1,m2,m3,L,[],C3,S3);
      % disp('Done general-bottom row')
    end
    G=(-1).^m1(:)'.*sqrt((2*l1(:)'+1).*(2*l2(:)'+1).*(2*l3(:)'+1)/4/pi)...
      .*wm.*w0;
    l=l1;
   otherwise
    error('Specify valid method')
  end
elseif strcmp(l1,'demo1')
  % For Table I in Sebilleau (1995)
  % Compare Guseinovs' approach with the Wigner 3j approach (no signs)
  difer(indeks(abs(gaunt(10,10,12,-9,3,-12,'recursive')),'end')-...
	indeks(gaunt(10,10,12,-9,3,-12,'guseinov'),'end')); 
  difer(indeks(abs(gaunt(12,15,5,-2,3,-5,'recursive')),'end')-...
	indeks(gaunt(12,15,5,-2,3,-5,'guseinov'),'end')); 
  difer(indeks(abs(gaunt(20,20,40,-1,-1,0,'recursive')),'end')-...
	indeks(gaunt(20,20,40,-1,-1,0,'guseinov'),'end')); 
  difer(indeks(abs(gaunt(29,29,34,-10,-5,-5,'recursive')),'end')-...
	indeks(gaunt(29,29,34,-10,-5,-5,'guseinov'),'end')); 
elseif strcmp(l1,'demo2')
  % For Table I in Guseinov (1995)
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
    % In Condon-Shortley convention; recursively; with sign
    t0=cputime;
    G1=indeks(gaunt(l1,l2,l3,m1,m2,m3,'recursive'),'end');
    tg1=cputime-t0; t0=cputime;
    % In Condon-Shortley convention; Gauss-Legendre; with sign
    G2=indeks(gaunt(l1,l2,l3,m1,m2,m3,'gl'),'end');
    tg2=cputime-t0; t0=cputime;
    % In Guseinov's stupid phase convention; no sign
    G3=min(abs(gaunt(l1,l2,l3,m1,m2,m3,'guseinov')));
    tg3=cputime-t0;
    if G3>1e3; G3=NaN; end
    comp(ind,:)=[G1 G2 G3];
    cpu(ind,:)=[tg1 tg2 tg3];
  end
  % Take the recursive calculation to be the reference... based on the
  % table it indeed appears to be the most accurate;
  % The first comparison should be down to the sign
  incro=max([L1 ; L2 ; L3]); [k,l]=sort(incro);
  figure(1); clf
  semilogy(k,abs(comp(l,1)-comp(l,2)),'-o')
  hold on
  % The second comparison is regardless of the sign
  semilogy(k,abs(abs(comp(l,1))-abs(comp(l,3))),'r-+')
  % And so is the third comparison
  semilogy(k,abs(abs(comp(l,2))-abs(comp(l,3))),'g-x')
  hold off
  ll=legend('Recursion vs Gauss-Legendre',...
	    'Recursion vs Guseinov',...
	    'Gauss-Legendre vs Guseinov');
  tt=title('Accuracy of Gaunt coefficient determination');
  grid on
  xl=xlabel('maximum degree');
  yl=ylabel('difference');
  
  incro=sum([L1 ; L2 ; L3]); [k,l]=sort(incro);
  figdisp([],1)
  figure(2); clf
  plot(k,cpu(l,1),'-o')
  hold on
  plot(k,cpu(l,2),'r-+')
  plot(k,cpu(l,3),'g-x')
  hold off
  ll=legend('Recursion','Gauss-Legendre','Guseinov');
  tt=title('CPU times of Gaunt coefficient determination');
  grid on
  xl=xlabel('maximum degree');
  yl=ylabel('time');
  figdisp([],2)
elseif strcmp(l1,'demo3')
  difer(indeks(gaunt(10,10,12,-9,3,-12,'recursive'),'end')-...
	gaunt(10,10,12,-9,3,-12,'table')); 
  difer(indeks(gaunt(12,15,5,-2,3,-5,'recursive'),'end')-...
	gaunt(12,15,5,-2,3,-5,'table')); 
  L1=[20 20 20 20];
  L2=[15 15 17  9];
  L3=[13 31 15 15];
  M1=[-3  3 -3 -3];
  M2=[2  -2 -5 -5];
  M3=M1-M2;
  for ind=1:length(L1)
    l1=L1(ind); l2=L2(ind); l3=L3(ind);
    m1=M1(ind); m2=M2(ind); m3=M3(ind);
    difer(indeks(gaunt(l1,l2,l2,m1,m2,m3,'recursive'),'end')-...
	  gaunt(l1,l2,l2,m1,m2,m3,'table')); 
    difer(gaunt(l1,l2,l2,m1,m2,m3,'gl')-...
	  gaunt(l1,l2,l2,m1,m2,m3,'table'));
  end
elseif strcmp(l1,'demo4')
  L1=[20 20 20 20];
  L2=[15 15 17  9];
  L3=[13 19 15 15];
  M1=[-3  3 -3 -3];
  M2=[2  -2 -5 -5];
  M3=M1-M2;
  Gt=gaunt(L1,L2,L3,M1,M2,M3,'table');
  for index=1:length(Gt)
    Gg(index)=gaunt(L1(index),L2(index),L3(index),...
	    M1(index),M2(index),M3(index),'gl');
  end
  difer(Gt-Gg)
else
  error('Specify valid option')
end
