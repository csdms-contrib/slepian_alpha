function [W,l,l2,l3]=wigner0j(L,l2,l3,xver)
% [W,l,l2,l3]=WIGNER0J(L,l2,l3,xver)
%
% Computes, recursively, the Wigner3j symbol with a bottom row of zero
%
% INPUT:
%
% L        Maximum degree, bandwidth
% l2       Second degree in the Wigner 3j symbol
% l3       Third degree in the Wigner 3j symbol
% xver     Uses orthogonality relation as a numerical check
%
% OUTPUT:
%
% W        Vector with all coefficients where the first degree ranges
%          from 0 up to the bandwidth L; most are zero
% l        The first degrees at which this symbol is evaluated
% l2       Second degree
% l3       Third degree
%
% In practice, only the even degrees between from l=abs(l2-l3)
% through l=l2+l3 need to be computed; all the others are zero. 
% The sign is correct. But note that it doesn't change when the degrees
% are permuted, since their sum is always even.
%
% See also ZEROJ, WIGNER3JM, and WIGNER3J, which this uses
%
% Last modified by fjsimons-at-alum.mit.edu, 05/31/2011

defval('l2',10)
defval('l3',10)
defval('L',l2+l3)
defval('xver',0)

% Check for saved calculations; note: order didn't matter
fname=fullfile(getenv('IFILES'),'WIGNER','0J',...
	       sprintf('W0J-L-%i-%i.mat',l2,l3));
fname2=fullfile(getenv('IFILES'),'WIGNER','0J',...
	       sprintf('W0J-L-%i-%i.mat',l3,l2));

if exist(fname)==2
  load(fname)
%  disp(sprintf('Loaded %s',fname))
elseif exist(fname2)==2
  load(fname2)
%  disp(sprintf('Loaded %s',fname))
else
  % Initialize
  W=repmat(0,1,l2+l3+1);
  
  % Provide the first nonzero term by Gauss-Legendre integration
  W(abs(l2-l3)+1)=(-1)^([abs(l2-l3)+l2+l3]/2)*...
      abs(wigner3j(abs(l2-l3),l2,l3,0,0,0,'gl',1));

  % Recursive evaluation of the other terms; compute them ALL
  % http://functions.wolfram.com/07.39.17.0009.01
  for l=abs(l2-l3)+2:2:l2+l3;
    W(l+1)=(sqrt(-l+l2+l3+2)*sqrt(l-l2+l3-1)*sqrt(l+l2-l3-1)*sqrt(l+l2+l3))/...
	   (sqrt(-l+l2+l3+1)*sqrt(l-l2+l3)  *sqrt(l+l2-l3)  *sqrt(l+l2+l3+1))*...
	   abs(W(l-1))*(-1)^([l+l2+l3]/2);
  end
  % Don't save anymore for small L - got tired of it
  if any([l2 l3]>=100)
    save(fname,'W','l2','l3')
  end
end

l=0:l2+l3;
if xver==1
  difer(sum(W.^2.*(2*l+1))-1,[],[],'WIGNER0J Normalization passed')
end

% Output
l=0:L;
if L<=l2+l3
  % Truncate
  W=W(1:L+1);
else
  % Append
  W=[W repmat(0,1,L-l2-l3)];
end

