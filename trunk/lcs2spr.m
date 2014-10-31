function A=lcs2spr(Lin,Cin,Sin,m,n)
% mat=LCS2SPR(L,C,S,m,n)
%
% Gives a sparse (mxn) matrix representation of the matrix
% represented by 
% 'L' with cumulative entries per row
% 'C' column number of entries
% 'S' element values
%
% Input may be variables or filenames. '.bin' will be appended.
%
% SEE ALSO: LCSD, LCSG
%
% Last modified by fjsimons-at-alum.mit.edu, 08/20/2007

if isstr(Lin); L=loadb([Lin '.bin'],'int32'); else;  L=Lin(:); end
if isstr(Cin); C=loadb([Cin '.bin'],'int32'); else; C=Cin(:); end
if isstr(Sin); S=loadb([Sin '.bin'],'float32'); else; S=Sin(:); end

% Now have local variables named L, S and C.
if ~all(size(S)==size(C))
  error('Sizes of S and C must be equal ')
end

if ~L(end)==size(C,1)
  error('No valid LCS system')
end

LL=cumsum(full(sparse([1;L(1:end-1)+1],1,ones(1,length(L)),L(end),1)));
A=sparse(LL,C,S,m,n,length(S));
