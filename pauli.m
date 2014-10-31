function P=pauli(v,f)
% P=PAULI(v,f)
%
% Returns a Pauli-type matrix useful to generate Toeplitz matrices 
%
% INPUT:
%
% v     A vector of numbers (if not a vector, is made into one)
% f     The column length of the new matrix [default: 2]
%
% OUTPUT:
%
% P     A matrix with the elements of v progressively repeated over the
%       rows of length f
%
% EXAMPLES:
% 
% imagesc(pauli(rand(1,100)*30,50)) % returns a wallpaper pattern
% difer(pauli([2:11],10)-[2:11]')   % is an identity operation
%
% Last modified by fjsimons-at-alum.mit.edu, 07/31/2008

defval('f',2)

v=v(:);

help1=ones(1,(length(v)-f+1)*f);
help1(f+1:f:length(help1))=-(f-2);
help1=cumsum(help1);
help1=reshape(help1,f,length(help1)/f)';

P=v(help1);
