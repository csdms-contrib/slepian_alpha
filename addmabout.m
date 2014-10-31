function a=addmabout(L,l,m)
% a=ADDMABOUT(L,l,m)
%
% Find an index in a block spherical harmonic array belonging to a
% particular degree l and order m and a bandwidth L. No input testing.
%
% INPUT:
%
% L         Spherical harmonic bandwidth  [may be vector]
% l         Spherical harmonic degree (0 -> L) [may be vector]
% m         Spherical harmonic order (-l -> l) [default: 0]
%
% OUTPUT:
%
% a         The running index in the degree/order block arrays
%
% EXAMPLE:
%
% L=round(rand*10);
% [EM,EL,mz,blkm,dblk]=addmout(L); EM=EM(blkm); EL=EL(blkm);
% for l=0:L; for m=-l:l; 
%   difer(EM(addmabout(L,l,m))-m);
%   difer(EL(addmabout(L,l,m))-l);
% end; end
%
% Last modified by fjsimons-at-alum.mit.edu, 12/28/2006

% WATCH OUT, THIS DOES NOT KNOW ANYTHING USEFUL ABOUT L AND M
% I.E. YOU ALWAYS GET AN ANSWER - SHOULD BUILD THIS IN!

defval('m',0)

a=(l-abs(m)+1)+(abs(m)>0)*(L+1+(abs(m)-1)*(2*L-abs(m)+2)+...
			 (m>0)*(L-m+1));
