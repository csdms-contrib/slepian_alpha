function ind=dindeks(i,d,sais)
% ind=DINDEKS(i,d,[n1 n2 n3])
%
% Returns the indices to the ith plane perpendicular to the dth dimension
% of a threedimensional array. Return everything as a linear array.
%
% EXAMPLE:
%
% vec=rand(28,13,78); chex=12;
% difer(indeks(rindeks(vec,chex),':')-vec(dindeks(chex,1,size(vec))))
% difer(indeks(squeeze(kindeks(vec,chex)),':')-vec(dindeks(chex,2,size(vec))))
% difer(indeks(squeeze(tindeks(vec,chex)),':')-vec(dindeks(chex,3,size(vec))))
%
% EXAMPLE:
%
% 
%
% SEE ALSO: RINDEKS, KINDEKS, TINDEKS
%
% Last modified by fjsimons-at-alum.mit.edu, 03/11/2010

defval('d',1)

% Put in the right dimensions
n1=sais(1);
n2=sais(2);
n3=sais(3);

% Check, in SV's code, with the statement
% return dnode->data[i + j*(dnode->len1) + k*(dnode->len1)*(dnode->len2)];
switch d
 case 1
  ind=repmat(i-n1+[1:n2]*n1,n3,1)'+repmat([0:n3-1]'*n2*n1,1,n2)';
 case 2
  ind=repmat(n1*(i-1)+[1:n1],n3,1)'+repmat([0:n3-1]'*n2*n1,1,n1)';
 case 3
  ind=repmat(n1*n2*(i-1)+[1:n1],n2,1)'+repmat([0:n2-1]'*n1,1,n1)';
end

% The test works for here:
% difer(squeeze(rindeks(vec,chex))-vec(dindeks(chex,1,size(vec))))
% difer(squeeze(kindeks(vec,chex))-vec(dindeks(chex,2,size(vec))))
% difer(squeeze(tindeks(vec,chex))-vec(dindeks(chex,3,size(vec))))

% Now the output
ind=ind(:);
