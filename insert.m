function varargout=insert(in,derin,pos)
% [out,outi]=INSERT(in,derin,pos)
% 
% INPUT:
%
% in       Original data matrix
% derin    Additional elements
% pos      These linear indices in old matrix replaced, rest shifted   
%
% OUTPUT:
%
% out      Row data vector with the new elements in it
% outi     Indexing output so you can be smart the next time
% 
% EXAMPLE:
%
% insert([1 2 3],[11 12],[2 3])
% insert([1 2 3],[11 13],[2 2]) 
%
% Last modified by fjsimons-at-alum.mit.edu, 03/18/2011

% What if nothing needs to happen?
if isempty(pos)
  out=in; outi=NaN;
  vars={out,outi};
  varargout=vars(1:nargout);
  return
end

defval('flag',[])

% Everything is a giant row vector, all indices simply linear
in=in(:)';
derin=derin(:)';
pos=pos(:)';

% If there is one unique last position
if pos(end)==(length(in)+1) && pos(end-1)~=(length(in)+1)
  lest=derin(end);
  derin=derin(1:end-1);
  pos=pos(1:end-1);
  flag=1;
elseif pos(end)==(length(in)+1) && pos(end-1)==(length(in)+1)
  error('Cannot add multiple elements past the range')
end

% Figure out the multiplicity of these position elements
[a,b]=degamini(pos);
slot=gamini(b,b)+1;

% Need to provide the right number of new slots
indices=ones(1,length(in));
indices(pos)=slot;
% Replicate...
out=gamini(in,indices);
% ...then replace the replicates with the additionals
outi=pos+(1:length(pos))-1;
out(outi)=derin; 

% Take care of the last element
if flag==1
  out=[out lest];
end

% Output
vars={out,outi};
varargout=vars(1:nargout);
