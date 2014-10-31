function varargout=serre(ah,resjo,ww)
% old=SERRE(ah,resjo,ww)
%
% INPUT:
% 
% ah        Vector of axis handles; if matrix goes row by row
% resjo     Fraction of their distance by which plots are moved [default: 1/2]
% ww        'down' or 'across' [default]
%
% OUTPUT:
%
% old       The original positions of the axis handles
%
% EXAMPLE:
%
% clf
% [ah,ha,H]=krijetem(subnum(3,3)); ho1=getpos(ah);
% o1=serre(H,1/2,'across'); difer(ho1-o1); ho2=getpos(ha);
% o2=serre(H',1/2,'down'); difer(ho2-o2);
% unserre(ah,o1); difer(ho1-getpos(ah));
%
% Last modified by fjsimons-at-alum.mit.edu, 04/02/2009

defval('resjo',1/2)
defval('ww','across')

if prod(size(ah))~=length(ah)
  % It's a matrix of handles
  for index=1:size(ah,1)
    o{index}=serre(ah(index,:),resjo,ww);
  end
  varns={cat(1,o{:})};
else
  % Get the axis positions of a vector of handles
  vdist=getpos(ah(end-1),2)...      
	-getpos(ah(end),2)...
	-getpos(ah(end),4);
  hdist=getpos(ah(end),1)...
	-getpos(ah(end-1),1)...
	-getpos(ah(end-1),3);
  for index=1:length(ah)-1
    switch ww
     case 'down'
      [n(index,:),o(index,:)]=...
	  movev(ah(index),-vdist*resjo*(length(ah)-index));
     case 'across'
      [n(index,:),o(index,:)]=...
	  moveh(ah(index),hdist*resjo*(length(ah)-index));
     otherwise
      error('Specify valid option')
    end
  end
  varns={[o ; getpos(ah(end))]};
end

% Provided output if so desired
varargout=varns(1:nargout);
