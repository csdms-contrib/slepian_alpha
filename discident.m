function [disc,varargout]=discident(depths)
% disc=discident(depths)
% [disc,remoef]=discident(depths)
% [disc,remoef,discval]=discident(depths)
% [disc,remoef,discval,avma,indi]=discident(depths)
%
% For 'depths' with twice repeated values, yields a vector
% with +1 and -1 for the discontinuities; +1 being the
% shallower value (so for depths, not radii)
%
% For 'depths' with 3X repeated values, yields a vector
% with +1, 0 and -1 for the discontinuities; +1 being
% at the lowest depth, and the rest NaN's.
%
% For 3p discontinuities, can return 'remoef' so that 
% 'depths'('remoef') is 'depths' without the middle one.
% For 2p, 'remoef' just includes everything.
% 
% For continuous functions, returns empties.
%
% For 2p discons, makes an averaging matrix (to get one
% average value over the discontinuous points)
%
% The values can be obtained from:
% depth(logical(abs(discident(depth)))) or
% depth(find(abs(discident(depth))))
%
% Pairs of indices are returned in indi
%
% Example
%
% depths=[1 1 1 2 2 2 3 4 5 5 5 6 7 8 8 8]'
% [depths discident(depths)]
% [a,b]=discident(depths);
% [depths(b) discident(depths(b))]
%
% depths=[1 1 2 2 3 4 5 5 6 7 8 8]'
% [depths discident(depths)]
%
% Last modified by fjsimons-at-alum.mit.edu, October 22th, 2002
%
% NEEDS TO BE MODIFIED AGAIN; NEW MATLAB DOES NOT 
% RETURN ZERO FOR ~NAN AS IN CASE 3

flag=0;
depths=depths(:);

if depths(1)>depths(end)
  depths=flipud(depths);
  flag=1;
end

% Two- or three point discontinuities?
[a,fold]=degamini(depths);
discval=a(~~(fold-1));
switch max(fold)
  case 1
    if nargout>=2
      disc=[];
      varargout{1}=[];
    end    
  case  2
    disc=[~diff(depths) ;  0]-[0 ; ~diff(depths)];
    if nargout>=2
      varargout{1}=logical(ones(size(depths)));
    end    
case 3
    disc=repmat(NaN,length(depths),1);
    disc(find(~diff(depths))+1)=-1;
    disc(find(~diff(depths)))=1;
    disc(find(~diff(disc))+1)=0;   
    if nargout>=2
      varargout{1}=~~disc;
    end    
end

if nargout>=3
  varargout{2}=discval;
end

if flag==1
  depths=flipud(depths);
  disc=flipud(disc);
end

if nargout>=4
  % Also make an averaging matrix in case they are 2p
  if ~any(isnan(disc)) & ~isempty(disc)
    avma=zeros(length(disc),length(discval));
    avma(find(disc==1)+(length(disc)*[0:length(discval)-1]'))=1;
    avma(find(disc==-1)+(length(disc)*[0:length(discval)-1]'))=-1;
  else
    avma=[];
  end
  if flag==1
    avma=flipud(avma);
  end
  varargout{3}=abs(avma);
  if nargout>=5
    [jk,indi]=intersect(depths,discval);
    if depths(indi)==depths(indi-1)
      varargout{4}=[indi(:)-1 indi(:)];
    end
  end
end

if ~exist('disc')
  disc=[];
end
