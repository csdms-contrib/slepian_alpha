function varargout=kelicol(wit)
% kcol=KELICOL(wit)
%
% Returns red and blue color map designed by H. Karason
%
% INPUT:
%
% wit   1 Make center range of this+1 many values really white
%
% Last modified by fjsimons-at-alum.mit.edu, 09/12/2012

defval('wit',0)

if ~exist('kelim')
  try
    load(fullfile(getenv('IFILES'),'COLORMAPS','kelim'))
  catch
    load('kelim')
  end
end

if wit>0
  lk=length(kelim);
  kelim(round((lk-wit)/2):round((lk+wit)/2),:)=...
      repmat(1,wit+1,3);
end

if ~nargout
  colormap(kelim)
end

fake=NaN;
varns={kelim,fake,fake};
varargout=varns(1:nargout);

