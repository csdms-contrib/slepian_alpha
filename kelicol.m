function varargout=kelicol(wit)
% kcol=KELICOL(wit)
%
% Returns, or applies red-white-blue color map designed by Hrafnkell Karason
%
% INPUT:
%
% wit     1 Make center range of this+1 many values really white [default: 0]
%
% OUTPUT:
%
% kcol    The RGB color scale. If no output, simply applies it
%
% Last modified by fjsimons-at-alum.mit.edu, 05/26/2021

defval('wit',0)

% Need to have a centrally located color file KELIM.MAT distributed with
% this code. Use environmental variables or stash locally (gasp!)
try
  load(fullfile(getenv('IFILES'),'COLORMAPS','kelim'))
catch
  load('kelim')
end

if wit>0
  lk=length(kelim);
  kelim(round((lk-wit)/2):round((lk+wit)/2),:)=...
      repmat(1,wit+1,3);
end

% Optional output and action if no output requested
if ~nargout
  colormap(kelim)
end

varns={kelim};
varargout=varns(1:nargout);
