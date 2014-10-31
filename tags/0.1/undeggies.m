function undeggies(ah,varargin)
% UNDEGGIES(ah)
%
% Gets rid of the degree sign - e.g. when you're changing the tick marks
%
% Last modified by fjsimons-at-alum.mit.edu, June 8th, 2004

for index=1:length(ah)
  set(ah(index),'xticklabel',get(ah(index),'xtick'))
  set(ah(index),'yticklabel',get(ah(index),'ytick'))
end
