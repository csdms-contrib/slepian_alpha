function [bh,th]=label(ah,pos,siz,offset,kees,hitmul,widmul,posxmul,posymul)
% [bh,th]=label(ah,pos,siz,offset,kees,hitmul,widmul,posxmul,posymul)
%
% Uses BOXTEX to put lettered boxes on axis handles
%
% INPUT:
%
% ah         Set of axis handles that you want to label
%            [default: the children of the current figure]
% pos        Position: 'll' lower left [default]
%                      'lr' lower right
%                      'ul' upper left
%                      'ur' upper right
% siz        Font size [default: 12]
% offset     Skip this many letters in the alphabet [default: 0]
% kees       1 Upper case [default]
%            0 or anything else: lower case
%            2 using actual numbers rather than letters
% hitmul     Multiplies height of the box by this factor [default: 1]
% widmul     Multiplies width of the box by this factor [default: 1]
% posxmul    Multiplies x-margin [default: set in BOXTEX]
% posymul    Multiplies y-margin [default: set in BOXTEX]
%
% OUTPUT:
%
% bh       Handles to the boxes
% th       Handles to the text in the boxes
%
% For logarithmic plots, see quick fix in EVENTS_ILL and ALLEN4,
% and more recently, in MERMAID05
%
% Last modified by fjsimons-at-mit.edu, 05/20/2009

defval('ah',flipud(getkids(gcf)))
defval('pos','ll')
defval('siz',12)
defval('offset',0)
defval('kees',0)
defval('hitmul',1)
defval('widmul',1)
% Don't supply defaults here, they come from BOXTEX
defval('posxmul',[])
defval('posymul',[])

for index=1:length(ah)
  [bh(index),th(index)]=...
      boxtex(pos,ah(index),index+offset,siz,kees,hitmul,widmul,...
	     posxmul,posymul);
end
