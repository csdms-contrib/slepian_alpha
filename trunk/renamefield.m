function S=renamefield(S,name1,name2)
% S=RENAMEFIELD(S,name1,name2)
%
% INPUT:
%
% S         The structure array
% name1     The old field name string
% name2     The new field name string
%
% OUTPUT:
%
% S         The new structure array
%
% SEE ALSO:
%
% STRUCT2STR
%
% Last modified by fjsimons-at-alum.mit.edu, 05/29/2013

[S.(name2)]=S.(name1);
S=rmfield(S,name1);


