function S=structempty(fields)
% S=STRUCTEMPTY(fields)
%
% Initializes a structure array with empties
%
% INPUT:
%
% fields     A cell array with field name strings
%
% OUTPUT:
%
% S      The initialized structure array
%
% SEE ALSO:
%
% CELLNAN
%
% Last modified by fjsimons-at-alum.mit.edu, 05/28/2010

% Defaults
defval('fields',{'nothing'})
if isstr(fields)
  fields={fields};
end

% Do it exactly like OPTIMSET
structinput=cell(2,length(fields));
% Fields go in first row
structinput(1,:)=fields';
% []'s go in second row
structinput(2,:)={[]};
% Turn it into correctly ordered comma separated list and call struct
S=struct(structinput{:});
