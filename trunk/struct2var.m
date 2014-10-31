function struct2var(S)
% STRUCT2VAR(S)
%
% Pulls all of the fieldnames out of a structure that only has one depth
% level and assigns them as named variables in to the calling workspace.
%
% OUTPUT:
%
% S       Structure name
% 
% SEE ALSO:
%
% DEFSTRUCT, GETFIELDR
%
% Last modified by fjsimons-at-alum.mit.edu, 10/07/2014

% Get all the field names  
fields=fieldnames(S);

% Extract the variables explicitly from this structure
for ind=1:length(fields)
  % Obsolete but save for instructional uses
  % eval(sprintf('%s=params.(fields{ind});',fields{ind}))
  assignin('caller',fields{ind},S.(fields{ind}))
end

