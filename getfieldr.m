function varargout = getfieldr(S,getFields,fname)
% GETFIELDR(S,getFields,fname)
% [field1,field2,field3,...] = GETFIELDR(S,getFields,fname)
% [allFields,allFieldpaths] = GETFIELDR(S,[],fname)
% 
% GETFIELDR is an improvement on MATLAB's native GETFIELD function in that
% it will retrieve fields of any depth in a structure.  The information 
% can be retrieved from the structure in several ways depending on the 
% nature of the "S" and "getFields" input variables and how many outputs
% are requested. Before retrieving fields, the function tests for
% uniqueness in the fieldnames of the structure because if one of the
% requested fields could refer to multiple fields within the structure,
% then the function cannot know which of the fields was actually intended.
%
% INPUT:
%
% S             Structure from which fields are to be retrieved
%                 1. If a structure, fields are extracted directly
%                 2. If a string, then a structure of that name is loaded
%                    from the provided savefile "fname"
% getFields     Fields to be retrieved
%                 1. If empty/not provided, no fields are retrieved. Rather,
%                    the "allFields" and "allPaths" are output 
%                 2. If it is the string ':', then all field values of "S"
%                    are requested and outputted
%                 3. For any other string or cell array of strings, those
%                    strings are treated as the requested fields of "S"
% fName         A file where the structure may be stored
% 
% OUTPUT:
%
% <none>        Field values assigned directly in the calling workspace with 
%               variable names given by their in-structure fieldnames
% field1,...    Field values passed out in the same order as "getFields"
% allFields     Cell array of all fields
% allPaths      Cell array of fieldpaths to the "allFields" fields
%
% See also EVAL, FIELDNAMESR, NONUNIQUE, UNIQUE
%
% Last modified by gleggers-at-alumni.princeton.edu, 05/20/2014
% Last modified by fjsimons-at-alum.mit.edu, 10/07/2014

% If necessary, load the requested structure from its savefile
if ischar(S)
    S=modload(fname,S);
end

% Get the fieldpaths to all fields (terminal or otherwise)
allPaths = fieldnamesr(S,'full','prefix');

% Truncate those fieldpaths to just the fieldnames
for i=1:length(allPaths)
    pl=strfind(allPaths{i},'.');
    allFields{i}=allPaths{i}(pl(end)+1:end);
end

% If a "getFields" is not passed, pass out the listing of
% all terminal fields and their fieldpaths
if ~exist('getFields','var') || isempty(getFields)
  varargout={allFields allPaths'};
  return
  
  % If "getFields" is a string...
elseif ischar(getFields)
  % If it is a ":", then all fields are being requested
  if strcmp(getFields,':')
    getFields=allFields;
    
    % Otherwise, assume just a single field is being requested,
    % but place it in a cell array for later convenience
  else
    getFields={getFields};
  end
end

% Get listing of all unique fieldnamess in the structure
uniqueFields=unique(unique(allFields));

% If the number of unique fieldnames is any less than the total number of
% fields, then the structure has duplicate fieldnames, to be addressed 
if length(uniqueFields) ~= length(allFields)
  % Get complementary list of nonunique fieldnames
    nonuniqueFields=nonunique(allFields);
    
    % If one of the requested fields is a nonunique field, the function
    % cannot distinguish between the same-named fields, so throw an error
    if any(ismember(getFields,nonuniqueFields))
        error(['Certain requested fields correspond to duplicate field '...
               'names in the structure. Please check and rename fields.'])
    
	% If none of the requested fields are the nonunique fieldnames, then
	% the function will be fine, but throw a warning so the user is aware
    else
        warning(['The structure from which fields are being requested'...
                 'has duplicate field names, but none requested']);
    end
end

% Depending on the number of outputs, the requested field values are 
% assigned differently
switch nargout
  % For no requested output, assign requested field values directly into
  % the calling workspace with fieldnames giving the new variable names
 case 0
  for i=1:length(getFields)
    assignin('caller',getFields{i},...
	     S.(getFields{i}))
%	     eval(allPaths{strcmp(getFields{i},allFields)}))
  end
  
  % For any other number of output, collect the field values into a 
  % cell array for output 
 otherwise
  for i=1:length(getFields) % To nargout?
    eval(sprintf('vars{i}=%s;',...
	     S.(getFields{i})));
    %		 allPaths{strcmp(getFields{i},allFields)}));
  end
  % Output
  varargout=vars(1:nargout);
end
