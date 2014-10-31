function out = modload(fname,vname)
% out = MODLOAD(fname,vname)
%
% A modified version of MATLAB's native LOAD function. This function takes 
% a filename and the name of a variable within that file and loads it,
% providing it as output. This permits the loading of a variable with a
% specific name to a general name suitable for large-scale processing.
%
% INPUT:
%
% fname   Name of the file where the requested variable is stored
% vname   Name of the requested variable
% 
% OUTPUT:
%
% out     The requested variable from file "fname"
%
% See also EVAL, LOAD, MODSAVE
%
% Last modified by gleggers-at-alumni.princeton.edu, 05/20/2014

% Get list of contents from the file
fNameContents = whos('-file',fname);

% Check if the requested variable is in the file and throw an error if not
if ~ismember(vname,{fNameContents.name})
    error('Requested variable "%s" is not in the storage file "%s"',...
        vname,fname)

% Otherwise, load the variable and assign it to the output 
else
    load(fname,vname);
    eval(sprintf('out = %s;',vname));
end
