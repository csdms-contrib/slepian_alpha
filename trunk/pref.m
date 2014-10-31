function prefix=pref(string,delim)
% prefix=PREF(string,delim)
%
% Finds the prefix or head in a string, say a filename.
% The delimiter is a '.' (period) by default, but may be set.
%
% INPUT:
%
% string        The string to be parsed
% delim         The string that separates prefix from suffix
%
% See also SUF
%
% Last modified by fjsimons-at-alum.mit.edu, 05/24/2010

defval('delim','.');

% Something fishy is going on here... I forget why
% Something to do with having more than one delimiter?
% See hsfigX
%[prefix,suffix]=strtok(fliplr(string),delim);

prefix=strtok(string,delim);

%prefix=fliplr(suffix);
%prefix=prefix(1:end-length(delim));


