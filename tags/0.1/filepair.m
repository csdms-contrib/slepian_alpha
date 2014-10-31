function [names,tits]=filepair(ddir)
% [names,tits]=FILEPAIR(ddir)
%
% Picks up a pair of filenames from a directory.
%
% INPUT:
%
% ddir        is the directory with filenames
%
% OUTPUT:
%
% names        a cell array with pairs of filenames
% tits         a cell array with pairs of filenames suitable for titles
%
% Last modified by fjsimons-at-alum.mit.edu, 05/24/2010

defval('ddir','/home/fjsimons/MERMAID/09102004/EVENTS/')

files=ls2cell(ddir);

for index=1:2:length(files)-1
  names{(index-1)/2+1}={files{index} files{index+1}};
  tits{(index-1)/2+1}=...
      nounder(sprintf('%s / %s',files{index},files{index+1}));
end

% If the number of files is odd must add last one
if mod(length(files),2)
  names{(index+3)/2}={files{index+2}};
  tits{(index+3)/2}=nounder(sprintf('%s',files{index+2}));
end




