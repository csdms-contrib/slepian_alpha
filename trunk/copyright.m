function str=copyright
% str=COPYRIGHT
%
% Produces a copyright string used by ID. Depending on the version, need
% to fix the PostScript later on using FIXCOPY.PL 
%
% Last modified by fjsimons-at-alum.mit.edu, 10/22/2012

deet=date;
deet(abs(deet)==45)= ' ';

% Sometimes it works, sometimes it doesn't
str=[char(169),' Frederik J Simons, Princeton University, ',deet];
str=[169,' Frederik J Simons, Princeton University, ',deet];
str=[' Frederik J Simons, Princeton University, ',deet];



