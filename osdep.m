function of=osdep
% of=OSDEP
%
% Returns the value of the read option to be used
% on the LOCAL operating system for files created
% on the LOCAL operating system.
%
% OUTPUT:
%
% of      Said value
%
% Last modified by fjsimons-at-alum.mit.edu, 06/16/2019

defval('of','l')

% Modified by jdsimon-at-princeton.edu on 11/15/2014 for MAC
if strcmp(getenv('OSTYPE'),'linux') || isunix==1
  of= 'l';  
end
if strcmp(getenv('OSTYPE'),'solaris')
  of= 'b';  
end
