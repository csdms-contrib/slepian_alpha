function x=loadb(file,frmt,of);
% x=LOADB(file,frmt,of)
%
% Loads a binary file.
%
% INPUT:
%
% file    Filename string
% frmt    Format (default: 'float32')
% of      File-opening flag 
%
% Last modified by fjsimons-at-alum.mit.edu, 04/30/2009

defval('of',osdep);
defval('frmt','float32')

fid=fopen(file,'r',of);
x=fread(fid,frmt);
fclose(fid);


