function filesize=fsize(filename)
% filesize=FSIZE(filename)
%
% Gets the size of a file, in bytes
%
% Last modified by fjsimons-at-alum.mit.edu, 07/31/2008

fid=fopen(filename,'r','l');
fseek(fid,0,1);
filesize=ftell(fid);
%fseek(fid,0,-1); % No need to rewind if file is closed!
fclose(fid);
