function writeb(x,name,varargin)
% WRITEB(x,name,tipe)
%
% Writes binary data to a new file.
%
% INPUT:
%
% x       A vector with the data you want to write out
% name    The name of the new file you want to create
% tipe    The binary format [default: float32]
%
% Last modified by fjsimons-at-alum.mit.edu, 09/12/2007

if nargin==3
  tipe=varargin{1};
else
  tipe= 'float32';
end

fid=fopen(name,'wb');
fwrite(fid,x,tipe);
fclose(fid);
