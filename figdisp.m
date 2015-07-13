function varargout=figdisp(name,ext,opt,act,form)
% FIGDISP(name,ext,opt,act,form)
% [fname,pstring]=FIGDISP(...)
%
% Suggests print command for figures.
% Assumes environment variable 'EPS' is set.
% Use PAINTERS to override, sometimes.
%
% INPUT:
%
% name          Filename (default: calling function)
% ext           Extension (number or string)
% opt           An option string, e.g. '-zbuffer'
% act           1 Actually print the figure
%               0 Don't print the figure [default]
%               2 Makes a PDF from an EPS which it removes
% form          A graphics format [pdf, png, etc: default: epsc]
%
% OUTPUT:
%
% fname         The full file name
% pstring       The plot string
%
% Last modified by fjsimons-at-alum.mit.edu, 07/13/2015

[p,n]=star69;

defval('name',n)
defval('ext',[])
defval('opt',[])
defval('act',0)
defval('form','epsc')

% Calls itself and cleans up afterward
if act==2
  [fname,pstring]=figdisp(name,ext,opt,1,'epsc');
  system(sprintf('degs %s.eps',fname));
  system(sprintf('epstopdf %s.eps',fname));
  system(sprintf('rm -f %s.eps',fname));
  varns={fname,pstring};
  varargout=varns(1:nargout);
  return
end

if ~isstr(ext)
  ext=num2str(ext);
end

if ~isempty(ext)
  fname=fullfile(getenv('EPS'),sprintf('%s_%s',name,ext));
else
  fname=fullfile(getenv('EPS'),sprintf('%s',name));
end

if ~isempty(opt)
  pstring=sprintf('print(''-d%s'',''%s'',''%s'')',form,opt,fname);
else
  pstring=sprintf('print(''-d%s'',''%s'')',form,fname);
end

disp(pstring)

if act>=1
  eval(pstring)
end

varns={fname,pstring};
varargout=varns(1:nargout);
