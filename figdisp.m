function varargout=figdisp(name,ext,opt,act,form,convo)
% FIGDISP(name,ext,opt,act,form,convo)
% [fname,pstring]=FIGDISP(...)
%
% Suggests print command for figures.
% Assumes environment variable 'EPS' is set.
% Use PAINTERS to override, sometimes.
%
% INPUT:
%
% name          Filename root, no extension (default: calling function)
% ext           An optional additional extension (number or string)
% opt           An option string, e.g. '-zbuffer'
% act           1 Actually print the figure
%               0 Don't print the figure [default]
%               2 Makes a PDF from an EPS which it removes
% form          A graphics format [pdf, png, etc: default: epsc]
% convo         A graphics converter string [default: 'ps2raster -Tf']
%               such as 'epstopdf' or 'ps2raster' or ...
%
% OUTPUT:
%
% fname         The full file name
% pstring       The plot string
%
% Last modified by fjsimons-at-alum.mit.edu, 07/01/2016

[p,n]=star69;

defval('name',n)
defval('ext',[])
defval('opt',[])
defval('act',0)
defval('convo','ps2raster -Tf')
defval('form','epsc')

% Calls itself and cleans up afterward
if act==2
  [fname,pstring]=figdisp(name,ext,opt,1,'epsc')
  system(sprintf('degs %s.eps',fname));
  system(sprintf('%s %s.eps',convo,fname));
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

% Optional output
varns={fname,pstring};
varargout=varns(1:nargout);
