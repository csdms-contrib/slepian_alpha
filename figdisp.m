function varargout=figdisp(name,ext,opt,act,form,convo)
% FIGDISP(name,ext,opt,act,form,convo)
% [fname,pstring]=FIGDISP(...)
%
% Suggests print command for figures. Assumes environment variable 'EPS' is set.
% Use PAINTERS to override, sometimes. Might need EPSTOPDF or PS2RASTER.
%
% INPUT:
%
% name          Filename, preferably without graphics extension or any
%               dots in the filename, really [default: calling function]
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
% fname         The full file name, but definitely no graphics extension
% pstring       The plot string
%
% SEE ALSO:
%
% DEGS
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)
% Last modified by fjsimons-at-alum.mit.edu, 03/20/2020

[p,n]=star69;

defval('name',n)
defval('ext',[])
defval('opt',[])
defval('act',0)
defval('form','epsc')
defval('convo','ps2raster -Tf')

% Calls itself and cleans up afterward
if act==2
  [fname,pstring]=figdisp(name,ext,opt,1,'epsc');
  % PRINT will keep ANY extension, defined as any period in the filename,
  % but also it will make a GRAPHICS extension if the filename had no
  % periods in it...  so if you GIVE it an extension coming in, careful!
  % What normally comes out of FIGDISP as 'fname' has no extension
  if ~any(fname==46) ; xtra='.eps'; else xtra=[]; end
  % You need to have DEGS available in your search path
  system(sprintf('degs %s%s',fname,xtra));
  % You need to have EPSTOPDF or PS2RASTER available in your search path
  system(sprintf('%s %s%s',convo,fname,xtra));
  % Now, depending on the behavior of the conversion, fix ITS extension;
  % epstopdf and ps2raster appear to STRIP ANY extension reading from the
  % BACK and substitute the right extension in (overwriting a file that
  % may have come in with the wrong extension!), so keep on fixing,
  % remove the original UNLESS it had the .pdf at the end at this point
  if ~all(fname(end-3:end)=='.pdf')
    system(sprintf('rm -f %s%s',fname,xtra));
  end
  % Shortcut and get out
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
