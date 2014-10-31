function plotoncube3(vname,tipe,kind,vscale,colmap,phil)
% PLOTONCUBE3(vname,tipe,kind,vscale,colmap,phil)
%
% Plots ALL DEPTHS of gridded whole-cubed-sphere data. 
% A wrapper for PLOTONCUBE2.
%
% INPUT:
%
% vname      A string identifying the vector [default: 'vecx.petsc'], or
%            A structure array with all of the information in there already
% tipe       See the option list in PLOTONCUBE
% kind       See the option list in PLOTONCUBE
% vscale     See the option list in PLOTONCUBE
% colmap     See the option list in PLOTONCUBE
% phil       See the option list in PLOTONCUBE
% 
% Last modified by fjsimons-at-alum.mit.edu, 04/21/2010

defval('vname','vecx.petsc')
defval('tipe','2D')
defval('kind',1)
defval('vscale',[-1 1])
defval('colmap','kelicol')
defval('phil',1)

% Hardwire the 37 in here, fix later, if it ain't right it will break
% anyway. 
more off
for index=1:37
  plotoncube2(vname,index,tipe,kind,vscale,colmap,phil);
  pause
end

% Now suggest the booklet command
if isstruct(vname); vname=inputname(1); end
if any(abs(vname)==46); vname=pref(vname); end
allname=sprintf('plotoncube2_%s_*','vname');
mergename=sprintf('plotoncube2_%s_all',vname);
disp(sprintf(...
    'gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=%s.pdf -dBATCH %s',...
    mergename,allname'))




