function varargout=rottick(ah,loc,only)
% t=ROTTICK(ah,loc,only)
% 
% Rotates tick mark labels on a certain axis by abs(90) degrees....
% they better be in cell arrays!! 
%
% INPUT:
%
% ah     Axis handle
% loc    't' or 'b' for the X-tick labels
%        'l' or 'r' for the Y-tick labels
% only   Indices into the tick label vector needing rotation
%
% OUTPUT:
%
% t      Vector with handles to the text objects
%
% Last modified by fjsimons-at-alum.mit.edu, 10/25/2010
 
defval('ah',gca)
defval('loc','b')

% First the specific part
switch loc
  case 'l'
    witch= 'Y'; val=0.9; rot=90; al='left'; posi=1; kaka=1;
  case 'r'
    witch= 'Y'; val=1.02; rot=-90; al='right'; posi=1; kaka=1;
  case 't'
    witch= 'X'; val=1.01; rot=90; al='left'; posi=2; kaka=2;
  case 'b'
    witch= 'X'; val=0.9; rot=-90; al='left'; posi=2; kaka=2;
end

% Now the generic part
axes(ah)
thelab=[witch 'Tick'];
thelbl=[witch 'TickLabel'];
popo=get(ah,thelab);
lab=get(ah,thelbl);
newl=lab;
% This used to be length not size; question is, is lab a string or a
% cell, and we should probably make allowance for this
defval('only',1:size(lab,1))

if length(only)==size(lab,1)
  % This used to be active
  newl=[];
  % only=1:length(lab);
  xval=0;
else
  for inx=only
    if iscell(lab)
      newl{inx}=[];
    else
      newl(inx,:)=' ';
    end
    % Should be an adjustment based on the width of the longest annotation
    xval=0;
  end
end
set(ah,thelbl,newl)

% Take the position from the relevant axis label (not tick; start early!)
pos=get(get(ah,[witch 'Label']),'Position');
y=val*pos(posi)+xval;

for index=only
  if ~iscell(lab)
    if kaka==1; t(index) = text(y,popo(index),num2str(lab(index,:))); end
    if kaka==2; t(index) = text(popo(index),y,num2str(lab(index,:))); end
  else
    if kaka==1; t(index) = text(y,popo(index),lab{index}); end
    if kaka==2; t(index) = text(popo(index),y,lab{index}); end
  end
end
t=t(~~t);
set(t,'Rotation',rot,'HorizontalAlignment',al)

% Generate output
varns={t};
varargout=varns(1:nargout);
