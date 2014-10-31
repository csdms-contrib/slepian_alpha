function fig2print(fignrs,oren)
% FIG2PRINT(fignrs,oren)
%
% Makes the figure window WYSIWYG.
% Also works on a vector of figure numbers
%
% Last modified by fjsimons-at-alum.mit.edu, Jan 5th, 2002

defval('fignrs',gcf)

oren=lower(oren);

for index=1:length(fignrs)
  fignr=fignrs(index);
  switch oren
   case { 'tall','landscape','portrait'}
    orient(fignr,oren)
   case { 'flandscape'};
    pu=get(fignr,'PaperUnits');
    set(fignr,'PaperUnits','normalized')
    set(fignr,'Paperorientation','landscape','paperposition',[0 0 1 1]);
    set(fignr,'PaperUnits',pu)
   case { 'fportrait','ftall'}
    pu=get(fignr,'Paperunits');
    set(fignr,'PaperUnits','normalized')
    set(fignr,'Paperorientation','portrait','paperposition',[0 0 1 1]);
    set(fignr,'PaperUnits',pu)
   otherwise
    error('Illegal option');
  end
  
  ppos = get(fignr,'PaperPosition');
  su = get(fignr,'Units');
  pu = get(fignr,'PaperUnits');  
  set(fignr,'Units',pu);
  spos = get(fignr,'Position');
  set(fignr,'Position',[spos(1)-(ppos(3)-spos(3)) spos(2)-(ppos(4)-spos(4)) ppos(3) ppos(4)])
  set(fignr,'Units',su)
end
