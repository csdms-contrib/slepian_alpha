function [p,ex,ey]=errorxy2(X,Y,DXp,DXm,DYp,DYm)
% [p,ex,ey]=ERRORXY(X,Y,DXp,DXm,DYp,DYm)
%
% Plots error bars on a plot
%
% Like ERRORXY but with asymmetric error bars.
%
% See also ERRORXY, ERRORBAR
%
% Last modified by fjsimons-at-alum.mit.edu, 04/05/2007

defval('DXp',[])
defval('DXm',[])
defval('DYp',[])
defval('DYm',[])

X=X(:)'; DXp=DXp(:)';  DXm=DXm(:)';
Y=Y(:)'; DYp=DYp(:)';  DYm=DYm(:)';

p=plot(X,Y,'ks');
hold on

if ~isempty(DYm)
  YDY=[Y-DYm ; Y+DYp];
  ey=plot([X ; X],YDY,'k');
else
  ey=[];
end

if ~isempty(DXm)
  XDX=[X-DXm ; X+DXp];
  ex=plot(XDX,[Y ; Y],'k');
else
  ex=[];
end

hold off
