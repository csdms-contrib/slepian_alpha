function [p,ex,ey]=errorxy(X,Y,DX,DY)
% [p,ex,ey]=ERRORXY(X,Y,DX,DY)
%
% Plots error bars on a plot
% 
% EXAMPLE:
%
% X=-10:10;
% Y=3*X.^2+4*X+rand(1,length(X));
% DX=1.1*rand(1,length(X));
% DY=20.3*rand(1,length(X));
% errorxy(X,Y,DX,DY)
%
% See also ERRORXY2, ERRORBAR
%
% Last modified by fjsimons-at-alum.mit.edu, 04/05/2007

defval('DX',[])
defval('DY',[])

X=X(:)'; DX=DX(:)';
Y=Y(:)'; DY=DY(:)';

if isempty(DX); DX=zeros(size(X)) ; end
if isempty(DY); DY=zeros(size(Y)) ; end

p=plot(X,Y,'ks');
hold on

if ~isempty(DY)
  YDY=[Y-DY ; Y+DY];
  ey=plot([X ; X],YDY,'k');
else
  ey=0;
end

if ~isempty(DX)
  XDX=[X-DX ; X+DX];
  ex=plot(XDX,[Y ; Y],'k');
else
  ex=[];
end

hold off
