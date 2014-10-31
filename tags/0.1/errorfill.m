function varargout=errorfill(X,Y,DY)
% ERRORFILL(X,Y,DY)
% [rp,e68,e95]=ERRORFILL(X,Y,DY)
% 
% Plots a function with shaded error zone and, optionally,
% returns handles to them.
%
% Last modified by fjsimons-at-alum.mit.edu, 25.2.2005

X=X(:);
Y=Y(:);
DY=DY(:);

s8=[0.7 0.7 0.7];
n5=[0.9 0.9 0.9];

e95=fill([X ; flipud(X)],[Y+2*DY ; flipud(Y-2*DY)],n5); hold on
e68=fill([X ; flipud(X)],[Y+DY ; flipud(Y-DY)],s8);
rp=plot(X,Y,'w');
hold off

varnam=[{'rp'},{'e68'},{'e95'}];
for index=1:nargout
  varargout{index}=eval(varnam{index});
end
