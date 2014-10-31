function Z=gauss2(X,Y,stdn,stdm)
% Z=GAUSS2(X,Y,stdn,stdm)
%
% Returns 2D gaussian curve with standard
% deviations 'stdn' and 'stdm' in 'X' ans 'Y' directions,
% respectively, when 'X' and 'Y' come out of meshgrid.
%
% EXAMPLE:
%
% [X,Y]=meshgrid(-10:0.1:10);
% subplot(211)
% surf(X,Y,gauss2(X,Y,2,4)) ; shading flat ; xlabel('x') ; ylabel('y')
% subplot(212)
% pcolor(X,Y,gauss2(X,Y,2,4)) ; shading flat
% axis ij ; axis image ;  xlabel('x') ; ylabel('y')

% Last modified by fjsimons-at-mit.edu, Dec 19th, 2001

defval('X',meshgrid(linspace(-3,3,256)))
defval('Y',meshgrid(linspace(-3,3,256))')

Z=gauss(X,stdn).*gauss(Y,stdm);
