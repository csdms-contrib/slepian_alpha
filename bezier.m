function XYb=bezier(XY,N)
% XYb=BEZIER(XY)
%
% An attempt at smoothing a coastline by B-spline fitting.
%
% INPUT:
%
% XY        The set of points, make sure XY(end,:)=XY(1,:)
% N         Number of times this needs to be smoothed
%
% http://www.me.cmu.edu/faculty1/shimada/gm98/project/ivan/project/
% Fujio Yamaguci "Curves and Surfaces in Computer Aided Geometric
% Design", Springer-Verlag, Berlin, 1988 
% http://mathworld.wolfram.com/CubicSpline.html
%
% Last modified by fjsimons-at-alum.mit.edu, June 4rd, 2004

n=size(XY,1);

% Calculate cumulative distance between these points
% This is the "knot vector"
t=cumsum([0 ; grcdist(XY(1:end-1,:),XY(2:end,:))]);
t=t/t(n);

% New, finer, linear distance between the points
tb=linspace(0,1,n*N);

XYb=spline(t,XY',tb)';


