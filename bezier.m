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
% Last modified by charig-at-princeton.edu, April 24th, 2015


if ~any(isnan(XY)) % There are no nans (i.e. one line)
n=size(XY,1);

% Calculate cumulative distance between these points
% This is the "knot vector"
t=cumsum([0 ; grcdist(XY(1:end-1,:),XY(2:end,:))]);
t=t/t(n);

% New, finer, linear distance between the points
tb=linspace(0,1,n*N);

XYb=spline(t,XY',tb)';

else % There was a NaN (i.e. we have a multi segment line)
    [latcell,loncell] = polysplit(XY(:,2),XY(:,1));
    XYb = [];
    for i=1:max(size(latcell))
        newxy = bezier([loncell{i} latcell{i}],N);
        XYb = [XYb; NaN NaN; newxy];
    end
end


