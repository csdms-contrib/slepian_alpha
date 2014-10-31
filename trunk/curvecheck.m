function curvecheck(lon,lat,pnum)
% CURVECHECK(lon,lat,pnum)
%
% Checks the sense of a closed curve by plotting it
%
% INPUT:
%
% lon,lat   The spherical coordinates of the curve
% pnum      The pause parameter in seconds [default: 0.1]
%
% SEE ALSO: POLY2CW
%
% Last modified by fjsimons-at-alum.mit.edu, 11/23/2011

defval('pnum',0.1)

XY=[lon(:) lat(:)];
hold on
for index=1:length(XY)
  pm(index)=plot(XY(index,1),XY(index,2),'o');
  set(pm(index),'MarkerE','k','MarkerF',[1 1 1]*(index-1)/length(XY))
  title(num2str(index))
  axis([minmax(XY(:,1)) minmax(XY(:,2))])
  pause(pnum)
end
hold off

