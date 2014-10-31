function varargout=greenland(res,buf)   
% XY=GREENLAND(res,buf)
% GREENLAND(...) % Only makes a plot
%
% Finds the coordinates of Greenland, potentially buffered by some
% amount, and perhaps diminished by Ellesmere and Baffin islands
%
% INPUT:
%
% res      0 The standard, default values
%          N Splined values at N times the resolution
% buf      Distance in degrees that the region outline will be enlarged
%          by BUFFERM, not necessarily integer, possibly negative
%          [default: 0]
% nearby   1 Subtract the glacial coordinates of the nearby islands of
%            Ellesmere and Baffin from your coordinates, or else:
%          0 don't [default]
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the continent
%
% Last modified by charig-at-princeton.edu, 07/10/2014
% Last modified by fjsimons-at-alum.mit.edu, 09/23/2014

defval('res',0)
defval('buf',0)
defval('nearby',1)

% Parameters that make this the region in question
regn='greenland';
c11=[286.7 83.75];
cmn=[349 59.75];
xunt=1:352;

% Do it! Make it, load it, save it
XY=regselect(regn,c11,cmn,xunt,res,buf);

% Subtract the nearby island of Ellsemere
if nearby
  XY2 = ellesmereg(10,buf);
  [x,y] = polybool('subtraction',XY(:,1),XY(:,2),XY2(:,1),XY2(:,2));
  XY = [x y];
  if buf>2
    XY2 = baffing(10,2);
    [x,y] = polybool('subtraction',XY(:,1),XY(:,2),XY2(:,1),XY2(:,2));
    XY = [x y];
  end
end

if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis image; grid on
end

% Prepare optional output
varns={XY};
varargout=varns(1:nargout);

