function varargout=greenland(res,buf,nearby)   
% XY=GREENLAND(res,buf,nearby)
% GREENLAND(...) % Only makes a plot
%
% Finds the coordinates of Greenland, potentially buffered by some amount.
%
% INPUT:
%
% res      0 The standard, default values
%          N Splined values at N times the resolution
% buf      Distance in degrees that the region outline will be enlarged
%          by BUFFERM, not necessarily integer, possibly negative
%          [default: 0]
% nearby   Subtract the nearby islands of Ellesmere and Baffin
%          from your coordinates [default: 1] or 0
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the continent
%
% Last modified by charig-at-princeton.edu, 07/10/2014
% Last modified by fjsimons-at-alum.mit.edu, 11/23/2011

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

if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis image; grid on
end

% Subtract the nearby island of Ellsemere
if nearby
    % NOTE: if you don't have the nearby stuff made at the desired buffer
    % then this will just end in a recursive error. i.e. You should make
    % Ellesmere first before trying to subtract it.
    XY2 = ellesmereg(10,buf);
    [x,y] = polybool('subtraction',XY(:,1),XY(:,2),XY2(:,1),XY2(:,2));
    XY = [x y];
    if buf>=2
        XY2 = baffing(10,1.5);
        [x,y] = polybool('subtraction',XY(:,1),XY(:,2),XY2(:,1),XY2(:,2));
        XY = [x y];
    end
end



% Prepare optional output
varns={XY};
varargout=varns(1:nargout);

